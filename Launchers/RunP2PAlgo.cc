#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <csv.h>
#include <routingkit/nested_dissection.h>

#include "Algorithms/CTL/BalancedTopologyCentricTreeHierarchy.h"
#include "Algorithms/CTL/TruncatedTreeLabelling.h"
#include "Algorithms/CTL/CTLMetric.h"
#include "Algorithms/CTL/CTLQuery.h"
#include "Algorithms/CCH/CCH.h"
#include "Algorithms/CCH/CCHMetric.h"
#include "Algorithms/CCH/EliminationTreeQuery.h"
#include "Algorithms/CTNR/CTNR.h"
#include "Algorithms/CH/CH.h"
#include "Algorithms/CH/CHQuery.h"
#include "Algorithms/Dijkstra/BiDijkstra.h"
#include "Algorithms/Dijkstra/Dijkstra.h"
#include "DataStructures/Graph/Attributes/LatLngAttribute.h"
#include "DataStructures/Graph/Attributes/LengthAttribute.h"
#include "DataStructures/Graph/Attributes/TravelTimeAttribute.h"
#include "DataStructures/Graph/Graph.h"
#include "DataStructures/Labels/BasicLabelSet.h"
#include "DataStructures/Labels/ParentInfo.h"
#include "DataStructures/Partitioning/SeparatorDecomposition.h"
#include "DataStructures/Partitioning/nested_strict_dissection.h"
#include "Tools/CommandLine/CommandLineParser.h"
#include "Tools/StringHelpers.h"
#include "Tools/Timer.h"
#include <ctlsa/road_network.h>
#include "Tools/CommandLine/ProgressBar.h"

inline void printUsage() {
    std::cout <<
              "Usage: RunP2PAlgo -a CH         -o <file> -g <file>\n"

              "       RunP2PAlgo -a CCH        -o <file> -g <file> [-b <balance>]\n"
              "       RunP2PAlgo -a CTL        -o <file> -g <file> [-b <balance>]\n"
              "       RunP2PAlgo -a CTNR       -o <file> -g <file> [-b <balance>]\n\n"

              "       RunP2PAlgo -a CCH-custom -o <file> -g <file> -s <file> [-n <num>]\n"
              "       RunP2PAlgo -a CTL-custom -o <file> -g <file> -s <file> [-n <num>]\n\n"

              "       RunP2PAlgo -a Dij        -o <file> -g <file> -d <file>\n"
              "       RunP2PAlgo -a Bi-Dij     -o <file> -g <file> -d <file>\n"
              "       RunP2PAlgo -a CH         -o <file> -h <file> -d <file>\n"
              "       RunP2PAlgo -a CCH-Dij    -o <file> -g <file> -d <file> -s <file>\n"
              "       RunP2PAlgo -a CCH-tree   -o <file> -g <file> -d <file> -s <file>\n"
              "       RunP2PAlgo -a CTL        -o <file> -g <file> -d <file> -s <file>\n"
              "       RunP2PAlgo -a CTNR       -o <file> -g <file> -d <file> -s <file>\n\n"

              "Runs the preprocessing, customization or query phase of various point-to-point\n"
              "shortest-path algorithms, such as Dijkstra, bidirectional search, CH, CCH, CTL, and CTNR.\n\n"

              "  -l                use physical lengths as metric (default: travel times)\n"
              "  -no-stall         do not use the stall-on-demand technique\n"
              "  -a <algo>         run algorithm <algo>\n"
              "  -b <balance>      balance parameter in % for nested dissection (default: 30)\n"
              "  -n <num>          run customization <num> times (default: 1000)\n"
              "  -g <file>         input graph in binary format\n"
              "  -s <file>         separator decomposition of input graph\n"
              "  -h <file>         weighted contraction hierarchy\n"
              "  -d <file>         file that contains OD pairs (queries)\n"
              "  -o <file>         place output in <file>\n"
              "  -help             display this help and exit\n";
}

// Some helper aliases.
using VertexAttributes = VertexAttrs<LatLngAttribute>;
using EdgeAttributes = EdgeAttrs<LengthAttribute, TravelTimeAttribute>;
using InputGraph = StaticGraph<VertexAttributes, EdgeAttributes>;
using LabelSet = BasicLabelSet<0, ParentInfo::NO_PARENT_INFO>;

// The query algorithms.
using Dij = Dijkstra<InputGraph, TravelTimeAttribute, LabelSet>;
using BiDij = BiDijkstra<Dij>;
template<bool useStalling>
using CCHDij = CHQuery<LabelSet, useStalling>;
using CCHTree = EliminationTreeQuery<LabelSet>;

// Writes the header line of the output CSV file.
template<typename AlgoT>
inline void writeHeaderLine(std::ofstream &out, AlgoT &) {
    out << "distance,query_time" << '\n';
}

template<>
inline void writeHeaderLine(std::ofstream &out, CTNRQuery<InputGraph> &) {
    out << "distance,query_time,mode" << '\n';
}

// Writes a record line of the output CSV file, containing statistics about a single query.
template<typename AlgoT>
inline void writeRecordLine(std::ofstream &out, AlgoT &algo, const int, const int64_t elapsed) {
    out << algo.getDistance() << ',' << elapsed << '\n';
}

template<>
inline void writeRecordLine(std::ofstream &out, Dij &algo, const int dst, const int64_t elapsed) {
    out << algo.getDistance(dst) << ',' << elapsed << '\n';
}

// Runs the specified P2P algorithm on the given OD pairs.
template<typename AlgoT, typename T>
inline void runQueries(AlgoT &algo, const std::string &demand, std::ofstream &out, T translate) {
    Timer timer;
    int src, dst, rank;
    using TrimPolicy = io::trim_chars<>;
    using QuotePolicy = io::no_quote_escape<','>;
    using OverflowPolicy = io::throw_on_overflow;
    using CommentPolicy = io::single_line_comment<'#'>;
    io::CSVReader<3, TrimPolicy, QuotePolicy, OverflowPolicy, CommentPolicy> demandFile(demand);
    const auto ignore = io::ignore_extra_column | io::ignore_missing_column;
    demandFile.read_header(ignore, "origin", "destination", "dijkstra_rank");
    const auto hasRanks = demandFile.has_column("dijkstra_rank");
    if (hasRanks) out << "dijkstra_rank,";
    writeHeaderLine(out, algo);
    while (demandFile.read_row(src, dst, rank)) {
        src = translate(src);
        dst = translate(dst);
        timer.restart();
        algo.run(src, dst);
        const auto elapsed = timer.elapsed<std::chrono::nanoseconds>();
        if (hasRanks) out << rank << ',';
        writeRecordLine(out, algo, dst, elapsed);
        if constexpr (std::is_same_v<AlgoT, CTNRQuery<InputGraph>>) {
            out.seekp(-1, std::ios_base::cur); // overwrite newline
            out << ',' << algo.getLastMode() << '\n';
        }
    }
}

// Invoked when the user wants to run the query phase of a P2P algorithm.
inline void runQueries(const CommandLineParser &clp) {
    const auto useLengths = clp.isSet("l");
    const auto noStalling = clp.isSet("no-stall");
    const auto algorithmName = clp.getValue<std::string>("a");
    const auto graphFileName = clp.getValue<std::string>("g");
    const auto sepFileName = clp.getValue<std::string>("s");
    const auto chFileName = clp.getValue<std::string>("h");
    const auto demandFileName = clp.getValue<std::string>("d");
    auto outputFileName = clp.getValue<std::string>("o");

    static constexpr uint64_t BYTES_PER_MB = 1 << 20;

    // Open the output CSV file.
    if (!endsWith(outputFileName, ".csv"))
        outputFileName += ".csv";
    std::ofstream outputFile(outputFileName);
    if (!outputFile.good())
        throw std::invalid_argument("file cannot be opened -- '" + outputFileName + ".csv'");

    if (algorithmName == "Dij") {

        // Run the query phase of Dijkstra's algorithm.
        std::ifstream graphFile(graphFileName, std::ios::binary);
        if (!graphFile.good())
            throw std::invalid_argument("file not found -- '" + graphFileName + "'");
        InputGraph graph(graphFile);
        graphFile.close();
        if (useLengths)
            FORALL_EDGES(graph, e)graph.travelTime(e) = graph.length(e);

        outputFile << "# Graph: " << graphFileName << '\n';
        outputFile << "# OD pairs: " << demandFileName << '\n';

        Dij algo(graph);
        runQueries(algo, demandFileName, outputFile, [](const int v) { return v; });

    } else if (algorithmName == "Bi-Dij") {

        // Run the query phase of bidirectional search.
        std::ifstream graphFile(graphFileName, std::ios::binary);
        if (!graphFile.good())
            throw std::invalid_argument("file not found -- '" + graphFileName + "'");
        InputGraph graph(graphFile);
        graphFile.close();
        if (useLengths)
            FORALL_EDGES(graph, e)graph.travelTime(e) = graph.length(e);

        outputFile << "# Graph: " << graphFileName << '\n';
        outputFile << "# OD pairs: " << demandFileName << '\n';

        InputGraph reverseGraph = graph.getReverseGraph();
        BiDij algo(graph, reverseGraph);
        runQueries(algo, demandFileName, outputFile, [](const int v) { return v; });

    } else if (algorithmName == "CH") {

        // Run the query phase of CH.
        std::ifstream chFile(chFileName, std::ios::binary);
        if (!chFile.good())
            throw std::invalid_argument("file not found -- '" + chFileName + "'");
        CH ch(chFile);
        chFile.close();

        outputFile << "# CH: " << chFileName << '\n';
        outputFile << "# OD pairs: " << demandFileName << '\n';

        if (noStalling) {
            CCHDij<false> algo(ch);
            runQueries(algo, demandFileName, outputFile, [&](const int v) { return ch.rank(v); });
        } else {
            CCHDij<true> algo(ch);
            runQueries(algo, demandFileName, outputFile, [&](const int v) { return ch.rank(v); });
        }

    } else if (algorithmName == "CCH-Dij") {

        // Run the Dijkstra-based query phase of CCH.
        std::ifstream graphFile(graphFileName, std::ios::binary);
        if (!graphFile.good())
            throw std::invalid_argument("file not found -- '" + graphFileName + "'");
        InputGraph graph(graphFile);
        graphFile.close();

        std::ifstream sepFile(sepFileName, std::ios::binary);
        if (!sepFile.good())
            throw std::invalid_argument("file not found -- '" + sepFileName + "'");
        SeparatorDecomposition sepDecomp;
        sepDecomp.readFrom(sepFile);
        sepFile.close();

        CCH cch;
        cch.preprocess(graph, sepDecomp);
        CCHMetric metric(cch, useLengths ? &graph.length(0) : &graph.travelTime(0));
        const auto minCH = metric.buildMinimumWeightedCH();

        outputFile << "# Graph: " << graphFileName << '\n';
        outputFile << "# Separator: " << sepFileName << '\n';
        outputFile << "# OD pairs: " << demandFileName << '\n';

        if (noStalling) {
            CCHDij<false> algo(minCH);
            runQueries(algo, demandFileName, outputFile, [&](const int v) { return minCH.rank(v); });
        } else {
            CCHDij<true> algo(minCH);
            runQueries(algo, demandFileName, outputFile, [&](const int v) { return minCH.rank(v); });
        }

    } else if (algorithmName == "CCH-tree") {

        // Run the elimination-tree-based query phase of CCH.
        std::ifstream graphFile(graphFileName, std::ios::binary);
        if (!graphFile.good())
            throw std::invalid_argument("file not found -- '" + graphFileName + "'");
        InputGraph graph(graphFile);
        graphFile.close();

        std::ifstream sepFile(sepFileName, std::ios::binary);
        if (!sepFile.good())
            throw std::invalid_argument("file not found -- '" + sepFileName + "'");
        SeparatorDecomposition sepDecomp;
        sepDecomp.readFrom(sepFile);
        sepFile.close();

        CCH cch;
        cch.preprocess(graph, sepDecomp);
        CCHMetric metric(cch, useLengths ? &graph.length(0) : &graph.travelTime(0));
        const auto minCH = metric.buildMinimumWeightedCH();

        outputFile << "# Graph: " << graphFileName << '\n';
        outputFile << "# Separator: " << sepFileName << '\n';
        outputFile << "# OD pairs: " << demandFileName << '\n';

        CCHTree algo(minCH, cch.getEliminationTree());
        outputFile << "# Memory usage CCH: " << (cch.sizeInBytes()) / BYTES_PER_MB << " MB" << '\n';
        outputFile << "# Memory usage CCHMetric: " << (metric.sizeInBytes()) / BYTES_PER_MB << " MB" << '\n';
        outputFile << "# Memory usage EliminationTreeQuery: " << (algo.sizeInBytes()) / BYTES_PER_MB << " MB" << '\n';
        outputFile << "# Memory usage total: "
                   << (cch.sizeInBytes() + metric.sizeInBytes() + algo.sizeInBytes()) / BYTES_PER_MB << " MB" << '\n';

        runQueries(algo, demandFileName, outputFile, [&](const int v) { return minCH.rank(v); });

    } else if (algorithmName == "CTL") {

        // Run truncated tree labelling (CTL) queries
        std::ifstream graphFile(graphFileName, std::ios::binary);
        if (!graphFile.good())
            throw std::invalid_argument("file not found -- '" + graphFileName + "'");
        InputGraph graph(graphFile);
        graphFile.close();

        std::ifstream sepFile(sepFileName, std::ios::binary);
        if (!sepFile.good())
            throw std::invalid_argument("file not found -- '" + sepFileName + "'");
        SeparatorDecomposition sepDecomp;
        sepDecomp.readFrom(sepFile);
        sepFile.close();

        CCH cch;
        cch.preprocess(graph, sepDecomp);

        BalancedTopologyCentricTreeHierarchy treeHierarchy;
        treeHierarchy.preprocess(graph, sepDecomp);

        using CTLLabelSet = std::conditional_t<CTL_SIMD_LOGK == 0,
                BasicLabelSet<0, ParentInfo::NO_PARENT_INFO>,
                SimdLabelSet<CTL_SIMD_LOGK, ParentInfo::NO_PARENT_INFO>>;
        using LabellingT = TruncatedTreeLabelling<CTLLabelSet::K, CTLLabelSet::KEEP_PARENT_EDGES>;
        LabellingT ctl(treeHierarchy);
        ctl.init();

        CTLMetric<LabellingT, CTLLabelSet, CTL_USE_PERFECT_CUSTOMIZATION> metric(treeHierarchy, cch, useLengths ? &graph.length(0) : &graph.travelTime(0));
        metric.buildCustomizedCTL(ctl);

        outputFile << "# Graph: " << graphFileName << '\n';
        outputFile << "# Separator: " << sepFileName << '\n';
        outputFile << "# OD pairs: " << demandFileName << '\n';

        CTLQuery<CTLMetric<LabellingT, CTLLabelSet, CTL_USE_PERFECT_CUSTOMIZATION>::SearchGraph, LabellingT, CTLLabelSet> algo(treeHierarchy, metric.upwardGraph(),
                                                                      metric.downwardGraph(), metric.upwardWeights(),
                                                                      metric.downwardWeights(), ctl);

        outputFile << "# Memory usage CCH: " << (cch.sizeInBytes()) / BYTES_PER_MB << " MB" << '\n';
        outputFile << "# Memory usage TreeHierarchy: " << (treeHierarchy.sizeInBytes()) / BYTES_PER_MB << " MB" << '\n';
        outputFile << "# Memory usage Labelling: " << (ctl.sizeInBytes()) / BYTES_PER_MB << " MB" << '\n';
        outputFile << "# Memory usage CTLMetric: " << (metric.sizeInBytes()) / BYTES_PER_MB << " MB" << '\n';
        outputFile << "# Memory usage CTLQuery: " << (algo.sizeInBytes()) / BYTES_PER_MB << " MB" << '\n';
        outputFile << "# Memory usage total: " <<
                   (cch.sizeInBytes() + treeHierarchy.sizeInBytes() + ctl.sizeInBytes() + metric.sizeInBytes() +
                    algo.sizeInBytes()) / BYTES_PER_MB << " MB" << '\n';
        runQueries(algo, demandFileName, outputFile, [&](const int v) { return cch.getRanks()[v]; });

    } else if (algorithmName == "CTNR") {

        // Run customizable transit node routing (CTNR) queries
        std::ifstream graphFile(graphFileName, std::ios::binary);
        if (!graphFile.good())
            throw std::invalid_argument("file not found -- '" + graphFileName + "'");
        InputGraph graph(graphFile);
        graphFile.close();

        std::ifstream sepFile(sepFileName, std::ios::binary);
        if (!sepFile.good())
            throw std::invalid_argument("file not found -- '" + sepFileName + "'");
        SeparatorDecomposition sepDecomp;
        sepDecomp.readFrom(sepFile);
        sepFile.close();

        // Build CTNR
        CTNR<InputGraph> ctnr(sepDecomp, 5); // Use top k levels as transit nodes
        ctnr.preprocess(graph);

        // Customize CTNR
        std::vector<int32_t> edgeWeights(graph.numEdges());
        FORALL_EDGES(graph, e) {
            edgeWeights[e] = useLengths ? graph.length(e) : graph.travelTime(e);
        }
        
        ctnr.customize(edgeWeights.data());
        
        outputFile << "# Graph: " << graphFileName << '\n';
        outputFile << "# OD pairs: " << demandFileName << '\n';
        // outputFile << "# Memory usage total: " << ctnr.sizeInBytes() / BYTES_PER_MB << " MB" << '\n';

        // Use generic runQueries with CTNRQuery; pass CCH rank IDs to the algo
        const auto &metric = ctnr.getMetric();
        CTNRQuery<InputGraph> algo(metric);
        runQueries(algo, demandFileName, outputFile, [&](const int v) { return ctnr.getCCH().getRanks()[v]; });

    } else {

        throw std::invalid_argument("invalid P2P algorithm -- '" + algorithmName + "'");

    }
}


// Invoked when the user wants to run the preprocessing or customization phase of a P2P algorithm.
inline void runPreprocessing(const CommandLineParser &clp) {
    const auto useLengths = clp.isSet("l");
    const auto imbalance = clp.getValue<int>("b", 30);
    const auto numCustomRuns = clp.getValue<int>("n", 1000);
    const auto algorithmName = clp.getValue<std::string>("a");
    const auto graphFileName = clp.getValue<std::string>("g");
    const auto sepFileName = clp.getValue<std::string>("s");
    auto outputFileName = clp.getValue<std::string>("o");

    // Read the input graph.
    std::ifstream graphFile(graphFileName, std::ios::binary);
    if (!graphFile.good())
        throw std::invalid_argument("file not found -- '" + graphFileName + "'");
    InputGraph graph(graphFile);
    graphFile.close();
    if (useLengths)
        FORALL_EDGES(graph, e)graph.travelTime(e) = graph.length(e);

    std::cout << "Graph has " << graph.numVertices() << " vertices and " << graph.numEdges() << " edges" << std::endl;

    if (algorithmName == "CH") {

        // Run the preprocessing phase of CH.
        if (!endsWith(outputFileName, ".ch.bin"))
            outputFileName += ".ch.bin";
        std::ofstream outputFile(outputFileName, std::ios::binary);
        if (!outputFile.good())
            throw std::invalid_argument("file cannot be opened -- '" + outputFileName);

        std::cout << "Constructing CH for " << graphFileName << "... " << std::flush;
        CH ch;
        Timer timer;
        ch.preprocess<TravelTimeAttribute>(graph);
        const auto preprocessTime = timer.elapsed<std::chrono::microseconds>();
        ch.writeTo(outputFile);

        std::cout << " finished (" << preprocessTime << " microseconds)." << std::endl;

    } else if (algorithmName == "CCH") {

        // Run the preprocessing phase of CCH.
        std::cout << "Constructing separator decomposition for CCH for " << graphFileName << "... " << std::flush;
        Timer timer;
        if (imbalance < 0)
            throw std::invalid_argument("invalid imbalance -- '" + std::to_string(imbalance) + "'");

        // Convert the input graph to RoutingKit's graph representation.
        std::vector<float> lats(graph.numVertices());
        std::vector<float> lngs(graph.numVertices());
        std::vector<unsigned int> tails(graph.numEdges());
        std::vector<unsigned int> heads(graph.numEdges());
        FORALL_VERTICES(graph, u) {
            lats[u] = graph.latLng(u).latInDeg();
            lngs[u] = graph.latLng(u).lngInDeg();
            FORALL_INCIDENT_EDGES(graph, u, e) {
                tails[e] = u;
                heads[e] = graph.edgeHead(e);
            }
        }

        // Compute a separator decomposition for the input graph.
        const auto fragment = RoutingKit::make_graph_fragment(graph.numVertices(), tails, heads);
        auto computeSep = [&](const RoutingKit::GraphFragment &fragment) {
            const auto cut = inertial_flow(fragment, imbalance, lats, lngs);
            return derive_separator_from_cut(fragment, cut.is_node_on_side);
        };
        const auto decomp = compute_separator_decomposition(fragment, computeSep);

        // Convert the separator decomposition to our representation.
        SeparatorDecomposition sepDecomp;
        for (const auto &n: decomp.tree) {
            SeparatorDecomposition::Node node;
            node.leftChild = n.left_child;
            node.rightSibling = n.right_sibling;
            node.firstSeparatorVertex = n.first_separator_vertex;
            node.lastSeparatorVertex = n.last_separator_vertex;
            sepDecomp.tree.push_back(node);
        }
        sepDecomp.order.assign(decomp.order.begin(), decomp.order.end());


        const auto preprocessTime = timer.elapsed<std::chrono::microseconds>();
        std::cout << " finished (" << preprocessTime << " microseconds)." << std::endl;

        if (!endsWith(outputFileName, ".sep.bin"))
            outputFileName += ".sep.bin";
        std::ofstream outputFile(outputFileName, std::ios::binary);
        if (!outputFile.good())
            throw std::invalid_argument("file cannot be opened -- '" + outputFileName);
        sepDecomp.writeTo(outputFile);
    } else if (algorithmName == "CCH-custom") {

        // Run the customization phase of CCH.
        std::ifstream sepFile(sepFileName, std::ios::binary);
        if (!sepFile.good())
            throw std::invalid_argument("file not found -- '" + sepFileName + "'");
        SeparatorDecomposition decomp;
        decomp.readFrom(sepFile);
        sepFile.close();

        if (!endsWith(outputFileName, ".csv"))
            outputFileName += ".csv";
        std::ofstream outputFile(outputFileName);
        if (!outputFile.good())
            throw std::invalid_argument("file cannot be opened -- '" + outputFileName + ".csv'");
        outputFile << "# Graph: " << graphFileName << '\n';
        outputFile << "# Separator: " << sepFileName << '\n';

        Timer timer;
        CCH cch;
        cch.preprocess(graph, decomp);
        const auto preprocessTime = timer.elapsed<std::chrono::microseconds>();
        outputFile << "# Preprocess time (for given sepdecomp): " << preprocessTime << " microseconds.\n";

        outputFile << "basic_customization,perfect_customization,construction,total_time\n";
        int64_t basicCustom, perfectCustom, construct, tot;
        for (auto i = 0; i < numCustomRuns; ++i) {
//            {
//                CCHMetric metric(cch, &graph.travelTime(0));
//                timer.restart();
//                metric.customize();
//                basicCustom = timer.elapsed<std::chrono::microseconds>();
//                timer.restart();
//                metric.runPerfectCustomization();
//                perfectCustom = timer.elapsed<std::chrono::microseconds>();
//            }
//            {
//                CCHMetric metric(cch, &graph.travelTime(0));
//                timer.restart();
//                metric.buildMinimumWeightedCH();
//                tot = timer.elapsed<std::chrono::microseconds>();
//            }

            timer.restart();
            CCHMetric metric(cch, &graph.travelTime(0));
            metric.buildMinimumWeightedCH<Timer>(basicCustom, perfectCustom, construct);
            tot = timer.elapsed<std::chrono::microseconds>();

            outputFile << basicCustom << ',' << perfectCustom << ',' << construct << ',' << tot << '\n';
        }

    } else if (algorithmName == "CTL") {
        // Run the preprocessing phase of CTL.
        std::cout << "Constructing separator decomposition with strict dissection for CTL for " << graphFileName
                  << "... " << std::flush;
        Timer timer;
        if (imbalance < 0)
            throw std::invalid_argument("invalid imbalance -- '" + std::to_string(imbalance) + "'");

        // Convert the input graph to RoutingKit's graph representation.
        std::vector<float> lats(graph.numVertices());
        std::vector<float> lngs(graph.numVertices());
        std::vector<unsigned int> tails(graph.numEdges());
        std::vector<unsigned int> heads(graph.numEdges());
        FORALL_VERTICES(graph, u) {
            lats[u] = graph.latLng(u).latInDeg();
            lngs[u] = graph.latLng(u).lngInDeg();
            FORALL_INCIDENT_EDGES(graph, u, e) {
                tails[e] = u;
                heads[e] = graph.edgeHead(e);
            }
        }

        // Compute a strict bisection separator decomposition for the input graph.
        auto fragment = RoutingKit::make_graph_fragment(graph.numVertices(), tails, heads);
        auto computeCut = [&](const RoutingKit::GraphFragment &fragment) {
            return inertial_flow(fragment, imbalance, lats, lngs);
        };
        RoutingKit::BitVector all(fragment.node_count(), true);
        RoutingKit::BitVector none = ~all;
        const auto decomp = compute_separator_decomposition_with_strict_dissection(std::move(fragment), computeCut,
                                                                                   std::move(none), std::move(all));

        // Convert the separator decomposition to our representation.
        SeparatorDecomposition sepDecomp;
        for (const auto &n: decomp.tree) {
            SeparatorDecomposition::Node node;
            node.leftChild = n.left_child;
            node.rightSibling = n.right_sibling;
            node.firstSeparatorVertex = n.first_separator_vertex;
            node.lastSeparatorVertex = n.last_separator_vertex;
            sepDecomp.tree.push_back(node);
        }
        sepDecomp.order.assign(decomp.order.begin(), decomp.order.end());


        const auto preprocessTime = timer.elapsed<std::chrono::microseconds>();
        std::cout << " finished (" << preprocessTime << " microseconds)." << std::endl;

        if (!endsWith(outputFileName, ".strict_bisep.bin"))
            outputFileName += ".strict_bisep.bin";
        std::ofstream outputFile(outputFileName, std::ios::binary);
        if (!outputFile.good())
            throw std::invalid_argument("file cannot be opened -- '" + outputFileName);
        sepDecomp.writeTo(outputFile);

    } else if (algorithmName == "CTNR") {
        std::cout << "CTNR preprocessing (using existing separator decomposition) for " << graphFileName
                  << "... " << std::flush;
        Timer timer;
        std::string sepFileName = graphFileName;
        size_t lastDot = sepFileName.find_last_of('.');
        if (lastDot != std::string::npos) {
            sepFileName = sepFileName.substr(0, lastDot);
        }
        sepFileName += ".strict_bisep.bin";
        
        // 读取现有的分隔分解文件
        std::ifstream sepFile(sepFileName, std::ios::binary);
        if (!sepFile.good()) {
            std::cout << "Separator decomposition file not found: " << sepFileName << std::endl;
            std::cout << "Please run CTL preprocessing first to generate separator decomposition." << std::endl;
            return;
        }
        
        SeparatorDecomposition sepDecomp;
        sepDecomp.readFrom(sepFile);
        sepFile.close();
        
        std::cout << "Loaded separator decomposition with " << sepDecomp.tree.size() << " nodes" << std::endl;

        // Build CTNR
        CTNR<InputGraph> ctnr(sepDecomp, 5); // Use top 5 levels as transit nodes
        ctnr.preprocess(graph);
        
        std::cout << "CTNR preprocessing completed with " << ctnr.getTransitNodes().size() 
                  << " transit nodes" << std::endl;

        if (!endsWith(outputFileName, ".ctnr.bin"))
            outputFileName += ".ctnr.bin";
        std::ofstream outputFile(outputFileName, std::ios::binary);
        if (!outputFile.good())
            throw std::invalid_argument("file cannot be opened -- '" + outputFileName);
        // CTNR currently has no serialization; keep placeholder to match other modes
        
        const auto preprocessTime = timer.elapsed<std::chrono::microseconds>();
        std::cout << " finished (" << preprocessTime << " microseconds)." << std::endl;

    } else if (algorithmName == "CTL-custom") {

        // Run the customization phase of CCH.
        std::ifstream sepFile(sepFileName, std::ios::binary);
        if (!sepFile.good())
            throw std::invalid_argument("file not found -- '" + sepFileName + "'");
        SeparatorDecomposition decomp;
        decomp.readFrom(sepFile);
        sepFile.close();

        if (!endsWith(outputFileName, ".csv"))
            outputFileName += ".csv";
        std::ofstream outputFile(outputFileName);
        if (!outputFile.good())
            throw std::invalid_argument("file cannot be opened -- '" + outputFileName + ".csv'");
        outputFile << "# Graph: " << graphFileName << '\n';
        outputFile << "# Separator: " << sepFileName << '\n';

        Timer timer;
        CCH cch;
        cch.preprocess(graph, decomp);
        BalancedTopologyCentricTreeHierarchy treeHierarchy;
        treeHierarchy.preprocess(graph, decomp);
        using CTLLabelSet = std::conditional_t<CTL_SIMD_LOGK == 0,
                BasicLabelSet<0, ParentInfo::NO_PARENT_INFO>,
                SimdLabelSet<CTL_SIMD_LOGK, ParentInfo::NO_PARENT_INFO>>;
        using LabellingT = TruncatedTreeLabelling<CTLLabelSet::K, CTLLabelSet::KEEP_PARENT_EDGES>;
        LabellingT ctl(treeHierarchy);
        ctl.init();
        const auto preprocessTime = timer.elapsed<std::chrono::microseconds>();
        outputFile << "# Preprocess time (for given sepdecomp): " << preprocessTime << " microseconds.\n";

        outputFile << "cch_customization,ctl_customization,total_time\n";
        timer.restart();
        int cchCustom, ctlCustom, tot;
        for (auto i = 0; i < numCustomRuns; ++i) {
            {
                CCHMetric metric(cch, &graph.travelTime(0));
                timer.restart();
                if constexpr (CTL_USE_PERFECT_CUSTOMIZATION) {
                    metric.buildMinimumWeightedCH();
                } else {
                    metric.customize();
                }
                cchCustom = timer.elapsed<std::chrono::microseconds>();
            }
            {
                CTLMetric<LabellingT, CTLLabelSet, CTL_USE_PERFECT_CUSTOMIZATION> metric(treeHierarchy, cch, &graph.travelTime(0));
                timer.restart();
                metric.buildCustomizedCTL(ctl);
                tot = timer.elapsed<std::chrono::microseconds>();
            }
            ctlCustom = tot - cchCustom;
            outputFile << cchCustom << ',' << ctlCustom << ',' << tot << '\n';
        }
    } else if (algorithmName == "CTLSACCH-custom") {

        // TODO: Allow using CCH from CTLSA again.
        throw std::invalid_argument("CTLSACCH-custom is not supported at the moment.");
//      // CTLSACCH can only deal with undirected graphs.
//      // Graph can be considered undirected if every edge exists both ways and has same travel time both ways.
//      FORALL_VALID_EDGES(graph, u, e) {
//              KASSERT(graph.template get<TravelTimeAttribute>(e) != TravelTimeAttribute::defaultValue());
//              const auto eBack = graph.uniqueEdgeBetween(graph.edgeHead(e), u);
//              KASSERT(eBack >= 0 && eBack < graph.numEdges());
//              KASSERT(graph.template get<TravelTimeAttribute>(e) == graph.template get<TravelTimeAttribute>(eBack));
//          }
//
//        // Convert the input graph to CTLSA representation.
//        ctlsa::road_network::Graph ctlsaGraph;
//      ctlsaGraph.resize(graph.numVertices());
//      FORALL_VALID_EDGES(graph, u, e) {
//              // CTLSA graph node IDs start at 1
//              ctlsaGraph.add_edge(u + 1, graph.edgeHead(e) + 1, graph.template get<TravelTimeAttribute>(e), false);
//          }
//
//
//      if (!endsWith(outputFileName, ".csv"))
//          outputFileName += ".csv";
//      std::ofstream outputFile(outputFileName);
//      if (!outputFile.good())
//          throw std::invalid_argument("file cannot be opened -- '" + outputFileName + ".csv'");
//      outputFile << "# Graph: " << graphFileName << '\n';
//      outputFile << "setup,customization,total\n";
//
//      // (Do not) contract degree 1 nodes
//      std::vector<ctlsa::road_network::Neighbor> closest;
//      ctlsaGraph.contract(closest, false);
//
//      // Build balanced tree hierarchy
//      std::vector<ctlsa::road_network::CutIndex> cutIndex;
//      static constexpr double CUT_BALANCE = 0.2;
//      static constexpr size_t LEAF_SIZE_THRESHOLD = 0; // CCH only works with THETA = 0.
//      ctlsaGraph.create_cut_index(cutIndex, CUT_BALANCE, LEAF_SIZE_THRESHOLD);
//
//      // Reset graph to original form before contractions
//      ctlsaGraph.reset();
//
//      // Initialize shortcut graph and labels
//      ctlsa::road_network::ContractionHierarchy ch;
//      ctlsaGraph.initialize(ch, cutIndex, closest);
//
//      Timer timer;
//      int64_t customTime, setupTime;
//      for (auto i = 0; i < numCustomRuns; ++i) {
//
//          timer.restart();
//
//          // Customize CCH and HL with new metric on edges:
//          ctlsaGraph.reset(ch);
//          std::vector<ctlsa::road_network::Edge> edges;
//          ctlsaGraph.get_edges(edges);
//          for (auto &e: edges) {
//              // CTLSA graph node IDs start at 1
//              const auto tail = e.a - 1;
//              const auto head = e.b - 1;
//              const auto eInInputGraph = graph.uniqueEdgeBetween(tail, head);
//              KASSERT(eInInputGraph >= 0 && eInInputGraph < graph.numEdges());
//              e.d = graph.travelTime(eInInputGraph);
//          }
//
//          setupTime = timer.elapsed<std::chrono::microseconds>();
//            timer.restart();
//
//          ctlsaGraph.customise_shortcut_graph(ch, edges);
//
//          customTime = timer.elapsed<std::chrono::microseconds>();
//          outputFile << setupTime << "," << customTime << "," << (setupTime + customTime) << '\n';
//      }
    } else if (algorithmName == "CTLSA-custom") {

        // CTLSA can only deal with undirected graphs.
        // Graph can be considered undirected if every edge exists both ways and has same travel time both ways.
        FORALL_VALID_EDGES(graph, u, e) {
                KASSERT(graph.template get<TravelTimeAttribute>(e) != TravelTimeAttribute::defaultValue());
                const auto eBack = graph.uniqueEdgeBetween(graph.edgeHead(e), u);
                KASSERT(eBack >= 0 && eBack < graph.numEdges());
                KASSERT(graph.template get<TravelTimeAttribute>(e) == graph.template get<TravelTimeAttribute>(eBack));
                if (eBack < 0 && eBack >= graph.numEdges())
                    throw std::invalid_argument(
                            "graph is not undirected -- '" + graphFileName + "': Reverse edge of (" +
                            std::to_string(u) + ", " + std::to_string(graph.edgeHead(e)) + ") does not exist.");
                if (graph.template get<TravelTimeAttribute>(e) != graph.template get<TravelTimeAttribute>(eBack))
                    throw std::invalid_argument(
                            "graph is not undirected -- '" + graphFileName + "': Travel time of edge (" +
                            std::to_string(u) + ", " + std::to_string(graph.edgeHead(e)) +
                            ") does not match reverse edge.");
            }

        // Convert the input graph to CTLSA representation.
        ctlsa::road_network::Graph ctlsaGraph;
        ctlsaGraph.resize(graph.numVertices());
        FORALL_VALID_EDGES(graph, u, e) {
                // CTLSA graph node IDs start at 1
                ctlsaGraph.add_edge(u + 1, graph.edgeHead(e) + 1, graph.template get<TravelTimeAttribute>(e), false);
            }

        if (!endsWith(outputFileName, ".csv"))
            outputFileName += ".csv";
        std::ofstream outputFile(outputFileName);
        if (!outputFile.good())
            throw std::invalid_argument("file cannot be opened -- '" + outputFileName + ".csv'");
        outputFile << "# Graph: " << graphFileName << '\n';
        outputFile << "setup,cch_customization,ctlsa_customization,total\n";

        std::cout << "Preprocessing..." << std::flush;
        // (Do not) contract degree 1 nodes
        std::vector<ctlsa::road_network::Neighbor> closest;
        ctlsaGraph.contract(closest, false);

        // Build balanced tree hierarchy
        std::vector<ctlsa::road_network::CutIndex> cutIndex;
        static constexpr double CUT_BALANCE = 0.2;
        static constexpr size_t LEAF_SIZE_THRESHOLD = 0;
        ctlsaGraph.create_cut_index(cutIndex, CUT_BALANCE, LEAF_SIZE_THRESHOLD);

        // Reset graph to original form before contractions
        ctlsaGraph.reset();

        // Initialize shortcut graph and labels
        ctlsa::road_network::ContractionHierarchy ch;
        ctlsaGraph.initialize(ch, cutIndex, closest);
        ctlsa::road_network::ContractionIndex ci(cutIndex, closest);

        std::cout << " done." << std::endl;

        std::cout << "Running customization " << numCustomRuns << " times... " << std::flush;
        Timer timer;
        int64_t setupTime, cchCustomTime, ctlsaCustomTime;
        ProgressBar progressBar(numCustomRuns);
        for (auto i = 0; i < numCustomRuns; ++i) {

            timer.restart();

            // Customize CCH and HL with new metric on edges:
            ctlsaGraph.reset(ch, ci);
            std::vector<ctlsa::road_network::Edge> edges;
            ctlsaGraph.get_edges(edges);
            for (auto &e: edges) {
                // CTLSA graph node IDs start at 1
                const auto tail = e.a - 1;
                const auto head = e.b - 1;
                const auto eInInputGraph = graph.uniqueEdgeBetween(tail, head);
                KASSERT(eInInputGraph >= 0 && eInInputGraph < graph.numEdges());
                e.d = graph.travelTime(eInInputGraph);
            }

            setupTime = timer.elapsed<std::chrono::microseconds>();
            timer.restart();

            ctlsaGraph.customise_shortcut_graph(ch, ci, edges);

            cchCustomTime = timer.elapsed<std::chrono::microseconds>();
            timer.restart();

            ctlsaGraph.customise_hub_labelling(ch, ci);

            ctlsaCustomTime = timer.elapsed<std::chrono::microseconds>();
            outputFile << setupTime << "," << cchCustomTime << "," << ctlsaCustomTime << ","
                       << (setupTime + cchCustomTime + ctlsaCustomTime) << '\n';
            ++progressBar;
        }
        progressBar.finish();
        std::cout << " done." << std::endl;
    } else {

        throw std::invalid_argument("invalid P2P algorithm -- '" + algorithmName + "'");

    }
}

int main(int argc, char *argv[]) {
    std::cout << "Using " << NUM_THREADS << " threads for parallel regions in RunP2PAlgo." << std::endl;
    omp_set_num_threads(NUM_THREADS);
    try {
        CommandLineParser clp(argc, argv);
        if (clp.isSet("help"))
            printUsage();
        else if (clp.isSet("d"))
            runQueries(clp);
        else
            runPreprocessing(clp);
    } catch (std::exception &e) {
        std::cerr << argv[0] << ": " << e.what() << std::endl;
        std::cerr << "Try '" << argv[0] << " -help' for more information." << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
