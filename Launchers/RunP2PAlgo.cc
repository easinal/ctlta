#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

#include <csv.h>
#include <routingkit/nested_dissection.h>

#include "Algorithms/CTNR/CTNRData.h"
#include "Algorithms/CTNR/CTNRMetric.h"
#include "Algorithms/CTNR/CTNRQuery.h"
#include "Algorithms/CTNR/CTNRPreprocessor.h"
#include "Algorithms/FHL/core/HubHierarchy.h"
#include "Algorithms/FHL/core/HubAccessPreprocessor.h"
#include "Algorithms/FHL/FHLData.h"
#include "Algorithms/FHL/FHLMetric.h"
#include "Algorithms/FHL/FHLQuery.h"
#include "Algorithms/CTL/BalancedTopologyCentricTreeHierarchy.h"
#include "Algorithms/CTL/TruncatedTreeLabelling.h"
#include "Algorithms/CTL/CTLMetric.h"
#include "Algorithms/CTL/CTLQuery.h"
#include "Algorithms/CCH/CCH.h"
#include "Algorithms/CCH/CCHMetric.h"
#include "Algorithms/CCH/EliminationTreeQuery.h"
#include "DataStructures/Partitioning/SeparatorDecompositionWalk.h"
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
              "       RunP2PAlgo -a FHL        -o <file> -g <file> [-b <balance>]\n\n"

              "       RunP2PAlgo -a CCH-custom -o <file> -g <file> -s <file> [-n <num>]\n"
              "       RunP2PAlgo -a CTL-custom -o <file> -g <file> -s <file> [-n <num>]\n\n"

              "       RunP2PAlgo -a Dij        -o <file> -g <file> -d <file>\n"
              "       RunP2PAlgo -a Bi-Dij     -o <file> -g <file> -d <file>\n"
              "       RunP2PAlgo -a CH         -o <file> -h <file> -d <file>\n"
              "       RunP2PAlgo -a CCH-Dij    -o <file> -g <file> -d <file> -s <file>\n"
              "       RunP2PAlgo -a CCH-tree   -o <file> -g <file> -d <file> -s <file>\n"
              "       RunP2PAlgo -a CTL        -o <file> -g <file> -d <file> -s <file>\n"
              "       RunP2PAlgo -a CTNR       -o <file> -g <file> -d <file> -s <file>\n\n"
              "       RunP2PAlgo -a FHL        -o <file> -g <file> -d <file> -s <file>\n\n"

              "Runs the preprocessing, customization or query phase of various point-to-point\n"
              "shortest-path algorithms, such as Dijkstra, bidirectional search, CH, CCH, CTL, and CTNR.\n\n"

              "  -l                use physical lengths as metric (default: travel times)\n"
              "  -no-stall         do not use the stall-on-demand technique\n"
              "  -a <algo>         run algorithm <algo>\n"
              "  -b <balance>      balance parameter in % for nested dissection (default: 30)\n"
              "  -n <num>          run customization <num> times (default: 1000 in the *-custom modes; 1 in query runs, which report every run and the median)\n"
              "  -g <file>         input graph in binary format\n"
              "  -s <file>         separator decomposition of input graph\n"
              "  -h <file>         weighted contraction hierarchy\n"
              "  -d <file>         file that contains OD pairs (queries)\n"
              "  -o <file>         place output in <file>\n"
              "  -ctnr-thresh <num|auto>  top decomposition levels that become transit nodes (CTNR, and FHL without -fhl-region-size); auto sizes the cut so T ~ 17.6 sqrt(n) (default: 5)\n"
              "  -ctnr-prune-thresh <num>  access node pruning level threshold. Set to 0 for no pruning, to 127 for full pruning or to a value > ctnr-thresh for partial pruning (default: 0)\n"
              "  -fhl-region-size <num|auto>  leaf regions are maximal subtrees of <= num vertices; auto = n/16384, min 64\n"
              "  -fhl-subtree-rebuild <0|1>  middle update gear: rebuild only the changed separator subtree instead of everything (keeps the span-1 tables resident)\n"
              "  -cch-cache <file>         cache the metric-independent CCH here; first run writes it, later runs load it (USA: 18.5 s -> well under 1 s of start-up)\n"
              "  -fhl-subregion-size <num>  SCHEME 2 (in-region hubs): cut each leaf region again at <num> vertices; same-region queries fold over the region's own hubs instead of searching\n"
              "  -fhl-top-labels <num>  SCHEME 3 (top labels): materialize per-vertex fused values for columns with level < num; far queries become one contiguous min\n"
              "  -fhl-overlay-every <num>  SCHEME 4 (overlay): region-restricted boundary-point cliques every <num> levels; slower queries, locally updatable\n"
              "  -fhl-overlay-cuts <l1,l2,...>  custom overlay cut depths\n"
              "  -overlay-bench-updates <num>  measure overlay partial re-customization over <num> regions\n"
              "  -fhl-bench-updates <num>  measure the three update gears over <num> perturbation rounds\n"
              "  -fhl-bench-updates-verify <num>  before timing, check <num> perturbed states against a from-scratch customization, array by array\n"
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
inline void writeHeaderLine(std::ofstream &out, CTNRQuery &) {
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

template<>
inline void writeRecordLine(std::ofstream &out, CTNRQuery &algo, const int, const int64_t elapsed) {
    out << algo.getDistance() << ',' << elapsed << ',' << algo.getLastMode() << '\n';
}

// T[L] = the transit nodes a cut at L levels selects: the separator vertices of every
// decomposition node at depth < L (root = depth 0). T.back() is the whole graph.
inline std::vector<int64_t> transitCountByLevels(const SeparatorDecomposition &sd) {
    std::vector<int64_t> perDepth;
    const auto add = [&](const int node, const size_t d) {
        if (perDepth.size() <= d)
            perDepth.resize(d + 1, 0);
        perDepth[d] += sd.lastSeparatorVertex(node) - sd.firstSeparatorVertex(node);
    };
    size_t cur = 0;
    add(0, 0);
    sepdecomp::forEachNodeInDfsOrder(sd, [&](const int, const int child) { add(child, ++cur); },
                                     [&](const int, const int) { --cur; });
    std::vector<int64_t> T(perDepth.size() + 1, 0);
    for (size_t d = 0; d < perDepth.size(); ++d)
        T[d + 1] = T[d] + perDepth[d];
    return T;
}

// Number of top decomposition levels whose vertices become transit nodes (-ctnr-thresh).
// `auto` picks the level count whose transit count T is closest, as a ratio, to 17.6 sqrt(n).
// That is classic TNR sizing: the table is T x T, so T ~ sqrt(n) keeps it proportional to n,
// the same cost per vertex on every graph. A fixed level count cuts relatively deeper into a
// small graph (at 14 levels T/sqrt(n) is 26 on FLA but 17.6 on USA); a fixed fraction of the
// depth over-corrects the other way, because the depth is only ~log n and 2^(a*depth) ~ n^a
// regions grows polynomially with the graph. 17.6 is USA's value at 14 levels, so USA keeps
// the setting the doc's tables were measured with.
inline int transitLevels(const CommandLineParser &clp, const SeparatorDecomposition &sd) {
    const auto arg = clp.getValue<std::string>("ctnr-thresh", "5");
    if (arg != "auto")
        return std::stoi(arg);
    const auto T = transitCountByLevels(sd);
    const double target = 17.6 * std::sqrt(static_cast<double>(T.back()));
    const auto miss = [&](const int L) {  // how far off target, as a ratio; empty cuts never win
        return T[L] > 0 ? std::abs(std::log(T[L] / target)) : HUGE_VAL;
    };
    int best = 1;
    for (int L = 2; L < static_cast<int>(T.size()); ++L)
        if (miss(L) < miss(best))
            best = L;
    std::cout << "Transit levels (auto): " << best << " of " << T.size() - 1 << ", T = "
              << T[best] << ", target 17.6 sqrt(n) = " << static_cast<int64_t>(target)
              << "\n  T by level count:";
    for (size_t L = 1; L < T.size(); ++L)
        std::cout << ' ' << L << ':' << T[L];
    std::cout << std::endl;
    return best;
}

// The CCH is metric-independent, so for a fixed (graph, separator decomposition) it only has
// to be built once. Building it dominates start-up on the large graphs (USA: 18.5 s of a 20 s
// start-up), which makes parameter sweeps painful. With -cch-cache <file> the first run writes
// the file and later runs load it. The file is NOT self-describing: pointing two different
// graphs at one cache path silently produces nonsense, so name the files after the graph.
// Where each run's time went, written to the csv header so preprocessing and customization
// are reported separately. Uniform rule for every scheme: "scheme preprocessing" is everything
// from a ready CCH up to the customization call; "customization" is that call. The CCH line
// says whether the CCH was built (a real preprocessing time) or merely loaded from the cache.
struct PhaseTimes {
    double cchMs = 0;
    bool cchLoaded = false;
    double schemeMs = 0;
    double customizationMs = 0;              // median of the runs below
    std::vector<double> customizationRuns;   // every run, in order (-n, as in the *-custom modes)
};

inline void writePhaseTimes(std::ofstream &out, const PhaseTimes &t) {
    out << "# Preprocessing CCH: " << t.cchMs << " ms ("
        << (t.cchLoaded ? "loaded from cache" : "built") << ")\n";
    out << "# Preprocessing scheme: " << t.schemeMs << " ms\n";
    out << "# Customization: " << t.customizationMs << " ms (median of "
        << t.customizationRuns.size() << ")\n";
    out << "# Customization runs (ms):";
    for (const double ms: t.customizationRuns)
        out << ' ' << ms;
    out << '\n';
}

inline double msSince(const Timer &t) { return t.elapsed<std::chrono::microseconds>() / 1000.0; }

// -n N: customize N times and report every run plus the median -- the same repetition format
// the *-custom modes (e.g. CTL-custom) use, so every scheme's customization is measured the same
// way; a single run varies by ~5%. The caller times the first run (firstMs) and `again` repeats
// exactly that call. The queries then run on the last customization, so the correctness check
// covers the repeats.
template<typename F>
void recordCustomization(const CommandLineParser &clp, PhaseTimes &times, const double firstMs,
                         F &&again) {
    const int n = std::max(1, clp.getValue<int>("n", 1));
    times.customizationRuns = {firstMs};
    for (int i = 1; i < n; ++i) {
        Timer t;
        again();
        times.customizationRuns.push_back(msSince(t));
    }
    std::vector<double> sorted = times.customizationRuns;
    std::sort(sorted.begin(), sorted.end());
    times.customizationMs = sorted[n / 2];
}

template<typename CchT, typename InputGraphT>
void buildOrLoadCch(const CommandLineParser &clp, const InputGraphT &graph,
                    const SeparatorDecomposition &sepDecomp, CchT &cch, PhaseTimes &times) {
    const auto cache = clp.getValue<std::string>("cch-cache", "");
    if (!cache.empty()) {
        if (std::ifstream in(cache, std::ios::binary); in.good()) {
            Timer t;
            cch.readFrom(in);
            times.cchMs = msSince(t);
            times.cchLoaded = true;
            std::cout << "Loaded CCH from " << cache << " in "
                      << t.elapsed<std::chrono::milliseconds>() << " ms." << std::endl;
            return;
        }
    }
    Timer t;
    cch.preprocess(graph, sepDecomp);
    times.cchMs = msSince(t);
    std::cout << "Built CCH in " << t.elapsed<std::chrono::milliseconds>() << " ms." << std::endl;
    if (!cache.empty()) {
        // Write to a temporary name and rename only once every byte is on disk. A plain
        // ofstream does not throw when the disk fills up: it would leave a TRUNCATED cache
        // under the real name, and every later run would load garbage from it.
        const auto tmp = cache + ".partial";
        {
            std::ofstream out(tmp, std::ios::binary);
            cch.writeTo(out);
            out.close();
            if (!out) {
                std::remove(tmp.c_str());
                throw std::runtime_error("failed writing the CCH cache (disk full?) -- '" +
                                         cache + "'");
            }
        }
        if (std::rename(tmp.c_str(), cache.c_str()) != 0)
            throw std::runtime_error("cannot rename '" + tmp + "' to '" + cache + "'");
        std::cout << "Wrote CCH to " << cache << "." << std::endl;
    }
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
    int count = 0;
    while (demandFile.read_row(src, dst, rank)) {
        src = translate(src);
        dst = translate(dst);
        timer.restart();
        algo.run(src, dst);
        const auto elapsed = timer.elapsed<std::chrono::nanoseconds>();
        if (hasRanks) out << rank << ',';
        writeRecordLine(out, algo, dst, elapsed);
        ++count;
    }
    if constexpr (requires { algo.getLocalPairs(); })
        std::cout << "Local (same-region) pairs: " << algo.getLocalPairs() << " of " << count
                  << ", answered by the in-region fold: " << algo.getLocalFolds()
                  << std::endl;
    if constexpr (requires { algo.getLabelOnlyPairs(); })
        std::cout << "Label-only slices (whole scan inside the label prefix): "
                  << algo.getLabelOnlyPairs() << " of " << count << std::endl;
    if constexpr (requires { algo.getEarlyStops(); }) {
        std::cout << "Overlay ladder cutoffs: " << algo.getEarlyStops() << " of " << count
                  << " queries; by level (mi at break):";
        for (size_t i = 0; i < algo.getStopLevels().size(); ++i)
            if (algo.getStopLevels()[i] > 0)
                std::cout << ' ' << i << ':' << algo.getStopLevels()[i];
        std::cout << std::endl;
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
        PhaseTimes times;
        buildOrLoadCch(clp, graph, sepDecomp, cch, times);
        Timer phase;
        CCHMetric metric(cch, useLengths ? &graph.length(0) : &graph.travelTime(0));
        times.schemeMs = msSince(phase);
        phase.restart();
        const auto minCH = metric.buildMinimumWeightedCH();
        recordCustomization(clp, times, msSince(phase), [&] { metric.buildMinimumWeightedCH(); });

        outputFile << "# Graph: " << graphFileName << '\n';
        outputFile << "# Separator: " << sepFileName << '\n';
        outputFile << "# OD pairs: " << demandFileName << '\n';
        writePhaseTimes(outputFile, times);

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
        PhaseTimes times;
        buildOrLoadCch(clp, graph, sepDecomp, cch, times);
        Timer phase;
        CCHMetric metric(cch, useLengths ? &graph.length(0) : &graph.travelTime(0));
        times.schemeMs = msSince(phase);
        phase.restart();
        const auto minCH = metric.buildMinimumWeightedCH();
        recordCustomization(clp, times, msSince(phase), [&] { metric.buildMinimumWeightedCH(); });

        outputFile << "# Graph: " << graphFileName << '\n';
        outputFile << "# Separator: " << sepFileName << '\n';
        outputFile << "# OD pairs: " << demandFileName << '\n';
        writePhaseTimes(outputFile, times);

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
        PhaseTimes times;
        buildOrLoadCch(clp, graph, sepDecomp, cch, times);
        Timer phase;

        BalancedTopologyCentricTreeHierarchy treeHierarchy;
        treeHierarchy.preprocess(graph, sepDecomp);

        using CTLLabelSet = std::conditional_t<CTL_SIMD_LOGK == 0,
                BasicLabelSet<0, ParentInfo::NO_PARENT_INFO>,
                SimdLabelSet<CTL_SIMD_LOGK, ParentInfo::NO_PARENT_INFO>>;
        using LabellingT = TruncatedTreeLabelling<CTLLabelSet::K, CTLLabelSet::KEEP_PARENT_EDGES>;
        LabellingT ctl(treeHierarchy);
        ctl.init();

        CTLMetric<LabellingT, CTLLabelSet, CTL_USE_PERFECT_CUSTOMIZATION> metric(treeHierarchy, cch,
                                                                                 useLengths ? &graph.length(0)
                                                                                            : &graph.travelTime(0));
        times.schemeMs = msSince(phase);
        phase.restart();
        metric.buildCustomizedCTL(ctl);
        recordCustomization(clp, times, msSince(phase), [&] { metric.buildCustomizedCTL(ctl); });

        outputFile << "# Graph: " << graphFileName << '\n';
        outputFile << "# Separator: " << sepFileName << '\n';
        outputFile << "# OD pairs: " << demandFileName << '\n';
        writePhaseTimes(outputFile, times);

        CTLQuery<CTLMetric<LabellingT, CTLLabelSet, CTL_USE_PERFECT_CUSTOMIZATION>::SearchGraph, LabellingT, CTLLabelSet> algo(
                treeHierarchy, metric.upwardGraph(),
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

        // Build CCH and tree hierarchy (non-layered: the layered variant serializes the
        // top layers of the triangle pass; SD index == rank here, so no vertex mapping)
        CCH cch;
        PhaseTimes times;
        buildOrLoadCch(clp, graph, sepDecomp, cch, times);
        Timer phase;

        const int levelThreshold = transitLevels(clp, sepDecomp);
        int pruneThreshold = clp.getValue<int>("ctnr-prune-thresh", 0);

        TransitNodeHierarchy hierarchy;
        hierarchy.preprocess(graph, levelThreshold, sepDecomp); // first levelThreshold levels are transit nodes
        // Build CTNR
        CTNRData data(hierarchy.numTransitNodes(), graph.numVertices());
        CTNRPreprocessor preprocessor;
        CTNRMetric<CCH> metric(hierarchy, cch, preprocessor.getAccessNodeEdges(),
                                       useLengths ? &graph.length(0) : &graph.travelTime(0), pruneThreshold);

        // Preprocess CTNR
        preprocessor.preprocess(hierarchy, cch, data);
        times.schemeMs = msSince(phase);
        std::cout << "Finished preprocessing" << std::endl;
        // Customize CTNR
        phase.restart();
        metric.customize(data);
        recordCustomization(clp, times, msSince(phase), [&] { metric.customize(data); });
        outputFile << "# Graph: " << graphFileName << '\n';
        outputFile << "# OD pairs: " << demandFileName << '\n';
        writePhaseTimes(outputFile, times);
        // outputFile << "# Memory usage total: " << ctnr.sizeInBytes() / BYTES_PER_MB << " MB" << '\n';

        // Use generic runQueries with CTNRQuery; pass CCH rank IDs to the algo
        CTNRQuery algo(hierarchy, data, metric.getLocalEliminationTree(), cch.getUpwardGraph(),
                                   metric.getCCHMetric().upwardWeights(), metric.getCCHMetric().downwardWeights());

        outputFile << "# Memory usage CCH: " << (cch.sizeInBytes()) / BYTES_PER_MB << " MB" << '\n';
        outputFile << "# Memory usage hierarchy: " << (hierarchy.sizeInBytes()) / BYTES_PER_MB << " MB" << '\n';
        outputFile << "# Memory usage distance table: " << (data.sizeDistanceTableInBytes()) / BYTES_PER_MB << " MB"
                   << '\n';
        outputFile << "# Memory usage access nodes: " << (data.sizeAccessNodesInBytes()) / BYTES_PER_MB << " MB"
                   << '\n';
        outputFile << "# Memory usage preprocessor: " << (preprocessor.sizeInBytes()) / BYTES_PER_MB << " MB" << '\n';
        outputFile << "# Memory usage query: " << (algo.sizeInBytes()) / BYTES_PER_MB << " MB" << '\n';
        outputFile << "# Memory usage total: " <<
                   (cch.sizeInBytes() + hierarchy.sizeInBytes() + data.sizeInBytes() + preprocessor.sizeInBytes() +
                    metric.sizeInBytes() +
                    algo.sizeInBytes()) / BYTES_PER_MB << " MB" << '\n';
        runQueries(algo, demandFileName, outputFile, [&](const int v) { return cch.getRanks()[v]; });

    } else if (algorithmName == "FHL") {

        // Run FHL queries (seed/rectangle family, or the overlay variant)
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

        CCH cch; // non-layered: the layered variant serializes the expensive top layers
                 // of the triangle pass (measured 10x slower customization)
        PhaseTimes times;
        buildOrLoadCch(clp, graph, sepDecomp, cch, times);
        Timer phase;

        const int levelThreshold = transitLevels(clp, sepDecomp);
        const int pruneThreshold = clp.getValue<int>("ctnr-prune-thresh", 0);
        const auto fhlRegionSizeStr = clp.getValue<std::string>("fhl-region-size", "0");
        // auto: ~n/16384 regions keeps the transit count at its ~45*sqrt(n) balance point
        const int32_t fhlRegionSize = fhlRegionSizeStr == "auto"
                ? std::max(64, graph.numVertices() / 16384)
                : std::stoi(fhlRegionSizeStr);
        const int32_t fhlTopLabels = clp.getValue<int>("fhl-top-labels", 0);
        const int32_t fhlBenchUpdates = clp.getValue<int>("fhl-bench-updates", 0);
        const bool fhlSubtreeRebuild = clp.getValue<int>("fhl-subtree-rebuild", 0) != 0;
        const int32_t fhlSubregionSize = clp.getValue<int>("fhl-subregion-size", 0);
        const int32_t fhlOverlayEvery = clp.getValue<int>("fhl-overlay-every", 0);
        const int32_t overlayBenchUpdates = clp.getValue<int>("overlay-bench-updates", 0);
        // 1 = seed/rectangle family (default), 2 = overlay
        const int32_t scheme = fhlOverlayEvery > 0 ? 2 : 1;
        std::vector<int32_t> fhlOverlayCuts;
        {
            const auto cutsStr = clp.getValue<std::string>("fhl-overlay-cuts", "");
            std::stringstream ss(cutsStr);
            std::string tok;
            while (std::getline(ss, tok, ','))
                if (!tok.empty())
                    fhlOverlayCuts.push_back(std::stoi(tok));
        }

        HubHierarchy hierarchy;
        if (fhlRegionSize > 0) {
            if (scheme != 1)
                throw std::invalid_argument("-fhl-region-size cannot be combined with "
                                            "-fhl-overlay");
            hierarchy.preprocessSizeCut(graph, fhlRegionSize, fhlSubregionSize, sepDecomp);
        } else {
            hierarchy.preprocess(graph, levelThreshold, sepDecomp); // SD index == rank for the non-layered CCH
        }

        fhl::Regions regions;
        {
            Timer regionTimer;
            regions.build(hierarchy, graph, cch.getRanks(), scheme == 2, fhlOverlayEvery,
                          fhlOverlayCuts, fhlRegionSize,
                          scheme == 1 && fhlSubregionSize > 0 ? 1 : 0);
            std::cout << "Built regions in "
                      << regionTimer.elapsed<std::chrono::milliseconds>() << " ms. ";
            regions.printStats(std::cout);
        }

        FHLData data(hierarchy.numHubs(), graph.numVertices());
        HubAccessPreprocessor preprocessor;
        // A plain customization frees the access-hub scaffolding at its end (it is only needed to
        // customize), which makes it one-shot. The update benchmarks and -n > 1 customize again,
        // so they keep it; for -n the scaffolding is released after the last run below, before
        // any memory is reported, so the reported footprint is the same as a single run's.
        const int customizationRuns = std::max(1, clp.getValue<int>("n", 1));
        const bool keepAccessHubs = fhlBenchUpdates > 0 || customizationRuns > 1;
        FHLMetric<CCH> metric(hierarchy, cch, preprocessor.getAccessHubEdges(),
                                                useLengths ? &graph.length(0) : &graph.travelTime(0),
                                                static_cast<fhl::Level>(pruneThreshold), scheme,
                                                &regions,
                                                fhlTopLabels, keepAccessHubs,
                                                fhlSubtreeRebuild && scheme == 1);

        preprocessor.preprocess(hierarchy, cch, data.getBaseData());
        times.schemeMs = msSince(phase);
        std::cout << "Finished preprocessing" << std::endl;
        Timer customizationTimer;
        metric.customize(data);
        const auto customizationTime = customizationTimer.elapsed<std::chrono::microseconds>();
        recordCustomization(clp, times, customizationTime / 1000.0, [&] { metric.customize(data); });
        if (customizationRuns > 1 && fhlBenchUpdates <= 0)
            FHLMetric<CCH>::releaseAccessHubs(data.getBaseData());
        std::cout << "Finished FHL customization in " << customizationTime
                  << " microseconds." << std::endl;
        const auto avgNonInfty = data.getBaseData().computeAverageNumberOfNonInftyAccessHubs();
        std::cout << "Average non-infty access nodes per vertex: forward " << avgNonInfty.first
                  << ", backward " << avgNonInfty.second << std::endl;

        if (overlayBenchUpdates > 0 && scheme == 2) {
            // Partial re-customization benchmark. Weights are unchanged, so the rebuilt
            // blocks must be identical; the query verification below runs on the
            // recustomized data and doubles as the correctness check.
            const int32_t numLeaf = regions.leaf().numRegions();
            Timer rt;
            for (int32_t i = 0; i < overlayBenchUpdates; ++i) {
                const int32_t r = static_cast<int32_t>((static_cast<int64_t>(i) * 7919 + 13) %
                                                       numLeaf);
                fhl::Overlay::recustomizeLeafRegion(
                        regions, hierarchy, cch, metric.getCCHMetric().upwardWeights(),
                        metric.getCCHMetric().downwardWeights(),
                        metric.getLocalEliminationTree(),
                        useLengths ? &graph.length(0) : &graph.travelTime(0), r, data.ladder);
            }
            const auto rcTime = rt.elapsed<std::chrono::microseconds>();
            std::cout << "Overlay partial re-customization: " << overlayBenchUpdates
                      << " regions, " << rcTime / std::max(1, overlayBenchUpdates)
                      << " us/region (full customization " << customizationTime << " us)."
                      << std::endl;
        }

        if (fhlBenchUpdates > 0 && scheme == 1) {
            // Incremental flat benchmark: perturb a few interior edges of a leaf region,
            // partially customize, restore the weights, partially customize back. The final
            // state equals the initial one, so the query verification below still applies.
            const int32_t numLeaf = regions.leaf().numRegions();
            auto *const weights = useLengths ? &graph.length(0) : &graph.travelTime(0);
            // Interior input edges per leaf region, with their CCH up-edge ids.
            // (input edge, reverse input edge, cch up-edge). Both directions are perturbed
            // together: a one-directional change makes the metric asymmetric, which the
            // symmetric structures cannot represent (the partial paths then escalate).
            struct RegionEdge {
                int32_t inEdge, revEdge, cchEdge;
            };
            std::vector<std::vector<RegionEdge>> regionEdges(numLeaf);
            const auto &upGraph = cch.getUpwardGraph();
            FORALL_VALID_EDGES(graph, u, e) {
                const int32_t ru = cch.getRanks()[u];
                const int32_t rv = cch.getRanks()[graph.edgeHead(e)];
                const int32_t r1 = regions.leafRegionOfVertex[ru];
                if (r1 < 0 || r1 != regions.leafRegionOfVertex[rv] || ru == rv)
                    continue;
                const int32_t t = std::min(ru, rv), h = std::max(ru, rv);
                int32_t cchEdge = -1;
                for (int ce = upGraph.firstEdge(t); ce < upGraph.lastEdge(t); ++ce)
                    if (upGraph.edgeHead(ce) == h) {
                        cchEdge = ce;
                        break;
                    }
                if (cchEdge < 0)
                    continue;
                int32_t rev = -1;
                const int32_t head = graph.edgeHead(e);
                for (int re = graph.firstEdge(head); re < graph.lastEdge(head); ++re)
                    if (graph.edgeHead(re) == u) {
                        rev = re;
                        break;
                    }
                regionEdges[r1].push_back({e, rev, cchEdge});
            }
            std::vector<int32_t> levelHist(32, 0);
            // Hard verification of the PERTURBED state (the timing loop only checks the
            // round trip): full-customize on perturbed weights = ground truth; compare the
            // partial result against it array by array.
            const int32_t verifyRounds = clp.getValue<int>("fhl-bench-updates-verify", 0);
            for (int32_t i = 0; i < verifyRounds; ++i) {
                const int32_t r = static_cast<int32_t>(
                        (static_cast<int64_t>(i) * 104729 + 7) % numLeaf);
                if (regionEdges[r].size() < 1)
                    continue;
                const auto re = regionEdges[r][regionEdges[r].size() / 2];
                const int32_t ie = re.inEdge, ce = re.cchEdge;
                const int32_t w0 = weights[ie];
                const int32_t wr0 = re.revEdge >= 0 ? weights[re.revEdge] : 0;
                weights[ie] = w0 + w0 / 20 + 1;
                if (re.revEdge >= 0)
                    weights[re.revEdge] = wr0 + wr0 / 20 + 1;
                metric.customize(data); // ground truth on perturbed weights
                std::vector<int32_t> gtSeed(data.ladder.seedFwd.data(),
                                            data.ladder.seedFwd.data() +
                                                    data.ladder.seedFwd.size());
                auto gtPtr = data.ladder.csrOut.ptr;
                auto gtEnd = data.ladder.csrOut.colEnd;
                auto gtRow = data.ladder.csrOut.row;
                auto gtVal = data.ladder.csrOut.val;
                auto gtMin = data.ladder.colMinOut;
                weights[ie] = w0;
                if (re.revEdge >= 0)
                    weights[re.revEdge] = wr0;
                metric.customize(data); // back to base state
                weights[ie] = w0 + w0 / 20 + 1;
                if (re.revEdge >= 0)
                    weights[re.revEdge] = wr0 + wr0 / 20 + 1;
                const auto st = metric.partialCustomizeFlat(data, {ce}, r);
                bool ok = true;
                if (st.outcome ==
                    FHLMetric<CCH>::PARTIAL_NEED_FULL) {
                    metric.customize(data);
                }
                const bool okSeed =
                        std::equal(gtSeed.begin(), gtSeed.end(), data.ladder.seedFwd.data());
                // Column-wise comparison: the incremental gears may rewrite a region's
                // segment in place or move it to the tail, so the arrays are equal only
                // per column, not element by element.
                const auto &csr = data.ladder.csrOut;
                bool okCsr = gtPtr.size() == csr.ptr.size();
                if (okCsr) {
                    const int64_t nCols = static_cast<int64_t>(gtPtr.size()) - 1;
                    for (int64_t j = 0; j < nCols && okCsr; ++j) {
                        const int64_t gs = gtPtr[j];
                        const int64_t ge = gtEnd.empty() ? gtPtr[j + 1] : gtEnd[j];
                        const int64_t ps = csr.ptr[j];
                        const int64_t pe = csr.endArray()[j];
                        if (ge - gs != pe - ps) {
                            okCsr = false;
                            break;
                        }
                        for (int64_t k = 0; k < ge - gs; ++k)
                            if (gtRow[gs + k] != csr.row[ps + k] ||
                                gtVal[gs + k] != csr.val[ps + k]) {
                                okCsr = false;
                                break;
                            }
                    }
                }
                const bool okMin = gtMin == data.ladder.colMinOut;
                ok = okSeed && okCsr && okMin;
                std::cout << "[verify] round " << i << " outcome " << st.outcome
                          << " changed=" << st.changed << " minLvl=" << st.minChangedLevel
                          << " seed:" << (okSeed ? "ok" : "BAD")
                          << " csr:" << (okCsr ? "ok" : "BAD")
                          << " colMin:" << (okMin ? "ok" : "BAD")
                          << (ok ? " MATCH" : " *** MISMATCH ***") << std::endl;
                weights[ie] = w0;
                if (re.revEdge >= 0)
                    weights[re.revEdge] = wr0;
                metric.customize(data); // restore base for the next round / final checks
            }
            int64_t fastUs = 0, fullUs = 0, superUs = 0;
            int32_t nNoop = 0, nFast = 0, nSuper = 0, nFull = 0;
            int64_t recomputed = 0, changedEdges = 0;
            int32_t rounds = 0;
            const auto applyUpdate = [&](const std::vector<int32_t> &cchEdges,
                                         const int32_t r) {
                Timer t;
                const auto st = metric.partialCustomizeFlat(data, cchEdges, r);
                ++levelHist[std::min(31, static_cast<int>(st.minChangedLevel))];
                recomputed += st.recomputed;
                changedEdges += st.changed;
                if (st.outcome == FHLMetric<CCH>::PARTIAL_NEED_FULL) {
                    metric.customize(data); // access structures were kept alive
                    fullUs += t.elapsed<std::chrono::microseconds>();
                    ++nFull;
                } else if (st.outcome == FHLMetric<CCH>::PARTIAL_SUPER) {
                    superUs += t.elapsed<std::chrono::microseconds>();
                    ++nSuper;
                } else {
                    fastUs += t.elapsed<std::chrono::microseconds>();
                    st.outcome == FHLMetric<CCH>::PARTIAL_FAST ? ++nFast
                                                                                 : ++nNoop;
                }
            };
            for (int32_t i = 0; rounds < fhlBenchUpdates && i < 4 * fhlBenchUpdates; ++i) {
                const int32_t r = static_cast<int32_t>(
                        (static_cast<int64_t>(i) * 7919 + 13) % numLeaf);
                if (regionEdges[r].size() < 4)
                    continue;
                ++rounds;
                std::vector<int32_t> cchEdges;
                std::vector<std::pair<int32_t, int32_t>> saved; // (input edge, old weight)
                const size_t stride = std::max<size_t>(1, regionEdges[r].size() / 4);
                for (size_t k = 0; k < regionEdges[r].size() && cchEdges.size() < 1;
                     k += stride) {
                    const auto re2 = regionEdges[r][k];
                    saved.emplace_back(re2.inEdge, weights[re2.inEdge]);
                    weights[re2.inEdge] += weights[re2.inEdge] / 20 + 1; // ~x1.05 increase
                    if (re2.revEdge >= 0) {
                        saved.emplace_back(re2.revEdge, weights[re2.revEdge]);
                        weights[re2.revEdge] += weights[re2.revEdge] / 20 + 1;
                    }
                    cchEdges.push_back(re2.cchEdge);
                }
                applyUpdate(cchEdges, r); // perturbed
                for (const auto &[ie, w] : saved)
                    weights[ie] = w;
                applyUpdate(cchEdges, r); // restored — final state == initial state
            }
            std::cout << "Escalation level histogram (min changed level per update):";
            for (int l = 0; l < 32; ++l)
                if (levelHist[l] > 0)
                    std::cout << ' ' << l << ':' << levelHist[l];
            std::cout << std::endl;
            const int32_t nPartial = nNoop + nFast;
            std::cout << "Incremental flat: " << 2 * rounds << " updates -> " << nNoop
                      << " no-op, " << nFast << " fast, " << nSuper << " super ("
                      << (nSuper > 0 ? superUs / nSuper : 0) << " us avg), " << nFull
                      << " full; avg partial "
                      << (nPartial > 0 ? fastUs / nPartial : 0) << " us"
                      << (nFull > 0 ? " (avg full " + std::to_string(fullUs / nFull) + " us)"
                                    : std::string())
                      << "; avg " << (2 * rounds > 0 ? recomputed / (2 * rounds) : 0)
                      << " edges rechecked / " << (2 * rounds > 0 ? changedEdges / (2 * rounds) : 0)
                      << " changed per update (full customization " << customizationTime
                      << " us)." << std::endl;
        }

        outputFile << "# Graph: " << graphFileName << '\n';
        outputFile << "# OD pairs: " << demandFileName << '\n';
        writePhaseTimes(outputFile, times);

        FHLQuery algo(hierarchy, data, metric.getLocalEliminationTree(),
                                   cch.getUpwardGraph(), metric.getCCHMetric().upwardWeights(),
                                   metric.getCCHMetric().downwardWeights(),
                                   &regions);

        outputFile << "# Memory usage CCH: " << (cch.sizeInBytes()) / BYTES_PER_MB << " MB" << '\n';
        outputFile << "# Memory usage hierarchy: " << (hierarchy.sizeInBytes()) / BYTES_PER_MB << " MB" << '\n';
        outputFile << "# Memory usage access nodes: " << (data.getBaseData().sizeAccessHubsInBytes()) / BYTES_PER_MB
                   << " MB" << '\n';
        outputFile << "# Memory usage hierarchical tables: " << (data.sizeTablesInBytes()) / BYTES_PER_MB << " MB"
                   << '\n';
        outputFile << "# Memory usage region data: " << (data.ladder.sizeInBytes()) / BYTES_PER_MB
                   << " MB" << '\n';
        outputFile << "# Memory usage regions: "
                   << regions.sizeInBytes() / BYTES_PER_MB
                   << " MB" << '\n';
        outputFile << "# Memory usage preprocessor: " << (preprocessor.sizeInBytes()) / BYTES_PER_MB << " MB" << '\n';
        outputFile << "# Memory usage query: " << (algo.sizeInBytes()) / BYTES_PER_MB << " MB" << '\n';
        // The Regions object is live query data (leafRegionOfVertex, flatLevelEnd, the
        // boundary-point lists), so it belongs in the total; it used to be reported on its
        // own line only.
        outputFile << "# Memory usage total: " <<
                   (cch.sizeInBytes() + hierarchy.sizeInBytes() + data.sizeInBytes() +
                    preprocessor.sizeInBytes() + metric.sizeInBytes() + algo.sizeInBytes() +
                    regions.sizeInBytes()) / BYTES_PER_MB << " MB"
                   << '\n';
        runQueries(algo, demandFileName, outputFile, [&](const int v) { return cch.getRanks()[v]; });

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

    } else if (algorithmName == "CTL" || algorithmName == "CTNR") {
        // Run the preprocessing phase of CTL.
        std::cout << "Constructing separator decomposition with strict dissection for CTL/CTNR for " << graphFileName
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

    }
//    else if (algorithmName == "CTNR") {
//        std::cout << "CTNR preprocessing (using existing separator decomposition) for " << graphFileName
//                  << "... " << std::flush;
//        Timer timer;
//        std::string sepFileName = graphFileName;
//        size_t lastDot = sepFileName.find_last_of('.');
//        if (lastDot != std::string::npos) {
//            sepFileName = sepFileName.substr(0, lastDot);
//        }
//        sepFileName += ".strict_bisep.bin";
//
//        // 读取现有的分隔分解文件
//        std::ifstream sepFile(sepFileName, std::ios::binary);
//        if (!sepFile.good()) {
//            std::cout << "Separator decomposition file not found: " << sepFileName << std::endl;
//            std::cout << "Please run CTL preprocessing first to generate separator decomposition." << std::endl;
//            return;
//        }
//
//        SeparatorDecomposition sepDecomp;
//        sepDecomp.readFrom(sepFile);
//        sepFile.close();
//
//        std::cout << "Loaded separator decomposition with " << sepDecomp.tree.size() << " nodes" << std::endl;
//
//        // Build CTNR
//        CTNR<InputGraph> ctnr(sepDecomp, 5); // Use top 5 levels as transit nodes
//        ctnr.preprocess(graph);
//
//        std::cout << "CTNR preprocessing completed with " << ctnr.getTransitNodes().size()
//                  << " transit nodes" << std::endl;
//
//        if (!endsWith(outputFileName, ".ctnr.bin"))
//            outputFileName += ".ctnr.bin";
//        std::ofstream outputFile(outputFileName, std::ios::binary);
//        if (!outputFile.good())
//            throw std::invalid_argument("file cannot be opened -- '" + outputFileName);
//        // CTNR currently has no serialization; keep placeholder to match other modes
//
//        const auto preprocessTime = timer.elapsed<std::chrono::microseconds>();
//        std::cout << " finished (" << preprocessTime << " microseconds)." << std::endl;
//
//    }
    else if (algorithmName == "CTNR-custom") {

        // Run the customization phase of CTNR.
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


        const int levelThreshold = transitLevels(clp, decomp);
        int pruneThreshold = clp.getValue<int>("ctnr-prune-thresh", 0);

        Timer timer;
        // Build CCH and tree hierarchy
        LayerCCH cch;
        cch.preprocess(graph, decomp);

        TransitNodeHierarchy hierarchy;
        // first levelThreshold levels of separator decomposition are transit nodes
        hierarchy.preprocess(graph, levelThreshold, decomp, cch.getSepDecompToCCHGraphVertexMapping());

        // Build CTNR
        CTNRData data(hierarchy.numTransitNodes(), graph.numVertices());
        CTNRPreprocessor preprocessor;
        preprocessor.preprocess(hierarchy, cch, data);

        const auto preprocessTime = timer.elapsed<std::chrono::microseconds>();
        outputFile << "# Preprocess time (for given sepdecomp): " << preprocessTime << " microseconds.\n";

        outputFile
                << "cch_basic_customization,distance_table_computation,access_node_computation,total_time,avg_noninfty_dists_forw,avg_noninfty_dists_backw\n";
        int64_t cchBasicCustom, distTableComp, accessNodeComp, tot;
        timer.restart();
        for (auto i = 0; i < numCustomRuns; ++i) {
            CTNRMetric<LayerCCH> metric(hierarchy, cch, preprocessor.getAccessNodeEdges(),
                                           useLengths ? &graph.length(0) : &graph.travelTime(0), pruneThreshold);
            timer.restart();
            metric.customizeWithMeasurements(data, cchBasicCustom, distTableComp, accessNodeComp);
            tot = timer.elapsed<std::chrono::microseconds>();
            const auto [avgNumNonInftyForward, avgNumNonInftyBackward] = data.computeAverageNumberOfNonInftyAccessNodes();
            std::cout << "CTNR: Average number of non-infinity forward/backward distances per vertex: "
                      << avgNumNonInftyForward << "/" << avgNumNonInftyBackward << std::endl;
            outputFile << cchBasicCustom << ',' << distTableComp << ',' << accessNodeComp << ',' << tot
                       << ',' << avgNumNonInftyForward << ',' << avgNumNonInftyBackward << '\n';

        }
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
                CCHMetric metric(cch, useLengths ? &graph.length(0) : &graph.travelTime(0));
                timer.restart();
                if constexpr (CTL_USE_PERFECT_CUSTOMIZATION) {
                    metric.buildMinimumWeightedCH();
                } else {
                    metric.customize();
                }
                cchCustom = timer.elapsed<std::chrono::microseconds>();
            }
            {
                CTLMetric<LabellingT, CTLLabelSet, CTL_USE_PERFECT_CUSTOMIZATION> metric(treeHierarchy, cch,
                                                                                         useLengths ? &graph.length(0)
                                                                                                    : &graph.travelTime(
                                                                                                 0));
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
    // Compile-time default, overridable at runtime for scaling measurements.
    int numThreads = NUM_THREADS;
    if (const char *const env = std::getenv("CTNR_THREADS"))
        numThreads = std::max(1, std::atoi(env));
    std::cout << "Using " << numThreads << " threads for parallel regions in RunP2PAlgo." << std::endl;
    omp_set_num_threads(numThreads);
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
