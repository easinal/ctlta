#pragma once

#include <cassert>
#include <cstdint>
#include <vector>

#include <omp.h>
#include <routingkit/constants.h>
#include <routingkit/customizable_contraction_hierarchy.h>

#include "DataStructures/Graph/Attributes/EdgeIdAttribute.h"
#include "DataStructures/Graph/Attributes/EdgeTailAttribute.h"
#include "DataStructures/Graph/Graph.h"
#include "DataStructures/Partitioning/SeparatorDecomposition.h"
#include "DataStructures/Utilities/Permutation.h"
#include "Tools/Constants.h"
#include "Tools/Workarounds.h"
#include "Tools/Timer.h"

// A metric-independent customizable contraction hierarchy. It uses a nested dissection order
// associated with a separator decomposition to order the vertices by importance.
class CCH {
public:
    using EliminationTree = std::vector<int32_t>;                             // The elimination tree.
    using UpGraph = StaticGraph<VertexAttrs<>, EdgeAttrs<EdgeTailAttribute>>; // The upward graph.
    using DownGraph = StaticGraph<VertexAttrs<>, EdgeAttrs<EdgeIdAttribute>>; // The downward graph.

    // Constructs an empty CCH.
    CCH() = default;

    // Constructs a CCH from the specified binary file.
    explicit CCH(std::ifstream &in) {
        readFrom(in);
    }

    // Builds the metric-independent CCH for the specified graph and separator decomposition.
    template<typename InputGraphT>
    void preprocess(const InputGraphT &inputGraph, const SeparatorDecomposition &sepDecomp) {
        assert(inputGraph.numVertices() == sepDecomp.order.size());
        std::vector<unsigned int> order(sepDecomp.order.begin(), sepDecomp.order.end());
        std::vector<unsigned int> tails(inputGraph.numEdges());
        std::vector<unsigned int> heads(inputGraph.numEdges());
        FORALL_VALID_EDGES(inputGraph, u, e) {
                tails[e] = u;
                heads[e] = inputGraph.edgeHead(e);
            }
        RoutingKit::CustomizableContractionHierarchy cch(order, tails, heads);

        decomp = sepDecomp;
        ranks.assign(cch.rank.begin(), cch.rank.end());
        eliminationTree.assign(cch.elimination_tree_parent.begin(), cch.elimination_tree_parent.end());
        eliminationTree.back() = INVALID_VERTEX;

        upGraph.reserve(inputGraph.numVertices(), cch.cch_arc_count());
        downGraph.reserve(inputGraph.numVertices(), cch.cch_arc_count());
        for (int v = 0; v != inputGraph.numVertices(); ++v) {
            upGraph.appendVertex();
            downGraph.appendVertex();
            for (int e = cch.up_first_out[v]; e != cch.up_first_out[v + 1]; ++e)
                upGraph.appendEdge(cch.up_head[e], cch.up_tail[e]);
            for (int e = cch.down_first_out[v]; e != cch.down_first_out[v + 1]; ++e)
                downGraph.appendEdge(cch.down_head[e], cch.down_to_up[e]);
        }

        firstUpInputEdge.resize(upGraph.numEdges() + 1);
        firstDownInputEdge.resize(upGraph.numEdges() + 1);
        FORALL_EDGES(upGraph, e) {
            firstUpInputEdge[e] = upInputEdges.size();
            firstDownInputEdge[e] = downInputEdges.size();
            if (cch.does_cch_arc_have_input_arc.is_set(e)) {
                const int i = cch.does_cch_arc_have_input_arc_mapper.to_local(e);
                if (cch.forward_input_arc_of_cch[i] != RoutingKit::invalid_id)
                    upInputEdges.push_back(cch.forward_input_arc_of_cch[i]);
                if (cch.backward_input_arc_of_cch[i] != RoutingKit::invalid_id)
                    downInputEdges.push_back(cch.backward_input_arc_of_cch[i]);
                if (cch.does_cch_arc_have_extra_input_arc.is_set(e)) {
                    const int j = cch.does_cch_arc_have_extra_input_arc_mapper.to_local(e);
                    const int firstExtraUpInputEdge = cch.first_extra_forward_input_arc_of_cch[j];
                    const int firstExtraDownInputEdge = cch.first_extra_backward_input_arc_of_cch[j];
                    const int lastExtraUpInputEdge = cch.first_extra_forward_input_arc_of_cch[j + 1];
                    const int lastExtraDownInputEdge = cch.first_extra_backward_input_arc_of_cch[j + 1];
                    for (int k = firstExtraUpInputEdge; k != lastExtraUpInputEdge; ++k)
                        upInputEdges.push_back(cch.extra_forward_input_arc_of_cch[k]);
                    for (int k = firstExtraDownInputEdge; k != lastExtraDownInputEdge; ++k)
                        downInputEdges.push_back(cch.extra_backward_input_arc_of_cch[k]);
                }
            }
        }
        firstUpInputEdge.back() = upInputEdges.size();
        firstDownInputEdge.back() = downInputEdges.size();

        buildLayering();
    }

    // Returns the separator decomposition used to build this CCH.
    const SeparatorDecomposition &getSeparatorDecomposition() const noexcept {
        return decomp;
    }

    // Returns the order in which vertices were contracted.
    const Permutation &getContractionOrder() const noexcept {
        return decomp.order;
    }

    // Returns the position of each vertex in the contraction order.
    const Permutation &getRanks() const noexcept {
        return ranks;
    }

    // Returns the elimination tree.
    const EliminationTree &getEliminationTree() const noexcept {
        return eliminationTree;
    }

    // Returns the upward graph.
    const UpGraph &getUpwardGraph() const noexcept {
        return upGraph;
    }

    // Applies func to each upward input edge mapping to the specified edge in the CCH.
    template<typename CallableT>
    bool forEachUpwardInputEdge(const int e, CallableT func) const {
        assert(e >= 0);
        assert(e < upGraph.numEdges());
        for (auto i = firstUpInputEdge[e]; i != firstUpInputEdge[e + 1]; ++i)
            if (!func(upInputEdges[i]))
                return false;
        return true;
    }

    // Applies func to each downward input edge mapping to the specified edge in the CCH.
    template<typename CallableT>
    bool forEachDownwardInputEdge(const int e, CallableT func) const {
        assert(e >= 0);
        assert(e < upGraph.numEdges());
        for (auto i = firstDownInputEdge[e]; i != firstDownInputEdge[e + 1]; ++i)
            if (!func(downInputEdges[i]))
                return false;
        return true;
    }
    
    template<typename CallableT>
    void forEachVertexBottomUp(CallableT func) const {
        forEachVertexBottomUp(0, upGraph.numVertices(), 0, func);
    }


    template <class SeqCallable, class ParaCallable>
    void forEachVertexBottomUpByLayer(SeqCallable Seqfunc, ParaCallable Parafunc) const {
        assert(numLayers != 0 || upGraph.numVertices() == 0);
        if (numLayers == 0)
            return;
        for (auto layer = 0; layer < numLayers; ++layer) {
            const auto layerSize = layerOffsets[layer + 1] - layerOffsets[layer];
            if (omp_get_num_threads()>1 && layerSize > 32 * omp_get_num_threads()) {
#pragma omp parallel for schedule(dynamic)
                for (auto i = layerOffsets[layer]; i < layerOffsets[layer + 1]; ++i){
                    Parafunc(layerVertices[i]);
                }
            } else {
                for (auto i = layerOffsets[layer]; i < layerOffsets[numLayers]; ++i){
                    Seqfunc(layerVertices[i]);
                }
                break;
            }
        }
    }

    template <class CallableT>
    void forEachVertexBottomUpByLayer(CallableT func) const {
        forEachVertexBottomUpByLayer(func, func);
    }

    template<typename CallableT>
    void forEachVertexTopDown(CallableT func) const {
        forEachVertexTopDown(0, upGraph.numVertices(), 0, func);
    }

    template <class SeqCallable, class ParaCallable>
    void forEachVertexTopDownByLayer(SeqCallable Seqfunc, ParaCallable Parafunc) const {
        assert(numLayers != 0 || upGraph.numVertices() == 0);
        for (auto layer = numLayers - 1; layer >= 0; --layer) {
            const auto layerSize = layerOffsets[layer + 1] - layerOffsets[layer];
            if (omp_get_num_threads()>1 && layerSize > 32 * omp_get_num_threads()) {
                #pragma omp parallel for schedule(dynamic)
                for (auto i = layerOffsets[layer+1]-1; i >= layerOffsets[layer]; --i){
                    Parafunc(layerVertices[i]);
                }
            } else {
                for (auto i = layerOffsets[layer+1]-1; i >= layerOffsets[0]; --i)
                    Seqfunc(layerVertices[i]);
                break;
            }
        }
    }
    template<typename CallableT>
    void forEachVertexTopDownByLayer(CallableT func) const {
        forEachVertexTopDownByLayer(func, func);
    }
    
    // Applies func to each lower triangle of the specified edge.
    template<typename CallableT>
    bool forEachLowerTriangle(const int tail, const int head, const int edge, CallableT func) const {
        unused(edge);
        assert(head == upGraph.edgeHead(edge));
        int edgeOnTail = downGraph.firstEdge(tail);
        int edgeOnHead = downGraph.firstEdge(head);
        const int lastEdgeOnTail = downGraph.lastEdge(tail);
        const int lastEdgeOnHead = downGraph.lastEdge(head);
        while (edgeOnTail != lastEdgeOnTail && edgeOnHead != lastEdgeOnHead) {
            const int neighborOfTail = downGraph.edgeHead(edgeOnTail);
            const int neighborOfHead = downGraph.edgeHead(edgeOnHead);
            if (neighborOfTail < neighborOfHead) {
                ++edgeOnTail;
            } else if (neighborOfTail > neighborOfHead) {
                ++edgeOnHead;
            } else {
                if (!func(neighborOfTail, downGraph.edgeId(edgeOnTail), downGraph.edgeId(edgeOnHead)))
                    return false;
                ++edgeOnTail;
                ++edgeOnHead;
            }
        }
        return true;
    }

    // Applies func to each upper triangle of the specified edge.
    template<typename CallableT>
    bool forEachUpperTriangle(const int tail, const int head, const int edge, CallableT func) const {
        assert(head == upGraph.edgeHead(edge));
        int edgeOnTail = edge + 1;
        int edgeOnHead = upGraph.firstEdge(head);
        const int lastEdgeOnTail = upGraph.lastEdge(tail);
        const int lastEdgeOnHead = upGraph.lastEdge(head);
        while (edgeOnTail != lastEdgeOnTail && edgeOnHead != lastEdgeOnHead) {
            const int neighborOfTail = upGraph.edgeHead(edgeOnTail);
            const int neighborOfHead = upGraph.edgeHead(edgeOnHead);
            if (neighborOfTail < neighborOfHead) {
                ++edgeOnTail;
            } else if (neighborOfTail > neighborOfHead) {
                ++edgeOnHead;
            } else {
                if (!func(neighborOfTail, edgeOnTail, edgeOnHead))
                    return false;
                ++edgeOnTail;
                ++edgeOnHead;
            }
        }
        return true;
    }

    // Reads the CCH from the specified binary file.
    void readFrom(std::ifstream &in) {
        decomp.readFrom(in);
        ranks.readFrom(in);
        bio::read(in, eliminationTree);
        upGraph.readFrom(in);
        downGraph.readFrom(in);
        bio::read(in, firstUpInputEdge);
        bio::read(in, firstDownInputEdge);
        bio::read(in, upInputEdges);
        bio::read(in, downInputEdges);

        buildLayering();
    }

    // Writes the CCH to the specified binary file.
    void writeTo(std::ofstream &out) const {
        decomp.writeTo(out);
        ranks.writeTo(out);
        bio::write(out, eliminationTree);
        upGraph.writeTo(out);
        downGraph.writeTo(out);
        bio::write(out, firstUpInputEdge);
        bio::write(out, firstDownInputEdge);
        bio::write(out, upInputEdges);
        bio::write(out, downInputEdges);
    }

    uint64_t sizeInBytes() const {
        return sizeof(*this)
               + decomp.sizeInBytes()
               + ranks.sizeInBytes()
               + eliminationTree.size() * sizeof(int32_t)
               + upGraph.sizeInBytes()
               + downGraph.sizeInBytes()
               + firstUpInputEdge.size() * sizeof(int32_t)
               + firstDownInputEdge.size() * sizeof(int32_t)
               + upInputEdges.size() * sizeof(int32_t)
               + downInputEdges.size() * sizeof(int32_t)
               + layerOffsets.size() * sizeof(int32_t)
               + layerVertices.size() * sizeof(int32_t);
    }

private:
    // Applies func to each vertex in bottom-up fashion, starting from from and proceeding to to - 1.
    // That is, func is applied to a vertex after it has been applied to each downward neighbor. If
    // this member function is called in a parallel region, the function calls are parallelized.
    template<typename CallableT>
    void forEachVertexBottomUp(int from, int to, const int node, CallableT func) const {
        assert(to == decomp.lastSeparatorVertex(node));
        assert(from >= 0);
        assert(from <= to);
        const auto threshold = upGraph.numVertices() / (32 * omp_get_num_threads());
        if (to - from <= threshold || omp_get_num_threads() == 1) {
            for (auto v = from; v < to; ++v)
                func(v);
        } else {
            for (auto child = decomp.leftChild(node); child != 0; child = decomp.rightSibling(child)) {
#pragma omp task
                forEachVertexBottomUp(from, decomp.lastSeparatorVertex(child), child, func);
                from = decomp.lastSeparatorVertex(child);
            }
#pragma omp taskwait
            for (auto v = from; v < to; ++v)
                func(v);
        }
    }

    // Applies func to each vertex in top-down fashion, starting from from and proceeding to to - 1.
    // That is, func is applied to a vertex after it has been applied to each upward neighbor. If
    // this member function is called in a parallel region, the function calls are parallelized.
    template<typename CallableT>
    void forEachVertexTopDown(int from, int to, const int node, CallableT func) const {
        assert(to == decomp.lastSeparatorVertex(node));
        assert(from >= 0);
        assert(from <= to);
        const auto threshold = upGraph.numVertices() / (32 * omp_get_num_threads());
        if (to - from <= threshold || omp_get_num_threads() == 1) {
            for (auto v = to - 1; v >= from; --v)
                func(v);
        } else {
            for (auto v = to - 1; v >= decomp.firstSeparatorVertex(node); --v)
                func(v);
            for (auto child = decomp.leftChild(node); child != 0; child = decomp.rightSibling(child)) {
#pragma omp task
                forEachVertexTopDown(from, decomp.lastSeparatorVertex(child), child, func);
                from = decomp.lastSeparatorVertex(child);
            }
        }
    }

    void buildLayering() {
        layerOffsets.clear();
        layerVertices.clear();
        numLayers = 0;

        const int n = upGraph.numVertices();
        if (n == 0)
            return;

        std::vector<int32_t> vertexLayer(n, 0);
        int32_t maxLayer = 0;
        for (int u = 0; u < n; ++u) {
            const int32_t baseLayer = vertexLayer[u];

            const auto propagate = [&](const int v) {
                const int32_t candidate = baseLayer + 1;
                if (candidate > vertexLayer[v]) {
                    vertexLayer[v] = candidate;
                    if (candidate > maxLayer)
                        maxLayer = candidate;
                }
            };

            for (int e = upGraph.firstEdge(u); e != upGraph.lastEdge(u); ++e)
                propagate(upGraph.edgeHead(e));
        }

        numLayers = maxLayer + 1;

        layerOffsets.assign(numLayers + 1, 0);
        for (const auto layer : vertexLayer)
            ++layerOffsets[layer + 1];
        for (int i = 0; i < numLayers; ++i){
            layerOffsets[i + 1] += layerOffsets[i];
        }

        layerVertices.resize(n);
        auto writePos = layerOffsets;
        for (int v = 0; v < n; ++v){
            layerVertices[writePos[vertexLayer[v]]] = v;
            writePos[vertexLayer[v]]++;
        }
    }

    SeparatorDecomposition decomp;   // The separator decomposition used to build this CCH.
    Permutation ranks;               // The position of each vertex in the contraction order.
    EliminationTree eliminationTree; // The associated elimination tree.

    UpGraph upGraph;     // The upward graph.
    DownGraph downGraph; // The downward graph.

    std::vector<int32_t> firstUpInputEdge;   // The idx of the 1st upward input edge for each edge.
    std::vector<int32_t> firstDownInputEdge; // The idx of the 1st downward input edge for each edge.
    std::vector<int32_t> upInputEdges;       // The upward input edges.
    std::vector<int32_t> downInputEdges;     // The downward input edges.

    std::vector<int32_t> layerOffsets;  // Prefix sums delimiting per-layer vertex ranges.
    std::vector<int32_t> layerVertices; // Vertices sorted by layer (depth in separator tree).
    int32_t numLayers = 0;              // Number of layers in the separator tree.
};
