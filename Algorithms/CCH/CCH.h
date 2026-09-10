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
template<bool ORDER_BY_LEVEL>
class CCHBase {
public:
    using EliminationTree = std::vector<int32_t>; // The elimination tree.
    using UpGraph = StaticGraph<VertexAttrs<>, EdgeAttrs<EdgeTailAttribute> >; // The upward graph.
    using DownGraph = StaticGraph<VertexAttrs<>, EdgeAttrs<EdgeIdAttribute> >; // The downward graph.

    // Constructs an empty CCH.
    CCHBase() = default;

    // Constructs a CCH from the specified binary file.
    explicit CCHBase(std::ifstream &in) {
        readFrom(in);
    }

    // Builds the metric-independent CCH for the specified graph and separator decomposition.
    template<typename InputGraphT>
    void preprocess(const InputGraphT &inputGraph, const SeparatorDecomposition &sepDecomp) {
        assert(inputGraph.numVertices() == sepDecomp.order.size());
        std::vector<unsigned int> tmpOrder(sepDecomp.order.begin(), sepDecomp.order.end());
        std::vector<unsigned int> tails(inputGraph.numEdges());
        std::vector<unsigned int> heads(inputGraph.numEdges());
        FORALL_VALID_EDGES(inputGraph, u, e) {
            tails[e] = u;
            heads[e] = inputGraph.edgeHead(e);
        }
        RoutingKit::CustomizableContractionHierarchy cch(tmpOrder, tails, heads);

        decomp = sepDecomp;
        order.assign(tmpOrder.begin(), tmpOrder.end());
        KASSERT(order.validate());
        ranks.assign(cch.rank.begin(), cch.rank.end());
        KASSERT(ranks.validate());
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

        if constexpr (ORDER_BY_LEVEL)
            reorderByLayers();
    }

    // Returns the separator decomposition used to build this CCH.
    const SeparatorDecomposition &getSeparatorDecomposition() const noexcept {
        return decomp;
    }

    const Permutation &getSepDecompToCCHGraphVertexMapping() const noexcept requires (ORDER_BY_LEVEL) {
        return decompRankToCchGraphId;
    }

    // Returns the order in which vertices were contracted.
    const Permutation &getContractionOrder() const noexcept {
        return order;
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
    void forEachVertexBottomUp(CallableT func) const requires (!ORDER_BY_LEVEL) {
        #pragma omp parallel
        #pragma omp single nowait
        forEachVertexBottomUpTree(0, upGraph.numVertices(), 0, func);
    }

    template<class SeqCallable, class ParaCallable>
    void forEachVertexBottomUp(SeqCallable /* ignored */, ParaCallable func) const requires (!ORDER_BY_LEVEL) {
        #pragma omp parallel
        #pragma omp single nowait
        forEachVertexBottomUpTree(0, upGraph.numVertices(), 0, func);
    }


    template<class SeqCallable, class ParaCallable>
    void forEachVertexBottomUp(SeqCallable Seqfunc, ParaCallable Parafunc) const requires (ORDER_BY_LEVEL) {
        assert(numLayers != 0 || upGraph.numVertices() == 0);
        if (numLayers == 0)
            return;
        const int startSeqLayer = NumberOfThreads == 1 ? 0 : firstSequentialLayer;
#pragma omp parallel
        {
            for (auto layer = 0; layer < startSeqLayer; ++layer) {
#pragma omp for schedule(dynamic, 32)
                for (auto v = layerOffsets[layer]; v < layerOffsets[layer + 1]; ++v) {
                    Parafunc(v);
                }
            }
        }
        // Process top layers sequentially
        for (auto v = layerOffsets[startSeqLayer]; v < layerOffsets[numLayers]; ++v) {
            Seqfunc(v);
        }
    }

    template<class CallableT>
    void forEachVertexBottomUp(CallableT func) const requires (ORDER_BY_LEVEL) {
        forEachVertexBottomUp(func, func);
    }

    template<typename CallableT>
    void forEachVertexTopDown(CallableT func) const requires (!ORDER_BY_LEVEL) {
#pragma omp parallel
#pragma omp single nowait
        forEachVertexTopDownTree(0, upGraph.numVertices(), 0, func);
    }

    template<class SeqCallable, class ParaCallable>
    void forEachVertexTopDown(SeqCallable /*ignored*/, ParaCallable func) const requires (!ORDER_BY_LEVEL) {
#pragma omp parallel
#pragma omp single nowait
        forEachVertexTopDownTree(0, upGraph.numVertices(), 0, func);
    }

    template<class SeqCallable, class ParaCallable>
    void forEachVertexTopDown(SeqCallable Seqfunc, ParaCallable Parafunc) const requires (ORDER_BY_LEVEL) {
        assert(numLayers != 0 || upGraph.numVertices() == 0);
        const int lastSeqLayer = NumberOfThreads == 1 ? 0 : firstSequentialLayer;
        for (auto v = layerOffsets[numLayers] - 1; v >= layerOffsets[lastSeqLayer]; --v) {
            Seqfunc(v);
        }
#pragma omp parallel
        {
            for (auto layer = lastSeqLayer - 1; layer >= 0; --layer) {
#pragma omp for schedule(dynamic, 32)
                for (auto v = layerOffsets[layer + 1] - 1; v >= layerOffsets[layer]; --v) {
                    Parafunc(v);
                }
            }
        }
    }

    template<typename CallableT>
    void forEachVertexTopDown(CallableT func) const requires (ORDER_BY_LEVEL) {
        forEachVertexTopDown(func, func);
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
        order.readFrom(in);
        bio::read(in, eliminationTree);
        upGraph.readFrom(in);
        downGraph.readFrom(in);
        bio::read(in, firstUpInputEdge);
        bio::read(in, firstDownInputEdge);
        bio::read(in, upInputEdges);
        bio::read(in, downInputEdges);

        // Mirror preprocess(): the layered variant permutes after building, the plain one
        // does not. Reordering unconditionally here would hand back a DIFFERENTLY NUMBERED
        // CCH than the one that was written, for CCH = CCHBase<false>.
        if constexpr (ORDER_BY_LEVEL)
            reorderByLayers();
    }

    // Writes the CCH to the specified binary file.
    void writeTo(std::ofstream &out) const {
        decomp.writeTo(out);
        ranks.writeTo(out);
        order.writeTo(out);
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
               + layerOffsets.size() * sizeof(int32_t);
    }

private:
    // Applies func to each vertex in bottom-up fashion, starting from from and proceeding to to - 1.
    // That is, func is applied to a vertex after it has been applied to each downward neighbor. If
    // this member function is called in a parallel region, the function calls are parallelized.
    template<typename CallableT>
    void forEachVertexBottomUpTree(int from, int to, const int node, CallableT func) const {
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
                forEachVertexBottomUpTree(from, decomp.lastSeparatorVertex(child), child, func);
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
    void forEachVertexTopDownTree(int from, int to, const int node, CallableT func) const {
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
                forEachVertexTopDownTree(from, decomp.lastSeparatorVertex(child), child, func);
                from = decomp.lastSeparatorVertex(child);
            }
        }
    }

    void reorderByLayers() {
        const auto n = upGraph.numVertices();

        // Compute layer of every vertex, and number of layers
        const auto vertexLayer = computeLayerOfVertices();

        // Compute number of vertices per layer and generate start offsets by prefix sum.
        layerOffsets.clear();
        layerOffsets.assign(numLayers + 1, 0);
        for (const auto layer: vertexLayer)
            ++layerOffsets[layer];

        int sum = 0;
        for (int i = 0; i < numLayers; ++i) {
            std::swap(sum, layerOffsets[i]);
            sum += layerOffsets[i];
        }
        layerOffsets[numLayers] = sum;
        KASSERT(sum == n);

        // TODO: better order within layer?
        // Find new ID of every vertex
        Permutation vertexPerm(n);
        for (int v = 0; v < n; ++v)
            vertexPerm[v] = layerOffsets[vertexLayer[v]]++;
        KASSERT(vertexPerm.validate());
        // Restore offsets
        int cur = 0;
        for (int i = 0; i < numLayers; ++i)
            std::swap(cur, layerOffsets[i]);
        KASSERT(cur == n);


        // Permute upGraph and downGraph, retrieve mapping of edges in upGraph
        Permutation edgePerm;
        upGraph.permuteVertices(vertexPerm, edgePerm);
        downGraph.permuteVertices(vertexPerm);

        // Update mapping from downGraph edges to upGraph edges
        FORALL_EDGES(downGraph, e) {
            downGraph.edgeId(e) = edgePerm[downGraph.edgeId(e)];
        }

        // Update vertex IDs in ranks
        for (int v = 0; v < n; ++v) {
            ranks[v] = vertexPerm[ranks[v]];
        }

        // Update vertex IDs in order
        order = ranks.getInversePermutation();

        // Update vertex IDs in elimination tree
        EliminationTree newElimTree(n);
        for (int v = 0; v < n; ++v) {
            const int oldChild = v;
            const int oldParent = eliminationTree[v];
            if (oldParent == INVALID_VERTEX) {
                newElimTree[vertexPerm[oldChild]] = INVALID_VERTEX;
                continue;
            }
            newElimTree[vertexPerm[oldChild]] = vertexPerm[oldParent];
        }
        eliminationTree = std::move(newElimTree);

        for (int v = 0; v < n; ++v) {
            KASSERT(eliminationTree[v] == INVALID_VERTEX || v < eliminationTree[v]);
        }


        // Reorder input edge information for new edge IDs
        const auto newEdgeToOldEdge = edgePerm.getInversePermutation();
        std::vector<int32_t> firstInputEdge;
        std::vector<int32_t> inputEdges;
        const auto m = upGraph.numEdges();
        KASSERT(firstUpInputEdge.size() == m + 1);
        firstInputEdge.resize(firstUpInputEdge.size());
        firstInputEdge[0] = 0;
        inputEdges.reserve(upInputEdges.size());
        for (int e = 0; e < m; ++e) {
            const int oldEdge = newEdgeToOldEdge[e];
            const auto first = firstUpInputEdge[oldEdge];
            const auto last = firstUpInputEdge[oldEdge + 1];
            firstInputEdge[e + 1] = inputEdges.size() + last - first;
            inputEdges.insert(inputEdges.end(), upInputEdges.begin() + first, upInputEdges.begin() + last);
        }
        KASSERT(inputEdges.size() == upInputEdges.size());
        KASSERT(firstInputEdge.back() == inputEdges.size());
        firstUpInputEdge = std::move(firstInputEdge);
        upInputEdges = std::move(inputEdges);

        inputEdges.clear();
        KASSERT(firstDownInputEdge.size() == m + 1);
        firstInputEdge.resize(firstDownInputEdge.size());
        firstInputEdge[0] = 0;
        inputEdges.reserve(downInputEdges.size());
        for (int e = 0; e < m; ++e) {
            const int oldEdge = newEdgeToOldEdge[e];
            const auto first = firstDownInputEdge[oldEdge];
            const auto last = firstDownInputEdge[oldEdge + 1];
            firstInputEdge[e + 1] = inputEdges.size() + last - first;
            inputEdges.insert(inputEdges.end(), downInputEdges.begin() + first, downInputEdges.begin() + last);
        }
        KASSERT(inputEdges.size() == downInputEdges.size());
        KASSERT(firstInputEdge.back() == inputEdges.size());
        firstDownInputEdge = std::move(firstInputEdge);
        downInputEdges = std::move(inputEdges);

        // Remember vertex permutation to allow mapping from separator decomposition to new graph ordering.
        decompRankToCchGraphId = std::move(vertexPerm);


        // Set first layer to be processed sequentially based on size
        firstSequentialLayer = 0;
        while (firstSequentialLayer < numLayers && layerOffsets[firstSequentialLayer + 1] - layerOffsets[
                   firstSequentialLayer] >= MIN_SIZE_PARALLEL_LAYER)
            ++firstSequentialLayer;
    }

    std::vector<int32_t> computeLayerOfVertices() {
        numLayers = 0;

        const int n = upGraph.numVertices();

        std::vector<int32_t> vertexLayer(n, 0);
        int32_t maxLayer = 0;
        for (int u = 0; u < n; ++u) {
            const int32_t baseLayer = vertexLayer[u];

            for (int e = upGraph.firstEdge(u); e != upGraph.lastEdge(u); ++e) {
                const auto v = upGraph.edgeHead(e);
                const int32_t candidate = baseLayer + 1;
                if (candidate > vertexLayer[v]) {
                    vertexLayer[v] = candidate;
                    if (candidate > maxLayer)
                        maxLayer = candidate;
                }
            }
        }

        numLayers = maxLayer + 1;
        return vertexLayer;
    }

    SeparatorDecomposition decomp; // The separator decomposition used to build this CCH.

    // Permutation from vertex IDs in the separator decomposition to vertex IDs in the CCH graph.
    // Only used if CCH graph is reordered according to layers.
    Permutation decompRankToCchGraphId;

    // The contraction order. If this is a layer CCH, this is not the order that vertices were contracted in but it
    // contains all vertices of the input graph ordered by layer (lower layers first).
    // Use this to map vertex IDs in the CCH graph to vertex IDs in the input graph.
    Permutation order;

    // The position of each vertex in the contraction order. If this is a layer CCH, this is not the position in the
    // contraction order but the position in an order by layer (lower layers first).
    // Use this to map vertex IDs in the input graph to vertex IDs in the CCH graph.
    Permutation ranks;

    EliminationTree eliminationTree; // The associated elimination tree.

    UpGraph upGraph; // The upward graph.
    DownGraph downGraph; // The downward graph.

    std::vector<int32_t> firstUpInputEdge; // The idx of the 1st upward input edge for each edge.
    std::vector<int32_t> firstDownInputEdge; // The idx of the 1st downward input edge for each edge.
    std::vector<int32_t> upInputEdges; // The upward input edges.
    std::vector<int32_t> downInputEdges; // The downward input edges.

    int32_t numLayers = 0; // Number of layers in the CCH graph. Only used if ordered by layers.
    std::vector<int32_t> layerOffsets;
    // Prefix sums delimiting per-layer vertex ranges. Only used if ordered by layers.

    static constexpr int NumberOfThreads = NUM_THREADS;
    static constexpr int MIN_SIZE_PARALLEL_LAYER = NumberOfThreads * 32;
    int32_t firstSequentialLayer;
};

using CCH = CCHBase<false>;
using LayerCCH = CCHBase<true>;
