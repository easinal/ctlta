#pragma once

#include <parallel/algorithm>

#include <limits>
#include <stdexcept>

// FHL's metric-independent access-node pass: for every vertex, which hubs its
// shortest paths can first reach, and where their distances live. FHL consumes this
// only as SCAFFOLDING -- the seed fold turns it into per-region seeds during
// customization and then releases it (see the [cust] releases phase).
// Forked from CTNR's copy on purpose: the classic table keeps these arrays alive as
// its query structure, FHL throws them away, so the two are free to diverge.
#include "Algorithms/FHL/core/HubHierarchy.h"
#include "Algorithms/FHL/core/HubAccessData.h"
#include "Algorithms/FHL/core/FHLConstants.h"
#include "Tools/Timer.h"

// Computes FHL's metric-independent preprocessing.
class HubAccessPreprocessor {

public:

    HubAccessPreprocessor() = default;

    const std::vector<HubAccessEdge>& getAccessHubEdges() const {
        return accessHubEdges;
    }

    // Determines access nodes of each vertex and allocates distance entries.
    template<typename CchT>
    void preprocess(const HubHierarchy &hierarchy, const CchT& cch, HubAccessData &data) {

        int numVertices = cch.getUpwardGraph().numVertices();
        std::vector<int32_t> elimTreeFirstChild;
        std::vector<int32_t> elimTreeChildren;
        convertInTreeToOutTree(cch.getEliminationTree(), elimTreeFirstChild, elimTreeChildren);


        const auto rankToIdx = [&](const int r) {
            return data.rankToIdx(r);
        };

        // Debug output:
        std::vector<int32_t> numHubsToRoot(cch.getUpwardGraph().numVertices());
        countHubsToRoot(numHubsToRoot, rankToIdx, hierarchy, elimTreeFirstChild, elimTreeChildren);
        std::cout << "FHL: average hubs to root: "
                  << static_cast<double>(std::accumulate(numHubsToRoot.begin(), numHubsToRoot.end(), 0))
                     / numVertices << std::endl;

        // Count number of access nodes per vertex
        data.pos.resize(numVertices + 1);
        countMetricIndependentAccessNodes(data.pos, cch.getUpwardGraph(), rankToIdx, hierarchy, elimTreeFirstChild, elimTreeChildren);

        if constexpr (HubAccessData::USE_SIMD) {
            // Pad counts to multiple of SIMD width
            for (int32_t &count: data.pos) {
                count = roundUpToMultiple(count, HubAccessData::K);
            }
        }

        // Prefix sum to get offsets. We can then be sure that the range for a vertex is able to
        // accommodate all access nodes of vertex v.
        // pos is int32-indexed by design, so the total entry count must fit. It does in
        // every size-cut configuration (20 access hubs per vertex on USA); a too-shallow
        // level cut instead leaves hundreds of hubs reachable from every vertex and blows
        // past the limit, which used to surface as an opaque resize failure.
        int64_t sum = 0;
        for (int32_t i = 0; i < numVertices; ++i) {
            const int32_t size = data.pos[i];
            data.pos[i] = static_cast<int32_t>(sum);
            sum += size;
            if (sum > std::numeric_limits<int32_t>::max())
                throw std::overflow_error(
                        "FHL: too many access hubs to index with 32 bits -- cut deeper "
                        "(-fhl-region-size) or raise the hub level threshold");
        }
        data.pos[numVertices] = static_cast<int32_t>(sum);

        // Allocate space for access nodes and distances
        data.accessHubs.resize(sum);
        data.forwardDistances.resize(sum);
        data.backwardDistances.resize(sum);

        writeMetricIndependentAccessNodes(data.pos, data.accessHubs, cch.getUpwardGraph(), rankToIdx, hierarchy, elimTreeFirstChild, elimTreeChildren);

        // Tens of millions of entries at continental scale: the largest single sort in the
        // module, and parallelising it is nearly free. Preprocessing only -- metric-
        // independent, so it never shows up in the customization numbers.
        __gnu_parallel::sort(accessHubEdges.begin(), accessHubEdges.end(), [](const HubAccessEdge &a, const HubAccessEdge &b) {
            return a.edge < b.edge;
        });

        std::cout << "FHL: average access hubs per vertex: "
                  << static_cast<double>(data.accessHubs.size()) / numVertices << std::endl;
    }

    uint64_t sizeInBytes() const {
        uint64_t size = sizeof(HubAccessPreprocessor);
        size += accessHubEdges.size() * sizeof(decltype(accessHubEdges)::value_type);
        return size;
    }


private:

    int roundUpToMultiple(int value, int multiple) const {
        return ((value + multiple - 1) / multiple) * multiple;
    }

    // An active vertex during a DFS, i.e., a vertex that has been reached but not finished.
    struct ActiveVertex {
        // Constructs an active vertex.
        ActiveVertex(const int id, const int nextUnexploredEdge)
                : id(id), nextUnexploredEdge(nextUnexploredEdge) {}

        int id = INVALID_ID;                 // The ID of the active vertex.
        int nextUnexploredEdge = INVALID_EDGE; // The next unexplored incident edge.
    };

    static void convertInTreeToOutTree(const std::vector<int> &parent,
                                       std::vector<int> &firstChild,
                                       std::vector<int> &children) {
        const auto numVertices = parent.size();
        // Build the elimination out-tree from the elimination in-tree.
        firstChild = std::vector<int>(numVertices + 1);
        children = std::vector<int>(numVertices - 1);
        for (auto v = 0; v < numVertices; ++v) {
            const auto p = parent[v];
            if (p == -1 || p == v) // Root of tree may be marked by -1 or edge to self
                continue;
            ++firstChild[parent[v]];
        }
        auto firstEdge = 0; // The index of the first edge out of the current/next vertex.
        for (auto v = 0; v <= numVertices; ++v) {
            std::swap(firstEdge, firstChild[v]);
            firstEdge += firstChild[v];
        }
        for (auto v = 0; v < numVertices; ++v) {
            const auto p = parent[v];
            if (p == -1 || p == v) // Root of tree may be marked by -1 or edge to self
                continue;
            children[firstChild[p]++] = v;
        }
        for (auto v = numVertices - 1; v > 0; --v)
            firstChild[v] = firstChild[v - 1];
        firstChild[0] = 0;
    }

    // Run a DFS for the given tree in out format.
    // Call callbacks when recursing or backtracking.
    template<typename RecurseCallBack,
            typename BacktrackCallBack>
    static void dfsOnTree(
            const std::vector<int> &firstChild,
            const std::vector<int> &children,
            RecurseCallBack recurse,
            BacktrackCallBack backtrack) {
        const int numVertices = static_cast<int>(firstChild.size()) - 1;
        std::stack<ActiveVertex, std::vector<ActiveVertex>> activeVertices;
        activeVertices.emplace(numVertices - 1, firstChild[numVertices - 1]); // add root
        while (!activeVertices.empty()) {
            auto &v = activeVertices.top();
            if (v.nextUnexploredEdge == firstChild[v.id + 1]) {
                const int vId = v.id; // save before pop
                activeVertices.pop();
                if (!activeVertices.empty())
                    backtrack(vId, activeVertices.top().id);
                continue;
            }
            // Advance to next child
            const auto child = children[v.nextUnexploredEdge];
            recurse(v.id, child);
            ++v.nextUnexploredEdge; // When backtracking from child later, look at next sibling
            activeVertices.emplace(child, firstChild[child]);
        }
    }

    template<typename GraphT, typename RankToIdxT>
    void countMetricIndependentAccessNodes(std::vector<int32_t> &outCounts, const GraphT &upGraph,
    const RankToIdxT &rankToIdx,
    const HubHierarchy &hierarchy,
    const std::vector<int32_t> &elimTreeFirstChild,
    const std::vector<int32_t> &elimTreeChildren) const {

        KASSERT(outCounts.size() == upGraph.numVertices() + 1);

        const int root = upGraph.numVertices() - 1;
        if (hierarchy.isHub(root)) {
            outCounts[rankToIdx(root)] = 0;
        }

        std::stack<int, std::vector<int>> numAdded;
        BitVector isActive(hierarchy.numHubs());
        std::stack<int, std::vector<int>> active;
        const auto recurse = [&](const int /*parent*/, const int child) {
            numAdded.push(0);
            if (hierarchy.isHub(child)) {
                outCounts[rankToIdx(child)] = 0;
                return;
            }
            FORALL_INCIDENT_EDGES(upGraph, child, e) {
                const int neighbor = upGraph.edgeHead(e);
                if (!hierarchy.isHub(neighbor))
                    continue;
                const int node = hierarchy.hubIdOfRank(neighbor);
                if (isActive[node])
                    continue;
                isActive[node] = true;
                active.push(node);
                ++numAdded.top();
            }
            outCounts[rankToIdx(child)] = static_cast<int32_t>(active.size());
        };

        const auto backtrack = [&](const int /*child*/, const int /*parent*/) {
            for (int i = 0; i < numAdded.top(); ++i) {
                const int node = active.top();
                isActive[node] = false;
                active.pop();
            }
            numAdded.pop();
        };

        dfsOnTree(elimTreeFirstChild, elimTreeChildren, recurse, backtrack);
    }

    template<typename GraphT, typename RankToIdxT>
    void writeMetricIndependentAccessNodes(const std::vector<int32_t> &pos, std::vector<fhl::HubId> &entries,
                                           const GraphT &upGraph,
                                           const RankToIdxT &rankToIdx,
    const HubHierarchy &hierarchy,
                                const std::vector<int32_t> &elimTreeFirstChild,
                                const std::vector<int32_t> &elimTreeChildren) {

        KASSERT(pos.size() == upGraph.numVertices() + 1);
        std::vector<int> curNumEntries(upGraph.numVertices(), 0);

        const auto recurse = [&](const int parent, const int child) {
            const int childIdx = rankToIdx(child);
            const int startChild = pos[childIdx];
            if (hierarchy.isHub(child)) {
                return;
            }

            // Add all hubs from parent to child
            const int parentIdx = rankToIdx(parent);
            KASSERT(curNumEntries[parentIdx] <= pos[parentIdx + 1] - pos[parentIdx] &&
                    curNumEntries[parentIdx] >= pos[parentIdx + 1] - pos[parentIdx] - HubAccessData::K);
            for (int i = 0; i < curNumEntries[parentIdx]; ++i) {
                entries[startChild + i] = entries[pos[parentIdx] + i];
            }
            curNumEntries[childIdx] = curNumEntries[parentIdx];

            // If any upward neighbor is a hub that has not been seen on this branch, add it to back of list
            FORALL_INCIDENT_EDGES(upGraph, child, e) {
                const int neighbor = upGraph.edgeHead(e);
                if (!hierarchy.isHub(neighbor))
                    continue;
                const fhl::HubId node = hierarchy.hubIdOfRank(neighbor);
                int i = 0;
                for (; i < curNumEntries[childIdx]; ++i) {
                    if (entries[startChild + i] == node) {
                        break;
                    }
                }
                // Mark access node edge:
                accessHubEdges.emplace_back(e, startChild + i);
                if (i < curNumEntries[childIdx])
                    continue; // Neighbor is already an access node
                // New access node
                entries[startChild + curNumEntries[childIdx]] = node;
                ++curNumEntries[childIdx];
            }
            KASSERT(curNumEntries[childIdx] <= pos[childIdx + 1] - startChild &&
                    curNumEntries[childIdx] >= pos[childIdx + 1] - startChild - HubAccessData::K);
        };

        const auto backtrack = [&](const int /*child*/, const int /*parent*/) {
            // no op
        };

        // Serial on purpose: one pass over the elimination tree, a negligible share here.
        dfsOnTree(elimTreeFirstChild, elimTreeChildren, recurse, backtrack);
    }


    template<typename RankToIdxT>
    void countHubsToRoot(std::vector<int32_t> &outCounts,
                                const RankToIdxT &rankToIdx,
    const HubHierarchy &hierarchy,
                                const std::vector<int32_t> &elimTreeFirstChild,
                                const std::vector<int32_t> &elimTreeChildren) const {

        int curCount = 0;
        const int root = outCounts.size() - 1;
        if (hierarchy.isHub(root)) {
            ++curCount;
        }
        outCounts[rankToIdx(root)] = curCount;

        const auto recurse = [&](const int /*parent*/, const int child) {
            if (hierarchy.isHub(child)) {
                ++curCount;
            }
            outCounts[rankToIdx(child)] = curCount;
        };

        const auto backtrack = [&](const int child, const int /*parent*/) {
            if (hierarchy.isHub(child)) {
                --curCount;
            }
        };

        dfsOnTree(elimTreeFirstChild, elimTreeChildren, recurse, backtrack);
    }

    // Information on edges leading from non-hubs directly to access nodes which is a special case during
    // customization.
    std::vector<HubAccessEdge> accessHubEdges;

};
