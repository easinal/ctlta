#pragma once

// Computes metric-independent preprocessing for CTNR.
class CTNRPreprocessor {

public:

    CTNRPreprocessor(const TransitNodeHierarchy &hierarchy, const CCH &cch)
            : hierarchy(hierarchy), cch(cch) {
        // Convert elimination tree to out-tree format
        convertInTreeToOutTree(cch.getEliminationTree(), elimTreeFirstChild, elimTreeChildren);



    }

    // Determines access nodes of each vertex and allocates distance entries.
    void preprocess(CTNRData &data) const {

        int numVertices = cch.getUpwardGraph().numVertices();

        const auto rankToIdx = [&](const int r) {
            return data.rankToIdx(r);
        };

        // Debug output:
        std::vector<int32_t> numTransitNodesToRoot(cch.getUpwardGraph().numVertices());
        countTransitNodesToRoot(numTransitNodesToRoot, rankToIdx);
        std::cout << "CTNR: Average number of transit nodes to root: "
                  << static_cast<double>(std::accumulate(numTransitNodesToRoot.begin(), numTransitNodesToRoot.end(), 0))
                     / numVertices << std::endl;

        // Count number of access nodes per vertex
        data.pos.resize(numVertices + 1);
        countMetricIndependentAccessNodes(data.pos, cch.getUpwardGraph(), rankToIdx);

        if constexpr (CTNRData::USE_SIMD) {
            // Pad counts to multiple of SIMD width
            for (int32_t &count: data.pos) {
                count = roundUpToMultiple(count, CTNRData::K);
            }
        }

        // Prefix sum to get offsets. We can then be sure that the range for a vertex is able to
        // accommodate all access nodes of vertex v.
        int32_t sum = 0;
        for (int32_t i = 0; i < numVertices; ++i) {
            const int32_t size = data.pos[i];
            data.pos[i] = sum;
            sum += size;
        }
        data.pos[numVertices] = sum;

        // Allocate space for access nodes and distances
        data.accessNodes.resize(sum);
        data.forwardDistances.resize(sum);
        data.backwardDistances.resize(sum);

        writeMetricIndependentAccessNodes(data.pos, data.accessNodes, cch.getUpwardGraph(), rankToIdx);

        // TODO: remove debug
        std::cout << "CTNR: Average access nodes per vertex: "
                  << static_cast<double>(data.accessNodes.size()) / numVertices << std::endl;
    }

    uint64_t sizeInBytes() const {
        uint64_t size = sizeof(CTNRPreprocessor);
        size += elimTreeFirstChild.size() * sizeof(int);
        size += elimTreeChildren.size() * sizeof(int);
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

        int id;                 // The ID of the active vertex.
        int nextUnexploredEdge; // The next unexplored incident edge.
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
                activeVertices.pop();
                if (!activeVertices.empty())
                    backtrack(v.id, activeVertices.top().id);
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
                                           const RankToIdxT &rankToIdx) const {

        KASSERT(outCounts.size() == upGraph.numVertices() + 1);

        const int root = upGraph.numVertices() - 1;
        if (hierarchy.isTransitNode(root)) {
            outCounts[rankToIdx(root)] = 1;
        }

        std::stack<int, std::vector<int>> numAdded;
        BitVector isActive(hierarchy.numTransitNodes());
        std::stack<int, std::vector<int>> active;
        const auto recurse = [&](const int /*parent*/, const int child) {
            numAdded.push(0);
            if (hierarchy.isTransitNode(child)) {
                outCounts[rankToIdx(child)] = 1;
                return;
            }
            FORALL_INCIDENT_EDGES(upGraph, child, e) {
                const int neighbor = upGraph.edgeHead(e);
                if (!hierarchy.isTransitNode(neighbor))
                    continue;
                const int node = hierarchy.getTransitNodeIndexOfRank(neighbor);
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
    void writeMetricIndependentAccessNodes(const std::vector<int32_t> &pos, std::vector<int32_t> &entries,
                                           const GraphT &upGraph,
                                           const RankToIdxT &rankToIdx) const {

        KASSERT(pos.size() == upGraph.numVertices() + 1);
        std::vector<int> curNumEntries(upGraph.numVertices(), 0);

        const int root = upGraph.numVertices() - 1;
        if (hierarchy.isTransitNode(root)) {
            const int idx = rankToIdx(root);
            entries[pos[idx]] = hierarchy.getTransitNodeIndexOfRank(root);
            ++curNumEntries[idx];
        }

        const auto recurse = [&](const int parent, const int child) {
            const int childIdx = rankToIdx(child);
            if (hierarchy.isTransitNode(child)) {
                entries[pos[childIdx]] = hierarchy.getTransitNodeIndexOfRank(child);
                ++curNumEntries[childIdx];
                return;
            }

            // Add all transit nodes from parent to child
            const int parentIdx = rankToIdx(parent);
            KASSERT(curNumEntries[parentIdx] <= pos[parentIdx + 1] - pos[parentIdx] &&
                    curNumEntries[parentIdx] >= pos[parentIdx + 1] - pos[parentIdx] - CTNRData::K);
            for (int i = 0; i < curNumEntries[parentIdx]; ++i) {
                entries[pos[childIdx] + i] = entries[pos[parentIdx] + i];
            }
            curNumEntries[childIdx] = curNumEntries[parentIdx];

            // If any upward neighbor is a transit node that has not been seen on this branch, add it to back of list
            FORALL_INCIDENT_EDGES(upGraph, child, e) {
                const int neighbor = upGraph.edgeHead(e);
                if (!hierarchy.isTransitNode(neighbor))
                    continue;
                const int node = hierarchy.getTransitNodeIndexOfRank(neighbor);
                if (contains(entries.begin() + pos[childIdx],
                             entries.begin() + pos[childIdx] + curNumEntries[childIdx], node))
                    continue;
                entries[pos[childIdx] + curNumEntries[childIdx]] = node;
                ++curNumEntries[childIdx];
            }
            KASSERT(curNumEntries[childIdx] <= pos[childIdx + 1] - pos[childIdx] &&
                    curNumEntries[childIdx] >= pos[childIdx + 1] - pos[childIdx] - CTNRData::K);
        };

        const auto backtrack = [&](const int /*child*/, const int /*parent*/) {
            // no op
        };

        // TODO: parallelize
        dfsOnTree(elimTreeFirstChild, elimTreeChildren, recurse, backtrack);
    }


    template<typename RankToIdxT>
    void countTransitNodesToRoot(std::vector<int32_t> &outCounts,
                                const RankToIdxT &rankToIdx) const {

        int curCount = 0;
        const int root = outCounts.size() - 1;
        if (hierarchy.isTransitNode(root)) {
            ++curCount;
        }
        outCounts[rankToIdx(root)] = curCount;

        const auto recurse = [&](const int /*parent*/, const int child) {
            if (hierarchy.isTransitNode(child)) {
                ++curCount;
            }
            outCounts[rankToIdx(child)] = curCount;
        };

        const auto backtrack = [&](const int child, const int /*parent*/) {
            if (hierarchy.isTransitNode(child)) {
                --curCount;
            }
        };

        dfsOnTree(elimTreeFirstChild, elimTreeChildren, recurse, backtrack);
    }


    const TransitNodeHierarchy &hierarchy;
    const CCH &cch;

    // CCH elimination as out-tree
    std::vector<int> elimTreeFirstChild;
    std::vector<int> elimTreeChildren;

};