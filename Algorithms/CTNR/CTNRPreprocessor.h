#pragma once

// Computes metric-independent preprocessing for CTNR.
class CTNRPreprocessor {

public:

    CTNRPreprocessor(const TransitNodeHierarchy &hierarchy, const CCH &cch)
            : hierarchy(hierarchy), cch(cch) {
        // Convert elimination tree to out-tree format
        convertInTreeToOutTree(cch.getEliminationTree(), elimTreeFirstChild, elimTreeChildren);
    }

    const std::vector<AccessNodeEdge>& getAccessNodeEdges() const {
        return accessNodeEdges;
    }

    // Determines access nodes of each vertex and allocates distance entries.
    void preprocess(CTNRData &data) {

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

        std::sort(accessNodeEdges.begin(), accessNodeEdges.end(), [](const AccessNodeEdge &a, const AccessNodeEdge &b) {
            return a.edge < b.edge;
        });

        // TODO: remove debug
        std::cout << "CTNR: Average access nodes per vertex: "
                  << static_cast<double>(data.accessNodes.size()) / numVertices << std::endl;
    }

    uint64_t sizeInBytes() const {
        uint64_t size = sizeof(CTNRPreprocessor);
        size += elimTreeFirstChild.size() * sizeof(decltype(elimTreeFirstChild)::value_type);
        size += elimTreeChildren.size() * sizeof(decltype(elimTreeChildren)::value_type);
        size += accessNodeEdges.size() * sizeof(decltype(accessNodeEdges)::value_type);
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
            outCounts[rankToIdx(root)] = 0;
        }

        std::stack<int, std::vector<int>> numAdded;
        BitVector isActive(hierarchy.numTransitNodes());
        std::stack<int, std::vector<int>> active;
        const auto recurse = [&](const int /*parent*/, const int child) {
            numAdded.push(0);
            if (hierarchy.isTransitNode(child)) {
                outCounts[rankToIdx(child)] = 0;
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
    void writeMetricIndependentAccessNodes(const std::vector<int32_t> &pos, std::vector<ctnr::TransitNodeId> &entries,
                                           const GraphT &upGraph,
                                           const RankToIdxT &rankToIdx) {

        KASSERT(pos.size() == upGraph.numVertices() + 1);
        std::vector<int> curNumEntries(upGraph.numVertices(), 0);

        const auto recurse = [&](const int parent, const int child) {
            const int childIdx = rankToIdx(child);
            const int startChild = pos[childIdx];
            if (hierarchy.isTransitNode(child)) {
                return;
            }

            // Add all transit nodes from parent to child
            const int parentIdx = rankToIdx(parent);
            KASSERT(curNumEntries[parentIdx] <= pos[parentIdx + 1] - pos[parentIdx] &&
                    curNumEntries[parentIdx] >= pos[parentIdx + 1] - pos[parentIdx] - CTNRData::K);
            for (int i = 0; i < curNumEntries[parentIdx]; ++i) {
                entries[startChild + i] = entries[pos[parentIdx] + i];
            }
            curNumEntries[childIdx] = curNumEntries[parentIdx];

            // If any upward neighbor is a transit node that has not been seen on this branch, add it to back of list
            FORALL_INCIDENT_EDGES(upGraph, child, e) {
                const int neighbor = upGraph.edgeHead(e);
                if (!hierarchy.isTransitNode(neighbor))
                    continue;
                const ctnr::TransitNodeId node = hierarchy.getTransitNodeIndexOfRank(neighbor);
                int i = 0;
                for (; i < curNumEntries[childIdx]; ++i) {
                    if (entries[startChild + i] == node) {
                        break;
                    }
                }
                // Mark access node edge:
                accessNodeEdges.emplace_back(e, startChild + i);
                if (i < curNumEntries[childIdx])
                    continue; // Neighbor is already an access node
                // New access node
                entries[startChild + curNumEntries[childIdx]] = node;
                ++curNumEntries[childIdx];
            }
            KASSERT(curNumEntries[childIdx] <= pos[childIdx + 1] - startChild &&
                    curNumEntries[childIdx] >= pos[childIdx + 1] - startChild - CTNRData::K);
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

    // Information on edges leading from non-transit nodes directly to access nodes which is a special case during
    // customization.
    std::vector<AccessNodeEdge> accessNodeEdges;

};