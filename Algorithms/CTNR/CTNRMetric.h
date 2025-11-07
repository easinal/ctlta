#pragma once

#include "Algorithms/CTL/BalancedTopologyCentricTreeHierarchy.h"
#include "Algorithms/CCH/CCH.h"
#include "Algorithms/CCH/CCHMetric.h"
#include "Algorithms/CH/CH.h"
#include "DataStructures/Partitioning/SeparatorDecomposition.h"
#include "Tools/Constants.h"
#include "Algorithms/CTNR/CTNRData.h"
#include "TransitNodeHierarchy.h"
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <memory>
#include "Algorithms/CCH/EliminationTreeQuery.h"
#include "Algorithms/CH/CHQuery.h"
#include "DataStructures/Labels/BasicLabelSet.h"
#include "DataStructures/Labels/ParentInfo.h"
#include <algorithm>
#include <iostream>

class CTNRMetric {

    struct IndexRange {
        int32_t start = INVALID_INDEX;
        int32_t end = INVALID_INDEX;
    };

    class AccessNodeUnifier {

    public:

        explicit AccessNodeUnifier(const int numTransitNodes) : nodes(), distances(numTransitNodes, INFTY) {
            nodes.reserve(numTransitNodes);
        }

        void addAccessNode(const int32_t nodeIndex, const int32_t distance) {
            int32_t &entry = distances[nodeIndex];
            if (entry == INFTY) {
                // New entry
                entry = distance;
                nodes.push_back(nodeIndex);
            } else {
                // Existing entry for this vertex
                if (distance < entry) {
                    entry = distance;
                }
            }
        }

        template<typename It>
        void flushAccessNodes(It outRange) {
            std::sort(nodes.begin(), nodes.end());
            int next = 0;
            for (const int &node: nodes) {
                int32_t &entry = distances[node];
                CTNRData::AccessNode an(node, entry);
                outRange[next++] = an;
                entry = INFTY; // reset for next use
            }
            nodes.clear();
        }

        size_t sizeOfUnion() const {
            return nodes.size();
        }

    private:
        std::vector<int> nodes; // list of transit node indices
        std::vector<int32_t> distances;
    };

public:

    // Constructor
    CTNRMetric(const TransitNodeHierarchy &hierarchy, const CCH &cch, const int32_t *const inputWeights)
            : hierarchy(hierarchy), cch(cch), cchMetric(cch, inputWeights)
//            ,
//              forwardRange(cch.getUpwardGraph().numVertices()), backwardRange(cch.getUpwardGraph().numVertices())
    {

        // Todo: remove debug
//        // Find depth in elimination tree of every transit node
//        std::vector<int> depth(hierarchy.numTransitNodes(), 0);
//        const auto &parent = cch.getEliminationTree();
//        for (const auto& n : hierarchy.getTransitNodes()) {
//            int c = n;
//            const auto idxN = hierarchy.getTransitNodeIndexOfRank(n);
//            while (c != INVALID_INDEX) {
//                c = parent[c];
//                depth[idxN]++;
//            }
//        }
//
//        // Print average and max depth of transit nodes
//        double avgDepth = 0.0;
//        int maxDepth = 0;
//        for (const auto& d : depth) {
//            avgDepth += d;
//            if (d > maxDepth) {
//                maxDepth = d;
//            }
//        }
//        avgDepth /= depth.size();
//        std::cout << "CTNR: Average elimination tree depth of transit nodes: " << avgDepth << std::endl;
//        std::cout << "CTNR: Maximum elimination tree depth of transit nodes: " << maxDepth << std::endl;


//        const auto& upGraph = cch.getUpwardGraph();
//        std::vector<int> cchLevel(upGraph.numVertices(), 0);
//        cch.forEachVertexBottomUp([&](int32_t v) {
//            const int levelV = cchLevel[v];
//            FORALL_INCIDENT_EDGES(upGraph, v, e) {
//                const int neighbor = upGraph.edgeHead(e);
//                if (levelV + 1 > cchLevel[neighbor]) {
//                    cchLevel[neighbor] = levelV + 1;
//                }
//            }
//        });
//        // Count number of vertices in each CCH level and print
//        std::map<int, size_t> levelCounts;
//        for (const auto &level: cchLevel) {
//            ++levelCounts[level];
//        }
//        for (const auto& [level, count] : levelCounts) {
//            std::cout << "CCH Level " << level << ": " << count << " vertices" << std::endl;
//        }

    }

    // Customization phase
    void customize(CTNRData &data) {
        int64_t dummy1, dummy2, dummy3;
        customizeWithMeasurements(data, dummy1, dummy2, dummy3);
    }

    // Sets measurement parameters to times for each step in microseconds.
    void customizeWithMeasurements(CTNRData &data, int64_t &cchCustomizationTime, int64_t &accessNodeComputationTime,
                                   int64_t &distanceTableComputationTime) {
        Timer timer;
        minCH = cchMetric.buildMinimumWeightedCH();
        cchCustomizationTime = timer.elapsed<std::chrono::microseconds>();
        timer.restart();
        computeDistanceTable(data);
        distanceTableComputationTime = timer.elapsed<std::chrono::microseconds>();
        timer.restart();
        computeAccessNodes(data);
        accessNodeComputationTime = timer.elapsed<std::chrono::microseconds>();
    }

    const CH &getMinCH() const { return minCH; }

    // Memory usage calculation including node levels
    uint64_t sizeInBytes() const {
        uint64_t size = sizeof(CTNRMetric);
        size += cchMetric.sizeInBytes();
        size += minCH.sizeInBytes();

        return size;
    }

private:

    void computeAccessNodes(CTNRData &data) {

        int numVertices = cch.getUpwardGraph().numVertices();

        // Compute upper bound on number of forward/backward access nodes per vertex
        const auto rankToIdx = [&](const int r) {
            return data.rankToIdx(r);
        };
        std::vector<int> maxNumForward(numVertices + 1, 0);
        countTransitNodesInSearchSpace(maxNumForward, minCH.upwardGraph(), rankToIdx);
        std::vector<int> maxNumBackward(numVertices + 1, 0);
        countTransitNodesInSearchSpace(maxNumBackward, minCH.downwardGraph(), rankToIdx);

        // Prefix sum over maxNumForward/maxNumBackward to get initial (suboptimal) offsets. We can be sure that the
        // range maxNumForward[rankToIdx(v)]..maxNumForward[rankToIdx(v)+1]-1 (similarly for backward) is able to
        // accommodate all access nodes of vertex v.
        int32_t forwardSum = 0;
        int32_t backwardSum = 0;
        for (int32_t i = 0; i < numVertices; ++i) {
            const int32_t fSize = maxNumForward[i];
            const int32_t bSize = maxNumBackward[i];
            maxNumForward[i] = forwardSum;
            maxNumBackward[i] = backwardSum;
            forwardSum += fSize;
            backwardSum += bSize;
        }
        maxNumForward[numVertices] = forwardSum;
        maxNumBackward[numVertices] = backwardSum;

        // Count the actual number of access nodes per vertex in data.forwardPos/data.backwardPos. Later, a prefix sum
        // in these vectors will give actual offsets into flat representation without gaps.
        data.forwardPos.clear();
        data.backwardPos.clear();
        data.forwardPos.resize(numVertices + 1, INVALID_INDEX);
        data.backwardPos.resize(numVertices + 1, INVALID_INDEX);
        data.forwardAccess.clear();
        data.backwardAccess.clear();
        data.forwardAccess.resize(forwardSum, CTNRData::AccessNode());
        data.backwardAccess.resize(backwardSum, CTNRData::AccessNode());

        // Collect ranges of access nodes into these temporary vectors first with arbitrary order of vertices.
//        const auto &upGraph = minCH.upwardGraph();
//        const auto &downGraph = minCH.downwardGraph();

#pragma omp parallel
#pragma omp single nowait
        cch.forEachVertexTopDown([&](int32_t rv) {
            const int idx = data.rankToIdx(rv);
            if (hierarchy.isTransitNode(rv)) {
                KASSERT(maxNumForward[idx + 1] - maxNumForward[idx] >= 1);
                KASSERT(maxNumBackward[idx + 1] - maxNumBackward[idx] >= 1);
                data.forwardAccess[maxNumForward[idx]] = {hierarchy.getTransitNodeIndexOfRank(rv), 0};
                data.forwardPos[idx] = 1;
                data.backwardAccess[maxNumBackward[idx]] = {hierarchy.getTransitNodeIndexOfRank(rv), 0};
                data.backwardPos[idx] = 1;
            } else {
                // Compute access nodes by unifying access nodes of upward neighbors
                AccessNodeUnifier unifier(hierarchy.numTransitNodes());
                computeAccessNodesForVertex(rv, maxNumForward, minCH.upwardGraph(), rankToIdx, data.forwardPos,
                                            data.forwardAccess, unifier);
                computeAccessNodesForVertex(rv, maxNumBackward, minCH.downwardGraph(), rankToIdx, data.backwardPos,
                                            data.backwardAccess, unifier);

                // Prune access nodes based on domination between each other
                pruneAccessNodesForVertex<true>(idx, maxNumForward, data.forwardPos, data.forwardAccess, data);
                pruneAccessNodesForVertex<false>(idx, maxNumBackward, data.backwardPos, data.backwardAccess, data);
            }
        });


        // Compute prefix sums over actual counts to get precise offsets
        forwardSum = 0;
        backwardSum = 0;
        for (int32_t i = 0; i < numVertices; ++i) {
            const int32_t fSize = data.forwardPos[i];
            const int32_t bSize = data.backwardPos[i];
            data.forwardPos[i] = forwardSum;
            data.backwardPos[i] = backwardSum;
            forwardSum += fSize;
            backwardSum += bSize;
        }
        data.forwardPos[numVertices] = forwardSum;
        data.backwardPos[numVertices] = backwardSum;

        // Move ranges to fill gaps
        for (int32_t i = 0; i < numVertices; ++i) {
            const int forwardCount = data.forwardPos[i + 1] - data.forwardPos[i];
            std::copy(data.forwardAccess.begin() + maxNumForward[i],
                      data.forwardAccess.begin() + maxNumForward[i] + forwardCount,
                      data.forwardAccess.begin() + data.forwardPos[i]);

            const int backwardCount = data.backwardPos[i + 1] - data.backwardPos[i];
            std::copy(data.backwardAccess.begin() + maxNumBackward[i],
                      data.backwardAccess.begin() + maxNumBackward[i] + backwardCount,
                      data.backwardAccess.begin() + data.backwardPos[i]);

            KASSERT(std::is_sorted(data.forwardAccess.begin() + data.forwardPos[i],
                                   data.forwardAccess.begin() + data.forwardPos[i + 1],
                                   [&](const CTNRData::AccessNode &a, const CTNRData::AccessNode &b) {
                                       return a.nodeIndex < b.nodeIndex;
                                   }));
            KASSERT(std::is_sorted(data.backwardAccess.begin() + data.backwardPos[i],
                                   data.backwardAccess.begin() + data.backwardPos[i + 1],
                                   [&](const CTNRData::AccessNode &a, const CTNRData::AccessNode &b) {
                                       return a.nodeIndex < b.nodeIndex;
                                   }));
        }
        data.forwardAccess.resize(forwardSum);
        data.backwardAccess.resize(backwardSum);

        // TODO: remove debug
        std::cout << "CTNR: Average forward access nodes per vertex: "
                  << static_cast<double>(data.forwardAccess.size()) / numVertices << std::endl;
        std::cout << "CTNR: Average backward access nodes per vertex: "
                  << static_cast<double>(data.backwardAccess.size()) / numVertices << std::endl;

        std::cout << "CTNR: Average forward upper bound in search space: "
                  << static_cast<double>(maxNumForward.back()) / numVertices << std::endl;
        std::cout << "CTNR: Average backward upper bound in search space: "
                  << static_cast<double>(maxNumBackward.back()) / numVertices << std::endl;
    }

    template<typename RankToIdx>
    void computeAccessNodesForVertex(const int rv,
                                     const std::vector<int32_t> &offset,
                                     const CH::SearchGraph &graph,
                                     const RankToIdx &rankToIdx,
                                     std::vector<int32_t> &dataPos,
                                     std::vector<CTNRData::AccessNode> &dataAccess,
                                     AccessNodeUnifier &unifier) const {
        const int idx = rankToIdx(rv);
        FORALL_INCIDENT_EDGES(graph, rv, e) {
            const int neighbor = graph.edgeHead(e);
            const int idxNeighbor = rankToIdx(neighbor);
            const int w = graph.traversalCost(e);
            KASSERT(w != INFTY);

            const int end = offset[idxNeighbor] + dataPos[idxNeighbor];
            for (auto i = offset[idxNeighbor]; i < end; ++i) {
                const auto &an = dataAccess[i];
                const int dist = an.distance + w;
                unifier.addAccessNode(an.nodeIndex, dist);
            }
        }
        const auto startOfRange = dataAccess.begin() + offset[idx];
        const int sizeOfRange = unifier.sizeOfUnion();
        KASSERT(sizeOfRange <= offset[idx + 1] - offset[idx]);
        dataPos[idx] = sizeOfRange;
        unifier.flushAccessNodes(startOfRange);
    }

    // Do not inline
    template<bool forward>
    [[gnu::noinline]] void pruneAccessNodesForVertex(const int idx,
                                                     const std::vector<int32_t> &offset,
                                                     std::vector<int32_t> &dataPos,
                                                     std::vector<CTNRData::AccessNode> &dataAccess,
                                                     CTNRData &data) const {
        int endOfNonDominated = offset[idx];
        int end = offset[idx] + dataPos[idx];
        for (int i = offset[idx]; i < end; ++i) {
            bool dominated = false;
            for (int j = offset[idx]; j < endOfNonDominated; ++j) {
                // TODO: These accesses to the distance table are all over the place, leading to about 20% of all cache
                //  misses during customization. We could reduce this by not using the actual distance table here but
                //  an optimized version instead. We can base this on the fact that these queries always have node j
                //  lower than node i, so we only care about upward/reverse downward distances between transit nodes
                //  here. Moreover, node i is always on the elim-tree branch of node j.
                //  Data structure: for each transit node, store upward distances to all higher transit nodes on its
                //  elim-tree branch, indexed by the depth on the branch from the root (root has 0). This depth is
                //  constant for every transit node so we can store it once.
                //  Let depth[t] be the depth of transit node t on its elim-tree branch.
                //  Let Dj[d] be the upward distance from transit node j to the transit node at depth d on j's elim-tree
                //  branch.
                //  Then, dTransit[j, i] = Dj[depth[i]].
                //  Can be done analogously but separately for reverse downward distances.
                //  Memory overhead should be okay since the elim-tree branch up to a transit node is short.
                //  This only improves cache behavior if we iterate over j (the lower node) in the outer loop and i
                //  (the higher node) in the inner loop. Though that should not be a problem (hopefully).
                const int dTransit = forward ?
                               data.getDistanceBetweenTransitNodes(dataAccess[j].nodeIndex, dataAccess[i].nodeIndex) :
                               data.getDistanceBetweenTransitNodes(dataAccess[i].nodeIndex, dataAccess[j].nodeIndex);
                KASSERT(dTransit != INFTY);
                if (dataAccess[j].distance + dTransit <= dataAccess[i].distance) {
                    dominated = true;
                    break;
                }
            }
            if (!dominated) {
                dataAccess[endOfNonDominated] = dataAccess[i];
                ++endOfNonDominated;
            }
        }
        dataPos[idx] = endOfNonDominated - offset[idx];
    }

//TODO: use PHAST to accelerate distance table computation
    void computeDistanceTable(CTNRData &data) {
        const int n = hierarchy.numTransitNodes();
        data.resetDistanceTable();
        using LabelSet = BasicLabelSet<0, ParentInfo::NO_PARENT_INFO>;
#pragma omp parallel
        {
            EliminationTreeQuery<LabelSet> chq(minCH, cch.getEliminationTree());
#pragma omp for
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    if (i == j) {
                        data.setDistanceBetweenTransitNodes(i, j, 0);
                        continue;
                    }
                    chq.run(hierarchy.getRankOfTransitNodeIndex(i), hierarchy.getRankOfTransitNodeIndex(j));
                    data.setDistanceBetweenTransitNodes(i, j, chq.getDistance());
                }
            }
        }
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
    void countTransitNodesInSearchSpace(std::vector<int32_t> &outCounts, const GraphT &upGraph,
                                        const RankToIdxT &rankToIdx) const {
        std::vector<int> firstChild;
        std::vector<int> children;
        convertInTreeToOutTree(cch.getEliminationTree(), firstChild, children);

//        const auto& upGraph = cch.getUpwardGraph();
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

        dfsOnTree(firstChild, children, recurse, backtrack);
    }


    const TransitNodeHierarchy &hierarchy;
    const CCH &cch;
    CCHMetric cchMetric;
    CH minCH;

    // Temporary data used during access node computation
//    std::vector<IndexRange> forwardRange;
//    std::vector<IndexRange> backwardRange;
//    std::vector<CTNRData::AccessNode> forwardAccessTemp;
//    std::vector<CTNRData::AccessNode> backwardAccessTemp;
};