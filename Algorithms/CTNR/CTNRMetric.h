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

        // Compute elimination tree in out format for DFS traversals
        convertInTreeToOutTree(cch.getEliminationTree(), elimTreeFirstChild, elimTreeChildren);
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

        // Compute partial distance table for better accesses to upward/downward distances between transit nodes during
        // pruning
        computeUpAndDownTransitDistances(data);

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
        data.forwardPos.resize(numVertices + 1);
        data.backwardPos.resize(numVertices + 1);
        data.forwardAccess.resize(forwardSum);
        data.backwardAccess.resize(backwardSum);

        // Collect ranges of access nodes into these temporary vectors first with arbitrary order of vertices.

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
                pruneAccessNodesForVertex<true>(idx, maxNumForward, data.forwardPos, data.forwardAccess);
                pruneAccessNodesForVertex<false>(idx, maxNumBackward, data.backwardPos, data.backwardAccess);
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
                                                     std::vector<CTNRData::AccessNode> &dataAccess) const {

        int& num = dataPos[idx];
        std::vector<bool> isDominated(num, false);

        const int start = offset[idx];
        KASSERT(std::is_sorted(dataAccess.begin() + start, dataAccess.begin() + start + num,
                               [&](const CTNRData::AccessNode &a, const CTNRData::AccessNode &b) {
                                   const int depthA = transitDistPos[a.nodeIndex + 1] - transitDistPos[a.nodeIndex];
                                   const int depthB = transitDistPos[b.nodeIndex + 1] - transitDistPos[b.nodeIndex];
                                   return depthA > depthB;
                               }));

        int endOfNonDominated = start;
        int end = start + num;
        for (int j = start; j < end; ++j) {
            if (isDominated[j - start])
                continue;
            const int nodeJ = dataAccess[j].nodeIndex;
            const int distJ = dataAccess[j].distance;
            const int transitDistStartJ = transitDistPos[nodeJ];
            int prevDepthI = -1;
            for (int i = end - 1; i > j; --i) {
                if (isDominated[i - start])
                    continue;
                const int nodeI = dataAccess[i].nodeIndex;
                KASSERT(nodeJ < nodeI);
                const int depthI = transitDistPos[nodeI + 1] - transitDistPos[nodeI];
                KASSERT(depthI > prevDepthI);
                const int dTransit = forward? upTransitDist[transitDistStartJ + depthI] : downTransitDist[transitDistStartJ + depthI];
                KASSERT(dTransit != INFTY);
                if (distJ + dTransit <= dataAccess[i].distance) {
                    isDominated[i - start] = true;
                }
                prevDepthI = depthI;
                unused(prevDepthI);
            }
            dataAccess[endOfNonDominated] = dataAccess[j];
            ++endOfNonDominated;
        }
        num = endOfNonDominated - start;
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

        KASSERT(std::is_sorted(firstChild.begin(), firstChild.end()));
    }

    // Run a DFS for the given tree in out format.
    // Call callbacks when recursing or backtracking.
    template<typename RecurseCallBack,
            typename BacktrackCallBack>
    static void dfsOnTree(
            const std::vector<int> &firstChild,
            const std::vector<int> &children,
            const RecurseCallBack& recurse,
            const BacktrackCallBack& backtrack) {
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
            const bool pruneAtChild = recurse(v.id, child);
            ++v.nextUnexploredEdge; // When backtracking from child later, look at next sibling
            if (pruneAtChild)
                continue;
            activeVertices.emplace(child, firstChild[child]);
        }
    }

    template<typename GraphT, typename RankToIdxT>
    void countTransitNodesInSearchSpace(std::vector<int32_t> &outCounts, const GraphT &upGraph,
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
                return false; // do not prune
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
            return false; // do not prune
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

    void computeElimTreeDepthOfTransitNodes(std::vector<int32_t>& depth) {

        const int numVertices = cch.getUpwardGraph().numVertices();
        depth[hierarchy.getTransitNodeIndexOfRank(numVertices - 1)] = 0; // root

        int curDepth = 0;

        const auto recurse = [&](const int /*parent*/, const int child) {
            if (!hierarchy.isTransitNode(child))
                return true; // prune
            depth[hierarchy.getTransitNodeIndexOfRank(child)] = ++curDepth;
//            ++curDepth;
//            if (hierarchy.isTransitNode(child))
//                depth[hierarchy.getTransitNodeIndexOfRank(child)] = curDepth;
            return false;
        };

        const auto backtrack = [&](const int child, const int /*parent*/) {
            unused(child);
            KASSERT(hierarchy.isTransitNode(child));
            --curDepth;
        };

        dfsOnTree(elimTreeFirstChild, elimTreeChildren, recurse, backtrack);
    }

    void computeUpAndDownTransitDistances(const CTNRData& data) {
        const int n = hierarchy.numTransitNodes();
        transitDistPos.resize(n + 1, 0);

        // Depths of transit nodes in elimination tree
        computeElimTreeDepthOfTransitNodes(transitDistPos);

        // Prefix sum to get offsets
        int32_t sum = 0;
        for (int32_t i = 0; i < n; ++i) {
            const int32_t d = transitDistPos[i];
            transitDistPos[i] = sum;
            sum += d;
        }
        transitDistPos[n] = sum;
        upTransitDist.resize(sum, INFTY);
        downTransitDist.resize(sum, INFTY);

        // fill upTransitDist and downTransitDist with distances to upper transit nodes on elim-tree branch
        // using DFS
        std::vector<int> onBranch; // stack of transit node indices on current elim-tree branch
        const int numVertices = cch.getUpwardGraph().numVertices();
        if (hierarchy.isTransitNode(numVertices - 1)) {
            onBranch.push_back(hierarchy.getTransitNodeIndexOfRank(numVertices - 1));
            KASSERT(transitDistPos[onBranch.back() + 1] - transitDistPos[onBranch.back()] == 0);
        }
        const auto recurse = [&](const int /*parent*/, const int child) {
            if (!hierarchy.isTransitNode(child))
                return true; // prune
            const int idxChild = hierarchy.getTransitNodeIndexOfRank(child);
            KASSERT(transitDistPos[idxChild + 1] - transitDistPos[idxChild] == static_cast<int>(onBranch.size()));
            const int start = transitDistPos[idxChild];
            for (auto u = 0; u < onBranch.size(); ++u) {
                KASSERT(onBranch[u] > idxChild);
                upTransitDist[start + u] = data.getDistanceBetweenTransitNodes(idxChild, onBranch[u]);
                downTransitDist[start + u] = data.getDistanceBetweenTransitNodes(onBranch[u], idxChild);
            }
            onBranch.push_back(idxChild);
            return false; // do not prune
        };

        const auto backtrack = [&](const int child, const int /*parent*/) {
            unused(child);
            KASSERT(hierarchy.isTransitNode(child));
            onBranch.pop_back();
        };

        dfsOnTree(elimTreeFirstChild, elimTreeChildren, recurse, backtrack);
    }


    const TransitNodeHierarchy &hierarchy;
    const CCH &cch;
    std::vector<int> elimTreeFirstChild;
    std::vector<int> elimTreeChildren;

    CCHMetric cchMetric;
    CH minCH;

    std::vector<int32_t> transitDistPos;
    std::vector<int32_t> upTransitDist;
    std::vector<int32_t> downTransitDist;

};