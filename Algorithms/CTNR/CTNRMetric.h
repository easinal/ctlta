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
#include "DataStructures/Labels/SimdLabelSet.h"
#include <algorithm>
#include <iostream>

class CTNRMetric {

public:

    // Constructor
    CTNRMetric(const TransitNodeHierarchy &hierarchy, const CCH &cch,
               const std::vector<AccessNodeEdge> &accessNodeEdges,
               const int32_t *const inputWeights)
            : hierarchy(hierarchy), cch(cch), accessNodeEdges(accessNodeEdges),
              cchMetric(cch, inputWeights),
              localEliminationTree(cch.getEliminationTree()) {

        // Build local elimination tree
        for (int &v: localEliminationTree) {
            if (v != INVALID_VERTEX && hierarchy.isTransitNode(v))
                v = INVALID_VERTEX;
        }
    }

    // Customization phase
    void customize(CTNRData &data) {
        int64_t dummy1, dummy2, dummy3, dummy4;
        customizeWithMeasurements(data, dummy1, dummy2, dummy3, dummy4);
    }

    // Sets measurement parameters to times for each step in microseconds.
    void customizeWithMeasurements(CTNRData &data, int64_t &cchCustomizationTime, int64_t &accessNodeComputationTime,
                                   int64_t &distanceTableComputationTime, int64_t &buildLocalMinCHTime) {
        Timer timer;
        minCH = cchMetric.buildMinimumWeightedCH();
        cchCustomizationTime = timer.elapsed<std::chrono::microseconds>();
        timer.restart();
        computeDistanceTable(data);
        distanceTableComputationTime = timer.elapsed<std::chrono::microseconds>();
        timer.restart();
        computeAccessNodes(data);
        accessNodeComputationTime = timer.elapsed<std::chrono::microseconds>();
        timer.restart();
        localMinCH = buildLocalMinCH();
        buildLocalMinCHTime = timer.elapsed<std::chrono::microseconds>();

//        // Debug information
//        int64_t sumNonInftyAfterPruningForward = 0;
//        int64_t sumNonInftyAfterPruningBackward = 0;
//        for (int32_t rv = 0; rv < cch.getUpwardGraph().numVertices(); ++rv) {
//            const int idx = data.rankToIdx(rv);
//            for (auto i = data.pos[idx]; i < data.pos[idx + 1]; ++i) {
//                if (data.forwardDistances[i] != CTNR_INFTY)
//                    sumNonInftyAfterPruningForward++;
//                if (data.backwardDistances[i] != CTNR_INFTY)
//                    sumNonInftyAfterPruningBackward++;
//            }
//        }
//
//        std::cout << "CTNR: Average number of non-infinity forward distances per vertex: "
//                  << static_cast<double>(sumNonInftyAfterPruningForward) / cch.getUpwardGraph().numVertices()
//                  << std::endl;
//        std::cout << "CTNR: Average number of non-infinity backward distances per vertex: "
//                  << static_cast<double>(sumNonInftyAfterPruningBackward) / cch.getUpwardGraph().numVertices()
//                  << std::endl;

        std::cout << "Finished CTNR customization in " << (cchCustomizationTime + accessNodeComputationTime +
                                                           distanceTableComputationTime + buildLocalMinCHTime)
                  << " microseconds." << std::endl;
    }

    const CH &getMinCH() const { return minCH; }

    const CH &getLocalMinCH() const { return localMinCH; }

    const std::vector<int32_t> &getLocalEliminationTree() const { return localEliminationTree; }

    // Memory usage calculation including node levels
    uint64_t sizeInBytes() const {
        uint64_t size = sizeof(CTNRMetric);
        size += cchMetric.sizeInBytes();
        size += minCH.sizeInBytes();

        return size;
    }

private:

    void computeAccessNodes(CTNRData &data) {

        const auto rankToIdx = [&](const int r) {
            return data.rankToIdx(r);
        };
        std::fill(data.forwardDistances.begin(), data.forwardDistances.end(), CTNR_INFTY);
        std::fill(data.backwardDistances.begin(), data.backwardDistances.end(), CTNR_INFTY);

        // Relax edges leading directly from non-transit nodes to access nodes first, then ignore them during sweep.
        const auto cchUpWeights = cchMetric.upwardWeights();
        const auto cchDownWeights = cchMetric.downwardWeights();
#pragma omp parallel for schedule(static)
        for (const auto& accessNodeEdge : accessNodeEdges) {
            const auto& e = accessNodeEdge.edge;
            auto& f = data.forwardDistances[accessNodeEdge.position];
            auto & b = data.backwardDistances[accessNodeEdge.position];
            f = std::min(f, cchUpWeights[e]);
            b = std::min(b, cchDownWeights[e]);
        }


        // TODO: Debug for USA network
#pragma omp parallel
#pragma omp single nowait
        cch.forEachVertexTopDown([&](int32_t rv) {
                // Compute access node distances by using access node distances of upward neighbors
                computeAccessNodeDistancesForVertex(rv, minCH.upwardGraph(), rankToIdx, data.pos,
                                                    data.accessNodes, data.forwardDistances);
                computeAccessNodeDistancesForVertex(rv, minCH.downwardGraph(), rankToIdx, data.pos,
                                                    data.accessNodes, data.backwardDistances);

                // Prune access nodes based on domination between each other
//                pruneAccessNodesForVertex<true>(idx, data.pos, data.accessNodes, data.forwardDistances, data);
//                pruneAccessNodesForVertex<false>(idx, data.pos, data.accessNodes, data.backwardDistances, data);
//            }
        });
    }

    template<typename RankToIdx>
    void computeAccessNodeDistancesForVertex(const int rv,
                                             const CH::SearchGraph &graph,
                                             const RankToIdx &rankToIdx,
                                             const std::vector<int32_t> &dataPos,
                                             const std::vector<int32_t> &dataNodes,
                                             CTNRData::DistanceVector<int32_t> &dataDistances) const {
        const int idx = rankToIdx(rv);
        const int startThis = dataPos[idx];
        FORALL_INCIDENT_EDGES(graph, rv, e) {
            const int neighbor = graph.edgeHead(e);
            const int idxNeighbor = rankToIdx(neighbor);
            const int w = graph.traversalCost(e);
            KASSERT(w < INFTY);

            // If neighbor is not a transit node, propagate all its access nodes. The list of access nodes of the
            // neighbor is a prefix of the list of access nodes of this vertex. Thus, we can just update by index.
            const auto startNeighbor = dataPos[idxNeighbor];
            const auto numNeighbor = dataPos[idxNeighbor + 1] - startNeighbor;
            KASSERT(numNeighbor <= dataPos[idx + 1] - startThis);
            if constexpr (!CTNRData::USE_SIMD) {
                for (auto i = 0; i < numNeighbor; ++i) {
                    KASSERT(dataNodes[startThis + i] == dataNodes[startNeighbor + i]);
                    const int newDist = dataDistances[startNeighbor + i] + w;
                    auto &entry = dataDistances[startThis + i];
                    if (newDist < entry)
                        entry = newDist;
                }
            } else {
                KASSERT(numNeighbor % CTNRData::K == 0);
                const auto numBatches = numNeighbor / CTNRData::K;
                CTNRData::DistanceLabel distancesThis;
                CTNRData::DistanceLabel distancesNeighbor;
                for (auto b = 0; b < numBatches; ++b) {
                    const auto offset = b * CTNRData::K;
                    distancesThis.load(&dataDistances[startThis + offset]);
                    distancesNeighbor.load(&dataDistances[startNeighbor + offset]);
                    const CTNRData::DistanceLabel wLabel = w;
                    const CTNRData::DistanceLabel newDistances = distancesNeighbor + wLabel;
                    distancesThis.min(newDistances);
                    distancesThis.store(&dataDistances[startThis + offset]);
                }
            }
        }
    }

    template<bool forward>
    DEBUG_NOINLINE
    void pruneAccessNodesForVertex(const int idx,
                                   const std::vector<int32_t> &dataPos,
                                   const std::vector<int32_t> &dataNodes,
                                   CTNRData::DistanceVector<int32_t> &dataDistances,
                                   const CTNRData &data) const {
        // Set distances for all dominated access nodes to CTNR_INFTY
        const auto start = dataPos[idx];
        const auto end = dataPos[idx + 1];
        for (auto i = start; i < end; ++i) {
            if (dataDistances[i] == CTNR_INFTY)
                continue;
            for (auto j = start; j < end; ++j) {
                if (i == j)
                    continue;
                if (dataDistances[j] == CTNR_INFTY)
                    continue;
                KASSERT(dataNodes[j] != dataNodes[i]);
                const int dTransit = forward ?
                                     data.getDistanceBetweenTransitNodes(dataNodes[j], dataNodes[i]) :
                                     data.getDistanceBetweenTransitNodes(dataNodes[i], dataNodes[j]);
                KASSERT(dTransit < CTNR_INFTY);
                if (dataDistances[j] + dTransit <= dataDistances[i]) {
                    dataDistances[i] = CTNR_INFTY;
                    break;
                }
            }
        }
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

    // Construct subgraph CH restricted to vertices below transit nodes, which is enough for local queries.
    CH buildLocalMinCH() const {
        CH::SearchGraph subUpGraph = minCH.upwardGraph();
        CH::SearchGraph subDownGraph = minCH.downwardGraph();
        const auto eraseEdgeToTransitNode = [&](const int, const int v) {
            return hierarchy.isTransitNode(v);
        };
        subUpGraph.eraseEdges(eraseEdgeToTransitNode);
        subDownGraph.eraseEdges(eraseEdgeToTransitNode);
        KASSERT(subUpGraph.isDefrag() && subUpGraph.validate());
        KASSERT(subDownGraph.isDefrag() && subDownGraph.validate());
        return {std::move(subUpGraph), std::move(subDownGraph), minCH.getOrderPermutation(),
                minCH.getRanksPermutation()};
    }


    const TransitNodeHierarchy &hierarchy;
    const CCH &cch;
    const std::vector<AccessNodeEdge> &accessNodeEdges;
    CCHMetric cchMetric;
    CH minCH;

    // Minimum CH restricted to non-transit nodes for local queries.
    CH localMinCH;

    // Elimination tree of the CCH restricted to non-transit nodes for local queries.
    std::vector<int32_t> localEliminationTree;
};