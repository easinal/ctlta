#pragma once

#include "Algorithms/CTL/BalancedTopologyCentricTreeHierarchy.h"
#include "Algorithms/CCH/CCH.h"
#include "Algorithms/CCH/CCHMetric.h"
#include "Algorithms/CH/CH.h"
#include "DataStructures/Partitioning/SeparatorDecomposition.h"
#include "Tools/Constants.h"
#include "Algorithms/CTNR/CTNRData.h"
#include "TransitNodeHierarchy.h"
#include "TransitDistanceTableBuilder.h"
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
#include <boost/dynamic_bitset.hpp>

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
        int64_t dummy1, dummy2, dummy3, dummy4, dummy5, dummy6;
        customizeWithMeasurements(data, dummy1, dummy2, dummy3, dummy4, dummy5, dummy6);
    }

    // Sets measurement parameters to times for each step in microseconds.
    void customizeWithMeasurements(CTNRData &data,
                                   int64_t &cchBasicCustomizationTime,
                                   int64_t &cchPerfectCustomizationTime,
                                   int64_t &cchConstructChTime,
                                   int64_t &accessNodeComputationTime,
                                   int64_t &distanceTableComputationTime,
                                   int64_t &buildLocalMinCHTime) {
        cchPerfectCustomizationTime = 0;
        cchConstructChTime = 0;
        buildLocalMinCHTime = 0;
        Timer timer;
        cchMetric.customize();
        cchBasicCustomizationTime = timer.elapsed<std::chrono::microseconds>();
//        localMinCH = customizeCCHAndBuildLocalMinCH(cchBasicCustomizationTime, cchPerfectCustomizationTime, buildLocalMinCHTime);
        timer.restart();
        computeDistanceTable(data);
        distanceTableComputationTime = timer.elapsed<std::chrono::microseconds>();
        timer.restart();
        computeAccessNodes(data);
        accessNodeComputationTime = timer.elapsed<std::chrono::microseconds>();
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

        std::cout << "Finished CTNR customization in " << (cchBasicCustomizationTime + cchPerfectCustomizationTime +
                                                           cchConstructChTime + accessNodeComputationTime +
                                                           distanceTableComputationTime + buildLocalMinCHTime)
                  << " microseconds." << std::endl;
        const auto [avgNumNonInftyForward, avgNumNonInftyBackward] = computeAverageNumberOfNonInftyAccessNodes(data);
        std::cout << "CTNR: Average number of non-infinity forward/backward distances per vertex: "
                  << avgNumNonInftyForward << "/" << avgNumNonInftyBackward << std::endl;

    }

//    const CH &getLocalMinCH() const { return localMinCH; }

    const std::vector<int32_t> &getLocalEliminationTree() const { return localEliminationTree; }

    const CCHMetric &getCCHMetric() const { return cchMetric; }

    // Memory usage calculation including node levels
    uint64_t sizeInBytes() const {
        uint64_t size = sizeof(CTNRMetric);
        size += cchMetric.sizeInBytes();
//        size += localMinCH.sizeInBytes();

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
        for (const auto &accessNodeEdge: accessNodeEdges) {
            const auto &e = accessNodeEdge.edge;
            auto &f = data.forwardDistances[accessNodeEdge.position];
            auto &b = data.backwardDistances[accessNodeEdge.position];
            f = std::min(f, cchUpWeights[e]);
            b = std::min(b, cchDownWeights[e]);
        }


        // TODO: Debug for USA network
        const auto &cchGraph = cch.getUpwardGraph();
#pragma omp parallel
#pragma omp single nowait
        cch.forEachVertexTopDown([&](int32_t rv) {
            // Compute access node distances by using access node distances of upward neighbors
            computeAccessNodeDistancesForVertex(rv, cchGraph, cchUpWeights, rankToIdx, data.pos,
                                                data.accessNodes, data.forwardDistances);
            computeAccessNodeDistancesForVertex(rv, cchGraph, cchDownWeights, rankToIdx, data.pos,
                                                data.accessNodes, data.backwardDistances);

            // Prune access nodes based on domination between each other
//                pruneAccessNodesForVertex<true>(idx, data.pos, data.accessNodes, data.forwardDistances, data);
//                pruneAccessNodesForVertex<false>(idx, data.pos, data.accessNodes, data.backwardDistances, data);
//            }
        });
    }

    template<typename GraphT, typename RankToIdx>
    void computeAccessNodeDistancesForVertex(const int rv,
                                             const GraphT &graph,
                                             int const *const weights,
                                             const RankToIdx &rankToIdx,
                                             const std::vector<int32_t> &dataPos,
                                             const std::vector<ctnr::TransitNodeId> &dataNodes,
                                             CTNRData::DistanceVector<int32_t> &dataDistances) const {
        const int idx = rankToIdx(rv);
        const int startThis = dataPos[idx];
        FORALL_INCIDENT_EDGES(graph, rv, e) {
            const int neighbor = graph.edgeHead(e);
            const int idxNeighbor = rankToIdx(neighbor);
            const int w = weights[e];
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

    void computeDistanceTable(CTNRData &data) {
        // TODO: implement minCH for only transit node subgraph to use in distance table computation instead of CCH graph?
        TransitDistanceTableBuilder builder(hierarchy, cch.getUpwardGraph(), cchMetric.upwardWeights(),
                                            cchMetric.downwardWeights());
        builder.buildDistanceTable(data);
    }


    // An active vertex during a DFS, i.e., a vertex that has been reached but not finished.
    struct ActiveVertex {
        // Constructs an active vertex.
        ActiveVertex(const int id, const int nextUnexploredEdge)
                : id(id), nextUnexploredEdge(nextUnexploredEdge) {}

        int id;                 // The ID of the active vertex.
        int nextUnexploredEdge; // The next unexplored incident edge.
    };

    // Customize CCH including perfect customization and construct CH that only has edges required for the given metric.
    // CH is restricted to vertices below transit nodes, which is enough for local queries.
    CH customizeCCHAndBuildLocalMinCH(int64_t &cchBasicCustomizationTime, int64_t &cchPerfectCustomizationTime,
                                      int64_t &buildLocalMinCHTime) {

        Timer timer;
        cchMetric.customize();
        cchBasicCustomizationTime = timer.elapsed<std::chrono::microseconds>();

        timer.restart();
        const auto &cchGraph = cch.getUpwardGraph();
        std::vector<int8_t> keepUpEdge;
        std::vector<int8_t> keepDownEdge;

#pragma omp parallel sections
        {
#pragma omp section
            keepUpEdge.resize(cchGraph.numEdges() + 1, true);
#pragma omp section
            keepDownEdge.resize(cchGraph.numEdges() + 1, true);
        }

        keepUpEdge.back() = false;
        keepDownEdge.back() = false;
        // Run perfect customization marking all edges for removal that are not needed for this metric
        cchMetric.runPerfectCustomization(
                [&](const int e) { keepUpEdge[e] = false; },
                [&](const int e) { keepDownEdge[e] = false; });

        // Also mark all edges leading to transit nodes for removal
        FORALL_EDGES(cchGraph, e) {
            const auto head = cchGraph.edgeHead(e);
            if (hierarchy.isTransitNode(head)) {
                keepUpEdge[e] = false;
                keepDownEdge[e] = false;
            }
        }

        cchPerfectCustomizationTime = timer.elapsed<std::chrono::microseconds>();

        timer.restart();
        CH ch = cchMetric.buildCHKeepingGivenEdges(keepUpEdge, keepDownEdge);
        buildLocalMinCHTime = timer.elapsed<std::chrono::microseconds>();

        return ch;
    }

    std::pair<double, double> computeAverageNumberOfNonInftyAccessNodes(const CTNRData &data) const {
        int64_t fSum = 0;
        int64_t bSum = 0;
        int64_t count = 0;
        const auto &cchGraph = cch.getUpwardGraph();
        FORALL_VERTICES(cchGraph, rv) {
            const int idx = data.rankToIdx(rv);
            int localFCount = 0;
            int localBCount = 0;
            for (auto i = data.pos[idx]; i < data.pos[idx + 1]; ++i) {
                if (data.forwardDistances[i] != CTNR_INFTY)
                    localFCount++;
                if (data.backwardDistances[i] != CTNR_INFTY)
                    localBCount++;
            }
            fSum += localFCount;
            bSum += localBCount;
            count++;
        }
        return {static_cast<double>(fSum) / count, static_cast<double>(bSum) / count};
    }

    const TransitNodeHierarchy &hierarchy;
    const CCH &cch;
    const std::vector<AccessNodeEdge> &accessNodeEdges;
    CCHMetric cchMetric;

    // Minimum CH restricted to non-transit nodes for local queries.
//    CH localMinCH;

    // Elimination tree of the CCH restricted to non-transit nodes for local queries.
    std::vector<int32_t> localEliminationTree;
};