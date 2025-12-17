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
               const int32_t *const inputWeights,
               const ctnr::Level pruneLevelThreshold = 0)
            : hierarchy(hierarchy), cch(cch), accessNodeEdges(accessNodeEdges),
              pruneLevelThreshold(pruneLevelThreshold),
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
        int64_t dummy1, dummy2, dummy3;
        customizeWithMeasurements(data, dummy1, dummy2, dummy3);
    }

    // Sets measurement parameters to times for each step in microseconds.
    void customizeWithMeasurements(CTNRData &data,
                                   int64_t &cchBasicCustomizationTime,
                                   int64_t &distanceTableComputationTime,
                                   int64_t &accessNodeComputationTime) {
        Timer timer;
        cchMetric.customize();
        cchBasicCustomizationTime = timer.elapsed<std::chrono::microseconds>();
        timer.restart();
        computeDistanceTable(data);
        distanceTableComputationTime = timer.elapsed<std::chrono::microseconds>();
        timer.restart();
        computeAccessNodes(data);
        accessNodeComputationTime = timer.elapsed<std::chrono::microseconds>();

        std::cout << "Finished CTNR customization in " << (cchBasicCustomizationTime + distanceTableComputationTime +
                                                           accessNodeComputationTime)
                  << " microseconds." << std::endl;
    }

    const std::vector<int32_t> &getLocalEliminationTree() const { return localEliminationTree; }

    const CCHMetric &getCCHMetric() const { return cchMetric; }

    // Memory usage calculation including node levels
    uint64_t sizeInBytes() const {
        uint64_t size = sizeof(CTNRMetric);
        size += cchMetric.sizeInBytes();
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


        std::cout << "Pruning access nodes for vertices of level < " << (int) pruneLevelThreshold << std::endl;
        // TODO: Debug for USA network
        const auto &cchGraph = cch.getUpwardGraph();
#pragma omp parallel
#pragma omp single nowait
        cch.forEachVertexTopDownByLayer([&](int32_t rv) {
            // Compute access node distances by using access node distances of upward neighbors
            computeAccessNodeDistancesForVertex(rv, cchGraph, cchUpWeights, rankToIdx, data.pos,
                                                data.accessNodes, data.forwardDistances);
            computeAccessNodeDistancesForVertex(rv, cchGraph, cchDownWeights, rankToIdx, data.pos,
                                                data.accessNodes, data.backwardDistances);

            // Prune access nodes for this metric based on domination between access nodes.
            // Important to do this in the top-down sweep since already pruned upper neighbors lead to fewer
            // relevant access nodes at this vertex, which decreases the amount of work for pruning at this vertex.
            // Set pruneLevelThreshold to balance number of non-infinity access nodes for better query time and
            // customization time.
            if (hierarchy.getVertexLevel(rv) >= pruneLevelThreshold)
                return;
            const int idx = data.rankToIdx(rv);
            fullPruneAccessNodesForVertex<true>(idx, data.pos, data.accessNodes, data.forwardDistances, data);
            fullPruneAccessNodesForVertex<false>(idx, data.pos, data.accessNodes, data.backwardDistances, data);
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
        CTNRData::DistanceLabel distancesThis;
        CTNRData::DistanceLabel distancesNeighbor;
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
                // TODO: Optimize SIMD for multi-threading (use lower-level SIMD library to avoid load/store)
                KASSERT(numNeighbor % CTNRData::K == 0);
                const auto numBatches = numNeighbor / CTNRData::K;
                for (auto b = 0; b < numBatches; ++b) {
                    const auto offset = b * CTNRData::K;
                    distancesThis.load(&dataDistances[startThis + offset]);
                    distancesNeighbor.load(&dataDistances[startNeighbor + offset]);
                    distancesThis.min(distancesNeighbor + w);
                    distancesThis.store(&dataDistances[startThis + offset]);
                }
            }
        }
    }

    template<bool forward>
    DEBUG_NOINLINE
    void fullPruneAccessNodesForVertex(const int idx,
                                       const std::vector<int32_t> &dataPos,
                                       const std::vector<ctnr::TransitNodeId> &dataNodes,
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
        TransitDistanceTableBuilder builder(hierarchy, cch.getUpwardGraph(), cch.getEliminationTree(),
                                            cchMetric.upwardWeights(), cchMetric.downwardWeights());
        builder.buildDistanceTable(data);
    }

    const TransitNodeHierarchy &hierarchy;
    const CCH &cch;
    const std::vector<AccessNodeEdge> &accessNodeEdges;
    const ctnr::Level pruneLevelThreshold;
    CCHMetric cchMetric;

    // Minimum CH restricted to non-transit nodes for local queries.
//    CH localMinCH;

    // Elimination tree of the CCH restricted to non-transit nodes for local queries.
    std::vector<int32_t> localEliminationTree;
};