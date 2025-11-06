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

        void flushAccessNodes(std::vector<CTNRData::AccessNode> &outAccessNodes) {
            std::sort(nodes.begin(), nodes.end());
            for (const int & node : nodes) {
                int32_t &entry = distances[node];
                CTNRData::AccessNode an(node, entry);
                outAccessNodes.push_back(an);
                entry = INFTY; // reset for next use
            }
            nodes.clear();
        }

    private:
        std::vector<int> nodes; // list of transit node indices
        std::vector<int32_t> distances;
    };

public:

    // Constructor
    CTNRMetric(const TransitNodeHierarchy &hierarchy, const CCH &cch, const int32_t *const inputWeights)
            : hierarchy(hierarchy), cch(cch), cchMetric(cch, inputWeights),
              forwardRange(cch.getUpwardGraph().numVertices()), backwardRange(cch.getUpwardGraph().numVertices()) {}

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

        // Count the number of access nodes per vertex in data.forwardPos/data.backwardPos. Later, a prefix sum in
        // these vectors will give offsets into flat representation.
        data.forwardPos.clear();
        data.backwardPos.clear();
        data.forwardPos.resize(numVertices + 1, INVALID_INDEX);
        data.backwardPos.resize(numVertices + 1, INVALID_INDEX);

        // Collect ranges of access nodes into these temporary vectors first with arbitrary order of vertices.
        KASSERT(forwardRange.size() == numVertices && backwardRange.size() == numVertices);
        forwardAccessTemp.clear();
        backwardAccessTemp.clear();
        AccessNodeUnifier unifier(hierarchy.numTransitNodes());

        const auto &upGraph = minCH.upwardGraph();
        const auto &downGraph = minCH.downwardGraph();

        cch.forEachVertexTopDown([&](int32_t rv) {
            forwardRange[rv].start = static_cast<int32_t>(forwardAccessTemp.size());
            backwardRange[rv].start = static_cast<int32_t>(backwardAccessTemp.size());
            if (hierarchy.isTransitNode(rv)) {
                forwardAccessTemp.push_back({hierarchy.getTransitNodeIndexOfRank(rv), 0});
                backwardAccessTemp.push_back({hierarchy.getTransitNodeIndexOfRank(rv), 0});
                forwardRange[rv].end = static_cast<int32_t>(forwardAccessTemp.size());
                backwardRange[rv].end = static_cast<int32_t>(backwardAccessTemp.size());
            } else {
                FORALL_INCIDENT_EDGES(upGraph, rv, e) {
                    const int neighbor = upGraph.edgeHead(e);
                    const int wUp = upGraph.traversalCost(e);
                    KASSERT(wUp != INFTY);

                    ConstantVectorRange<CTNRData::AccessNode> neighborAccess(
                            forwardAccessTemp.begin() + forwardRange[neighbor].start,
                            forwardAccessTemp.begin() + forwardRange[neighbor].end);
                    for (const auto &an: neighborAccess) {
                        const int dist = wUp + an.distance;
                        unifier.addAccessNode(an.nodeIndex, dist);
                    }
                }
                unifier.flushAccessNodes(forwardAccessTemp);
                forwardRange[rv].end = static_cast<int32_t>(forwardAccessTemp.size());

                FORALL_INCIDENT_EDGES(downGraph, rv, e) {
                    const int neighbor = downGraph.edgeHead(e);
                    const int wDown = downGraph.traversalCost(e);
                    KASSERT(wDown != INFTY);

                    ConstantVectorRange<CTNRData::AccessNode> neighborAccess(
                            backwardAccessTemp.begin() + backwardRange[neighbor].start,
                            backwardAccessTemp.begin() + backwardRange[neighbor].end);
                    for (const auto &an: neighborAccess) {
                        const int dist = an.distance + wDown;
                        unifier.addAccessNode(an.nodeIndex, dist);
                    }
                }
                unifier.flushAccessNodes(backwardAccessTemp);
                backwardRange[rv].end = static_cast<int32_t>(backwardAccessTemp.size());

                // Prune access nodes based on domination between each other
                int endOfNonDominated = forwardRange[rv].start;
                for (int i = forwardRange[rv].start; i < forwardRange[rv].end; ++i) {
                    bool dominated = false;
                    for (int j = forwardRange[rv].start; j < endOfNonDominated; ++j) {
                        int dji = data.getDistanceBetweenTransitNodes(forwardAccessTemp[j].nodeIndex,
                                                                      forwardAccessTemp[i].nodeIndex);
                        if (dji != INFTY && forwardAccessTemp[j].distance + dji <= forwardAccessTemp[i].distance) {
                            dominated = true;
                            break;
                        }
                    }
                    if (!dominated) {
                        forwardAccessTemp[endOfNonDominated] = forwardAccessTemp[i];
                        ++endOfNonDominated;
                    }
                }
                forwardRange[rv].end = endOfNonDominated;
                forwardAccessTemp.erase(forwardAccessTemp.begin() + endOfNonDominated, forwardAccessTemp.end());

                endOfNonDominated = backwardRange[rv].start;
                for (int i = backwardRange[rv].start; i < backwardRange[rv].end; ++i) {
                    bool dominated = false;
                    for (int j = backwardRange[rv].start; j < endOfNonDominated; ++j) {
                        int dij = data.getDistanceBetweenTransitNodes(backwardAccessTemp[i].nodeIndex,
                                                                      backwardAccessTemp[j].nodeIndex);
                        if (dij != INFTY && dij + backwardAccessTemp[j].distance <= backwardAccessTemp[i].distance) {
                            dominated = true;
                            break;
                        }
                    }
                    if (!dominated) {
                        backwardAccessTemp[endOfNonDominated] = backwardAccessTemp[i];
                        ++endOfNonDominated;
                    }
                }
                backwardRange[rv].end = endOfNonDominated;
                backwardAccessTemp.erase(backwardAccessTemp.begin() + endOfNonDominated, backwardAccessTemp.end());
            }

            // Enter counts into data.forwardPos/data.backwardPos
            data.forwardPos[rv] = forwardRange[rv].end - forwardRange[rv].start;
            data.backwardPos[rv] = backwardRange[rv].end - backwardRange[rv].start;
        });


        // Compute prefix sums over counts to get offsets
        int32_t forwardSum = 0;
        int32_t backwardSum = 0;
        for (int32_t v = 0; v < numVertices; ++v) {
            const int32_t fSize = data.forwardPos[v];
            const int32_t bSize = data.backwardPos[v];
            data.forwardPos[v] = forwardSum;
            data.backwardPos[v] = backwardSum;
            forwardSum += fSize;
            backwardSum += bSize;
        }
        data.forwardPos[numVertices] = forwardSum;
        data.backwardPos[numVertices] = backwardSum;

        // Move ranges from unordered temporary vector to flat representation
        data.forwardAccess.resize(forwardAccessTemp.size());
        data.backwardAccess.resize(backwardAccessTemp.size());
        for (int32_t v = 0; v < numVertices; ++v) {
            KASSERT(data.forwardPos[v + 1] - data.forwardPos[v] == forwardRange[v].end - forwardRange[v].start);
            KASSERT(data.backwardPos[v + 1] - data.backwardPos[v] == backwardRange[v].end - backwardRange[v].start);
            std::copy(forwardAccessTemp.begin() + forwardRange[v].start,
                      forwardAccessTemp.begin() + forwardRange[v].end,
                      data.forwardAccess.begin() + data.forwardPos[v]);
            std::copy(backwardAccessTemp.begin() + backwardRange[v].start,
                      backwardAccessTemp.begin() + backwardRange[v].end,
                      data.backwardAccess.begin() + data.backwardPos[v]);
            KASSERT(std::is_sorted(data.forwardAccess.begin() + data.forwardPos[v],
                                   data.forwardAccess.begin() + data.forwardPos[v + 1],
                                   [&](const CTNRData::AccessNode &a, const CTNRData::AccessNode &b) {
                                       return a.nodeIndex < b.nodeIndex;
                                   }));
            KASSERT(std::is_sorted(data.backwardAccess.begin() + data.backwardPos[v],
                                   data.backwardAccess.begin() + data.backwardPos[v + 1],
                                   [&](const CTNRData::AccessNode &a, const CTNRData::AccessNode &b) {
                                       return a.nodeIndex < b.nodeIndex;
                                   }));
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

    const TransitNodeHierarchy &hierarchy;
    const CCH &cch;
    CCHMetric cchMetric;
    CH minCH;

    // Temporary data used during access node computation
    std::vector<IndexRange> forwardRange;
    std::vector<IndexRange> backwardRange;
    std::vector<CTNRData::AccessNode> forwardAccessTemp;
    std::vector<CTNRData::AccessNode> backwardAccessTemp;
};