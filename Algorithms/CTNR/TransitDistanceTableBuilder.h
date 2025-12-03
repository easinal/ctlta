#pragma once

#include "Algorithms/CH/CH.h"
#include "Algorithms/CTNR/CTNRData.h"
#include "TransitNodeHierarchy.h"
#include "Tools/Constants.h"

#include <algorithm>
#include <queue>
#include <vector>

class TransitDistanceTableBuilder {
public:
    TransitDistanceTableBuilder(const TransitNodeHierarchy &hierarchy,
                                const CCH::UpGraph &cchGraph,
                                const std::vector<int32_t>& eliminationTree,
                                int32_t const * const upWeights,
                                int32_t const * const downWeights)
            : hierarchy(hierarchy),
              cchGraph(cchGraph),
              eliminationTree(eliminationTree),
              upWeights(upWeights),
              downWeights(downWeights) {}

    void buildDistanceTable(CTNRData &data) const {
        const int numTransit = hierarchy.numTransitNodes();
        data.resetDistanceTable();
        if (numTransit == 0)
            return;

        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < numTransit; ++i) {
            runTransitSSSP(i, data);
        }
    }

private:
    void runTransitSSSP(const int sourceTransitIndex,
                        CTNRData &data ) const {
        const auto &transitIdToVertexId = hierarchy.getTransitNodes();
        const auto &vertexIdToTransitId = hierarchy.getTransitNodeIndexOfRankVector();
        const int sourceVertexId = transitIdToVertexId[sourceTransitIndex];
        int* distanceTableRow = data.getDistanceTableRow(sourceTransitIndex);
        distanceTableRow[sourceTransitIndex] = 0;
        int u = sourceVertexId;
        while (u != INVALID_VERTEX) {
            const int32_t d = distanceTableRow[vertexIdToTransitId[u]];
            FORALL_INCIDENT_EDGES(cchGraph, u, e) {
                const int v = cchGraph.edgeHead(e);
                const int neighborTransitIndex = vertexIdToTransitId[v];
                const int32_t w = upWeights[e];
                const int32_t nd = d + w;
                if (nd < distanceTableRow[neighborTransitIndex]) {
                    distanceTableRow[neighborTransitIndex] = nd;
                }
            }
            u = eliminationTree[u];
        }

        for (int i = hierarchy.numTransitNodes() - 1; i >= 0; --i) {
            const int u = transitIdToVertexId[i];
            FORALL_INCIDENT_EDGES(cchGraph, u, e) {
                const int neighborTransitIndex = vertexIdToTransitId[cchGraph.edgeHead(e)];
                const int32_t w = downWeights[e];
                const int32_t vd = distanceTableRow[neighborTransitIndex];
                const int32_t nd = vd + w;
                if (nd < distanceTableRow[i]) {
                    distanceTableRow[i] = nd;
                }
            }
        }
    }

    const TransitNodeHierarchy &hierarchy;
    const CCH::UpGraph &cchGraph;
    const std::vector<int32_t> &eliminationTree;
    int32_t const * const upWeights;
    int32_t const * const downWeights;
};


