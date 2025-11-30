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
                                int32_t const * const upWeights,
                                int32_t const * const downWeights)
            : hierarchy(hierarchy),
              cchGraph(cchGraph),
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
//        const auto &downGraph = minCH.downwardGraph();
        using PQEntry = std::pair<int32_t, int32_t>;
        std::priority_queue<PQEntry, std::vector<PQEntry>, std::greater<>> pq;
        const int sourceVertexId = transitIdToVertexId[sourceTransitIndex];
        data.setDistanceBetweenTransitNodes(sourceTransitIndex, sourceTransitIndex, 0);
        pq.emplace(0, sourceVertexId);
        while (!pq.empty()) {
            const auto [d, u] = pq.top();
            pq.pop();
            if (d != data.getDistanceBetweenTransitNodes(sourceTransitIndex, vertexIdToTransitId[u]))
                continue;
            FORALL_INCIDENT_EDGES(cchGraph, u, e) {
                const int v = cchGraph.edgeHead(e);
                const int neighborTransitIndex = vertexIdToTransitId[v];
                const int32_t w = upWeights[e];
                if (w == INFTY)
                    continue;
                const int32_t nd = d + w;
                if (nd < data.getDistanceBetweenTransitNodes(sourceTransitIndex, neighborTransitIndex)) {
                    data.setDistanceBetweenTransitNodes(sourceTransitIndex, neighborTransitIndex, nd);
                    pq.emplace(nd, v);
                }
            }
        }
        for (int i = hierarchy.numTransitNodes() - 1; i >= 0; --i) {
            const int u = transitIdToVertexId[i];
            FORALL_INCIDENT_EDGES(cchGraph, u, e) {
                const int neighborTransitIndex = vertexIdToTransitId[cchGraph.edgeHead(e)];
                const int32_t w = downWeights[e];
                const int32_t vd = data.getDistanceBetweenTransitNodes(sourceTransitIndex, neighborTransitIndex);
                if (w == INFTY || vd == INFTY)
                    continue;
                const int32_t nd = vd + w;
                if (nd < data.getDistanceBetweenTransitNodes(sourceTransitIndex, i))
                    data.setDistanceBetweenTransitNodes(sourceTransitIndex, i, nd);
            }
        }
    }

    const TransitNodeHierarchy &hierarchy;
    const CCH::UpGraph &cchGraph;
    int32_t const * const upWeights;
    int32_t const * const downWeights;
};


