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
                                CH transitCH,
                                std::vector<int32_t> vertexIdToTransitId,
                                std::vector<int32_t> transitIdToVertexId)
            : hierarchy(hierarchy),
              transitCH(std::move(transitCH)),
              vertexIdToTransitId(std::move(vertexIdToTransitId)),
              transitIdToVertexId(std::move(transitIdToVertexId)) {}

    void buildDistanceTable(CTNRData &data) const {
        const int numTransit = hierarchy.numTransitNodes();
        data.resetDistanceTable();
        if (numTransit == 0)
            return;

        // #pragma omp parallel for schedule(dynamic)
            for (int i = 0; i < numTransit; ++i) {
                runTransitSSSP(i, data);
            }
    }

private:
    void runTransitSSSP(const int sourceTransitIndex,
                        CTNRData &data) const {
        const auto &upGraph = transitCH.upwardGraph();
        const auto &downGraph = transitCH.downwardGraph();
        using PQEntry = std::pair<int32_t, int32_t>;
        std::priority_queue<PQEntry, std::vector<PQEntry>, std::greater<PQEntry>> pq;
        const int sourceVertexId = transitIdToVertexId[sourceTransitIndex];
        data.setDistanceBetweenTransitNodes(sourceTransitIndex, sourceTransitIndex, 0);
        pq.emplace(0, sourceVertexId);

        while (!pq.empty()) {
            const auto [d, u] = pq.top();
            pq.pop();
            if (d != data.getDistanceBetweenTransitNodes(sourceTransitIndex, vertexIdToTransitId[u]))
                continue;
            FORALL_INCIDENT_EDGES(upGraph, u, e) {
                const int v = upGraph.edgeHead(e);
                const int neighborTransitIndex = vertexIdToTransitId[v];
                const int32_t w = upGraph.traversalCost(e);
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
            const int32_t du = data.getDistanceBetweenTransitNodes(sourceTransitIndex, vertexIdToTransitId[u]);
            if (du == INFTY)
                continue;
            FORALL_INCIDENT_EDGES(downGraph, u, e) {
                
                const int v = downGraph.edgeHead(e);
                const int neighborTransitIndex = vertexIdToTransitId[v];
                const int32_t w = downGraph.traversalCost(e);

                if (w == INFTY)
                    continue;
                const int32_t nd = du + w;
                if (nd < data.getDistanceBetweenTransitNodes(sourceTransitIndex, neighborTransitIndex))
                    data.setDistanceBetweenTransitNodes(sourceTransitIndex, neighborTransitIndex, nd);
            }
        }
    }

    const TransitNodeHierarchy &hierarchy;
    CH transitCH;
    std::vector<int32_t> vertexIdToTransitId;
    std::vector<int32_t> transitIdToVertexId;
};


