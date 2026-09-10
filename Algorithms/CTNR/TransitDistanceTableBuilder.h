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
        if (numTransit == 0)
            return;
        constexpr int K = 16;
        const auto &transitIdToVertexId = hierarchy.getTransitNodes();
        const auto &vertexIdToTransitId = hierarchy.getTransitNodeIndexOfRankVector();

        const int numChunks = (numTransit + K - 1) / K;

        // Transit-subgraph CSR with heads already as transit indices and weights inlined,
        // so the down sweep never resolves a vertex id again.
        std::vector<int64_t> tFirst(numTransit + 1, 0);
        for (int i = 0; i < numTransit; ++i) {
            const int u = transitIdToVertexId[i];
            tFirst[i + 1] = tFirst[i] + (cchGraph.lastEdge(u) - cchGraph.firstEdge(u));
        }
        std::vector<int32_t> tHead(tFirst[numTransit]);
        std::vector<int32_t> tUpW(tFirst[numTransit]);
        std::vector<int32_t> tDownW(tFirst[numTransit]);
#pragma omp parallel for schedule(static)
        for (int i = 0; i < numTransit; ++i) {
            const int u = transitIdToVertexId[i];
            int64_t k = tFirst[i];
            FORALL_INCIDENT_EDGES(cchGraph, u, e) {
                tHead[k] = vertexIdToTransitId[cchGraph.edgeHead(e)];
                tUpW[k] = upWeights[e];
                tDownW[k] = downWeights[e];
                ++k;
            }
        }
#pragma omp parallel
        {
            std::vector<int32_t> buf(static_cast<size_t>(numTransit) * K);
            std::vector<int32_t *> rows(K);
#pragma omp for schedule(dynamic, 1)
            for (int chunk = 0; chunk < numChunks; ++chunk) {
                const int s0 = chunk * K;
                const int kn = std::min(K, numTransit - s0);
                std::fill(buf.begin(), buf.end(), CTNR_INFTY);

                // Up phase: one elimination-path sweep per source (scalar, cheap).
                for (int k = 0; k < kn; ++k) {
                    buf[static_cast<size_t>(s0 + k) * K + k] = 0;
                    for (int u = transitIdToVertexId[s0 + k]; u != INVALID_VERTEX;
                         u = eliminationTree[u]) {
                        const int i = vertexIdToTransitId[u];
                        const int32_t d = buf[static_cast<size_t>(i) * K + k];
                        if (d >= CTNR_INFTY)
                            continue;
                        for (int64_t e = tFirst[i]; e < tFirst[i + 1]; ++e) {
                            const int32_t nd = d + tUpW[e];
                            auto &slot = buf[static_cast<size_t>(tHead[e]) * K + k];
                            if (nd < slot)
                                slot = nd;
                        }
                    }
                }

                // Down phase: one shared sweep, K sources wide.
                for (int i = numTransit - 1; i >= 0; --i) {
                    int32_t *const __restrict dst = buf.data() + static_cast<size_t>(i) * K;
                    for (int64_t e = tFirst[i]; e < tFirst[i + 1]; ++e) {
                        const int32_t w = tDownW[e];
                        const int32_t *const __restrict srcRow =
                                buf.data() + static_cast<size_t>(tHead[e]) * K;
                        for (int k = 0; k < K; ++k) {
                            const int32_t nd = srcRow[k] + w;
                            if (nd < dst[k])
                                dst[k] = nd;
                        }
                    }
                }

                // Scatter the transposed batch into the table rows (K sequential streams).
                for (int k = 0; k < kn; ++k)
                    rows[k] = data.getDistanceTableRow(s0 + k);
                for (int j = 0; j < numTransit; ++j) {
                    const int32_t *const v = buf.data() + static_cast<size_t>(j) * K;
                    for (int k = 0; k < kn; ++k)
                        rows[k][j] = v[k];
                }
            }
        }
    }

private:
    const TransitNodeHierarchy &hierarchy;
    const CCH::UpGraph &cchGraph;
    const std::vector<int32_t> &eliminationTree;
    int32_t const * const upWeights;
    int32_t const * const downWeights;
};


