#pragma once

#include <algorithm>
#include <cstdint>
#include <memory>
#include <numeric>
#include <ostream>
#include <queue>
#include <stdexcept>
#include <vector>

#include <kassert/kassert.hpp>

#include "Algorithms/FHL/core/FHLConstants.h"
#include "Algorithms/FHL/core/HubHierarchy.h"
#include "Tools/Constants.h"
#include "Tools/Timer.h"
#include "Algorithms/FHL/core/FHLRegions.h"
#include "Algorithms/FHL/core/FHLStore.h"

namespace fhl {

// SCHEME 2 -- FHL + in-region hubs. The leaf region is cut a second time at S' vertices;
// the vertices left above that inner cut are the region's own hub set and every vertex
// stores its region-restricted distances to them. Same-region pairs in different
// sub-regions then fold two aligned rows instead of running the in-region CCH search.
struct InRegionHubs {
    // In-region rows for a list of regions (nullptr = all): d_R(v, col) for the region's own
    // top-K separator vertices, by the same descending-rank DP as the seeds — an in-region
    // up-neighbour is an elimination ancestor with a higher rank, so its row is final, and a
    // column's own row starts at 0.
    template<typename CchT>
    static void buildLocalRows(const Regions &R,
                               const HubHierarchy &hierarchy [[maybe_unused]],
                               const CchT &cch, const int32_t *const upW,
                               const int32_t *const downW, const bool symmetric,
                               const int32_t *const regionList, const int32_t numInList,
                               Values &out) {
        if (R.locLevels <= 0)
            return;
        const auto &graph = cch.getUpwardGraph();
        out.locLevels = R.locLevels;
        if (regionList == nullptr) {
            out.locPos.assign(R.numVertices + 1, 0);
            for (int32_t v = 0; v < R.numVertices; ++v) {
                const int32_t r = R.leafRegionOfVertex[v];
                out.locPos[v + 1] = out.locPos[v] +
                        (r < 0 ? 0 : R.locColFirst[r + 1] - R.locColFirst[r]);
            }
            out.locFwd.allocateForOverwrite(out.locPos[R.numVertices]);
            if (symmetric)
                out.locBwd.clear();
            else
                out.locBwd.allocateForOverwrite(out.locPos[R.numVertices]);
        }
        const auto doRegion = [&](const int32_t r) {
            const int64_t base = R.locColFirst[r];
            const int32_t n = static_cast<int32_t>(R.locColFirst[r + 1] - base);
            if (n == 0)
                return;
            const int32_t *const cols = R.locColIds.data() + base; // ascending ranks
            const auto colOf = [&](const int32_t x) -> int32_t {
                if (R.innerRegionOfVertex[x] >= 0)
                    return -1; // inside a sub-region, not a column
                const auto it = std::lower_bound(cols, cols + n, x);
                return (it != cols + n && *it == x) ? static_cast<int32_t>(it - cols) : -1;
            };
            for (int64_t k = R.leafVertFirst[r + 1] - 1; k >= R.leafVertFirst[r]; --k) {
                const int32_t v = R.leafVertIds[k];
                int32_t *const fw = out.locFwd.data() + out.locPos[v];
                int32_t *const bw = symmetric ? nullptr : out.locBwd.data() + out.locPos[v];
                std::fill(fw, fw + n, FHL_INFTY);
                if (bw)
                    std::fill(bw, bw + n, FHL_INFTY);
                const int32_t own = colOf(v);
                if (own >= 0) {
                    fw[own] = 0;
                    if (bw)
                        bw[own] = 0;
                }
                for (int e = graph.firstEdge(v); e < graph.lastEdge(v); ++e) {
                    const int32_t h = graph.edgeHead(e);
                    if (R.leafRegionOfVertex[h] != r) // leaves the region: not our business
                        continue;
                    const int32_t *const rowF = out.locFwd.data() + out.locPos[h];
                    const int32_t wu = upW[e];
                    for (int32_t j = 0; j < n; ++j) {
                        const int32_t dd = wu + rowF[j];
                        if (rowF[j] < FHL_INFTY && dd < fw[j])
                            fw[j] = dd;
                    }
                    if (bw) {
                        const int32_t *const rowB = out.locBwd.data() + out.locPos[h];
                        const int32_t wd = downW[e];
                        for (int32_t j = 0; j < n; ++j) {
                            const int32_t dd = wd + rowB[j];
                            if (rowB[j] < FHL_INFTY && dd < bw[j])
                                bw[j] = dd;
                        }
                    }
                }
            }
        };
        if (regionList != nullptr) {
#pragma omp parallel for schedule(dynamic, 4)
            for (int32_t i = 0; i < numInList; ++i)
                doRegion(regionList[i]);
        } else {
#pragma omp parallel for schedule(dynamic, 8)
            for (int32_t r = 0; r < R.leaf().numRegions(); ++r)
                doRegion(r);
        }
    }

};

} // namespace fhl
