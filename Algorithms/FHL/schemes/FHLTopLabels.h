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
#include "Algorithms/FHL/schemes/FHLSeedRectangle.h"

namespace fhl {

// SCHEME 3 -- FHL + top labels. For the columns above level K the fused value
// min_w(seed[w] + rect[w][j]) is vertex-independent work that a query would redo every
// time, so it is materialized per vertex. Far queries (LCA level < K) then reduce to one
// contiguous min over two label arrays -- no seeds, no CSR, no skip bounds.
struct TopLabels {
    // Optional per-vertex top labels: the fused outputs for columns with level < KL,
    // computed from the CSR survivors (valid for every seed vector, hence still exact).
    static void buildTopLabels(const Regions &R, const bool symmetric, const int32_t KL,
                               Values &out) {
        out.topLabLevels = KL;
        out.topLabPos.assign(R.numVertices + 1, 0);
        for (int32_t v = 0; v < R.numVertices; ++v) {
            const int32_t r = R.leafRegionOfVertex[v];
            out.topLabPos[v + 1] = out.topLabPos[v] +
                    (r < 0 ? 0
                           : R.flatLevelEnd[static_cast<int64_t>(r) * R.levelStride + KL - 1]);
        }
        out.topLabFwd.allocateForOverwrite(out.topLabPos[R.numVertices]);
        if (symmetric)
            out.topLabBwd.clear();
        else
            out.topLabBwd.allocateForOverwrite(out.topLabPos[R.numVertices]);
#pragma omp parallel for schedule(dynamic, 8)
        for (int32_t r = 0; r < R.leaf().numRegions(); ++r) {
            const int64_t base = R.flatColFirst[r];
            for (int64_t k = R.leafVertFirst[r]; k < R.leafVertFirst[r + 1]; ++k) {
                const int32_t v = R.leafVertIds[k];
                const int64_t nLab = out.topLabPos[v + 1] - out.topLabPos[v];
                if (nLab == 0)
                    continue;
                const auto fuse = [&](const Values::FlatCsr &csr, const int32_t *const seeds,
                                      int32_t *const lab) {
                    for (int64_t j = 0; j < nLab; ++j) {
                        int32_t best = FHL_INFTY;
                        const int64_t *const ce = csr.endArray();
                        for (int64_t kk = csr.ptr[base + j]; kk < ce[base + j]; ++kk) {
                            const int32_t d = seeds[csr.row[kk]] + csr.val[kk];
                            if (d < best)
                                best = d;
                        }
                        lab[j] = best > FHL_INFTY ? FHL_INFTY : best;
                    }
                };
                fuse(out.csrOut, out.seedFwd.data() + out.seedPos[v],
                     out.topLabFwd.data() + out.topLabPos[v]);
                if (!symmetric)
                    fuse(out.csrIn, out.seedBwd.data() + out.seedPos[v],
                         out.topLabBwd.data() + out.topLabPos[v]);
            }
        }
        std::cout << "Flat top labels: levels < " << KL << ", avg "
                  << static_cast<double>(out.topLabPos[R.numVertices]) /
                     std::max(1, R.numVertices)
                  << " entries/vertex." << std::endl;
    }

};

} // namespace fhl
