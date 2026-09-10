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

// SCHEME 1 -- plain FHL. Per-vertex seeds to the boundary points of the own leaf region, times one
// per-region [boundary points x path-ancestor hubs] rectangle, canonically pruned to a
// per-column CSR
// with column minima as admissible skip bounds. A query folds seeds against the column
// slice shared by both endpoints. This is the base every other scheme extends.
struct SeedRectangle {
    template<typename HData>
    static void buildSeedsFromRect(const Regions &R, const HubHierarchy &hierarchy,
                                   HData &hdata, const bool symmetric,
                                   const std::vector<int64_t> &valOff,
                                   const std::vector<int32_t> &denseOut,
                                   const std::vector<int32_t> &denseIn, Values &out) {
        const auto &base = hdata.getBaseData();
        const auto &c = R.leaf();
        out.seedPos.assign(R.numVertices + 1, 0);
        for (int32_t v = 0; v < R.numVertices; ++v) {
            const int32_t r = R.leafRegionOfVertex[v];
            out.seedPos[v + 1] = out.seedPos[v] + (r < 0 ? 0 : c.bSize(r));
        }
        out.seedFwd.allocateForOverwrite(out.seedPos[R.numVertices]);
        if (symmetric)
            out.seedBwd.clear();
        else
            out.seedBwd.allocateForOverwrite(out.seedPos[R.numVertices]);

        std::vector<int8_t> levelOfId(R.numHubs);
        for (int32_t id = 0; id < R.numHubs; ++id)
            levelOfId[id] = hierarchy.getVertexLevel(hierarchy.rankOfHub(id));

#pragma omp parallel
        {
#pragma omp for schedule(dynamic, 8)
        for (int32_t r = 0; r < c.numRegions(); ++r) {
            const int32_t n = static_cast<int32_t>(c.bSize(r));
            if (n == 0)
                continue;
            const int32_t *const lend = R.flatLevelRow(r);
            const int32_t *const cstart = R.flatColRow(r);
            // The rectangle is already column-major, so column j (all |B| boundary points against
            // ancestor j) is contiguous -- exactly what the per-access-hub fold below reads.
            // dOutCM column j = d(b_i -> col_j); dInCM column j = d(col_j -> b_i).
            const int32_t *const dOutCM = denseOut.data() + valOff[r];
            const int32_t *const dInCM =
                    symmetric ? dOutCM : denseIn.data() + valOff[r];
            for (int64_t k = R.leafVertFirst[r]; k < R.leafVertFirst[r + 1]; ++k) {
                const int32_t v = R.leafVertIds[k];
                int32_t *const fw = out.seedFwd.data() + out.seedPos[v];
                int32_t *const bw =
                        symmetric ? nullptr : out.seedBwd.data() + out.seedPos[v];
                std::fill(fw, fw + n, FHL_INFTY);
                if (bw)
                    std::fill(bw, bw + n, FHL_INFTY);
                const int32_t idx = base.rankToIdx(v);
                for (int32_t e = base.pos[idx]; e < base.pos[idx + 1]; ++e) {
                    const HubId a = base.accessHubs[e];
                    const int32_t fa = base.forwardDistances[e];
                    const int32_t ba =
                            base.backwardDistances.empty() ? fa : base.backwardDistances[e];
                    if (fa >= FHL_INFTY && ba >= FHL_INFTY)
                        continue;
                    // Locate a within its level segment of the (level, id)-ordered columns.
                    // The segment IS a contiguous run of bucketIds, so binary-search there
                    // instead of materializing the region's column list.
                    const int32_t la = levelOfId[a];
                    const int32_t segBeg = la == 0 ? 0 : lend[la - 1];
                    const HubId *const beg = R.bucketIds.data() + cstart[la];
                    const auto it = std::lower_bound(beg, beg + (lend[la] - segBeg), a);
                    KASSERT(it != beg + (lend[la] - segBeg) && *it == a);
                    const int64_t ca = segBeg + (it - beg);
                    if (fa < FHL_INFTY) {
                        const int32_t *const colIn = dInCM + ca * n; // d(a -> boundary_*)
                        for (int32_t j = 0; j < n; ++j) {
                            const int32_t d = fa + colIn[j];
                            if (d < fw[j])
                                fw[j] = d;
                        }
                    }
                    if (bw && ba < FHL_INFTY) {
                        const int32_t *const colOut = dOutCM + ca * n; // d(boundary_* -> a)
                        for (int32_t j = 0; j < n; ++j) {
                            const int32_t d = colOut[j] + ba;
                            if (d < bw[j])
                                bw[j] = d;
                        }
                    }
                }
            }
        }
        }
        (void) hdata;
    }


    // ---- the three per-region steps both the full and the single-region path need. They live
    // here once so that a layout change (the rectangle is COLUMN-major: entry
    // (boundary point i, column
    // j) is at j * n + i) cannot be applied to one path and missed in the other.

    // Boundary-to-boundary clique of one region, used as the domination test's hop table.
    template<typename HData>
    static void fillRegionClique(HData &hdata, const HubId *const rows, const int32_t n,
                                 int32_t *const M) {
        for (int32_t i = 0; i < n; ++i)
            for (int32_t j = 0; j < n; ++j)
                M[static_cast<int64_t>(i) * n + j] =
                        i == j ? 0 : hdata.getComparableDistance(rows[i], rows[j]);
    }

    // [boundary points x columns] rectangle of one region, column-major, both directions.
    template<typename HData>
    static void fillRegionRectangle(HData &hdata, const bool symmetric,
                                    const HubId *const rows, const int32_t n,
                                    const HubId *const cols, const int64_t nCols,
                                    int32_t *const rectOut, int32_t *const rectIn) {
        for (int64_t j = 0; j < nCols; ++j)
            for (int32_t i = 0; i < n; ++i) {
                rectOut[j * n + i] = hdata.getComparableDistance(rows[i], cols[j]);
                if (!symmetric)
                    rectIn[j * n + i] = hdata.getComparableDistance(cols[j], rows[i]);
            }
    }

    // Marks the survivors of the canonical prune, per column. perColumn may be null.
    static int64_t markSurvivors(const int32_t *const rect, const int32_t *const leafM,
                                 const int32_t n, const int64_t nCols, const bool inDir,
                                 uint8_t *const keep, int64_t *const perColumn) {
        int64_t total = 0;
        for (int64_t j = 0; j < nCols; ++j) {
            uint8_t *const kpc = keep + j * n;
            int32_t cnt = 0;
            for (int32_t w = 0; w < n; ++w) {
                const bool sv = survivesEntry(rect, leafM, n, w, j, inDir);
                kpc[w] = sv;
                cnt += sv;
            }
            if (perColumn != nullptr)
                perColumn[j] = cnt;
            total += cnt;
        }
        return total;
    }

    template<typename HData>
    static void buildLeafClique(const Regions &R, HData &hdata, Values &out) {
        const auto &c = R.leaf();
        out.leafClique.off.assign(c.numRegions() + 1, 0);
        for (int32_t r = 0; r < c.numRegions(); ++r)
            out.leafClique.off[r + 1] = out.leafClique.off[r] + c.bSize(r) * c.bSize(r);
        out.leafClique.M.resize(out.leafClique.off[c.numRegions()]);
#pragma omp parallel for schedule(dynamic, 8)
        for (int32_t r = 0; r < c.numRegions(); ++r) {
            const HubId *const ids = c.bIds.data() + c.bFirst[r];
            const int32_t n = static_cast<int32_t>(c.bSize(r));
            fillRegionClique(hdata, ids, n, out.leafClique.M.data() + out.leafClique.off[r]);
        }
    }


    // Hub value rows: the degenerate "seeds = {self: 0}" case; every column lies on the
    // node's own root path or in its own separator, so single lookups suffice.
    // ids == nullptr refreshes every hub and (re)allocates; otherwise only the listed hubs
    // are refreshed in place, which is what the incremental gears call.
    template<typename HData>
    static void buildHubRows(const Regions &R, HData &hdata, const bool symmetric, Values &out,
                             const HubId *const ids = nullptr, const int32_t numIds = 0) {
        if (ids == nullptr) {
            out.tRowFwd.assign(R.tColFirst[R.numHubs], FHL_INFTY);
            out.tRowBwd.assign(symmetric ? 0 : R.tColFirst[R.numHubs], FHL_INFTY);
        }
        const int32_t count = ids == nullptr ? R.numHubs : numIds;
#pragma omp parallel for schedule(dynamic, 64)
        for (int32_t i = 0; i < count; ++i) {
            const HubId x = ids == nullptr ? static_cast<HubId>(i) : ids[i];
            const int64_t base = R.tColFirst[x];
            R.forEachCol(R.tColRow(x), R.tLevelRow(x),
                         [&](const int32_t o, const HubId u) {
                             out.tRowFwd[base + o] =
                                     u == x ? 0 : hdata.getComparableDistance(x, u);
                             if (!symmetric)
                                 out.tRowBwd[base + o] =
                                         u == x ? 0 : hdata.getComparableDistance(u, x);
                         });
        }
    }


    // Canonically pruned per-column CSR of the [boundary x path ancestors] rectangles: entry
    // (w, u) is dropped when another boundary member w' gives d(w,w') + d(w',u) <= d(w,u)
    // (ties toward the smaller row index; the single-hop test subsumes chains). Sound for
    // every vertex because seeds are exact: seed(w') <= seed(w) + d(w,w'). Column minima
    // over the survivors equal those over all rows (a dropped entry's cover is as good).
    template<typename HData>
    static void buildCsr(const Regions &R, HData &hdata, const bool symmetric, Values &out,
                         std::vector<int64_t> &valOff, std::vector<int32_t> &denseOut,
                         std::vector<int32_t> &denseIn, const bool relocatable = false) {
        const auto &c = R.leaf();
        const int32_t numLeaf = c.numRegions();
        for (int32_t r = 0; r < numLeaf; ++r)
            if (static_cast<int32_t>(c.bSize(r)) > Values::FlatCsr::MAX_ROWS)
                throw std::invalid_argument(
                        "leaf region " + std::to_string(r) + " has " +
                        std::to_string(c.bSize(r)) + " boundary points, more than the " +
                        std::to_string(Values::FlatCsr::MAX_ROWS) +
                        " the CSR row index holds; lower -fhl-region-size.");
        const int64_t totalCols = R.flatColFirst[numLeaf];

        // Dense rectangles (row-major over the full column lists); handed back to the caller
        // for the seed fold, freed there.
        valOff.assign(numLeaf + 1, 0);
        for (int32_t r = 0; r < numLeaf; ++r)
            valOff[r + 1] = valOff[r] +
                    c.bSize(r) * (R.flatColFirst[r + 1] - R.flatColFirst[r]);
        denseOut.resize(valOff[numLeaf]);
        denseIn.resize(symmetric ? 0 : valOff[numLeaf]);
        std::vector<HubId> colBuf;
#pragma omp parallel for schedule(dynamic, 8) firstprivate(colBuf)
        for (int32_t r = 0; r < numLeaf; ++r) {
            const HubId *const rows = c.bIds.data() + c.bFirst[r];
            const int32_t numRows = static_cast<int32_t>(c.bSize(r));
            R.gatherFlatCols(r, colBuf);
            const HubId *const cols = colBuf.data();
            const int64_t numCols = R.flatColFirst[r + 1] - R.flatColFirst[r];
            fillRegionRectangle(hdata, symmetric, rows, numRows, cols, numCols,
                                denseOut.data() + valOff[r],
                                symmetric ? nullptr : denseIn.data() + valOff[r]);
        }

        const auto buildOne = [&](const std::vector<int32_t> &dense, const bool inDir,
                                  Values::FlatCsr &csr, std::vector<int32_t> &mins) {
            csr.ptr.assign(totalCols + 1, 0);
            // Verdicts cached (the domination test is O(|B|) per entry; running it in both
            // the counting and the filling pass doubled the dominant cost).
            std::vector<uint8_t> keep(valOff[numLeaf]);
#pragma omp parallel for schedule(dynamic, 8)
            for (int32_t r = 0; r < numLeaf; ++r) {
                const int32_t n = static_cast<int32_t>(c.bSize(r));
                const int64_t nCols = R.flatColFirst[r + 1] - R.flatColFirst[r];
                const int32_t *const rect = dense.data() + valOff[r];
                const int32_t *const leafM = out.leafClique.M.data() + out.leafClique.off[r];
                uint8_t *const kp = keep.data() + valOff[r];
                markSurvivors(rect, leafM, n, nCols, inDir, kp,
                              csr.ptr.data() + R.flatColFirst[r] + 1);
            }
            for (int64_t k = 0; k < totalCols; ++k)
                csr.ptr[k + 1] += csr.ptr[k];
            csr.row.resize(csr.ptr[totalCols]);
            // Gather mode: the values live in the hub rows already (every rectangle
            // entry is d(boundary point, ancestor) between two ancestor-related hubs), so
            // only the survivor list is stored and the query reads through.
            csr.val.resize(csr.ptr[totalCols]);
            mins.assign(totalCols, FHL_INFTY);
#pragma omp parallel for schedule(dynamic, 8)
            for (int32_t r = 0; r < numLeaf; ++r) {
                const int32_t n = static_cast<int32_t>(c.bSize(r));
                const int64_t nCols = R.flatColFirst[r + 1] - R.flatColFirst[r];
                const int32_t *const rect = dense.data() + valOff[r];
                const uint8_t *const kp = keep.data() + valOff[r];
                for (int64_t j = 0; j < nCols; ++j) {
                    int64_t k = csr.ptr[R.flatColFirst[r] + j];
                    auto &mn = mins[R.flatColFirst[r] + j];
                    const int32_t *const col = rect + j * n;
                    const uint8_t *const kpc = kp + j * n;
                    for (int32_t w = 0; w < n; ++w)
                        if (kpc[w]) {
                            csr.row[k] = static_cast<Values::FlatCsr::RowIdx>(w);
                            csr.val[k] = col[w];
                            mn = std::min(mn, col[w]);
                            ++k;
                        }
                    KASSERT(k == csr.ptr[R.flatColFirst[r] + j + 1]);
                }
            }
        };
        buildOne(denseOut, false, out.csrOut, out.colMinOut);
        if (!symmetric)
            buildOne(denseIn, true, out.csrIn, out.colMinIn);
        // Block minima, per region, over its own column range.
        out.colBlockFirst.assign(numLeaf + 1, 0);
        for (int32_t r = 0; r < numLeaf; ++r) {
            const int64_t nCols = R.flatColFirst[r + 1] - R.flatColFirst[r];
            out.colBlockFirst[r + 1] = out.colBlockFirst[r] +
                    (nCols + Values::COL_BLOCK - 1) / Values::COL_BLOCK;
        }
        const auto blocks = [&](const std::vector<int32_t> &mins,
                                std::vector<int32_t> &blockMin) {
            blockMin.assign(out.colBlockFirst[numLeaf], FHL_INFTY);
#pragma omp parallel for schedule(dynamic, 16)
            for (int32_t r = 0; r < numLeaf; ++r) {
                const int64_t colBase = R.flatColFirst[r];
                const int64_t nCols = R.flatColFirst[r + 1] - colBase;
                for (int64_t j = 0; j < nCols; ++j) {
                    auto &m = blockMin[out.colBlockFirst[r] + j / Values::COL_BLOCK];
                    m = std::min(m, mins[colBase + j]);
                }
            }
        };
        blocks(out.colMinOut, out.colBlockMinOut);
        if (!symmetric)
            blocks(out.colMinIn, out.colBlockMinIn);
        if (relocatable) {
            const auto setup = [&](Values::FlatCsr &csr) {
                if (!csr.active())
                    return;
                csr.colEnd.assign(totalCols, 0);
                for (int64_t k = 0; k < totalCols; ++k)
                    csr.colEnd[k] = csr.ptr[k + 1];
                csr.segStart.assign(numLeaf, 0);
                csr.segCap.assign(numLeaf, 0);
                for (int32_t r = 0; r < numLeaf; ++r) {
                    csr.segStart[r] = csr.ptr[R.flatColFirst[r]];
                    csr.segCap[r] = csr.ptr[R.flatColFirst[r + 1]] - csr.segStart[r];
                }
            };
            setup(out.csrOut);
            setup(out.csrIn);
        }
        std::cout << "Flat CSR: avg "
                  << (totalCols > 0
                              ? static_cast<double>(out.csrOut.ptr[totalCols]) / totalCols
                              : 0.0)
                  << " surviving rows per column." << std::endl;
    }



    // ---------------------------------------------------------- super-region rebuild parts
    // Rebuilds one leaf region's clique and pruned CSR segment from the (already refreshed)
    // span-1 tables. The segment is rewritten in place when the new survivor count fits its
    // original capacity, otherwise the whole segment moves to the tail of the arrays; per
    // column ends make either layout valid for the query. Requires the relocatable layout.
    template<typename HData>
    static void rebuildLeafRegionFromTables(const Regions &R, HData &hdata,
                                            const bool symmetric, const int32_t r,
                                            Values &out) {
        const auto &c = R.leaf();
        const int32_t n = static_cast<int32_t>(c.bSize(r));
        const HubId *const rows = c.bIds.data() + c.bFirst[r];
        int32_t *const M = out.leafClique.M.data() + out.leafClique.off[r];
        fillRegionClique(hdata, rows, n, M);
        if (n == 0)
            return;

        const int64_t colBase = R.flatColFirst[r];
        const int64_t nCols = R.flatColFirst[r + 1] - colBase;
        std::vector<HubId> colBuf;
        R.gatherFlatCols(r, colBuf);
        const HubId *const cols = colBuf.data();
        std::vector<int32_t> rect(static_cast<size_t>(n) * nCols);
        std::vector<int32_t> rectIn(symmetric ? 0 : static_cast<size_t>(n) * nCols);
        fillRegionRectangle(hdata, symmetric, rows, n, cols, nCols, rect.data(),
                            symmetric ? nullptr : rectIn.data());

        const auto writeSeg = [&](const std::vector<int32_t> &dense, const bool inDir,
                                  Values::FlatCsr &csr, std::vector<int32_t> &mins) {
            if (!csr.active())
                return;
            std::vector<uint8_t> keep(static_cast<size_t>(n) * nCols);
            const int64_t total =
                    markSurvivors(dense.data(), M, n, nCols, inDir, keep.data(), nullptr);
            int64_t at = csr.segStart[r];
            if (total > csr.segCap[r]) { // grew past its slot: relocate to the tail
                at = static_cast<int64_t>(csr.row.size());
                csr.row.resize(at + total);
                if (!csr.val.empty())
                    csr.val.resize(at + total);
                csr.segStart[r] = at;
                csr.segCap[r] = total;
            }
            for (int64_t j = 0; j < nCols; ++j) {
                csr.ptr[colBase + j] = at;
                int32_t mn = FHL_INFTY;
                for (int32_t w = 0; w < n; ++w)
                    if (keep[j * n + w]) {
                        csr.row[at] = static_cast<Values::FlatCsr::RowIdx>(w);
                        csr.val[at] = dense[j * n + w];
                        mn = std::min(mn, dense[j * n + w]);
                        ++at;
                    }
                csr.colEnd[colBase + j] = at;
                mins[colBase + j] = mn;
            }
        };
        writeSeg(rect, false, out.csrOut, out.colMinOut);
        if (!symmetric)
            writeSeg(rectIn, true, out.csrIn, out.colMinIn);
    }


    // The domination test of buildCsr, shared with the incremental region rebuild.
    // rect is COLUMN-major: entry (boundary point w, column j) is rect[j * n + w], so all of a
    // column's candidates sit on one cache line.
    static bool survivesEntry(const int32_t *const rect, const int32_t *const leafM,
                              const int32_t n, const int32_t w,
                              const int64_t j, const bool inDir) {
        const int32_t *const col = rect + j * n;
        const int32_t dw = col[w];
        if (dw >= FHL_INFTY)
            return false;
        for (int32_t w2 = 0; w2 < n; ++w2) {
            if (w2 == w)
                continue;
            const int32_t dw2 = col[w2];
            if (dw2 >= FHL_INFTY)
                continue;
            const int32_t hop = inDir ? leafM[static_cast<int64_t>(w2) * n + w]
                                      : leafM[static_cast<int64_t>(w) * n + w2];
            if (hop >= FHL_INFTY)
                continue;
            const int32_t via = hop + dw2;
            if (via < dw || (via == dw && w2 < w))
                return false;
        }
        return true;
    }

};

} // namespace fhl
