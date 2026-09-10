#pragma once

#include <algorithm>
#include <climits>
#include <vector>

#include "Algorithms/CCH/EliminationTreeQuery.h"
#include "Algorithms/FHL/FHLData.h"
#include "DataStructures/Labels/BasicLabelSet.h"
#include "DataStructures/Labels/ParentInfo.h"
#include "Tools/Constants.h"
#include "Algorithms/FHL/core/FHLConstants.h"

// The query side of all four FHL schemes.
//
// A shortest s-t path must cross the separator shared by s and t, so with l = LCA level,
//     d(s,t) = min over hubs u with level(u) <= l of  d(s,u) + d(u,t).
// Those candidate u ARE the first n columns of both endpoints' column lists, in the same
// order (FHLRegions guarantees the alignment), so the combine is index-parallel: no id
// comparisons, no search, no priority queue.
//
// Dispatch (run):
//   same leaf region   -> boundary-point clique bounded by the in-region answer, which is
//                         the in-region hub fold when -fhl-subregion-size is on and the two
//                         sit in different sub-regions, otherwise a CCH search
//   overlay            -> climb cut by cut, meeting at every common region
//   otherwise          -> one aligned column scan (runFlatQuery), deepest column first under
//                         an admissible skip bound; when -fhl-top-labels covers the whole slice
//                         it collapses to a contiguous min over two label arrays
class FHLQuery {
public:
    using LabelSet = BasicLabelSet<0, ParentInfo::NO_PARENT_INFO>;

    FHLQuery(const HubHierarchy &hierarchy,
                          const FHLData &data,
                          const std::vector<int> &localEliminationTree,
                          const LayerCCH::UpGraph &cchGraph,
                          int const *const cchUpWeights,
                          int const *const cchDownWeights,
                          const fhl::Regions *const ladderRegions = nullptr)
            : hierarchy(hierarchy),
              data(data),
              ladderRegions(ladderRegions),
              localQuery(cchGraph, cchUpWeights, cchDownWeights, localEliminationTree) {
        if (ladderRegions != nullptr) {
            // LCA level -> index of the smallest cut depth above it.
            ciMeetOfLevel.assign(ladderRegions->thresh, 0);
            for (int32_t l = 0; l < ladderRegions->thresh; ++l) {
                int32_t ci = 0;
                while (ladderRegions->cutDepths[ci] <= l)
                    ++ci;
                ciMeetOfLevel[l] = ci;
            }
        }
    }

    int32_t run(int32_t s, int32_t t) {
        const fhl::Level lcaLevel = hierarchy.getLevelOfLowestCommonAncestor(s, t);
        KASSERT(data.ladder.active());
        int32_t dist;
        {
            lastModeIsLocal = !hierarchy.isHub(s) && !hierarchy.isHub(t) &&
                              ladderRegions->leafRegionOfVertex[s] ==
                                      ladderRegions->leafRegionOfVertex[t];
            // FHL_LOCAL_VIA_COLUMNS: answer same-region pairs with the ordinary column
            // scan instead of the boundary-point clique. Sound: a path that leaves the region
            // meets at a common path ancestor, i.e. one of the region's columns, and both
            // endpoints share the column list, so the aligned fold covers it.
            // Same-region pairs in the seed/rectangle family: the in-region answer is
            // computed FIRST and then bounds the boundary-pair fold, which is almost always the
            // loser (a path that leaves the region and comes back is rarely shorter).
            const bool boundaryPairOnly =
                    data.ladder.mode != 2 && lastModeIsLocal && !localViaColumns;
            if (data.ladder.mode == 2)
                dist = runOverlayQuery(s, t, lcaLevel);
            else if (!boundaryPairOnly)
                dist = runFlatQuery(s, t, lcaLevel);
            else
                dist = FHL_INFTY;
            if (lastModeIsLocal) {
                ++localPairs;
                // In-region hub fold: if the pair's LCA lies in the stored top levels of
                // the region's own hierarchy, the crossing vertex is one of the region's own
                // hubs, so the fold is exact and the CCH fallback can be skipped entirely.
                bool covered = false;
                int32_t inRegion = FHL_INFTY;
                const auto &L = data.ladder;
                const auto &Rg = *ladderRegions;
                if (L.locLevels > 0) {
                    const int32_t r = Rg.leafRegionOfVertex[s];
                    // Different sub-regions => every common ancestor sits above the inner
                    // cut, i.e. among the stored columns, so the fold is exact.
                    if (Rg.innerRegionOfVertex[s] != Rg.innerRegionOfVertex[t] ||
                        Rg.innerRegionOfVertex[s] < 0) {
                        const int32_t n = static_cast<int32_t>(Rg.locColFirst[r + 1] -
                                                               Rg.locColFirst[r]);
                        const int32_t *const a = L.locFwd.data() + L.locPos[s];
                        const int32_t *const b =
                                (L.locBwd.empty() ? L.locFwd.data() : L.locBwd.data()) +
                                L.locPos[t];
                        int32_t best = FHL_INFTY;
                        for (int32_t j = 0; j < n; ++j) {
                            const int32_t d2 = a[j] + b[j];
                            if (d2 < best)
                                best = d2;
                        }
                        inRegion = best;
                        covered = true;
                        ++localFolds;
                    }
                }
                // Measurement switch: FHL_SKIP_LOCAL=1 drops the in-region CCH fallback
                // (WRONG answers for same-region pairs) to expose its cost share.
                if (!covered && !skipLocalFallback) {
                    localQuery.run(s, t);
                    inRegion = localQuery.getDistance();
                }
                dist = boundaryPairOnly ? runSeedLeafPairQuery(s, t, inRegion)
                                    : std::min(dist, inRegion);
            }
            lastDistance = dist;
            return dist;
        }
    }

    int32_t getDistance() const { return lastDistance; }
    int64_t getEarlyStops() const { return earlyStops; }
    int64_t getLocalPairs() const { return localPairs; }
    int64_t getLabelOnlyPairs() const { return labelOnlyPairs; }
    int64_t getLocalFolds() const { return localFolds; }
    inline static const bool useBlockBound = std::getenv("FHL_NO_BLOCK") == nullptr;
    // Ablation toggle for the contiguous label-prefix merge.
    inline static const bool useLabelFast = std::getenv("FHL_NO_LABFAST") == nullptr;
    inline static const bool localViaColumns =
            std::getenv("FHL_LOCAL_VIA_COLUMNS") != nullptr;
    int64_t localFolds = 0;
    inline static const bool skipLocalFallback =
            std::getenv("FHL_SKIP_LOCAL") != nullptr;
    int64_t localPairs = 0;
    mutable int64_t labelOnlyPairs = 0;
    const std::array<int64_t, 8> &getStopLevels() const { return stopAtLevel; }
    std::array<int64_t, 8> stopAtLevel{};
    // Experiment switch: FHL_NO_EARLY=1 disables the admissible ladder cutoff.
    inline static const bool overlayEarlyStop = std::getenv("FHL_NO_EARLY") == nullptr;
    int64_t earlyStops = 0;

    const char *getLastMode() const { return lastModeIsLocal ? "local" : "hub"; }

    uint64_t sizeInBytes() const {
        uint64_t size = sizeof(FHLQuery);
        size += (middleIds.capacity() + distToMiddle.capacity() + distFromMiddle.capacity()) *
                sizeof(int32_t);
        size += localQuery.sizeInBytes();
        return size;
    }

private:
    // ------------------------------------------------------------------ ladder-scheme queries

    struct LadderFrontier {
        int32_t region = -1;   // cut region at the meet cut; frontier values aligned to its B
        int32_t n = 0;
        const int32_t *vals = nullptr;
        int32_t selfId = -1;   // >= 0: hub endpoint inside the meet region's sep range,
        int32_t selfSPos = -1; //       participating as itself at distance 0
    };

    // Combines two frontiers on cut mi through the meet clique of their common region.
    int32_t ladderCliqueMeet(const LadderFrontier &fs, const LadderFrontier &ft,
                             const int32_t mi) {
        const auto &R = *ladderRegions;
        const auto &L = data.ladder;
        const auto &c = R.cuts[mi];
        const auto &ms = R.meets[mi];
        int32_t mr;
        if (fs.region >= 0)
            mr = c.parentRegion[fs.region];
        else if (ft.region >= 0)
            mr = c.parentRegion[ft.region];
        else
            mr = mi == 0 ? 0 : R.initRegionOfHub[fs.selfId];
        KASSERT(fs.region < 0 || ft.region < 0 ||
                c.parentRegion[fs.region] == c.parentRegion[ft.region]);
        const int32_t Sn = static_cast<int32_t>(ms.sSize(mr));
        if (Sn == 0)
            return FHL_INFTY;
        const int32_t *const M = L.meetClique[mi].M.data() + L.meetClique[mi].off[mr];

        const int32_t *sPosS, *valsS, *sPosT, *valsT;
        int32_t nS, nT;
        if (fs.region >= 0) {
            sPosS = c.sPosOfB.data() + c.bFirst[fs.region];
            valsS = fs.vals;
            nS = fs.n;
        } else {
            sPosS = &fs.selfSPos;
            valsS = &ladderZero;
            nS = 1;
        }
        if (ft.region >= 0) {
            sPosT = c.sPosOfB.data() + c.bFirst[ft.region];
            valsT = ft.vals;
            nT = ft.n;
        } else {
            sPosT = &ft.selfSPos;
            valsT = &ladderZero;
            nT = 1;
        }
        int32_t best = FHL_INFTY;
        for (int32_t i = 0; i < nS; ++i) {
            const int32_t v = valsS[i];
            if (v >= best) // clique and target values are non-negative
                continue;
            const int32_t *const row = M + static_cast<int64_t>(sPosS[i]) * Sn;
            for (int32_t j = 0; j < nT; ++j) {
                const int32_t d = v + row[sPosT[j]] + valsT[j];
                if (d < best)
                    best = d;
            }
        }
        return best;
    }

    // Local pair query over the leaf-region boundary clique (both endpoints share the leaf
    // region when the LCA level reaches the hub level threshold; both are non-hubs).
    // `best` comes in as the in-region answer: a path through a boundary point only
    // matters if it beats it, and that bound makes the outer skip fire from the first
    // boundary point instead of never.
    int32_t runSeedLeafPairQuery(const int32_t s, const int32_t t, int32_t best) {
        const auto &R = *ladderRegions;
        const auto &L = data.ladder;
        const int32_t r = R.leafRegionOfVertex[s];
        KASSERT(r >= 0 && r == R.leafRegionOfVertex[t]);
        const int32_t n = static_cast<int32_t>(R.leaf().bSize(r));
        if (n == 0)
            return FHL_INFTY;
        const int32_t *const fsv = L.seedFwd.data() + L.seedPos[s];
        const int32_t *const ftv =
                (L.seedBwd.empty() ? L.seedFwd.data() : L.seedBwd.data()) + L.seedPos[t];
        const int32_t *const M = L.leafClique.M.data() + L.leafClique.off[r];
        for (int32_t i = 0; i < n; ++i) {
            const int32_t v = fsv[i];
            if (v >= FHL_INFTY || v >= best)
                continue;
            const int32_t *const row = M + static_cast<int64_t>(i) * n;
            for (int32_t j = 0; j < n; ++j) {
                const int32_t d = v + row[j] + ftv[j];
                if (d < best)
                    best = d;
            }
        }
        return best;
    }

    // Flat query. Both endpoints provide exact values for the shared columns (all path
    // ancestors of the common prefix with level <= LCA level; both column lists agree there,
    // position by position). A non-hub endpoint fuses its seeds with the region's pruned
    // CSR — or reads its materialized top labels below the label cut; a hub endpoint
    // reads its stored value row. Columns are evaluated deepest-first under an admissible
    // lower-bound skip (minSeed + column minimum per side). Exact by the shared-path
    // middle-node argument: the max-rank vertex of a shortest path is a path ancestor of
    // both branches with level <= lambda, sits at the same position in both slices, and both
    // sides are tight there.
    int32_t runFlatQuery(const int32_t s, const int32_t t, const int32_t lcaLevel) {
        const auto &R = *ladderRegions;
        const auto &L = data.ladder;

        struct Side {
            const int32_t *rowVals = nullptr; // hub endpoint: stored value row
            const int32_t *seeds = nullptr;   // non-hub: seed row and CSR base
            int32_t region = -1;              // non-hub: its leaf region
            const fhl::Values::FlatCsr *csr = nullptr;
            const int64_t *csrEnd = nullptr;
            int64_t base = 0;
            const int32_t *lab = nullptr;     // non-hub: top label row (nullable)
            int32_t nLab = 0;
            const int32_t *cm = nullptr;      // admissible per-column minima
            const int32_t *bm = nullptr;      // admissible per-block minima (32 columns)
            int32_t ms = 0;                   // minimum seed value (0 for a hub)
            int32_t n = -1;                   // columns with level <= lcaLevel
        };
        const auto makeSide = [&](const int32_t v, const bool forward, Side &sd) {
            if (hierarchy.isHub(v)) {
                const int32_t x = hierarchy.hubIdOfRank(v);
                sd.rowVals = (forward || L.tRowBwd.empty() ? L.tRowFwd.data()
                                                           : L.tRowBwd.data()) +
                             R.tColFirst[x];
                sd.cm = sd.rowVals; // exact values are their own lower bounds
                sd.bm = nullptr;    // hub rows have no block index
                sd.n = R.tLevelEnd[static_cast<int64_t>(x) * R.levelStride + lcaLevel];
                return;
            }
            const int32_t r = R.leafRegionOfVertex[v];
            sd.seeds = (forward || L.seedBwd.empty() ? L.seedFwd.data() : L.seedBwd.data()) +
                       L.seedPos[v];
            sd.csr = forward || !L.csrIn.active() ? &L.csrOut : &L.csrIn;
            sd.csrEnd = sd.csr->endArray();
            sd.region = r;
            sd.base = R.flatColFirst[r];
            sd.cm = (forward || L.colMinIn.empty() ? L.colMinOut.data()
                                                   : L.colMinIn.data()) +
                    sd.base;
            if (!L.colBlockFirst.empty())
                sd.bm = (forward || L.colBlockMinIn.empty() ? L.colBlockMinOut.data()
                                                            : L.colBlockMinIn.data()) +
                        L.colBlockFirst[r];
            sd.n = R.flatLevelEnd[static_cast<int64_t>(r) * R.levelStride + lcaLevel];
            if (L.topLabLevels > 0) {
                sd.lab = (forward || L.topLabBwd.empty() ? L.topLabFwd.data()
                                                         : L.topLabBwd.data()) +
                         L.topLabPos[v];
                sd.nLab = std::min(
                        sd.n, R.flatLevelEnd[static_cast<int64_t>(r) * R.levelStride +
                                             L.topLabLevels - 1]);
            }
            sd.ms = FHL_INFTY;
            const int32_t nSeeds = static_cast<int32_t>(R.leaf().bSize(r));
            for (int32_t i = 0; i < nSeeds; ++i)
                sd.ms = std::min(sd.ms, sd.seeds[i]);
        };
        Side S, T;
        makeSide(s, true, S);
        makeSide(t, false, T);
        KASSERT(S.n == T.n);
        if (S.ms >= FHL_INFTY || T.ms >= FHL_INFTY)
            return FHL_INFTY;

        const auto eval = [&, this](const Side &sd, const int32_t j) {
            if (sd.rowVals != nullptr)
                return sd.rowVals[j];
            if (j < sd.nLab)
                return sd.lab[j];
            int32_t best = FHL_INFTY;
            const int64_t end = sd.csrEnd[sd.base + j];
            for (int64_t k = sd.csr->ptr[sd.base + j]; k < end; ++k) {
                const int32_t d = sd.seeds[sd.csr->row[k]] + sd.csr->val[k];
                if (d < best)
                    best = d;
            }
            return best;
        };
        // Whole slice inside the materialized label prefix (lcaLevel < topLabLevels): both
        // sides hold the EXACT fused value at every column, so the skip bounds cannot help --
        // they would only add two loads and a compare per column to prune nothing. Straight
        // contiguous min instead, the shape CTL's merge has (CTLQuery::computeMinDistance-
        // InLabels). No meeting-hub index to track, so it stays a pure horizontal min.
        if (useLabelFast && S.lab != nullptr && T.lab != nullptr &&
            S.nLab == S.n && T.nLab == T.n) {
            ++labelOnlyPairs;
            const int32_t *const a = S.lab;
            const int32_t *const b = T.lab;
            int32_t bestLab = FHL_INFTY;
            for (int32_t j = 0; j < S.n; ++j)
                bestLab = std::min(bestLab, a[j] + b[j]);
            return bestLab > FHL_INFTY ? FHL_INFTY : bestLab;
        }
        const int32_t msSum = S.ms + T.ms;
        int32_t best = FHL_INFTY;
        const auto scanRange = [&](const int32_t hi, const int32_t lo) {
            for (int32_t j = hi; j >= lo; --j) {
                if (msSum + S.cm[j] + T.cm[j] >= best)
                    continue;
                const int32_t bS = eval(S, j);
                if (bS >= best)
                    continue;
                const int32_t d = bS + eval(T, j);
                if (d < best)
                    best = d;
            }
        };
        // The block layer costs two extra cache lines (both sides' block minima), so it only
        // pays off once the slice is long enough for the saved per-column tests to exceed
        // them. On short slices it measured as a net loss, hence the threshold below.
        constexpr int32_t BLKMIN = fhl::Values::COL_BLOCK;
        if (S.bm != nullptr && T.bm != nullptr && useBlockBound && S.n >= 8 * BLKMIN) {
            // Block bound first: otherwise almost every column is tested individually.
            constexpr int32_t BLK = fhl::Values::COL_BLOCK;
            for (int32_t b = (S.n - 1) / BLK; b >= 0; --b) {
                if (msSum + S.bm[b] + T.bm[b] >= best)
                    continue;
                scanRange(std::min(S.n - 1, b * BLK + BLK - 1), b * BLK);
            }
        } else {
            scanRange(S.n - 1, 0);
        }
        return best > FHL_INFTY ? FHL_INFTY : best;
    }

    // One overlay climb step: frontier on cut ci -> cut ci-1, through the restricted clique
    // of the meet region above cut ci (climb rectangles are submatrices of these cliques in
    // overlay semantics, so no separate tables exist).
    void ladderOverlayClimb(const bool forward, int32_t &ci, int32_t &region,
                            const int32_t *&vals, int32_t &n, std::vector<int32_t> &buf) {
        const auto &R = *ladderRegions;
        const auto &L = data.ladder;
        const auto &c = R.cuts[ci];
        const int32_t mr = c.parentRegion[region];
        const auto &ms = R.meets[ci];
        const int32_t Sn = static_cast<int32_t>(ms.sSize(mr));
        const int32_t *const M = L.meetClique[ci].M.data() + L.meetClique[ci].off[mr];
        const int32_t *const rowPos = c.sPosOfB.data() + c.bFirst[region];
        const auto &pc = R.cuts[ci - 1];
        const int32_t nCols = static_cast<int32_t>(pc.bSize(mr));
        // The parent boundary points occupy S positions [0, nCols), so output column j IS
        // column j of the clique and the inner loop reads one contiguous row.
        KASSERT(ms.numParentBoundary[mr] == nCols);
        const bool symmetricClique = data.ladder.symmetric;
        buf.assign(std::max(nCols, 1), FHL_INFTY);
        for (int32_t i = 0; i < n; ++i) {
            const int32_t v = vals[i];
            if (v >= FHL_INFTY)
                continue;
            if (forward || symmetricClique) { // symmetric metric: clique is its own transpose
                const int32_t *const row = M + static_cast<int64_t>(rowPos[i]) * Sn;
                for (int32_t j = 0; j < nCols; ++j) {
                    const int32_t d = v + row[j];
                    if (d < buf[j])
                        buf[j] = d;
                }
            } else {
                for (int32_t j = 0; j < nCols; ++j) {
                    const int32_t d = v + M[static_cast<int64_t>(j) * Sn + rowPos[i]];
                    if (d < buf[j])
                        buf[j] = d;
                }
            }
        }
        vals = buf.data();
        n = nCols;
        region = mr;
        --ci;
    }

    // Overlay (scheme 1) query: meets at every common region from the deepest one up to the
    // root; each meet uses the pre-climb frontiers of both sides and the restricted clique.
    int32_t runOverlayQuery(const int32_t s, const int32_t t, const int32_t lcaLevel) {
        const auto &R = *ladderRegions;
        const auto &L = data.ladder;
        const int32_t m = R.numCuts();
        const bool localPair = lcaLevel >= R.thresh;
        const int32_t startMi = localPair ? m - 1 : ciMeetOfLevel[lcaLevel];
        int32_t best = localPair ? runSeedLeafPairQuery(s, t, FHL_INFTY) : FHL_INFTY;

        // Frontier minima: a climb only adds non-negative weights, so they never decrease.
        // Once minS + minT reaches the best distance found so far, no higher meet can beat
        // it and the ladder can stop — the admissible early termination that a Dijkstra-based
        // overlay search (CRP) gets for free, and which makes short trips cheap.
        const auto frontierMin = [](const int32_t *const v, const int32_t n) {
            int32_t mn = FHL_INFTY;
            for (int32_t i = 0; i < n; ++i)
                mn = std::min(mn, v[i]);
            return mn;
        };
        int32_t ciS = 0, regS = -1, nS = 0, selfS = -1, selfSPosS = -1;
        int32_t ciT = 0, regT = -1, nT = 0, selfT = -1, selfSPosT = -1;
        const int32_t *valsS = nullptr, *valsT = nullptr;
        const auto init = [&](const int32_t v, const bool fwd, int32_t &ci, int32_t &reg,
                              const int32_t *&vals, int32_t &n, int32_t &selfId,
                              int32_t &selfSPos) {
            if (hierarchy.isHub(v)) {
                const int32_t x = hierarchy.hubIdOfRank(v);
                ci = R.initCutOfHub[x];
                if (ci < startMi) {
                    selfId = x;
                    selfSPos = R.selfSPosOfHub[x];
                    reg = -1;
                    return;
                }
                reg = R.initRegionOfHub[x];
                vals = (fwd || L.tSeedBwd.empty() ? L.tSeedFwd.data() : L.tSeedBwd.data()) +
                       L.tSeedPos[x];
                n = static_cast<int32_t>(R.cuts[ci].bSize(reg));
            } else {
                ci = m - 1;
                reg = R.leafRegionOfVertex[v];
                vals = (fwd || L.seedBwd.empty() ? L.seedFwd.data() : L.seedBwd.data()) +
                       L.seedPos[v];
                n = static_cast<int32_t>(R.cuts[m - 1].bSize(reg));
            }
        };
        init(s, true, ciS, regS, valsS, nS, selfS, selfSPosS);
        init(t, false, ciT, regT, valsT, nT, selfT, selfSPosT);
        bool useA1 = false, useB1 = false;
        while (regS >= 0 && ciS > startMi) {
            useA1 = !useA1;
            ladderOverlayClimb(true, ciS, regS, valsS, nS, useA1 ? ladBufA1 : ladBufA2);
        }
        while (regT >= 0 && ciT > startMi) {
            useB1 = !useB1;
            ladderOverlayClimb(false, ciT, regT, valsT, nT, useB1 ? ladBufB1 : ladBufB2);
        }

        for (int32_t mi = startMi;; --mi) {
            LadderFrontier fs, ft;
            if (regS >= 0) {
                fs.region = regS;
                fs.n = nS;
                fs.vals = valsS;
            } else {
                fs.selfId = selfS;
                fs.selfSPos = selfSPosS;
            }
            if (regT >= 0) {
                ft.region = regT;
                ft.n = nT;
                ft.vals = valsT;
            } else {
                ft.selfId = selfT;
                ft.selfSPos = selfSPosT;
            }
            const int32_t d = ladderCliqueMeet(fs, ft, mi);
            if (d < best)
                best = d;
            if (mi == 0)
                break;

            // Climb both sides into cut mi-1; an above-cut endpoint enters with its seed row.
            if (regS < 0) {
                regS = R.initRegionOfHub[selfS];
                ciS = mi - 1;
                valsS = L.tSeedFwd.data() + L.tSeedPos[selfS];
                nS = static_cast<int32_t>(R.cuts[mi - 1].bSize(regS));
                selfS = -1;
            } else {
                useA1 = !useA1;
                ladderOverlayClimb(true, ciS, regS, valsS, nS, useA1 ? ladBufA1 : ladBufA2);
            }
            if (regT < 0) {
                regT = R.initRegionOfHub[selfT];
                ciT = mi - 1;
                valsT = (L.tSeedBwd.empty() ? L.tSeedFwd.data() : L.tSeedBwd.data()) +
                        L.tSeedPos[selfT];
                nT = static_cast<int32_t>(R.cuts[mi - 1].bSize(regT));
                selfT = -1;
            } else {
                useB1 = !useB1;
                ladderOverlayClimb(false, ciT, regT, valsT, nT, useB1 ? ladBufB1 : ladBufB2);
            }
            // Post-climb cutoff: every remaining meet happens at cut mi-1 or above, so it
            // costs at least the two climbed frontiers' minima. Admissible because climbs
            // only add non-negative weights (the minima never decrease).
            if (overlayEarlyStop) {
                const int32_t mnS = regS >= 0 ? frontierMin(valsS, nS) : 0;
                const int32_t mnT = regT >= 0 ? frontierMin(valsT, nT) : 0;
                if (mnS >= FHL_INFTY || mnT >= FHL_INFTY || mnS + mnT >= best) {
                    ++earlyStops;
                    ++stopAtLevel[std::min<size_t>(stopAtLevel.size() - 1, mi)];
                    break;
                }
            }
        }
        return best;
    }


    const HubHierarchy &hierarchy;
    const FHLData &data;
    int32_t lastDistance = FHL_INFTY;
    bool lastModeIsLocal = true;

    // Precomputed per hub index / per table (see constructor).

    // Ladder-scheme structures and scratch.
    const fhl::Regions *ladderRegions = nullptr;
    std::vector<int32_t> ciMeetOfLevel; // LCA level -> meet cut index
    std::vector<int32_t> ladBufA1, ladBufA2, ladBufB1, ladBufB2;
    int32_t ladderZero = 0;

    // Scratch space for the middle-node query.
    std::vector<int32_t> middleIds;
    std::vector<int32_t> distToMiddle;
    std::vector<int32_t> distFromMiddle;

    EliminationTreeQuery<LabelSet, false> localQuery;
};
