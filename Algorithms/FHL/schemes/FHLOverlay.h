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
#include "Algorithms/FHL/schemes/FHLInRegionHubs.h"
#include "Algorithms/FHL/schemes/FHLSeedRectangle.h"

namespace fhl {

// SCHEME 4 -- overlay. CRP-style boundary-to-boundary cliques per cut level, each computed under
// the metric RESTRICTED to its own region, built bottom-up from original edges only. A
// query climbs cut by cut and meets at every common region: much slower than the schemes
// above, but a metric change confined to one region rebuilds only that region's clique and
// its root chain, which is what makes partial re-customization possible.
struct Overlay {
    // -------------------------------------------------------------------- overlay mode build
    // Region-restricted semantics built from the original edges only: leaf cliques by
    // Dijkstra inside each leaf region, meet cliques bottom-up on small H-graphs (child
    // cliques + original edges among S), seeds by elimination-path sweeps restricted to the
    // own region (ascending CCH paths stay inside the region and their shortcut unpackings
    // stay below their tails, so the sweep is exact for the restricted metric).
    template<typename CchT>
    static void buildOverlay(const Regions &R, const HubHierarchy &hierarchy,
                             const CchT &cch, const int32_t *const upW,
                             const int32_t *const downW,
                             const std::vector<int32_t> &localElimTree,
                             const int32_t *const inputWeights, const bool symmetric,
                             Values &out) {
        KASSERT(R.overlay && !R.csrFirst.empty());
        out.mode = 2;
        out.symmetric = symmetric;
        Timer ot;
        buildOverlaySeeds(R, hierarchy, cch, upW, downW, localElimTree, symmetric,
                          /*regionList=*/nullptr, 0, out);
        const auto tSeeds = ot.elapsed<std::chrono::milliseconds>();
        ot.restart();
        buildOverlayLeafCliques(R, hierarchy, inputWeights, out); // composes the seed rows
        const auto tLeaf = ot.elapsed<std::chrono::milliseconds>();
        ot.restart();
        buildOverlayMeets(R, hierarchy, inputWeights, /*onlyRegion=*/-1, out);
        std::cout << "[cust] overlay: seeds " << tSeeds << " ms, leaf cliques " << tLeaf
                  << " ms, meet cliques " << ot.elapsed<std::chrono::milliseconds>()
                  << " ms\n";
        const auto mb = [](const uint64_t bytes) { return bytes / 1048576.0; };
        uint64_t meetBytes = 0;
        for (const auto &mc : out.meetClique)
            meetBytes += mc.M.size() * 4ull + mc.off.size() * 8ull;
        std::cout << "Overlay storage: leaf cliques " << mb(out.leafClique.M.size() * 4ull)
                  << " MB, meet cliques " << mb(meetBytes) << " MB, vertex seeds "
                  << mb(out.seedFwd.size() * 4ull + out.seedBwd.size() * 4ull +
                        out.seedPos.size() * 8ull)
                  << " MB, hub seeds "
                  << mb(out.tSeedFwd.size() * 4ull + out.tSeedBwd.size() * 4ull +
                        out.tSeedPos.size() * 8ull)
                  << " MB" << std::endl;
    }


    // Incremental flat (fast path): after a metric change whose CCH propagation stayed below
    // the hub zone, only the region's seeds (and top labels) can have changed — every
    // other flat structure is a function of hub-zone shortcuts alone. The new global
    // seeds are composed exactly as  d(v,u) = min_w d_R(v,w) + d(w,u):  a region-restricted
    // elimination sweep (no access nodes, no oracle tables needed) chained with the leaf
    // clique's global boundary-to-boundary distances (first-exit decomposition).
    template<typename CchT>
    static void refreshRegionAfterInteriorChange(const Regions &R,
                                                 const HubHierarchy &hierarchy,
                                                 const CchT &cch, const int32_t *const upW,
                                                 const int32_t *const downW,
                                                 const std::vector<int32_t> &localElimTree,
                                                 const int32_t r, Values &out) {
        refreshRegionsFromCliques(R, hierarchy, cch, upW, downW, localElimTree, &r, 1, out);
    }


    // Seeds (restricted sweep composed with the region's global boundary-to-boundary
    // clique) and top
    // labels for a list of leaf regions. The cliques and CSR must already be current.
    template<typename CchT>
    static void refreshRegionsFromCliques(const Regions &R,
                                          const HubHierarchy &hierarchy,
                                          const CchT &cch, const int32_t *const upW,
                                          const int32_t *const downW,
                                          const std::vector<int32_t> &localElimTree,
                                          const int32_t *const regionList,
                                          const int32_t numRegionsInList, Values &out) {
        const bool symmetric = out.seedBwd.empty();
        buildOverlaySeeds(R, hierarchy, cch, upW, downW, localElimTree, symmetric,
                          regionList, numRegionsInList, out);
#pragma omp parallel for schedule(dynamic, 4)
        for (int32_t i = 0; i < numRegionsInList; ++i)
            composeRegionSeeds(R, regionList[i], symmetric, out);
        // In-region rows depend only on this region's own weights.
        InRegionHubs::buildLocalRows(R, hierarchy, cch, upW, downW, symmetric, regionList,
                       numRegionsInList, out);
    }


    static void composeRegionSeeds(const Regions &R, const int32_t r, const bool symmetric,
                                   Values &out) {
        const auto &leaf = R.leaf();
        const int32_t n = static_cast<int32_t>(leaf.bSize(r));
        if (n == 0)
            return;
        // Compose with the region's global leaf clique (first-exit decomposition).
        const int32_t *const M = out.leafClique.M.data() + out.leafClique.off[r];
        // The forward fold needs column j of M (d(boundary_w -> boundary_j) over w). M is stored by
        // rows, so take its transpose: for a symmetric metric M IS its own transpose, so row
        // j serves; otherwise transpose once per region (n is the boundary-point count, so
        // n^2 is tiny).
        std::vector<int32_t> mt;
        const int32_t *mCol = M;
        if (!symmetric) {
            mt.resize(static_cast<size_t>(n) * n);
            for (int32_t a = 0; a < n; ++a)
                for (int32_t b = 0; b < n; ++b)
                    mt[static_cast<int64_t>(a) * n + b] = M[static_cast<int64_t>(b) * n + a];
            mCol = mt.data();
        }
        std::vector<int32_t> buf(n);
        for (int64_t k = R.leafVertFirst[r]; k < R.leafVertFirst[r + 1]; ++k) {
            const int32_t v = R.leafVertIds[k];
            int32_t *const fw = out.seedFwd.data() + out.seedPos[v];
            std::copy(fw, fw + n, buf.begin());
            for (int32_t j = 0; j < n; ++j) {
                const int32_t *const col = mCol + static_cast<int64_t>(j) * n;
                int32_t best = FHL_INFTY;
                for (int32_t w = 0; w < n; ++w) {
                    if (buf[w] >= FHL_INFTY)
                        continue;
                    const int32_t d = buf[w] + col[w];
                    if (d < best)
                        best = d;
                }
                fw[j] = best > FHL_INFTY ? FHL_INFTY : best;
            }
            if (symmetric)
                continue;
            int32_t *const bw = out.seedBwd.data() + out.seedPos[v];
            std::copy(bw, bw + n, buf.begin());
            for (int32_t j = 0; j < n; ++j) {
                int32_t best = FHL_INFTY;
                for (int32_t w = 0; w < n; ++w) {
                    if (buf[w] >= FHL_INFTY)
                        continue;
                    const int32_t d = M[static_cast<int64_t>(j) * n + w] + buf[w];
                    if (d < best)
                        best = d;
                }
                bw[j] = best > FHL_INFTY ? FHL_INFTY : best;
            }
        }
        // Refresh the affected vertices' top labels from the (unchanged) CSR.
        if (out.topLabLevels > 0) {
            const int64_t base = R.flatColFirst[r];
            const auto fuse = [&](const Values::FlatCsr &csr, const int32_t *const seeds,
                                  int32_t *const lab, const int64_t nLab) {
                const int64_t *const ce = csr.endArray();
                for (int64_t j = 0; j < nLab; ++j) {
                    int32_t best = FHL_INFTY;
                    for (int64_t kk = csr.ptr[base + j]; kk < ce[base + j]; ++kk) {
                        const int32_t d = seeds[csr.row[kk]] + csr.val[kk];
                        if (d < best)
                            best = d;
                    }
                    lab[j] = best > FHL_INFTY ? FHL_INFTY : best;
                }
            };
            for (int64_t k = R.leafVertFirst[r]; k < R.leafVertFirst[r + 1]; ++k) {
                const int32_t v = R.leafVertIds[k];
                const int64_t nLab = out.topLabPos[v + 1] - out.topLabPos[v];
                if (nLab == 0)
                    continue;
                fuse(out.csrOut, out.seedFwd.data() + out.seedPos[v],
                     out.topLabFwd.data() + out.topLabPos[v], nLab);
                if (!symmetric)
                    fuse(out.csrIn, out.seedBwd.data() + out.seedPos[v],
                         out.topLabBwd.data() + out.topLabPos[v], nLab);
            }
        }
    }


    // After a metric change confined to one leaf region, rebuild that leaf clique, the meet
    // cliques on its root chain, and the seeds of its vertices.
    template<typename CchT>
    static void recustomizeLeafRegion(const Regions &R, const HubHierarchy &hierarchy,
                                      const CchT &cch, const int32_t *const upW,
                                      const int32_t *const downW,
                                      const std::vector<int32_t> &localElimTree,
                                      const int32_t *const inputWeights,
                                      const int32_t leafRegion, Values &out) {
        const int32_t m = R.numCuts();
        // Seeds first: the leaf clique is now composed from them.
        buildOverlaySeeds(R, hierarchy, cch, upW, downW, localElimTree, out.seedBwd.empty(),
                          &leafRegion, 1, out);
        buildOverlayLeafClique(R, hierarchy, inputWeights, leafRegion, out);
        int32_t mr = R.cuts[m - 1].parentRegion[leafRegion];
        for (int32_t mi = m - 1; mi >= 0; --mi) {
            buildOverlayMeets(R, hierarchy, inputWeights, mr, out, mi);
            if (mi > 0)
                mr = R.cuts[mi - 1].parentRegion[mr];
        }
    }


    // Dijkstra rows of one leaf clique: restricted distances between boundary members with
    // all intermediates inside the region (local indices via the sorted vertex list).
    // Boundary-to-boundary restricted distances composed from the seed rows instead of one Dijkstra
    // per boundary point: any path from boundary point i through the region has a FIRST
    // interior vertex v, and
    // the rest of it is exactly seed_v (the restricted distance from v to every boundary point). So
    //     M[i][j] = min( direct edge i->j , min over v in N(i) inside R of w(i,v) + seed_v[j] )
    // — the entry-side mirror of the first-exit decomposition the seeds already encode.
    // Requires the seeds of this region to be current (buildOverlay runs seeds first).
    static void buildOverlayLeafClique(const Regions &R, const HubHierarchy &hierarchy,
                                       const int32_t *const inputWeights, const int32_t r,
                                       Values &out) {
        const auto &c = R.leaf();
        const int32_t n = static_cast<int32_t>(c.bSize(r));
        if (n == 0)
            return;
        const HubId *const ids = c.bIds.data() + c.bFirst[r];
        int32_t *const M = out.leafClique.M.data() + out.leafClique.off[r];
        std::fill(M, M + static_cast<int64_t>(n) * n, FHL_INFTY);
        for (int32_t i = 0; i < n; ++i) {
            const int32_t src = hierarchy.rankOfHub(ids[i]);
            int32_t *const row = M + static_cast<int64_t>(i) * n;
            for (int64_t e = R.csrFirst[src]; e < R.csrFirst[src + 1]; ++e) {
                const int32_t h = R.csrHead[e];
                const int32_t w = inputWeights[R.csrEdge[e]];
                if (R.leafRegionOfVertex[h] == r) { // enters the region here
                    const int32_t *const seed = out.seedFwd.data() + out.seedPos[h];
                    for (int32_t j = 0; j < n; ++j) {
                        const int32_t d = w + seed[j];
                        if (seed[j] < FHL_INFTY && d < row[j])
                            row[j] = d;
                    }
                    continue;
                }
                // Another boundary point of the same region: the direct edge is a valid connection.
                if (hierarchy.getVertexLevel(h) >= c.depth)
                    continue;
                const HubId hid = hierarchy.hubIdOfRank(h);
                const auto it = std::lower_bound(ids, ids + n, hid);
                if (it != ids + n && *it == hid && w < row[it - ids])
                    row[it - ids] = w;
            }
            row[i] = 0;
        }
    }


    static void buildOverlayLeafCliques(const Regions &R,
                                        const HubHierarchy &hierarchy,
                                        const int32_t *const inputWeights, Values &out) {
        const auto &c = R.leaf();
        out.leafClique.off.assign(c.numRegions() + 1, 0);
        for (int32_t r = 0; r < c.numRegions(); ++r)
            out.leafClique.off[r + 1] = out.leafClique.off[r] + c.bSize(r) * c.bSize(r);
        out.leafClique.M.resize(out.leafClique.off[c.numRegions()]);
#pragma omp parallel for schedule(dynamic, 4)
        for (int32_t r = 0; r < c.numRegions(); ++r)
            buildOverlayLeafClique(R, hierarchy, inputWeights, r, out);
    }


    // Bottom-up meet cliques. onlyRegion >= 0 (with onlyMi) rebuilds a single region block.
    static void buildOverlayMeets(const Regions &R, const HubHierarchy &hierarchy,
                                  const int32_t *const inputWeights, const int32_t onlyRegion,
                                  Values &out, const int32_t onlyMi = -1) {
        const int32_t m = R.numCuts();
        if (onlyRegion < 0) {
            out.meetClique.assign(m, {});
            for (int32_t mi = 0; mi < m; ++mi) {
                const auto &ms = R.meets[mi];
                auto &cl = out.meetClique[mi];
                cl.off.assign(ms.numRegions() + 1, 0);
                for (int32_t r = 0; r < ms.numRegions(); ++r)
                    cl.off[r + 1] = cl.off[r] + ms.sSize(r) * ms.sSize(r);
                cl.M.resize(cl.off[ms.numRegions()]);
            }
        }
        for (int32_t mi = m - 1; mi >= 0; --mi) {
            if (onlyMi >= 0 && mi != onlyMi)
                continue;
            const auto &ms = R.meets[mi];
#pragma omp parallel for schedule(dynamic, 2)
            for (int32_t r = 0; r < ms.numRegions(); ++r) {
                if (onlyRegion >= 0 && r != onlyRegion)
                    continue;
                buildOverlayMeetRegion(R, hierarchy, inputWeights, mi, r, out);
            }
        }
    }


    static void buildOverlayMeetRegion(const Regions &R, const HubHierarchy &hierarchy,
                                       const int32_t *const inputWeights, const int32_t mi,
                                       const int32_t r, Values &out) {
        const int32_t m = R.numCuts();
        const auto &ms = R.meets[mi];
        const int32_t n = static_cast<int32_t>(ms.sSize(r));
        if (n == 0)
            return;
        const HubId *const ids = ms.sIds.data() + ms.sFirst[r];
        int32_t *const M = out.meetClique[mi].M.data() + out.meetClique[mi].off[r];

        // H-graph over the S list: original edges among S plus child clique entries.
        std::vector<std::vector<std::pair<int32_t, int32_t>>> adj(n); // (target idx, weight)
        // The S list is two sorted runs (parent boundary points, then the band), so lookups go
        // through the meet level's segmented search.
        const auto idxOf = [&](const HubId id) -> int32_t { return ms.posOf(r, id); };
        for (int32_t i = 0; i < n; ++i) {
            const int32_t v = hierarchy.rankOfHub(ids[i]);
            for (int64_t e = R.csrFirst[v]; e < R.csrFirst[v + 1]; ++e) {
                const int32_t h = R.csrHead[e];
                if (hierarchy.getVertexLevel(h) >= R.thresh)
                    continue;
                const int32_t j = idxOf(hierarchy.hubIdOfRank(h));
                if (j >= 0)
                    adj[i].emplace_back(j, inputWeights[R.csrEdge[e]]);
            }
        }
        // Children: cut regions at cut mi (their B lists live inside this S list via
        // sPosOfB; their cliques are the next meet level for mi < m-1, the leaf cliques for
        // mi == m-1).
        const auto &cc = R.cuts[mi];
        for (int32_t crr = 0; crr < cc.numRegions(); ++crr) {
            if (cc.parentRegion[crr] != (mi == 0 ? 0 : r))
                continue;
            const int32_t bn = static_cast<int32_t>(cc.bSize(crr));
            if (bn == 0)
                continue;
            const int64_t bOff = cc.bFirst[crr];
            const int32_t *childM;
            int32_t childN;
            const int32_t *bPosInChild;
            std::vector<int32_t> scratch;
            if (mi == m - 1) { // child = leaf clique over exactly the B list
                childM = out.leafClique.M.data() + out.leafClique.off[crr];
                childN = bn;
                scratch.resize(bn);
                std::iota(scratch.begin(), scratch.end(), 0);
                bPosInChild = scratch.data();
            } else { // child = meet region at level mi+1 (indexed like cut region crr)
                childM = out.meetClique[mi + 1].M.data() + out.meetClique[mi + 1].off[crr];
                childN = static_cast<int32_t>(R.meets[mi + 1].sSize(crr));
                bPosInChild = R.meets[mi + 1].bPos.data() + bOff;
            }
            for (int32_t a = 0; a < bn; ++a) {
                const int32_t ia = cc.sPosOfB[bOff + a];
                for (int32_t b = 0; b < bn; ++b) {
                    if (a == b)
                        continue;
                    const int32_t w = childM[static_cast<int64_t>(bPosInChild[a]) * childN +
                                             bPosInChild[b]];
                    if (w < FHL_INFTY)
                        adj[ia].emplace_back(cc.sPosOfB[bOff + b], w);
                }
            }
        }

        // Dijkstra from every S member over the H-graph.
        std::vector<int32_t> dist(n);
        using QE = std::pair<int32_t, int32_t>;
        for (int32_t s = 0; s < n; ++s) {
            std::fill(dist.begin(), dist.end(), FHL_INFTY);
            std::priority_queue<QE, std::vector<QE>, std::greater<QE>> pq;
            dist[s] = 0;
            pq.emplace(0, s);
            while (!pq.empty()) {
                const auto [d, v] = pq.top();
                pq.pop();
                if (d > dist[v])
                    continue;
                for (const auto &[h, w] : adj[v])
                    if (d + w < dist[h]) {
                        dist[h] = d + w;
                        pq.emplace(dist[h], h);
                    }
            }
            for (int32_t j = 0; j < n; ++j)
                M[static_cast<int64_t>(s) * n + j] = dist[j];
        }
    }


    template<typename CchT>
    static void buildOverlaySeeds(const Regions &R, const HubHierarchy &hierarchy,
                                  const CchT &cch, const int32_t *const upW,
                                  const int32_t *const downW,
                                  const std::vector<int32_t> &localElimTree
                                          [[maybe_unused]],
                                  const bool symmetric, const int32_t *const regionList,
                                  const int32_t numRegionsInList, Values &out) {
        const auto &leaf = R.leaf();
        const auto &graph = cch.getUpwardGraph();

        if (regionList == nullptr) {
            out.seedPos.assign(R.numVertices + 1, 0);
            for (int32_t v = 0; v < R.numVertices; ++v) {
                const int32_t r = R.leafRegionOfVertex[v];
                out.seedPos[v + 1] = out.seedPos[v] + (r < 0 ? 0 : leaf.bSize(r));
            }
            out.seedFwd.assignFill(out.seedPos[R.numVertices], FHL_INFTY);
            if (symmetric)
                out.seedBwd.clear();
            else
                out.seedBwd.assignFill(out.seedPos[R.numVertices], FHL_INFTY);
        }

        const auto sweep = [&](const int32_t start, const auto &parent,
                               const int32_t minLevel, const HubId *const ids,
                               const int32_t n, int32_t *const fw, int32_t *const bw,
                               std::vector<int32_t> &path, std::vector<int32_t> &pUp,
                               std::vector<int32_t> &pDown) {
            std::fill(fw, fw + n, FHL_INFTY);
            if (bw)
                std::fill(bw, bw + n, FHL_INFTY);
            if (n == 0)
                return;
            path.clear();
            for (int32_t x = start;
                 x != INVALID_VERTEX && hierarchy.getVertexLevel(x) >= minLevel;
                 x = parent[x])
                path.push_back(x); // ascending ranks, all inside the region
            pUp.assign(path.size(), FHL_INFTY);
            pDown.assign(path.size(), FHL_INFTY);
            pUp[0] = pDown[0] = 0;
            for (size_t k = 0; k < path.size(); ++k) {
                const int32_t x = path[k];
                const int32_t du = pUp[k];
                const int32_t dd = pDown[k];
                if (du >= FHL_INFTY && dd >= FHL_INFTY)
                    continue;
                for (int e = graph.firstEdge(x); e < graph.lastEdge(x); ++e) {
                    const int32_t h = graph.edgeHead(e);
                    if (hierarchy.getVertexLevel(h) < minLevel) { // boundary hop
                        const HubId hid = hierarchy.hubIdOfRank(h);
                        const auto it = std::lower_bound(ids, ids + n, hid);
                        if (it == ids + n || *it != hid)
                            continue;
                        const auto j = it - ids;
                        if (du + upW[e] < fw[j])
                            fw[j] = du + upW[e];
                        if (bw && dd + downW[e] < bw[j])
                            bw[j] = dd + downW[e];
                        continue;
                    }
                    const auto pit = std::lower_bound(path.begin() + k + 1, path.end(), h);
                    if (pit == path.end() || *pit != h)
                        continue;
                    const auto pk = pit - path.begin();
                    if (du + upW[e] < pUp[pk])
                        pUp[pk] = du + upW[e];
                    if (dd + downW[e] < pDown[pk])
                        pDown[pk] = dd + downW[e];
                }
            }
        };

        // Region-restricted seeds by a descending-rank dynamic program instead of one
        // elimination-path sweep per vertex. Every in-region up-neighbour of v is an
        // elimination ancestor of v with a higher rank, so its own seed row is already
        // final when v is processed:
        //     row_v[j] = min( w(v,x) for a boundary point x = j , min over in-region (v,x) of
        //                     w(v,x) + row_x[j] )
        // Same "first exit" semantics as the sweep, but the region's shared upper paths are
        // walked once instead of once per vertex, and each row is a contiguous |B|-int run.
        const auto sweepRegion = [&](const int32_t r, std::vector<int32_t> &,
                                     std::vector<int32_t> &, std::vector<int32_t> &) {
            const int32_t n = static_cast<int32_t>(leaf.bSize(r));
            const HubId *const ids = leaf.bIds.data() + leaf.bFirst[r];
            const int32_t minLevel = R.regionDepth[r];
            for (int64_t k = R.leafVertFirst[r + 1] - 1; k >= R.leafVertFirst[r]; --k) {
                const int32_t v = R.leafVertIds[k];
                int32_t *const fw = out.seedFwd.data() + out.seedPos[v];
                int32_t *const bw =
                        symmetric ? nullptr : out.seedBwd.data() + out.seedPos[v];
                std::fill(fw, fw + n, FHL_INFTY);
                if (bw)
                    std::fill(bw, bw + n, FHL_INFTY);
                for (int e = graph.firstEdge(v); e < graph.lastEdge(v); ++e) {
                    const int32_t h = graph.edgeHead(e);
                    if (hierarchy.getVertexLevel(h) < minLevel) { // exits here
                        const HubId hid = hierarchy.hubIdOfRank(h);
                        const auto it = std::lower_bound(ids, ids + n, hid);
                        if (it == ids + n || *it != hid)
                            continue;
                        const auto j = it - ids;
                        if (upW[e] < fw[j])
                            fw[j] = upW[e];
                        if (bw && downW[e] < bw[j])
                            bw[j] = downW[e];
                        continue;
                    }
                    // In-region ancestor: fold its finished row.
                    const int32_t *const rowF = out.seedFwd.data() + out.seedPos[h];
                    const int32_t wu = upW[e];
                    for (int32_t j = 0; j < n; ++j) {
                        const int32_t d = wu + rowF[j];
                        if (rowF[j] < FHL_INFTY && d < fw[j])
                            fw[j] = d;
                    }
                    if (bw) {
                        const int32_t *const rowB = out.seedBwd.data() + out.seedPos[h];
                        const int32_t wd = downW[e];
                        for (int32_t j = 0; j < n; ++j) {
                            const int32_t d = wd + rowB[j];
                            if (rowB[j] < FHL_INFTY && d < bw[j])
                                bw[j] = d;
                        }
                    }
                }
            }
        };
#pragma omp parallel
        {
            std::vector<int32_t> path, pUp, pDown;
            if (regionList != nullptr) {
#pragma omp for schedule(dynamic, 4)
                for (int32_t i = 0; i < numRegionsInList; ++i)
                    sweepRegion(regionList[i], path, pUp, pDown);
            } else {
#pragma omp for schedule(dynamic, 8)
                for (int32_t r = 0; r < leaf.numRegions(); ++r)
                    sweepRegion(r, path, pUp, pDown);
            }
        }

        if (regionList != nullptr)
            return;
        const auto &elim = cch.getEliminationTree();
        out.tSeedPos.assign(R.numHubs + 1, 0);
        for (int32_t x = 0; x < R.numHubs; ++x) {
            const int32_t ci = R.initCutOfHub[x];
            out.tSeedPos[x + 1] = out.tSeedPos[x] +
                    (ci < 0 ? 0 : R.cuts[ci].bSize(R.initRegionOfHub[x]));
        }
        out.tSeedFwd.assign(out.tSeedPos[R.numHubs], FHL_INFTY);
        out.tSeedBwd.assign(symmetric ? 0 : out.tSeedPos[R.numHubs], FHL_INFTY);
#pragma omp parallel
        {
            std::vector<int32_t> path, pUp, pDown;
#pragma omp for schedule(dynamic, 256)
            for (int32_t x = 0; x < R.numHubs; ++x) {
                const int32_t ci = R.initCutOfHub[x];
                if (ci < 0)
                    continue;
                const auto &c = R.cuts[ci];
                const int32_t r = R.initRegionOfHub[x];
                sweep(hierarchy.rankOfHub(x), elim, c.depth,
                      c.bIds.data() + c.bFirst[r], static_cast<int32_t>(c.bSize(r)),
                      out.tSeedFwd.data() + out.tSeedPos[x],
                      symmetric ? nullptr : out.tSeedBwd.data() + out.tSeedPos[x], path, pUp,
                      pDown);
            }
        }
    }
};

} // namespace fhl
