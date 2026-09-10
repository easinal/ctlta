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

// METRIC-INDEPENDENT skeleton shared by all four FHL schemes: the cut ladder, the leaf
// regions and their boundary points, and the per-region column runs (which hubs a region's rows are
// aligned to). Built once from the separator decomposition; a metric change never touches it.
namespace fhl {

using HubId = fhl::HubId;

struct Regions {
    int32_t thresh = 0;
    int32_t numHubs = 0;
    int32_t numVertices = 0;
    bool overlay = false;
    std::vector<int32_t> cutDepths; // ascending, back() == thresh ({thresh} in flat mode)

    // Regions at cut depth cutDepths[ci]; boundary lists sorted by hub id. Flat mode
    // builds only the leaf level; overlay mode builds the whole ladder.
    struct CutLevel {
        int32_t depth = 0;
        std::vector<uint64_t> prefixOfRegion;
        std::vector<int32_t> regionOfPrefix; // size 1 << depth, -1 if absent
        std::vector<int64_t> bFirst;         // CSR region -> [bFirst[r], bFirst[r+1])
        std::vector<HubId> bIds;
        std::vector<int32_t> parentRegion;   // region index at cut ci-1 (0 for ci == 0:
                                             // the single root meet region)
        std::vector<int32_t> sPosOfB;        // parallel to bIds: position within the S
                                             // list of the meet region meets[ci]
        int32_t numRegions() const { return static_cast<int32_t>(prefixOfRegion.size()); }
        int64_t bSize(const int32_t r) const { return bFirst[r + 1] - bFirst[r]; }
    };
    std::vector<CutLevel> cuts;

    const CutLevel &leaf() const { return cuts.back(); }

    // Overlay meet level mi sits ABOVE cut mi: depth (mi == 0 ? 0 : cutDepths[mi-1]); its
    // separator levels reach down to cutDepths[mi]. For mi >= 1 the meet regions coincide
    // 1:1 with the cut regions at ci = mi-1 and share their indexing; mi == 0 is the
    // single root region. S = own separator levels ∪ B.
    struct MeetLevel {
        int32_t depth = 0;
        int32_t sepEnd = 0;
        std::vector<int64_t> sFirst;
        std::vector<HubId> sIds;
        std::vector<int32_t> bPos; // aligned to cuts[mi-1].bIds: position of each
                                   // boundary vertex within its meet region's S list
        std::vector<int32_t> numParentBoundary; // per meet region: sorted prefix length
        int32_t numRegions() const { return static_cast<int32_t>(sFirst.size()) - 1; }
        // Position of hub id within meet region r's S list (-1 if absent). The list
        // is two sorted runs: the parent boundary points, then the band's own vertices.
        int32_t posOf(const int32_t r, const HubId id) const {
            const auto beg = sIds.begin() + sFirst[r];
            const auto mid = beg + numParentBoundary[r];
            const auto end = sIds.begin() + sFirst[r + 1];
            auto it = std::lower_bound(beg, mid, id);
            if (it != mid && *it == id)
                return static_cast<int32_t>(it - beg);
            it = std::lower_bound(mid, end, id);
            return (it != end && *it == id) ? static_cast<int32_t>(it - beg) : -1;
        }
        int64_t sSize(const int32_t r) const { return sFirst[r + 1] - sFirst[r]; }
    };
    std::vector<MeetLevel> meets; // overlay only

    // Per-region cut metadata (uniform mode: constant thresh / the prefix table;
    // size-cut mode: per-region root depth and prefix). levelStride is the row stride
    // of flatLevelEnd/tLevelEnd (= thresh in uniform mode, max cut depth + 1 in size
    // mode); far-pair LCA levels are always below both regions' depths.
    int32_t levelStride = 0;
    std::vector<int32_t> regionDepth;
    std::vector<uint64_t> regionPrefix;

    // Optional in-region columns: the top `locLevels` levels of each leaf region's OWN
    // separator hierarchy. A same-region pair whose LCA falls in those levels can be
    // answered by a fold over these columns (the in-region middle-node argument), so the
    // CCH fallback is skipped entirely.
    int32_t locLevels = 0;             // != 0 marks "in-region columns enabled"
    std::vector<int64_t> locColFirst;  // per leaf region
    std::vector<int32_t> locColIds;    // ranks (ascending within a region)
    std::vector<int32_t> locLevelEnd;  // unused with the size-based inner cut
    std::vector<int32_t> innerRegionOfVertex; // by rank; -1 = the region's local hubs

    std::vector<int32_t> leafRegionOfVertex; // by rank; -1 for hubs
    std::vector<int64_t> leafVertFirst;      // CSR: leaf region -> its vertices (ranks,
    std::vector<int32_t> leafVertIds;        // ascending); overlay Dijkstra indexing

    // Flat columns: per leaf region, all hub path ancestors of its branch, ordered by
    // (level asc, id asc). Two leaf regions agree on every column with level <= their LCA
    // level, at the same positions — the flat combine exploits this. Hubs get
    // the same lists along their own prefix chain (including their own separator).
    std::vector<int64_t> flatColFirst; // by leaf region
    std::vector<int32_t> flatLevelEnd; // [r * thresh + l] = #cols of r with level <= l
    // A column list is a concatenation of contiguous per-level bucket runs, so the ids
    // are never materialized: bucketIds holds every hub exactly once, grouped
    // by level and sorted by side prefix within a level, and each (owner, level) pair
    // records only where its run starts -- the run's length is the levelEnd difference.
    // Materializing the ids instead costs two orders of magnitude more memory.
    std::vector<HubId> bucketIds; // numHubs, grouped by level
    std::vector<int32_t> flatColStart;    // [r * thresh + l] -> index into bucketIds
    // Per boundary point (parallel to leaf().bIds): its index within its region's column
    // list. Both lists come from the same (level, prefix) bucket walk and a boundary point
    // lies on the region's root path, so its own column list is a POSITION-FOR-POSITION
    // prefix of the region's, so its rectangle row can also be located inside the
    // hub rows at a computed offset.
    std::vector<int32_t> boundaryColPos;
    std::vector<int64_t> tColFirst;    // by hub id
    std::vector<int32_t> tLevelEnd;    // [x * thresh + l]
    std::vector<int32_t> tColStart;    // [x * thresh + l] -> index into bucketIds

    // Overlay hub helpers: deepest cut containing the node, the region there, and
    // the node's position within the S list of its own meet region.
    std::vector<int32_t> initCutOfHub;
    std::vector<int32_t> initRegionOfHub;
    std::vector<int32_t> selfSPosOfHub;

    // Original-graph adjacency in rank space (heads + input edge ids); overlay only.
    std::vector<int64_t> csrFirst;
    std::vector<int32_t> csrHead;
    std::vector<int32_t> csrEdge;

    int32_t numCuts() const { return static_cast<int32_t>(cutDepths.size()); }

    // Walk one owner's column list in order without materializing it; f(localIndex, id).
    template<typename F>
    void forEachCol(const int32_t *const start, const int32_t *const lend, F &&f) const {
        int32_t prev = 0;
        for (int32_t l = 0; l < levelStride; ++l) {
            const int32_t end = lend[l];
            for (int32_t o = prev; o < end; ++o)
                f(o, bucketIds[start[l] + (o - prev)]);
            prev = end;
        }
    }
    // Same, into a caller-owned buffer, for the loops that index columns out of order.
    void gatherCols(const int32_t *const start, const int32_t *const lend,
                    std::vector<HubId> &buf) const {
        buf.resize(lend[levelStride - 1]);
        int32_t prev = 0;
        for (int32_t l = 0; l < levelStride; ++l) {
            const int32_t end = lend[l];
            if (end > prev)
                std::copy_n(bucketIds.data() + start[l], end - prev, buf.data() + prev);
            prev = end;
        }
    }
    const int32_t *flatColRow(const int32_t r) const {
        return flatColStart.data() + static_cast<int64_t>(r) * levelStride;
    }
    const int32_t *flatLevelRow(const int32_t r) const {
        return flatLevelEnd.data() + static_cast<int64_t>(r) * levelStride;
    }
    const int32_t *tColRow(const int32_t x) const {
        return tColStart.data() + static_cast<int64_t>(x) * levelStride;
    }
    const int32_t *tLevelRow(const int32_t x) const {
        return tLevelEnd.data() + static_cast<int64_t>(x) * levelStride;
    }
    void gatherFlatCols(const int32_t r, std::vector<HubId> &buf) const {
        gatherCols(flatColRow(r), flatLevelRow(r), buf);
    }
    // One column of a region without materializing the whole list.
    HubId flatColId(const int32_t r, const int32_t j) const {
        const int32_t *const start = flatColRow(r);
        const int32_t *const lend = flatLevelRow(r);
        int32_t prev = 0;
        for (int32_t l = 0; l < levelStride; ++l) {
            if (j < lend[l])
                return bucketIds[start[l] + (j - prev)];
            prev = lend[l];
        }
        return HubId(-1);
    }

    uint64_t sizeInBytes() const {
        uint64_t s = sizeof(Regions);
        for (const auto &c : cuts)
            s += c.prefixOfRegion.size() * 8 + c.regionOfPrefix.size() * 4 +
                 c.bFirst.size() * 8 + c.bIds.size() * 4 + c.parentRegion.size() * 4 +
                 c.sPosOfB.size() * 4;
        for (const auto &ms : meets)
            s += ms.sFirst.size() * 8 +
                 (ms.sIds.size() + ms.bPos.size() + ms.numParentBoundary.size()) * 4;
        s += (leafRegionOfVertex.size() + leafVertIds.size()) * 4 +
             leafVertFirst.size() * 8;
        s += flatColFirst.size() * 8 +
             (bucketIds.size() + flatColStart.size() + flatLevelEnd.size()) * 4;
        s += regionDepth.size() * 4 + regionPrefix.size() * 8 + boundaryColPos.size() * 4;
        s += locColFirst.size() * 8 +
             (locColIds.size() + locLevelEnd.size() + innerRegionOfVertex.size()) * 4;
        s += tColFirst.size() * 8 + (tColStart.size() + tLevelEnd.size()) * 4;
        s += (initCutOfHub.size() + initRegionOfHub.size() +
              selfSPosOfHub.size()) * 4;
        s += csrFirst.size() * 8 + (csrHead.size() + csrEdge.size()) * 4;
        return s;
    }

    template<typename InputGraphT, typename RanksT>
    void build(const HubHierarchy &hierarchy, const InputGraphT &graph,
               const RanksT &ranks, const bool overlayMode, const int32_t blockHeight,
               const std::vector<int32_t> &customCuts,
               const int32_t targetRegionSize = 0, const int32_t localLevels = 0) {
        locLevels = localLevels;
        thresh = hierarchy.getHubLevelThreshold();
        numHubs = hierarchy.numHubs();
        numVertices = static_cast<int32_t>(hierarchy.numVertices());
        overlay = overlayMode;
        if (targetRegionSize > 0) {
            KASSERT(!overlayMode && hierarchy.isSizeCut());
            buildSizeCut(hierarchy, graph, ranks);
            return;
        }
        if (thresh <= 0 || thresh > 25)
            throw std::invalid_argument("regions: hub level threshold out of range");

        cutDepths.clear();
        if (overlay) {
            if (!customCuts.empty()) {
                for (const auto d : customCuts)
                    if (d > 0 && d < thresh)
                        cutDepths.push_back(d);
                std::sort(cutDepths.begin(), cutDepths.end());
                cutDepths.erase(std::unique(cutDepths.begin(), cutDepths.end()),
                                cutDepths.end());
            } else {
                for (int32_t d = blockHeight; d < thresh; d += blockHeight)
                    cutDepths.push_back(d);
            }
        }
        cutDepths.push_back(thresh);
        const int32_t m = numCuts();

        // Region enumeration per cut depth.
        cuts.assign(m, {});
        for (int32_t ci = 0; ci < m; ++ci) {
            auto &c = cuts[ci];
            c.depth = cutDepths[ci];
            c.regionOfPrefix.assign(1ULL << c.depth, -1);
            const uint64_t mask = (1ULL << c.depth) - 1;
            for (int32_t v = 0; v < numVertices; ++v) {
                if (hierarchy.getVertexLevel(v) < c.depth)
                    continue;
                const uint64_t p = hierarchy.getPackedSideId(v) & mask;
                if (c.regionOfPrefix[p] < 0) {
                    c.regionOfPrefix[p] = static_cast<int32_t>(c.prefixOfRegion.size());
                    c.prefixOfRegion.push_back(p);
                }
            }
        }

        // Boundary discovery from the original edges: for an edge (a, b) with level(b) <
        // level(a), b bounds every region containing a whose depth d satisfies
        // level(b) < d <= level(a).
        std::vector<std::vector<std::vector<HubId>>> bSets(m);
        for (int32_t ci = 0; ci < m; ++ci)
            bSets[ci].resize(cuts[ci].numRegions());
        const auto addBoundary = [&](const int32_t deep, const int32_t shallow) {
            const int32_t ld = hierarchy.getVertexLevel(deep);
            const int32_t ls = hierarchy.getVertexLevel(shallow);
            for (int32_t ci = 0; ci < m; ++ci) {
                const int32_t d = cutDepths[ci];
                if (d <= ls)
                    continue;
                if (d > ld)
                    break;
                const uint64_t mask = (1ULL << d) - 1;
                const int32_t r =
                        cuts[ci].regionOfPrefix[hierarchy.getPackedSideId(deep) & mask];
                KASSERT(r >= 0);
                bSets[ci][r].push_back(hierarchy.hubIdOfRank(shallow));
            }
        };
        if (overlay) {
            csrFirst.assign(numVertices + 1, 0);
            csrHead.clear();
            csrEdge.clear();
        }
        {
            std::vector<int64_t> outDeg(overlay ? numVertices : 0, 0);
            for (int u = 0; u < graph.numVertices(); ++u)
                for (int e = graph.firstEdge(u); e < graph.lastEdge(u); ++e) {
                    const int32_t a = ranks[u];
                    const int32_t b = ranks[graph.edgeHead(e)];
                    const int32_t la = hierarchy.getVertexLevel(a);
                    const int32_t lb = hierarchy.getVertexLevel(b);
                    if (la > lb)
                        addBoundary(a, b);
                    else if (lb > la)
                        addBoundary(b, a);
                    if (overlay)
                        ++outDeg[a];
                }
            if (overlay) {
                for (int32_t v = 0; v < numVertices; ++v)
                    csrFirst[v + 1] = csrFirst[v] + outDeg[v];
                csrHead.resize(csrFirst[numVertices]);
                csrEdge.resize(csrFirst[numVertices]);
                std::vector<int64_t> fill(csrFirst.begin(), csrFirst.end() - 1);
                for (int u = 0; u < graph.numVertices(); ++u)
                    for (int e = graph.firstEdge(u); e < graph.lastEdge(u); ++e) {
                        const int32_t a = ranks[u];
                        csrHead[fill[a]] = ranks[graph.edgeHead(e)];
                        csrEdge[fill[a]] = e;
                        ++fill[a];
                    }
            }
        }
        for (int32_t ci = 0; ci < m; ++ci) {
            auto &c = cuts[ci];
            c.bFirst.assign(c.numRegions() + 1, 0);
            for (int32_t r = 0; r < c.numRegions(); ++r) {
                auto &set = bSets[ci][r];
                std::sort(set.begin(), set.end());
                set.erase(std::unique(set.begin(), set.end()), set.end());
                c.bFirst[r + 1] = c.bFirst[r] + static_cast<int64_t>(set.size());
            }
            c.bIds.resize(c.bFirst[c.numRegions()]);
            for (int32_t r = 0; r < c.numRegions(); ++r)
                std::copy(bSets[ci][r].begin(), bSets[ci][r].end(),
                          c.bIds.begin() + c.bFirst[r]);
            bSets[ci].clear();
            bSets[ci].shrink_to_fit();

            c.parentRegion.assign(c.numRegions(), 0);
            if (ci > 0) {
                const uint64_t pmask = (1ULL << cutDepths[ci - 1]) - 1;
                for (int32_t r = 0; r < c.numRegions(); ++r) {
                    c.parentRegion[r] =
                            cuts[ci - 1].regionOfPrefix[c.prefixOfRegion[r] & pmask];
                    KASSERT(c.parentRegion[r] >= 0);
                }
            }
        }

        buildVertexMaps(hierarchy);
        levelStride = thresh;
        regionDepth.assign(leaf().numRegions(), thresh);
        regionPrefix = leaf().prefixOfRegion;
        buildFlatColumns(hierarchy);
        buildLocalColumns(hierarchy);
        if (overlay)
            buildMeets(hierarchy);
    }

    // Size-cut build: leaf regions come from the hierarchy walk (variable depths);
    // boundaries via one hub-mask edge scan; columns from the shared bucket walk.
    template<typename InputGraphT, typename RanksT>
    void buildSizeCut(const HubHierarchy &hierarchy, const InputGraphT &graph,
                      const RanksT &ranks) {
        const int32_t numLeaf = hierarchy.numLeafRegions();
        cutDepths.assign(1, hierarchy.getMaxCutDepth() + 1);
        levelStride = hierarchy.getMaxCutDepth() + 1;

        // Vertex maps straight from the hierarchy.
        leafRegionOfVertex = hierarchy.getLeafRegionRaw();
        leafVertFirst.assign(numLeaf + 1, 0);
        for (int32_t v = 0; v < numVertices; ++v)
            if (leafRegionOfVertex[v] >= 0)
                ++leafVertFirst[leafRegionOfVertex[v] + 1];
        for (int32_t r = 0; r < numLeaf; ++r)
            leafVertFirst[r + 1] += leafVertFirst[r];
        leafVertIds.resize(leafVertFirst[numLeaf]);
        {
            std::vector<int64_t> fill(leafVertFirst.begin(), leafVertFirst.end() - 1);
            for (int32_t v = 0; v < numVertices; ++v)
                if (leafRegionOfVertex[v] >= 0)
                    leafVertIds[fill[leafRegionOfVertex[v]]++] = v;
        }
        // Region root depth = min member level; prefix = side bits up to that depth.
        regionDepth.assign(numLeaf, 127);
        regionPrefix.assign(numLeaf, 0);
        for (int32_t v = 0; v < numVertices; ++v) {
            const int32_t r = leafRegionOfVertex[v];
            if (r >= 0)
                regionDepth[r] = std::min(regionDepth[r],
                                          static_cast<int32_t>(
                                                  hierarchy.getVertexLevel(v)));
        }
        for (int32_t r = 0; r < numLeaf; ++r)
            KASSERT(regionDepth[r] < 127);
        for (int64_t k = 0; k < static_cast<int64_t>(leafVertIds.size()); ++k) {
            const int32_t v = leafVertIds[k];
            const int32_t r = leafRegionOfVertex[v];
            const int32_t d = regionDepth[r];
            regionPrefix[r] = hierarchy.getPackedSideId(v) &
                              (d > 0 ? (1ULL << d) - 1 : 0);
        }

        // Boundary discovery: an original edge from a region interior to a hub
        // vertex makes that hub a boundary point of the region.
        std::vector<std::vector<HubId>> bSets(numLeaf);
        for (int u = 0; u < graph.numVertices(); ++u)
            for (int e = graph.firstEdge(u); e < graph.lastEdge(u); ++e) {
                const int32_t a = ranks[u];
                const int32_t b = ranks[graph.edgeHead(e)];
                const bool ta = hierarchy.isHub(a);
                const bool tb = hierarchy.isHub(b);
                if (!ta && tb)
                    bSets[leafRegionOfVertex[a]].push_back(
                            hierarchy.hubIdOfRank(b));
                else if (ta && !tb)
                    bSets[leafRegionOfVertex[b]].push_back(
                            hierarchy.hubIdOfRank(a));
            }
        // Invariant check (cheap, once per preprocessing): a genuine separator
        // decomposition admits no edge between the interiors of two different leaf
        // regions — every out-of-region neighbour must be a hub on the
        // region's root path. A violation would silently drop a boundary point.
        {
            int64_t crossing = 0;
            for (int u = 0; u < graph.numVertices(); ++u)
                for (int e = graph.firstEdge(u); e < graph.lastEdge(u); ++e) {
                    const int32_t a = ranks[u], b = ranks[graph.edgeHead(e)];
                    if (!hierarchy.isHub(a) && !hierarchy.isHub(b) &&
                        leafRegionOfVertex[a] != leafRegionOfVertex[b])
                        ++crossing;
                }
            if (crossing > 0)
                std::cout << "FHL WARNING: " << crossing
                          << " edges join two different leaf regions (separator "
                             "decomposition is not a true vertex cut)" << std::endl;
        }

        cuts.assign(1, {});
        auto &c = cuts.back();
        c.depth = levelStride; // informational
        c.prefixOfRegion = regionPrefix;
        c.bFirst.assign(numLeaf + 1, 0);
        for (int32_t r = 0; r < numLeaf; ++r) {
            auto &set = bSets[r];
            std::sort(set.begin(), set.end());
            set.erase(std::unique(set.begin(), set.end()), set.end());
            c.bFirst[r + 1] = c.bFirst[r] + static_cast<int64_t>(set.size());
        }
        c.bIds.resize(c.bFirst[numLeaf]);
        for (int32_t r = 0; r < numLeaf; ++r)
            std::copy(bSets[r].begin(), bSets[r].end(), c.bIds.begin() + c.bFirst[r]);

        buildFlatColumns(hierarchy);
        buildLocalColumns(hierarchy);
    }

    // In-region columns = each leaf region's OWN local hub set: the vertices left
    // above a second, finer size cut inside the region (same rule as the outer cut, one
    // scale down). Two vertices of the region that sit in different sub-regions must
    // meet at one of these, so their query is an aligned fold; two vertices of the SAME
    // sub-region may meet deeper and still fall back.
    void buildLocalColumns(const HubHierarchy &hierarchy) {
        const int32_t numLeaf = leaf().numRegions();
        locColFirst.assign(numLeaf + 1, 0);
        locColIds.clear();
        innerRegionOfVertex.clear();
        if (!hierarchy.hasInnerCut()) {
            locLevelEnd.clear();
            return;
        }
        innerRegionOfVertex = hierarchy.getInnerRegionRawVec();
        for (int32_t r = 0; r < numLeaf; ++r) {
            for (int64_t k = leafVertFirst[r]; k < leafVertFirst[r + 1]; ++k) {
                const int32_t v = leafVertIds[k];
                if (innerRegionOfVertex[v] < 0) // above the inner cut
                    locColIds.push_back(v);
            }
            locColFirst[r + 1] = static_cast<int64_t>(locColIds.size());
        }
        locLevelEnd.clear();
    }

    void printStats(std::ostream &out) const {
        if (overlay) { // per-cut boundary statistics: how wide the climb frontiers get
            for (int32_t ci = 0; ci < numCuts(); ++ci) {
                const auto &c = cuts[ci];
                int64_t sum = 0, mx = 0;
                for (int32_t r = 0; r < c.numRegions(); ++r) {
                    sum += c.bSize(r);
                    mx = std::max<int64_t>(mx, c.bSize(r));
                }
                out << "cut depth " << c.depth << ": " << c.numRegions()
                    << " regions, |B| avg "
                    << (c.numRegions() ? sum / c.numRegions() : 0) << " max " << mx
                    << '\n';
            }
        }
        out << (overlay ? "Overlay cuts:" : "Flat regions, thresh") << ' ';
        for (const auto d : cutDepths)
            out << d << ' ';
        const auto &c = leaf();
        std::vector<int64_t> bs(c.numRegions());
        for (int32_t r = 0; r < c.numRegions(); ++r)
            bs[r] = c.bSize(r);
        std::sort(bs.begin(), bs.end());
        out << "| leaf regions=" << c.numRegions()
            << " |B| med=" << (bs.empty() ? 0 : bs[bs.size() / 2])
            << " max=" << (bs.empty() ? 0 : bs.back())
            << " flat cols=" << flatColFirst[c.numRegions()] << "\n";
    }

    void buildVertexMaps(const HubHierarchy &hierarchy) {
        leafRegionOfVertex.assign(numVertices, -1);
        const uint64_t lmask = (1ULL << thresh) - 1;
        for (int32_t v = 0; v < numVertices; ++v) {
            if (hierarchy.getVertexLevel(v) < thresh)
                continue;
            leafRegionOfVertex[v] =
                    leaf().regionOfPrefix[hierarchy.getPackedSideId(v) & lmask];
            KASSERT(leafRegionOfVertex[v] >= 0);
        }
        const int32_t numLeaf = leaf().numRegions();
        leafVertFirst.assign(numLeaf + 1, 0);
        for (int32_t v = 0; v < numVertices; ++v)
            if (leafRegionOfVertex[v] >= 0)
                ++leafVertFirst[leafRegionOfVertex[v] + 1];
        for (int32_t r = 0; r < numLeaf; ++r)
            leafVertFirst[r + 1] += leafVertFirst[r];
        leafVertIds.resize(leafVertFirst[numLeaf]);
        std::vector<int64_t> fill(leafVertFirst.begin(), leafVertFirst.end() - 1);
        for (int32_t v = 0; v < numVertices; ++v)
            if (leafRegionOfVertex[v] >= 0)
                leafVertIds[fill[leafRegionOfVertex[v]]++] = v;
    }

    // Column lists from (level, prefix) buckets, for leaf regions and hubs.
    void buildFlatColumns(const HubHierarchy &hierarchy) {
        std::vector<std::vector<std::pair<uint64_t, HubId>>> byLevel(levelStride);
        for (int32_t id = 0; id < numHubs; ++id) {
            const int32_t v = hierarchy.rankOfHub(id);
            const int32_t lv = hierarchy.getVertexLevel(v);
            const uint64_t mask = lv > 0 ? (1ULL << lv) - 1 : 0;
            byLevel[lv].emplace_back(hierarchy.getPackedSideId(v) & mask, id);
        }
        for (auto &bucket : byLevel)
            std::sort(bucket.begin(), bucket.end());
        // Flatten the buckets once; every column run is a slice of this array.
        std::vector<int32_t> bucketFirst(levelStride + 1, 0);
        for (int32_t l = 0; l < levelStride; ++l)
            bucketFirst[l + 1] = bucketFirst[l] + static_cast<int32_t>(byLevel[l].size());
        bucketIds.resize(bucketFirst[levelStride]);
        for (int32_t l = 0; l < levelStride; ++l)
            for (size_t i = 0; i < byLevel[l].size(); ++i)
                bucketIds[bucketFirst[l] + i] = byLevel[l][i].second;
        // (start index into bucketIds, run length) for one owner's level-l columns.
        const auto bucketRun = [&](const int32_t l, const uint64_t prefix) {
            const uint64_t lp = l > 0 ? (prefix & ((1ULL << l) - 1)) : 0;
            const auto &bucket = byLevel[l];
            const auto lo = std::lower_bound(bucket.begin(), bucket.end(),
                                             std::make_pair(lp, HubId(-1)));
            auto hi = lo;
            while (hi != bucket.end() && hi->first == lp)
                ++hi;
            return std::make_pair(
                    bucketFirst[l] + static_cast<int32_t>(lo - bucket.begin()),
                    static_cast<int32_t>(hi - lo));
        };

        const int32_t numLeaf = leaf().numRegions();
        flatColFirst.assign(numLeaf + 1, 0);
        flatLevelEnd.assign(static_cast<int64_t>(numLeaf) * levelStride, 0);
        flatColStart.assign(static_cast<int64_t>(numLeaf) * levelStride, 0);
        for (int32_t r = 0; r < numLeaf; ++r) {
            const int64_t base = static_cast<int64_t>(r) * levelStride;
            int32_t total = 0;
            for (int32_t l = 0; l < levelStride; ++l) {
                if (l < regionDepth[r]) { // columns end at the region's own root depth
                    const auto run = bucketRun(l, regionPrefix[r]);
                    flatColStart[base + l] = run.first;
                    total += run.second;
                }
                flatLevelEnd[base + l] = total;
            }
            flatColFirst[r + 1] = flatColFirst[r] + total;
        }

        const auto &lc = leaf();
        boundaryColPos.assign(lc.bIds.size(), -1);
        int64_t boundaryOffPath = 0;
        for (int32_t r = 0; r < numLeaf; ++r) {
            const HubId *const bpts = lc.bIds.data() + lc.bFirst[r];
            const int32_t n = static_cast<int32_t>(lc.bSize(r));
            forEachCol(flatColRow(r), flatLevelRow(r),
                       [&](const int32_t j, const HubId id) {
                           const auto it = std::lower_bound(bpts, bpts + n, id);
                           if (it != bpts + n && *it == id)
                               boundaryColPos[lc.bFirst[r] + (it - bpts)] = j;
                       });
        }
        for (size_t k = 0; k < boundaryColPos.size(); ++k)
            if (boundaryColPos[k] < 0)
                ++boundaryOffPath;
        if (boundaryOffPath > 0)
            std::cout << "FHL WARNING: " << boundaryOffPath << " of "
                      << boundaryColPos.size()
                      << " boundary points are not path ancestors of their region" << std::endl;

        tColFirst.assign(numHubs + 1, 0);
        tLevelEnd.assign(static_cast<int64_t>(numHubs) * levelStride, 0);
        tColStart.assign(static_cast<int64_t>(numHubs) * levelStride, 0);
        for (int32_t x = 0; x < numHubs; ++x) {
            const int32_t v = hierarchy.rankOfHub(x);
            const int32_t lv = hierarchy.getVertexLevel(v);
            const int64_t base = static_cast<int64_t>(x) * levelStride;
            int32_t total = 0;
            for (int32_t l = 0; l < levelStride; ++l) {
                if (l <= lv) { // own separator included: same-level siblings can cross
                    const auto run = bucketRun(l, hierarchy.getPackedSideId(v));
                    tColStart[base + l] = run.first;
                    total += run.second;
                }
                tLevelEnd[base + l] = total;
            }
            tColFirst[x + 1] = tColFirst[x] + total;
        }
    }

    void buildMeets(const HubHierarchy &hierarchy) {
        const int32_t m = numCuts();
        meets.assign(m, {});
        for (int32_t mi = 0; mi < m; ++mi) {
            auto &ms = meets[mi];
            ms.depth = mi == 0 ? 0 : cutDepths[mi - 1];
            ms.sepEnd = cutDepths[mi];
            const int32_t numR = mi == 0 ? 1 : cuts[mi - 1].numRegions();
            std::vector<std::vector<HubId>> sSets(numR);
            const uint64_t mask = ms.depth > 0 ? (1ULL << ms.depth) - 1 : 0;
            for (int32_t id = 0; id < numHubs; ++id) {
                const int32_t v = hierarchy.rankOfHub(id);
                const int32_t lv = hierarchy.getVertexLevel(v);
                if (lv < ms.depth || lv >= ms.sepEnd)
                    continue;
                const int32_t r = mi == 0 ? 0
                                          : cuts[mi - 1].regionOfPrefix[
                                                    hierarchy.getPackedSideId(v) & mask];
                KASSERT(r >= 0);
                sSets[r].push_back(id);
            }
            // Layout: the PARENT cut's boundary points occupy the first slots of each meet
            // region's S list (they are the climb's output columns), the band's own
            // separator vertices follow. Both halves stay sorted, so lookups are two
            // binary searches, and the climb's inner loop reads M row-major with no
            // column indirection (the parent boundary point j sits at S position j).
            for (auto &set : sSets)
                std::sort(set.begin(), set.end());
            ms.numParentBoundary.assign(numR, 0);
            ms.sFirst.assign(numR + 1, 0);
            for (int32_t r = 0; r < numR; ++r) {
                const int64_t nPar =
                        mi == 0 ? 0
                                : cuts[mi - 1].bFirst[r + 1] - cuts[mi - 1].bFirst[r];
                ms.numParentBoundary[r] = static_cast<int32_t>(nPar);
                ms.sFirst[r + 1] = ms.sFirst[r] + nPar +
                                   static_cast<int64_t>(sSets[r].size());
            }
            ms.sIds.resize(ms.sFirst[numR]);
            for (int32_t r = 0; r < numR; ++r) {
                auto at = ms.sIds.begin() + ms.sFirst[r];
                if (mi > 0)
                    at = std::copy(cuts[mi - 1].bIds.begin() + cuts[mi - 1].bFirst[r],
                                   cuts[mi - 1].bIds.begin() + cuts[mi - 1].bFirst[r + 1],
                                   at);
                std::copy(sSets[r].begin(), sSets[r].end(), at);
            }
            if (mi > 0) {
                const auto &pc = cuts[mi - 1];
                ms.bPos.resize(pc.bIds.size());
                for (int32_t r = 0; r < numR; ++r)
                    for (int64_t k = pc.bFirst[r]; k < pc.bFirst[r + 1]; ++k)
                        ms.bPos[k] = static_cast<int32_t>(k - pc.bFirst[r]);
            }
        }

        for (int32_t ci = 0; ci < m; ++ci) {
            auto &c = cuts[ci];
            const auto &ms = meets[ci];
            c.sPosOfB.resize(c.bIds.size());
            for (int32_t r = 0; r < c.numRegions(); ++r) {
                const int32_t mr = c.parentRegion[r];
                for (int64_t k = c.bFirst[r]; k < c.bFirst[r + 1]; ++k) {
                    c.sPosOfB[k] = ms.posOf(mr, c.bIds[k]);
                    KASSERT(c.sPosOfB[k] >= 0);
                }
            }
        }

        initCutOfHub.assign(numHubs, -1);
        initRegionOfHub.assign(numHubs, -1);
        selfSPosOfHub.assign(numHubs, -1);
        for (int32_t id = 0; id < numHubs; ++id) {
            const int32_t v = hierarchy.rankOfHub(id);
            const int32_t lv = hierarchy.getVertexLevel(v);
            int32_t ci = -1;
            while (ci + 1 < m && cutDepths[ci + 1] <= lv)
                ++ci;
            initCutOfHub[id] = ci;
            const int32_t mi = ci + 1; // own meet level: depth <= lv < sepEnd
            KASSERT(mi < m && lv >= meets[mi].depth && lv < meets[mi].sepEnd);
            int32_t mr = 0;
            if (ci >= 0) {
                const uint64_t mask = (1ULL << cutDepths[ci]) - 1;
                const int32_t r =
                        cuts[ci].regionOfPrefix[hierarchy.getPackedSideId(v) & mask];
                KASSERT(r >= 0);
                initRegionOfHub[id] = r;
                mr = r;
            }
            selfSPosOfHub[id] = meets[mi].posOf(mr, id);
            KASSERT(selfSPosOfHub[id] >= 0);
        }
    }
};


} // namespace fhl
