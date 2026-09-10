#pragma once

#include <algorithm>
#include <cstdint>
#include <queue>
#include <vector>

#include "Algorithms/CCH/CCHMetric.h"
#include "Algorithms/FHL/core/HubAccessData.h"
#include "Algorithms/FHL/FHLData.h"
#include "Algorithms/FHL/core/HubTableBuilder.h"
#include "Algorithms/FHL/core/FHLConstants.h"
#include "Tools/Constants.h"
#include "Tools/Timer.h"

// Customization: refills every metric-dependent value while the skeleton (hierarchy,
// regions, column runs) stays untouched -- that is what makes the index customizable.
//
// Full pass for the seed/rectangle family, in order (the [cust] log lines follow it):
//   cch customize -> access hubs -> span-1 hub tables (transient scaffolding) ->
//   rectangles + canonical prune + CSR -> seeds -> optional in-region rows / top labels ->
//   release the scaffolding
// The overlay family replaces the middle of that with region-restricted cliques.
//
// Three update gears (-fhl-subtree-rebuild turns the middle one on): a change confined to one
// leaf region refreshes that region only; a change reaching the hub layer rebuilds the
// smallest separator subtree containing it; only a change at the top needs the full pass.
template<typename CchT>
class FHLMetric {
public:
    FHLMetric(const HubHierarchy &hierarchy, const CchT &cch,
                           const std::vector<HubAccessEdge> &accessNodeEdges,
                           const int32_t *const inputWeights,
                           const fhl::Level pruneLevelThreshold = 0,
                           // 1 = seed/rectangle family, 2 = overlay
                           const int32_t regionMode = 1,
                           const fhl::Regions *const regions = nullptr,
                           const int32_t flatTopLab = 0,
                           const bool keepForIncremental = false,
                           const bool superRebuild = false)
            : hierarchy(hierarchy),
              cch(cch),
              accessNodeEdges(accessNodeEdges),
              inputWeights(inputWeights),
              pruneLevelThreshold(pruneLevelThreshold),
              regionMode(regionMode),
              regions(regions),
              flatTopLab(flatTopLab),
              keepForIncremental(keepForIncremental),
              superRebuild(superRebuild),
              cchMetric(cch, inputWeights),
              localEliminationTree(cch.getEliminationTree()) {
        for (int &v : localEliminationTree) {
            if (v != INVALID_VERTEX && hierarchy.isHub(v))
                v = INVALID_VERTEX;
        }
    }

    void customize(FHLData &data) {
        symCache = -1; // weights change: recompute the symmetry flag on demand
        Timer cchTimer;
        cchMetric.customize();
        std::cout << "[cust] cch customize: "
                  << cchTimer.elapsed<std::chrono::milliseconds>() << " ms\n";
        if (regionMode == 2) { // overlay: cliques come from the original edges only
            fhl::Overlay::buildOverlay(*regions, hierarchy, cch, cchMetric.upwardWeights(),
                                     cchMetric.downwardWeights(), localEliminationTree,
                                     inputWeights, metricIsSymmetric(), data.ladder);
            releaseAccessHubs(data.getBaseData()); // preallocated, unused in this mode
            return;
        }
        Timer stageTimer;
        computeAccessNodes(data.getBaseData());
        if (regionMode == 1) { // seed/rectangle family: read off the span-1 tables
            std::cout << "[cust] access nodes: "
                      << stageTimer.elapsed<std::chrono::milliseconds>() << " ms\n";
            stageTimer.restart();
            HubTableBuilder tableBuilder;
            std::vector<HubTableBuilder::TableSpec> specs;
            buildSpan1Specs(specs);
            std::vector<HubTableBuilder::TableData> tables;
            tableBuilder.buildDistanceTables(specs, cch.getUpwardGraph(),
                                             cchMetric.upwardWeights(),
                                             cchMetric.downwardWeights(),
                                             cch.getEliminationTree(), hierarchy, tables);
            std::cout << "[cust] span-1 oracle build: "
                      << stageTimer.elapsed<std::chrono::milliseconds>() << " ms\n";
            stageTimer.restart();
            data.setTables(std::move(tables));
            std::cout << "[cust] setTables (index maps): "
                      << stageTimer.elapsed<std::chrono::milliseconds>() << " ms\n";
            buildFlat(*regions, hierarchy, data, metricIsSymmetric(),
                                  flatTopLab, data.ladder, superRebuild);
            if (regions->locLevels > 0) {
                Timer lt;
                fhl::InRegionHubs::buildLocalRows(*regions, hierarchy, cch,
                                           cchMetric.upwardWeights(),
                                           cchMetric.downwardWeights(), metricIsSymmetric(),
                                           /*regionList=*/nullptr, 0, data.ladder);
                std::cout << "[cust] in-region rows: "
                          << lt.elapsed<std::chrono::milliseconds>() << " ms\n";
            }
            stageTimer.restart();
            if (superRebuild) {
                buildTableSuperIndex(data); // tables stay alive: the middle gear rewrites them
            } else {
                data.releaseTables();
            }
            if (!keepForIncremental) // incremental mode reuses them on full rebuilds
                releaseAccessHubs(data.getBaseData());
            std::cout << "[cust] releases: " << stageTimer.elapsed<std::chrono::milliseconds>()
                      << " ms\n";
            return;
        }
    }

    // ------------------------------------------------------ incremental flat customization
    // Exact queue-driven partial CCH customization followed by a flat fast-path refresh.
    // Seeds: the CCH edges whose input weights changed (all inside leaf region `region`).
    // Each popped edge is recomputed FROM ITS BASE (respecting weight + all lower
    // triangles) in ascending (tail, edge) order — the same relaxation order as full
    // customization restricted to the affected closure, so increases and decreases are both
    // handled exactly. If no recomputed edge with a hub-zone tail changed, only the
    // region's seeds/top labels can differ (everything else is a function of hub-zone
    // shortcuts) and they are refreshed locally. Otherwise the caller must run a full
    // customization (NEED_FULL).
    // PARTIAL_SUPER: shortcuts between hubs changed, but every changed shortcut sits inside one
    // super region (a coarser size cut). Only that super region's tables, rectangles and
    // seeds are rebuilt — the middle gear between the seed-only fast path and a full rebuild.
    enum PartialOutcome {
        PARTIAL_NOOP = 0,
        PARTIAL_FAST = 1,
        PARTIAL_SUPER = 2,
        PARTIAL_NEED_FULL = 3
    };

    struct PartialStats {
        PartialOutcome outcome = PARTIAL_NOOP;
        int64_t recomputed = 0;
        int64_t changed = 0;
        int32_t minChangedLevel = 127; // shallowest level a shortcut actually changed at
        int32_t superRegion = -1;      // depth of the rebuilt scope (PARTIAL_SUPER) or -1
    };

    PartialStats partialCustomizeFlat(FHLData &data,
                                      const std::vector<int32_t> &seedCchEdges,
                                      const int32_t region) {
        KASSERT(regionMode == 1 && regions != nullptr);
        const auto &graph = cch.getUpwardGraph();
        auto &up = cchMetric.upWeights;
        auto &down = cchMetric.downWeights;
        if (pcInQueue.size() != static_cast<size_t>(graph.numEdges()))
            pcInQueue.assign(graph.numEdges(), 0);
        using QE = std::pair<int64_t, int32_t>; // ((tail << 32) | edge, edge)
        std::priority_queue<QE, std::vector<QE>, std::greater<QE>> heap;
        const auto push = [&](const int32_t e) {
            if (pcInQueue[e])
                return;
            pcInQueue[e] = 1;
            pcTouched.push_back(e);
            heap.emplace((static_cast<int64_t>(graph.edgeTail(e)) << 32) | e, e);
        };
        for (const auto e : seedCchEdges)
            push(e);

        PartialStats st;
        bool capped = false;
        bool anyHubChanged = false;
        bool becameAsymmetric = false; // a recomputed shortcut now differs by direction
        // Scope of the hub-zone changes: the deepest separator subtree containing every
        // changed shortcut, as (depth, side-id prefix). Depth 0 = the whole graph.
        int32_t scopeDepth = 127;
        uint64_t scopePrefix = 0;
        bool haveScope = false;
        int32_t minChangedLevel = 127;
        while (!heap.empty()) {
            const int32_t e = heap.top().second;
            heap.pop();
            const int32_t t = graph.edgeTail(e);
            const int32_t h = graph.edgeHead(e);
            ++st.recomputed;
            int32_t newUp = INFTY, newDown = INFTY;
            cch.forEachUpwardInputEdge(e, [&](const int ie) {
                newUp = std::min(newUp, inputWeights[ie]);
                return true;
            });
            cch.forEachDownwardInputEdge(e, [&](const int ie) {
                newDown = std::min(newDown, inputWeights[ie]);
                return true;
            });
            cch.forEachLowerTriangle(t, h, e, [&](int, const int ezt, const int ezh) {
                if (down[ezt] + up[ezh] < newUp)
                    newUp = down[ezt] + up[ezh];
                if (down[ezh] + up[ezt] < newDown)
                    newDown = down[ezh] + up[ezt];
                return true;
            });
            if (newUp == up[e] && newDown == down[e])
                continue;
            up[e] = newUp;
            down[e] = newDown;
            becameAsymmetric |= newUp != newDown;
            ++st.changed;
            minChangedLevel = std::min(minChangedLevel,
                                       static_cast<int32_t>(hierarchy.getVertexLevel(t)));
            if (hierarchy.isHub(t)) {
                anyHubChanged = true;
                // A changed shortcut only moves the labels of the tail's descendants, so the
                // affected subtree is the one containing every changed tail: intersect the
                // side-id prefixes and cap the depth at each tail's own level.
                const uint64_t side = hierarchy.getPackedSideId(t);
                const int32_t lv = hierarchy.getVertexLevel(t);
                if (!haveScope) {
                    scopePrefix = side;
                    scopeDepth = lv;
                    haveScope = true;
                } else {
                    const uint64_t diff = scopePrefix ^ side;
                    if (diff != 0)
                        scopeDepth = std::min(scopeDepth, lowestOneBit(diff));
                    scopeDepth = std::min(scopeDepth, lv);
                }
            }
            if (st.recomputed > 500000) {
                capped = true; // safety cap: force a full rebuild, state is only partial
                break;
            }
            // Dependents: (t,h) participates in a triangle at t with EVERY other up-edge
            // (t,x) — as the lower leg (x > h) or as the inter leg (x < h); the dependent
            // is the edge between the two heads whenever it exists. (Enumerating only the
            // x > h half — forEachUpperTriangle — under-propagates: caught by the
            // perturbed-state verification.)
            for (int e2 = graph.firstEdge(t); e2 < graph.lastEdge(t); ++e2) {
                if (e2 == e)
                    continue;
                const int32_t x = graph.edgeHead(e2);
                const int32_t lo = std::min(x, h), hi2 = std::max(x, h);
                int fe = graph.firstEdge(lo), fl = graph.lastEdge(lo);
                while (fe < fl && graph.edgeHead(fe) < hi2)
                    ++fe;
                if (fe < fl && graph.edgeHead(fe) == hi2)
                    push(fe);
            }
        }
        for (const auto e : pcTouched)
            pcInQueue[e] = 0;
        pcTouched.clear();

        st.minChangedLevel = minChangedLevel;
        st.superRegion = haveScope ? scopeDepth : -1;
        // Structures built for a symmetric metric store one direction only; an update that
        // makes the metric asymmetric cannot be represented incrementally.
        const bool structuresSymmetric = data.ladder.seedBwd.empty();
        if (st.changed == 0) {
            st.outcome = PARTIAL_NOOP;
        } else if (capped || (structuresSymmetric && becameAsymmetric)) {
            st.outcome = PARTIAL_NEED_FULL; // partial state: must rebuild fully
        } else if (anyHubChanged && superRebuild && haveScope && scopeDepth > 0 &&
                   rebuildScope(data, scopeDepth, scopePrefix)) {
            st.outcome = PARTIAL_SUPER;
        } else if (!anyHubChanged) {
            fhl::Overlay::refreshRegionAfterInteriorChange(
                    *regions, hierarchy, cch, cchMetric.upwardWeights(),
                    cchMetric.downwardWeights(), localEliminationTree, region, data.ladder);
            st.outcome = PARTIAL_FAST;
        } else {
            st.outcome = PARTIAL_NEED_FULL; // caller reruns the full flat build
        }
        return st;
    }

    // Table (level, prefix) coordinates, kept alongside the live tables so the middle gear
    // can select the ones inside a scope subtree in one linear pass.
    void buildTableSuperIndex(FHLData &data) {
        tableLevel.assign(data.numTables(), -1);
        tablePrefix.assign(data.numTables(), 0);
        for (int32_t ti = 0; ti < data.numTables(); ++ti) {
            const auto &t = data.getTable(ti);
            tableLevel[ti] = t.depthBase;
            tablePrefix[ti] = t.prefix;
        }
    }

    // Middle gear: the metric change reached the hub zone but is confined to the
    // separator subtree (scopeDepth, scopePrefix). Only that subtree's tables can move (a
    // changed shortcut alters the elimination-path labels of its descendants only, and no
    // table or leaf region outside the subtree has a row, column or boundary point
    // inside it), so the
    // rebuild is: subset tables -> their leaf cliques and CSR segments -> their hub rows
    // -> their seeds and top labels. Returns false when the scope is too large to pay off,
    // leaving the caller to run a full customization.
    bool rebuildScope(FHLData &data, const int32_t scopeDepth,
                      const uint64_t scopePrefix) {
        const uint64_t mask = (1ULL << scopeDepth) - 1;
        const uint64_t want = scopePrefix & mask;
        const auto inScope = [&](const int32_t depth, const uint64_t prefix) {
            return depth >= scopeDepth && (prefix & mask) == want;
        };
        std::vector<int32_t> regionList;
        for (int32_t r = 0; r < regions->leaf().numRegions(); ++r)
            if (inScope(regions->regionDepth[r], regions->regionPrefix[r]))
                regionList.push_back(r);
        // A scope covering most of the graph costs more than a plain full rebuild (which
        // also parallelizes better), so hand those back to the caller.
        if (regionList.size() * 3 > static_cast<size_t>(regions->leaf().numRegions()))
            return false;
        std::vector<int32_t> tableIdx;
        for (int32_t ti = 0; ti < static_cast<int32_t>(tableLevel.size()); ++ti)
            if (inScope(tableLevel[ti], tablePrefix[ti]))
                tableIdx.push_back(ti);
        if (tableIdx.empty())
            return false;

        HubTableBuilder tableBuilder;
        tableBuilder.rebuildTablesSubset(tableIdx, cch.getUpwardGraph(),
                                         cchMetric.upwardWeights(),
                                         cchMetric.downwardWeights(),
                                         cch.getEliminationTree(), hierarchy,
                                         data.mutableTables());
        const bool symmetric = data.ladder.seedBwd.empty();
        const int32_t numRegionsInList = static_cast<int32_t>(regionList.size());
#pragma omp parallel for schedule(dynamic, 4)
        for (int32_t i = 0; i < numRegionsInList; ++i)
            fhl::SeedRectangle::rebuildLeafRegionFromTables(*regions, data, symmetric, regionList[i],
                                                    data.ladder);
        // Hub rows of the scope: their columns are ancestors, whose distances the
        // refreshed tables now hold.
        std::vector<fhl::HubId> scopeHubs;
        for (int32_t x = 0; x < hierarchy.numHubs(); ++x) {
            const int32_t v = hierarchy.rankOfHub(x);
            if (inScope(hierarchy.getVertexLevel(v), hierarchy.getPackedSideId(v)))
                scopeHubs.push_back(x);
        }
        fhl::SeedRectangle::buildHubRows(*regions, data, symmetric, data.ladder,
                                         scopeHubs.data(),
                                         static_cast<int32_t>(scopeHubs.size()));
        fhl::Overlay::refreshRegionsFromCliques(*regions, hierarchy, cch,
                                              cchMetric.upwardWeights(),
                                              cchMetric.downwardWeights(),
                                              localEliminationTree, regionList.data(),
                                              numRegionsInList, data.ladder);
        return true;
    }

    // Span-1 table specs built from (level, prefix) buckets in output-linear time — the
    // generic recursive builder scans all hubs once PER SPEC.
    // Reproduces its exact layout: rows = own separator (ids ascending), ancestor columns in
    // globally ascending id order (= deepest ancestor level first).
    void buildSpan1Specs(std::vector<HubTableBuilder::TableSpec> &specs) {
        const int32_t numHubs = hierarchy.numHubs();
        const int32_t thresh = hierarchy.getHubLevelThreshold();
        struct Key {
            uint64_t prefix;
            int32_t level;
            int32_t id;
        };
        std::vector<Key> keys(numHubs);
        for (int32_t id = 0; id < numHubs; ++id) {
            const int32_t v = hierarchy.rankOfHub(id);
            const int32_t lv = hierarchy.getVertexLevel(v);
            const uint64_t mask = lv > 0 ? (1ULL << lv) - 1 : 0;
            keys[id] = {hierarchy.getPackedSideId(v) & mask, lv, id};
        }
        // Bucket start indices per (level, prefix): ids are already sorted by (level desc,
        // rank asc), so within one level the keys are contiguous and prefix-sortable.
        std::vector<int32_t> order(numHubs);
        for (int32_t i = 0; i < numHubs; ++i)
            order[i] = i;
        std::sort(order.begin(), order.end(), [&](const int32_t a, const int32_t b) {
            if (keys[a].level != keys[b].level)
                return keys[a].level > keys[b].level; // deepest first = ascending ids
            if (keys[a].prefix != keys[b].prefix)
                return keys[a].prefix < keys[b].prefix;
            return keys[a].id < keys[b].id;
        });
        // Per (level, prefix) group: one spec; ancestors = walk the chain buckets from the
        // deepest ancestor level up (globally ascending ids).
        struct Group {
            int32_t level;
            uint64_t prefix;
            int32_t first, last;
        };
        std::vector<Group> groups;
        for (int32_t i = 0; i < numHubs;) {
            int32_t j = i;
            while (j < numHubs && keys[order[j]].level == keys[order[i]].level &&
                   keys[order[j]].prefix == keys[order[i]].prefix)
                ++j;
            groups.push_back({keys[order[i]].level, keys[order[i]].prefix, i, j});
            i = j;
        }
        // Group lookup by (level, prefix) for the ancestor walks.
        std::vector<std::pair<uint64_t, int32_t>> groupAt; // key = (level << 40) | prefix
        groupAt.reserve(groups.size());
        for (int32_t g = 0; g < static_cast<int32_t>(groups.size()); ++g)
            groupAt.emplace_back((static_cast<uint64_t>(groups[g].level) << 40) |
                                         groups[g].prefix,
                                 g);
        std::sort(groupAt.begin(), groupAt.end());
        const auto findGroup = [&](const int32_t level, const uint64_t prefix) -> int32_t {
            const uint64_t key = (static_cast<uint64_t>(level) << 40) | prefix;
            const auto it = std::lower_bound(groupAt.begin(), groupAt.end(),
                                             std::make_pair(key, INT32_MIN));
            return (it != groupAt.end() && it->first == key) ? it->second : -1;
        };
        specs.clear();
        specs.resize(groups.size());
#pragma omp parallel for schedule(dynamic, 64)
        for (int32_t g = 0; g < static_cast<int32_t>(groups.size()); ++g) {
            const auto &grp = groups[g];
            auto &spec = specs[g];
            spec.sepNode = 0;
            spec.depthBase = grp.level;
            spec.span = 1;
            spec.prefixLevel = grp.level;
            spec.prefix = grp.prefix;
            spec.rowRanks.reserve(grp.last - grp.first);
            for (int32_t i = grp.first; i < grp.last; ++i)
                spec.rowRanks.push_back(
                        hierarchy.rankOfHub(keys[order[i]].id));
            for (int32_t l = grp.level - 1; l >= 0; --l) { // deepest ancestors first
                const uint64_t lp = l > 0 ? (grp.prefix & ((1ULL << l) - 1)) : 0;
                const int32_t ag = findGroup(l, lp);
                if (ag < 0)
                    continue;
                for (int32_t i = groups[ag].first; i < groups[ag].last; ++i)
                    spec.ancestorRanks.push_back(
                            hierarchy.rankOfHub(keys[order[i]].id));
            }
        }
        (void) thresh;
    }

    bool metricIsSymmetric() const {
        if (symCache < 0) {
            const int32_t numEdges = cch.getUpwardGraph().numEdges();
            const int32_t *const up = cchMetric.upwardWeights();
            const int32_t *const down = cchMetric.downwardWeights();
            bool equal = true;
#pragma omp parallel for schedule(static) reduction(&& : equal)
            for (int32_t e = 0; e < numEdges; ++e)
                equal = equal && up[e] == down[e];
            symCache = equal ? 1 : 0;
        }
        return symCache == 1;
    }

    static void releaseAccessHubs(HubAccessData &data) {
        data.pos.clear();
        data.pos.shrink_to_fit();
        data.accessHubs.clear();
        data.accessHubs.shrink_to_fit();
        data.forwardDistances.clear();
        data.forwardDistances.shrink_to_fit();
        data.backwardDistances.clear();
        data.backwardDistances.shrink_to_fit();
    }

    const std::vector<int32_t> &getLocalEliminationTree() const { return localEliminationTree; }

    const BaseCCHMetric<CchT> &getCCHMetric() const { return cchMetric; }

    uint64_t sizeInBytes() const {
        uint64_t size = sizeof(FHLMetric);
        size += cchMetric.sizeInBytes();
        return size;
    }

private:
    // Orchestrates one full customization of the seed/rectangle family: the shared
    // primitives first, then the optional per-scheme extras.
    // hdata must hold the transient span-1 subtree tables and the (unpruned) access nodes;
    // the caller releases both afterwards.
    template<typename HData>
    static void buildFlat(const fhl::Regions &R, const HubHierarchy &hierarchy,
                          HData &hdata, const bool symmetric, const int32_t topLabLevels,
                          fhl::Values &out, const bool relocatableCsr = false) {
        out.mode = 1;
        out.symmetric = symmetric;
        Timer t;
        fhl::SeedRectangle::buildLeafClique(R, hdata, out);
        std::cout << "[cust] leaf cliques: " << t.elapsed<std::chrono::milliseconds>()
                  << " ms\n";
        t.restart();
        fhl::SeedRectangle::buildHubRows(R, hdata, symmetric, out);
        std::cout << "[cust] hub rows: " << t.elapsed<std::chrono::milliseconds>()
                  << " ms\n";
        t.restart();
        // Dense rectangles live through the CSR build AND the seed fold (the access nodes
        // are columns of their own region, so seeds read rectangle columns instead of
        // probing the span-1 tables vertex by vertex).
        std::vector<int64_t> valOff;
        std::vector<int32_t> denseOut, denseIn;
        fhl::SeedRectangle::buildCsr(R, hdata, symmetric, out, valOff, denseOut, denseIn,
                                     relocatableCsr);
        std::cout << "[cust] rect+prune+CSR: " << t.elapsed<std::chrono::milliseconds>()
                  << " ms\n";
        t.restart();
        fhl::SeedRectangle::buildSeedsFromRect(R, hierarchy, hdata, symmetric, valOff, denseOut, denseIn, out);
        std::cout << "[cust] seeds: " << t.elapsed<std::chrono::milliseconds>() << " ms\n";
        denseOut.clear();
        denseOut.shrink_to_fit();
        denseIn.clear();
        denseIn.shrink_to_fit();
        if (topLabLevels > 0) {
            t.restart();
            fhl::TopLabels::buildTopLabels(R, symmetric, std::min(topLabLevels, R.levelStride), out);
            std::cout << "[cust] top labels: " << t.elapsed<std::chrono::milliseconds>()
                      << " ms\n";
        }
    }


    // Seeds via the TNR access-node theorem, folded from the dense rectangles: every
    // access node is a path ancestor of the branch and hence a COLUMN of its own region, so
    // d(a, b_j) is a rectangle read (dIn[j][col(a)]) instead of a span-1 table probe.
    // Iterating region by region keeps each rectangle cache-resident.


    void computeAccessNodes(HubAccessData &data) {
        const auto rankToIdx = [&](const int r) {
            return data.rankToIdx(r);
        };
        const auto cchUpWeights = cchMetric.upwardWeights();
        const auto cchDownWeights = cchMetric.downwardWeights();

        // Symmetric metric: forward and backward access distances coincide; store and sweep
        // one copy only (the accessors fall back to forward when backward is empty).
        const bool symmetric = metricIsSymmetric();
        if (symmetric) {
            std::cout << "Symmetric metric: storing one access node direction." << std::endl;
            data.backwardDistances.clear();
            data.backwardDistances.shrink_to_fit();
        } else if (data.backwardDistances.size() != data.forwardDistances.size()) {
            data.backwardDistances.resize(data.forwardDistances.size());
        }

        // Nearly 2 GB on USA: a single-threaded fill is ~270 ms regardless of the thread count.
        const auto fillInfinity = [](auto &vec) {
            auto *const p = vec.data();
            const int64_t n = static_cast<int64_t>(vec.size());
#pragma omp parallel for schedule(static)
            for (int64_t i = 0; i < n; ++i)
                p[i] = FHL_INFTY;
        };
        fillInfinity(data.forwardDistances);
        fillInfinity(data.backwardDistances);

#pragma omp parallel for schedule(static)
        for (const auto &accessNodeEdge : accessNodeEdges) {
            const auto &e = accessNodeEdge.edge;
            auto &f = data.forwardDistances[accessNodeEdge.position];
            f = std::min(f, cchUpWeights[e]);
            if (!symmetric) {
                auto &b = data.backwardDistances[accessNodeEdge.position];
                b = std::min(b, cchDownWeights[e]);
            }
        }

        const auto &cchGraph = cch.getUpwardGraph();
        cch.forEachVertexTopDown([&](int32_t rv) {
            computeAccessNodeDistancesForVertex(rv, cchGraph, cchUpWeights, rankToIdx, data.pos,
                                                data.accessHubs, data.forwardDistances);
            if (!symmetric)
                computeAccessNodeDistancesForVertex(rv, cchGraph, cchDownWeights, rankToIdx,
                                                    data.pos, data.accessHubs,
                                                    data.backwardDistances);
        });
    }

    template<typename GraphT, typename RankToIdx>
    void computeAccessNodeDistancesForVertex(const int rv,
                                             const GraphT &graph,
                                             int const *const weights,
                                             const RankToIdx &rankToIdx,
                                             const std::vector<int32_t> &dataPos,
                                             const std::vector<fhl::HubId> &dataNodes,
                                             HubAccessData::DistanceVector<int32_t> &dataDistances) const {
        const int idx = rankToIdx(rv);
        const int startThis = dataPos[idx];
        HubAccessData::DistanceLabel distancesThis;
        HubAccessData::DistanceLabel distancesNeighbor;
        FORALL_INCIDENT_EDGES(graph, rv, e) {
            const int neighbor = graph.edgeHead(e);
            const int idxNeighbor = rankToIdx(neighbor);
            const int w = weights[e];
            KASSERT(w < INFTY);

            const auto startNeighbor = dataPos[idxNeighbor];
            const auto numNeighbor = dataPos[idxNeighbor + 1] - startNeighbor;
            KASSERT(numNeighbor <= dataPos[idx + 1] - startThis);
            if constexpr (!HubAccessData::USE_SIMD) {
                for (auto i = 0; i < numNeighbor; ++i) {
                    KASSERT(dataNodes[startThis + i] == dataNodes[startNeighbor + i]);
                    const int newDist = dataDistances[startNeighbor + i] + w;
                    auto &entry = dataDistances[startThis + i];
                    if (newDist < entry)
                        entry = newDist;
                }
            } else {
                KASSERT(numNeighbor % HubAccessData::K == 0);
                const auto numBatches = numNeighbor / HubAccessData::K;
                for (auto b = 0; b < numBatches; ++b) {
                    const auto offset = b * HubAccessData::K;
                    distancesThis.load(&dataDistances[startThis + offset]);
                    distancesNeighbor.load(&dataDistances[startNeighbor + offset]);
                    distancesThis.min(distancesNeighbor + w);
                    distancesThis.store(&dataDistances[startThis + offset]);
                }
            }
        }
    }

    const HubHierarchy &hierarchy;
    const CchT &cch;
    const std::vector<HubAccessEdge> &accessNodeEdges;
    const int32_t *const inputWeights;
    const fhl::Level pruneLevelThreshold;
    const int32_t regionMode;
    const fhl::Regions *const regions;
    const int32_t flatTopLab;
    const bool keepForIncremental;
    const bool superRebuild;
    std::vector<int32_t> tableLevel;   // per live table: separator level and side prefix
    std::vector<uint64_t> tablePrefix;
    BaseCCHMetric<CchT> cchMetric;
    std::vector<int32_t> localEliminationTree;

    mutable int8_t symCache = -1; // -1 unknown; symmetry of the current weights

    // Scratch for partial customization (lazily sized, reset via the touched list).
    std::vector<uint8_t> pcInQueue;
    std::vector<int32_t> pcTouched;
};
