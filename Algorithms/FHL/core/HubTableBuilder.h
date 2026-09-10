#pragma once

#include <algorithm>
#include <cstdint>
#include <type_traits>
#include <map>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "Algorithms/CCH/EliminationTreeQuery.h"
#include "Algorithms/FHL/core/FHLConstants.h"
#include "DataStructures/Labels/BasicLabelSet.h"
#include "DataStructures/Labels/ParentInfo.h"
#include "DataStructures/Partitioning/SeparatorDecomposition.h"
#include "Tools/Timer.h"
#include "Algorithms/FHL/core/HubHierarchy.h"

// Builds distance tables for FHL. Hubs are partitioned into tables by
// bands of separator decomposition levels (and, in subtree mode, by subtree of the separator
// decomposition). Each table is rectangular: rows are the hubs owned by the table,
// columns are the rows plus all hubs on the separator tree path above the table's
// subtree ("ancestor columns"). Both directions are stored (row->col and col->row), so that
// for any pair (a, b) of hubs where one lies on the separator path of the other, the
// deeper node's table holds the exact distance in both directions. Together with the fact that
// TNR query candidates always lie on the shared root path of s and t, this makes every query
// resolvable by table lookups alone, without fallback graph searches.
// Hub-to-hub distance tables over bands of the separator decomposition. Shared
// infrastructure, NOT a scheme: FHL's customization builds span-1 tables with it as
// transient scaffolding (the [cust] span-1 oracle stage, its largest phase) and
// then discards them. Wider bands are supported but nothing uses them today.
class HubTableBuilder {
public:
    struct TableSpec {
        int32_t sepNode = 0;
        int32_t depthBase = 0;
        int32_t span = 0;
        int32_t prefixLevel = 0;  // depth of the subtree root; the prefix has this many bits
        uint64_t prefix = 0;      // packed side id of the subtree root
        std::vector<int32_t> rowRanks;      // CCH ranks of hubs owned by this table
        std::vector<int32_t> ancestorRanks; // CCH ranks of hubs above the band that are
                                            // ancestors of at least one row (path of the subtree
                                            // root plus subtree-interior nodes below prefixLevel)
    };

    // A table's slice of one allocation shared by every table of a build: all four arrays of
    // all ~24K tables (USA) live in a single pool. As ~96K separate allocations they cost
    // ~115 ms serially and contend when issued in parallel; one allocation, first touched in
    // parallel, costs neither. Ownership is shared, so the pool lives exactly as long as
    // some table refers to it -- the tables move into FHLData and die in the releases phase
    // as before. Incremental rebuilds (rebuildTablesSubset) rewrite slices in place.
    class PooledArray {
    public:
        PooledArray() = default;
        PooledArray(std::shared_ptr<int32_t[]> owner, const size_t offset, const size_t size)
            : pool(std::move(owner)), ptr(pool.get() + offset), len(size) {}
        int32_t &operator[](const size_t i) { return ptr[i]; }
        const int32_t &operator[](const size_t i) const { return ptr[i]; }
        int32_t *data() { return ptr; }
        const int32_t *data() const { return ptr; }
        const int32_t *begin() const { return ptr; }
        const int32_t *end() const { return ptr + len; }
        size_t size() const { return len; }
        bool empty() const { return len == 0; }

    private:
        std::shared_ptr<int32_t[]> pool;
        int32_t *ptr = nullptr;
        size_t len = 0;
    };

    struct TableData {
        int32_t sepNode = 0;
        int32_t depthBase = 0;
        int32_t span = 0;
        int32_t prefixLevel = 0;
        uint64_t prefix = 0;
        int32_t numRows = 0;
        int32_t numCols = 0; // == numRows + number of ancestor columns
        PooledArray rowIds; // hub ids of rows (fhl::HubId is int32_t)
        PooledArray colIds; // hub ids of columns: rows first, then ancestors
        PooledArray distOut; // numRows x numCols, distOut[i * numCols + j] = d(row i -> col j)
        PooledArray distIn;  // numRows x numCols, distIn[i * numCols + j] = d(col j -> row i);
                                      // empty on symmetric metrics (distOut serves both)

        // Backward-direction array: on symmetric metrics distIn is not stored and distOut has
        // identical values at identical offsets.
        const int32_t *inArr() const { return distIn.empty() ? distOut.data() : distIn.data(); }

        uint64_t sizeInBytes() const {
            return sizeof(TableData) +
                   rowIds.size() * sizeof(fhl::HubId) +
                   colIds.size() * sizeof(fhl::HubId) +
                   distOut.size() * sizeof(int32_t) +
                   distIn.size() * sizeof(int32_t);
        }
    };

    // Fills the distance entries of each table with exact shortest-path distances. Instead of
    // one elimination tree query per entry (which re-relaxes the dense top-of-hierarchy edges
    // for every pair), the two halves of the elimination tree query are decoupled: every column
    // does one upward and one downward sweep along its root path (labels stored sparsely), and
    // every row does the same into thread-local dense arrays. An entry is then just a merge of
    // the row's dense labels with the column's sparse labels over the column's root path:
    //   d(row -> col) = min over path(col) vertices u of  dup(row, u) + ddown(u, col)
    //   d(col -> row) = min over path(col) vertices u of  dup(col, u) + ddown(u, row)
    // Rows are processed in parallel.
    // Rebuilds the values of the tables in `tableIdx` in place, leaving every other table
    // untouched. Correct because a metric change confined to one super region can only move
    // the elimination-path labels of vertices INSIDE that region (a changed shortcut's tail
    // is on the path of its descendants only), and no table outside the region has a row or
    // column there. Phase A therefore sweeps only the ancestor closure of the affected rows,
    // and Phase B re-merges only those rows. Layout (rows, columns, offsets) never changes.
    template<typename GraphT>
    void rebuildTablesSubset(const std::vector<int32_t> &tableIdx,
                             const GraphT &cchGraph,
                             int const *const upWeights,
                             int const *const downWeights,
                             const std::vector<int32_t> &eliminationTree,
                             const HubHierarchy &hierarchy,
                             std::vector<TableData> &tables) const {
        if (tableIdx.empty())
            return;
        const int32_t numHubs = hierarchy.numHubs();
        // The tables' own layout decides the direction count: a symmetric build stores no
        // distIn. (The caller escalates to a full rebuild if the metric's symmetry flips.)
        const bool symmetric = tables[tableIdx.front()].distIn.empty();

        std::vector<int32_t> parentIdx(numHubs);
#pragma omp parallel for schedule(static)
        for (int32_t i = 0; i < numHubs; ++i) {
            const int32_t pu = eliminationTree[hierarchy.rankOfHub(i)];
            parentIdx[i] =
                    pu == INVALID_VERTEX ? -1 : hierarchy.hubIdOfRank(pu);
        }

        // Ancestor closure of every row AND column of the affected tables: the merge reads
        // the labels of both, and a column need not lie on any row's elimination path (it is
        // a separator ancestor, which may be unreachable — the merge then yields infinity).
        std::vector<int32_t> slotOf(numHubs, -1);
        std::vector<int32_t> nodes;
        const auto close = [&](const int32_t id) {
            for (int32_t u = id; u != -1 && slotOf[u] < 0; u = parentIdx[u]) {
                slotOf[u] = 0;
                nodes.push_back(u);
            }
        };
        for (const int32_t ti : tableIdx) {
            for (const auto id : tables[ti].rowIds)
                close(id);
            for (const auto id : tables[ti].colIds)
                close(id);
        }
        std::sort(nodes.begin(), nodes.end());
        for (size_t k = 0; k < nodes.size(); ++k)
            slotOf[nodes[k]] = static_cast<int32_t>(k);
        const int32_t numNodes = static_cast<int32_t>(nodes.size());

        // Hub-subgraph rows for the swept nodes only.
        std::vector<int64_t> tFirst(numNodes + 1, 0);
        for (int32_t k = 0; k < numNodes; ++k) {
            const int32_t u = hierarchy.rankOfHub(nodes[k]);
            tFirst[k + 1] = tFirst[k] + (cchGraph.lastEdge(u) - cchGraph.firstEdge(u));
        }
        std::vector<int32_t> tHead(tFirst[numNodes]), tUpW(tFirst[numNodes]),
                tDownW(tFirst[numNodes]);
#pragma omp parallel for schedule(static)
        for (int32_t k = 0; k < numNodes; ++k) {
            const int32_t u = hierarchy.rankOfHub(nodes[k]);
            int64_t at = tFirst[k];
            FORALL_INCIDENT_EDGES(cchGraph, u, e) {
                const int32_t h = hierarchy.hubIdOfRank(cchGraph.edgeHead(e));
                // Heads outside the closure are never on a swept node's path.
                tHead[at] = h < 0 ? -1 : slotOf[h];
                tUpW[at] = upWeights[e];
                tDownW[at] = downWeights[e];
                ++at;
            }
        }

        // Phase A on the closure (slot-indexed).
        std::vector<int32_t> labelOffset(numNodes + 1, 0);
        for (int32_t k = 0; k < numNodes; ++k) {
            int32_t len = 0;
            for (int32_t u = nodes[k]; u != -1; u = parentIdx[u])
                ++len;
            labelOffset[k + 1] = labelOffset[k] + len;
        }
        std::vector<int32_t> labelSlots(labelOffset[numNodes]);
        std::vector<int32_t> labelUp(labelOffset[numNodes]);
        std::vector<int32_t> labelDown(labelOffset[numNodes]);
#pragma omp parallel
        {
            std::vector<int32_t> distUp(numNodes, FHL_INFTY);
            std::vector<int32_t> distDown(numNodes, FHL_INFTY);
#pragma omp for schedule(dynamic, 32)
            for (int32_t k = 0; k < numNodes; ++k) {
                distUp[k] = 0;
                distDown[k] = 0;
                int32_t at = labelOffset[k];
                for (int32_t u = nodes[k]; u != -1; u = parentIdx[u]) {
                    const int32_t su = slotOf[u];
                    KASSERT(su >= 0);
                    const int32_t du = distUp[su], dd = distDown[su];
                    for (int64_t e = tFirst[su]; e < tFirst[su + 1]; ++e) {
                        const int32_t h = tHead[e];
                        if (h < 0)
                            continue;
                        if (du + tUpW[e] < distUp[h])
                            distUp[h] = du + tUpW[e];
                        if (dd + tDownW[e] < distDown[h])
                            distDown[h] = dd + tDownW[e];
                    }
                    labelSlots[at] = su;
                    labelUp[at] = du;
                    labelDown[at] = dd;
                    ++at;
                }
                for (at = labelOffset[k]; at < labelOffset[k + 1]; ++at) {
                    distUp[labelSlots[at]] = FHL_INFTY;
                    distDown[labelSlots[at]] = FHL_INFTY;
                }
            }
        }

        // Phase B for the affected rows only.
        std::vector<std::pair<int32_t, int32_t>> work;
        for (const int32_t ti : tableIdx)
            for (int32_t i = 0; i < tables[ti].numRows; ++i)
                work.emplace_back(ti, i);
#pragma omp parallel
        {
            std::vector<int32_t> denseUp(numNodes, FHL_INFTY);
            std::vector<int32_t> denseDown(numNodes, FHL_INFTY);
#pragma omp for schedule(dynamic, 16)
            for (size_t w = 0; w < work.size(); ++w) {
                auto &table = tables[work[w].first];
                const int32_t i = work[w].second;
                const int32_t rowSlot = slotOf[table.rowIds[i]];
                for (int32_t k = labelOffset[rowSlot]; k < labelOffset[rowSlot + 1]; ++k) {
                    denseUp[labelSlots[k]] = labelUp[k];
                    denseDown[labelSlots[k]] = labelDown[k];
                }
                const size_t rowOffset = static_cast<size_t>(i) * table.numCols;
                for (int32_t j = 0; j < table.numCols; ++j) {
                    const int32_t colSlot = slotOf[table.colIds[j]];
                    KASSERT(colSlot >= 0);
                    int32_t bestOut = FHL_INFTY, bestIn = FHL_INFTY;
                    for (int32_t k = labelOffset[colSlot]; k < labelOffset[colSlot + 1];
                         ++k) {
                        const int32_t u = labelSlots[k];
                        if (denseUp[u] + labelDown[k] < bestOut)
                            bestOut = denseUp[u] + labelDown[k];
                        if (!symmetric && labelUp[k] + denseDown[u] < bestIn)
                            bestIn = labelUp[k] + denseDown[u];
                    }
                    table.distOut[rowOffset + j] = bestOut;
                    if (!symmetric)
                        table.distIn[rowOffset + j] = bestIn;
                }
                for (int32_t k = labelOffset[rowSlot]; k < labelOffset[rowSlot + 1]; ++k) {
                    denseUp[labelSlots[k]] = FHL_INFTY;
                    denseDown[labelSlots[k]] = FHL_INFTY;
                }
            }
        }
    }

    template<typename GraphT>
    void buildDistanceTables(const std::vector<TableSpec> &specs,
                             const GraphT &cchGraph,
                             int const *const upWeights,
                             int const *const downWeights,
                             const std::vector<int32_t> &eliminationTree,
                             const HubHierarchy &hierarchy,
                             std::vector<TableData> &outTables) const {
        const int32_t numHubs = hierarchy.numHubs();
        outTables.clear();
        outTables.resize(specs.size());

        // Symmetric metric: d(r -> c) == d(c -> r), so the distIn arrays are redundant and
        // are not stored (TableData::inArr() serves reads from distOut instead).
        const bool symmetric =
                std::equal(upWeights, upWeights + cchGraph.numEdges(), downWeights);

        // Phase A: the elimination tree path of a hub consists of hubs only
        // (its elimination ancestors lie in ancestor separators), so one up/down sweep per
        // hub — indexed by hub index, not by rank — yields all labels any
        // table entry can need: labelUp[k] = dup(v, u_k) and labelDown[k] = d(u_k -> v) for
        // the k-th vertex on v's path. Sweeping every node exactly once removes the per-table
        // duplication of ancestor-column sweeps, and the dense scratch shrinks from
        // numVertices to numHubs entries (cache resident).
        Timer phaseTimer;
        // Hub-subgraph CSR with heads as hub indices (the per-edge rank->index
        // lookup is a numVertices-sized random access otherwise).
        std::vector<int64_t> tFirst(numHubs + 1, 0);
        for (int32_t i = 0; i < numHubs; ++i) {
            const int32_t u = hierarchy.rankOfHub(i);
            tFirst[i + 1] = tFirst[i] + (cchGraph.lastEdge(u) - cchGraph.firstEdge(u));
        }
        std::vector<int32_t> tHead(tFirst[numHubs]);
        std::vector<int32_t> tUpW(tFirst[numHubs]);
        std::vector<int32_t> tDownW(tFirst[numHubs]);
#pragma omp parallel for schedule(static)
        for (int32_t i = 0; i < numHubs; ++i) {
            const int32_t u = hierarchy.rankOfHub(i);
            int64_t k = tFirst[i];
            FORALL_INCIDENT_EDGES(cchGraph, u, e) {
                tHead[k] = hierarchy.hubIdOfRank(cchGraph.edgeHead(e));
                tUpW[k] = upWeights[e];
                tDownW[k] = downWeights[e];
                ++k;
            }
        }
        std::vector<int32_t> parentIdx(numHubs); // elimination parent as hub id
#pragma omp parallel for schedule(static)
        for (int32_t i = 0; i < numHubs; ++i) {
            const int32_t pu = eliminationTree[hierarchy.rankOfHub(i)];
            parentIdx[i] = pu == INVALID_VERTEX
                                   ? -1
                                   : hierarchy.hubIdOfRank(pu);
        }
        std::vector<int32_t> labelOffset(numHubs + 1, 0);
        for (int32_t id = 0; id < numHubs; ++id) {
            int32_t len = 0;
            for (int32_t u = id; u != -1; u = parentIdx[u])
                ++len;
            labelOffset[id + 1] = labelOffset[id] + len;
        }
        std::vector<int32_t> labelIds(labelOffset[numHubs]);
        std::vector<int32_t> labelUp(labelOffset[numHubs]);
        // Symmetric metric: the down sweep relaxes the SAME weights as the up sweep, so
        // labelDown would be a byte-for-byte copy of labelUp. Skip it (half the relaxations,
        // half the scratch) and let readers alias labelUp instead.
        std::vector<int32_t> labelDown(symmetric ? 0 : labelOffset[numHubs]);

        const auto sweepAll = [&](auto symTag) {
            constexpr bool SYM = decltype(symTag)::value;
#pragma omp parallel
            {
                std::vector<int32_t> distUp(numHubs, FHL_INFTY);
                std::vector<int32_t> distDown(SYM ? 0 : numHubs, FHL_INFTY);
#pragma omp for schedule(dynamic, 64)
                for (int32_t id = 0; id < numHubs; ++id) {
                    distUp[id] = 0;
                    if constexpr (!SYM)
                        distDown[id] = 0;
                    int32_t k = labelOffset[id];
                    for (int32_t uId = id; uId != -1; uId = parentIdx[uId]) {
                        const int32_t du = distUp[uId];
                        for (int64_t e = tFirst[uId]; e < tFirst[uId + 1]; ++e) {
                            const int32_t hId = tHead[e];
                            if (du + tUpW[e] < distUp[hId])
                                distUp[hId] = du + tUpW[e];
                        }
                        if constexpr (!SYM) {
                            const int32_t dd = distDown[uId];
                            for (int64_t e = tFirst[uId]; e < tFirst[uId + 1]; ++e) {
                                const int32_t hId = tHead[e];
                                if (dd + tDownW[e] < distDown[hId])
                                    distDown[hId] = dd + tDownW[e];
                            }
                            labelDown[k] = dd;
                        }
                        labelIds[k] = uId;
                        labelUp[k] = du;
                        ++k;
                    }
                    KASSERT(k == labelOffset[id + 1]);
                    for (k = labelOffset[id]; k < labelOffset[id + 1]; ++k) {
                        distUp[labelIds[k]] = FHL_INFTY;
                        if constexpr (!SYM)
                            distDown[labelIds[k]] = FHL_INFTY;
                    }
                }
            }
        };
        if (symmetric)
            sweepAll(std::true_type{});
        else
            sweepAll(std::false_type{});
        // Readers below take the backward label from here; on a symmetric metric it is labelUp.
        const int32_t *const labelDownRead = symmetric ? labelUp.data() : labelDown.data();

        std::cout << "[cust]   span1 phase A (sweeps): "
                  << phaseTimer.elapsed<std::chrono::milliseconds>() << " ms\n";
        phaseTimer.restart();
        // Table metadata and storage, plus a flattened (table, row) work list so that small
        // tables do not limit parallelism in Phase B. The serial pass is bookkeeping only: it
        // lays every array of every table out in one pool,
        //   [ all distOut | all distIn (asymmetric only) | all rowIds | all colIds ],
        // via prefix sums. The distance part is then filled in one contiguous parallel sweep,
        // and each table's hub-id lookups run in parallel.
        const size_t numTables = specs.size();
        std::vector<size_t> cellOff(numTables + 1, 0), rowOff(numTables + 1, 0),
                colOff(numTables + 1, 0);
        for (size_t si = 0; si < numTables; ++si) {
            const auto &spec = specs[si];
            auto &table = outTables[si];
            table.sepNode = spec.sepNode;
            table.depthBase = spec.depthBase;
            table.span = spec.span;
            table.prefixLevel = spec.prefixLevel;
            table.prefix = spec.prefix;
            table.numRows = static_cast<int32_t>(spec.rowRanks.size());
            table.numCols = table.numRows + static_cast<int32_t>(spec.ancestorRanks.size());
            cellOff[si + 1] = cellOff[si] + static_cast<size_t>(table.numRows) * table.numCols;
            rowOff[si + 1] = rowOff[si] + table.numRows;
            colOff[si + 1] = colOff[si] + table.numCols;
        }
        const size_t cells = cellOff.back();
        const size_t distSize = symmetric ? cells : 2 * cells; // distIn follows all distOut
        const size_t rowBase = distSize, colBase = distSize + rowOff.back();
        const auto pool = std::make_shared_for_overwrite<int32_t[]>(colBase + colOff.back());
        {
            int32_t *const p = pool.get();
#pragma omp parallel for schedule(static)
            for (int64_t i = 0; i < static_cast<int64_t>(distSize); ++i)
                p[i] = FHL_INFTY;
        }
#pragma omp parallel for schedule(dynamic, 64)
        for (int64_t si = 0; si < static_cast<int64_t>(numTables); ++si) {
            const auto &spec = specs[si];
            auto &table = outTables[si];
            const int32_t numRows = table.numRows, numCols = table.numCols;
            table.rowIds = PooledArray(pool, rowBase + rowOff[si], numRows);
            table.colIds = PooledArray(pool, colBase + colOff[si], numCols);
            for (int32_t i = 0; i < numRows; ++i)
                table.rowIds[i] = hierarchy.hubIdOfRank(spec.rowRanks[i]);
            for (int32_t j = 0; j < numRows; ++j)
                table.colIds[j] = table.rowIds[j];
            for (int32_t j = numRows; j < numCols; ++j)
                table.colIds[j] = hierarchy.hubIdOfRank(spec.ancestorRanks[j - numRows]);
            const size_t n = cellOff[si + 1] - cellOff[si];
            table.distOut = PooledArray(pool, cellOff[si], n);
            if (!symmetric)
                table.distIn = PooledArray(pool, cells + cellOff[si], n);
        }
        std::vector<std::pair<int32_t, int32_t>> work;
        for (size_t si = 0; si < specs.size(); ++si) {   // same (table, row) order as before
            for (int32_t i = 0; i < outTables[si].numRows; ++i)
                work.emplace_back(static_cast<int32_t>(si), i);
        }

        // Phase B: per row, scatter the row's labels into a dense (by hub index) array
        // and merge every column's sparse labels against it:
        //   distOut(r, c) = min over u on path(c) of  dup(r, u) + d(u -> c)
        //   distIn(r, c)  = min over u on path(c) of  dup(c, u) + d(u -> r)
        // Vertices u outside path(r) contribute infinity and never win.
        std::cout << "[cust]   span1 alloc+worklist: "
                  << phaseTimer.elapsed<std::chrono::milliseconds>() << " ms\n";
        phaseTimer.restart();
        // Rows and columns of a span-1 table live in the same band, so a column is almost
        // always an ancestor of the row and its label list is then a SUFFIX of the row's
        // (both are bottom-up root paths). One comparison detects that and gives the
        // alignment arithmetically; only the rare crossing case needs the search. Either way
        // the merge is a single sequential pass -- no dense scatter, no clear pass.
        const auto mergeAll = [&](auto symTag) {
            constexpr bool SYM = decltype(symTag)::value;
#pragma omp parallel for schedule(dynamic, 16)
            for (size_t w = 0; w < work.size(); ++w) {
                auto &table = outTables[work[w].first];
                const int32_t i = work[w].second;
                const int32_t rowId = table.rowIds[i];
                const int32_t rs = labelOffset[rowId], re = labelOffset[rowId + 1];
                const size_t rowOffset = static_cast<size_t>(i) * table.numCols;
                for (int32_t j = 0; j < table.numCols; ++j) {
                    const int32_t colId = table.colIds[j];
                    const int32_t cs = labelOffset[colId], ce = labelOffset[colId + 1];
                    const int32_t clen = ce - cs;
                    int32_t ra, cb, n;
                    if (clen <= re - rs && labelIds[re - clen] == labelIds[cs]) {
                        ra = re - clen;
                        cb = cs;
                        n = clen;
                    } else { // paths cross instead of nest: find where they diverge
                        int32_t a = re - 1, b = ce - 1;
                        while (a >= rs && b >= cs && labelIds[a] == labelIds[b]) {
                            --a;
                            --b;
                        }
                        ra = a + 1;
                        cb = b + 1;
                        n = re - ra;
                    }
                    int32_t bestOut = FHL_INFTY;
                    int32_t bestIn = FHL_INFTY;
                    for (int32_t k = 0; k < n; ++k) {
                        const int32_t out = labelUp[ra + k] + labelDownRead[cb + k];
                        if (out < bestOut)
                            bestOut = out;
                        if constexpr (!SYM) {
                            const int32_t in = labelUp[cb + k] + labelDown[ra + k];
                            if (in < bestIn)
                                bestIn = in;
                        }
                    }
                    table.distOut[rowOffset + j] = bestOut;
                    if constexpr (!SYM)
                        table.distIn[rowOffset + j] = bestIn;
                }
            }
        };
        if (symmetric)
            mergeAll(std::true_type{});
        else
            mergeAll(std::false_type{});
        std::cout << "[cust]   span1 phase B (merges): "
                  << phaseTimer.elapsed<std::chrono::milliseconds>() << " ms\n";
    }

private:
    struct SepNodeInfo {
        int32_t level = 0;
        uint64_t prefix = 0;
    };

    static void buildSepNodeInfo(const SeparatorDecomposition &sd,
                                 const int node,
                                 const int depth,
                                 const uint64_t prefix,
                                 std::vector<SepNodeInfo> &info) {
        // depth starts at 1 for root; level is depth - 1 (root level = 0).
        info[node].level = depth - 1;
        info[node].prefix = prefix;

        const int left = sd.leftChild(node);
        if (left == 0)
            return;

        buildSepNodeInfo(sd, left, depth + 1, prefix, info);
        const int right = sd.rightSibling(left);
        if (right != 0) {
            const uint64_t rightPrefix = prefix | (1ULL << (depth - 1));
            buildSepNodeInfo(sd, right, depth + 1, rightPrefix, info);
        }
    }

    static bool inSubtree(const HubHierarchy &hierarchy,
                          const int32_t v,
                          const int32_t rootLevel,
                          const uint64_t rootPrefix) {
        if (rootLevel <= 0)
            return true;
        const uint64_t mask = (1ULL << rootLevel) - 1;
        return (hierarchy.getPackedSideId(v) & mask) == rootPrefix;
    }

    static void collectHubsInSubtree(const HubHierarchy &hierarchy,
                                        const int32_t rootLevel,
                                        const uint64_t rootPrefix,
                                        const int32_t depthBase,
                                        const int32_t span,
                                        const int32_t thresh,
                                        std::vector<int32_t> &out) {
        const int32_t upper = std::min(depthBase + span, thresh);
        const int32_t numHubs = hierarchy.numHubs();
        for (int32_t id = 0; id < numHubs; ++id) {
            const int32_t v = hierarchy.rankOfHub(id);
            const int32_t level = hierarchy.getVertexLevel(v);
            if (level < depthBase || level >= upper)
                continue;
            if (!inSubtree(hierarchy, v, rootLevel, rootPrefix))
                continue;
            out.push_back(v);
        }
    }

    // Collects every hub above the band (level < depthBase) that is an ancestor of at
    // least one row of the table: nodes on the path from the root to the subtree root (level <
    // prefixLevel) as well as nodes inside the subtree at intermediate levels [prefixLevel,
    // depthBase). Both are characterized by the side id prefix matching the subtree root's
    // prefix up to min(level(v), prefixLevel) bits.
    static void collectAncestorColumns(const HubHierarchy &hierarchy,
                                       const uint64_t rootPrefix,
                                       const int32_t prefixLevel,
                                       const int32_t depthBase,
                                       std::vector<int32_t> &out) {
        const int32_t numHubs = hierarchy.numHubs();
        for (int32_t id = 0; id < numHubs; ++id) {
            const int32_t v = hierarchy.rankOfHub(id);
            const int32_t level = hierarchy.getVertexLevel(v);
            if (level >= depthBase)
                continue;
            const int32_t bits = std::min(level, prefixLevel);
            const uint64_t mask = bits > 0 ? (1ULL << bits) - 1 : 0;
            if (((hierarchy.getPackedSideId(v) ^ rootPrefix) & mask) != 0)
                continue;
            out.push_back(v);
        }
    }

    static void buildTablesForSubtree(const SeparatorDecomposition &sd,
                                      const int32_t sepNode,
                                      const int32_t depthBase,
                                      const int32_t span,
                                      const bool fixedSpan,
                                      const int32_t thresh,
                                      const HubHierarchy &hierarchy,
                                      const std::vector<SepNodeInfo> &info,
                                      std::vector<TableSpec> &outSpecs) {
        if (depthBase >= thresh)
            return;

        TableSpec spec;
        spec.sepNode = sepNode;
        spec.depthBase = depthBase;
        spec.span = span;
        spec.prefixLevel = info[sepNode].level;
        spec.prefix = info[sepNode].prefix;
        collectHubsInSubtree(hierarchy, info[sepNode].level, info[sepNode].prefix,
                                depthBase, span, thresh, spec.rowRanks);
        collectAncestorColumns(hierarchy, info[sepNode].prefix, info[sepNode].level,
                               depthBase, spec.ancestorRanks);
        outSpecs.push_back(std::move(spec));

        const int32_t nextDepthBase = depthBase + span;
        if (nextDepthBase >= thresh)
            return;

        const int32_t nextSpan = fixedSpan ? span : span * 2;

        const int left = sd.leftChild(sepNode);
        if (left == 0)
            return;
        buildTablesForSubtree(sd, left, nextDepthBase, nextSpan, fixedSpan, thresh, hierarchy,
                              info, outSpecs);
        const int right = sd.rightSibling(left);
        if (right != 0) {
            buildTablesForSubtree(sd, right, nextDepthBase, nextSpan, fixedSpan, thresh,
                                  hierarchy, info, outSpecs);
        }
    }

public:
};
