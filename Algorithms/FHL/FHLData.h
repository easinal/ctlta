#pragma once

#include <algorithm>
#include <cstdint>
#include <vector>

#ifdef __AVX2__
#include <immintrin.h>
#endif

#include "Algorithms/FHL/core/HubAccessData.h"
#include "Algorithms/FHL/core/HubTableBuilder.h"
#include "Algorithms/FHL/schemes/FHLInRegionHubs.h"
#include "Algorithms/FHL/schemes/FHLOverlay.h"
#include "Algorithms/FHL/schemes/FHLSeedRectangle.h"
#include "Algorithms/FHL/schemes/FHLTopLabels.h"

// Holds the access-hub scaffolding plus the hub distance tables.
// The base HubAccessData is created without the full hub distance table; the hierarchical
// tables replace it entirely. All index structures are flat arrays over hub indices,
// so every lookup is O(1) without hashing.
class FHLData {
public:
    explicit FHLData(const int numHubs, const int numVertices)
            : baseData(numVertices), // no T x T table: FHL answers from seeds and rectangles
              numHubs(numHubs) {}

    HubAccessData &getBaseData() { return baseData; }
    const HubAccessData &getBaseData() const { return baseData; }

    void setTables(std::vector<HubTableBuilder::TableData> &&newTables) {
        tables = std::move(newTables);
        buildIndexMaps();
        buildBandIndex();
    }

    int32_t numTables() const { return static_cast<int32_t>(tables.size()); }

    const HubTableBuilder::TableData &getTable(const int32_t ti) const {
        return tables[ti];
    }

    // Returns the column of hub index `id` in table `ti`, or -1 if absent. O(1) via
    // the flat maps; binary search over the sorted column ids when the flat maps would be too
    // large (many small tables, e.g. span=1 subtree mode).
    int32_t colIndexIn(const int32_t ti, const int32_t id) const {
        if (flatCols)
            return colOf[ti][id];
        const auto &cols = tables[ti].colIds;
        const auto it = std::lower_bound(cols.begin(), cols.end(), id);
        if (it == cols.end() || *it != id)
            return -1;
        return static_cast<int32_t>(it - cols.begin());
    }

    // Returns the exact distance from hub `from` to hub `to` (both given as
    // hub indices) using table lookups only. The typical case is a single lookup: if
    // one node lies on the separator path of the other, the deeper node's table stores the
    // distance directly. Otherwise (both endpoints are hubs in unrelated subtrees),
    // the distance is assembled over the common path ancestors of the two tables, which always
    // contain the maximum-rank vertex of a shortest path.
    int32_t getDirectedDistance(const int32_t from, const int32_t to) const {
        if (from == to)
            return 0;

        const int32_t tf = tableOf[from];
        const auto &fromTable = tables[tf];
        const size_t fromRow = static_cast<size_t>(rowOf[from]) * fromTable.numCols;
        const int32_t cTo = colIndexIn(tf, to);
        if (cTo >= 0)
            return fromTable.distOut[fromRow + cTo];

        const int32_t tt = tableOf[to];
        const auto &toTable = tables[tt];
        const size_t toRow = static_cast<size_t>(rowOf[to]) * toTable.numCols;
        const int32_t cFrom = colIndexIn(tt, from);
        if (cFrom >= 0)
            return toTable.inArr()[toRow + cFrom];

        int32_t best = FHL_INFTY;
        for (int32_t j = 0; j < fromTable.numCols; ++j) {
            const int32_t cu = colIndexIn(tt, fromTable.colIds[j]);
            if (cu < 0)
                continue;
            const int32_t d = fromTable.distOut[fromRow + j] + toTable.inArr()[toRow + cu];
            if (d < best)
                best = d;
        }
        return best;
    }

    // Distance between two ancestor-related hubs (one lies on the root path of the
    // other; smaller index means deeper). Single unconditional lookup — use
    // getDirectedDistance when the relation is unknown.
    int32_t getComparableDistance(const int32_t from, const int32_t to) const {
        if (from == to)
            return 0;
        if (from < to) { // from is deeper; to is a column of from's table
            const int32_t tf = tableOf[from];
            const auto &t = tables[tf];
            const int32_t c = colIndexIn(tf, to);
            KASSERT(c >= 0);
            return t.distOut[static_cast<size_t>(rowOf[from]) * t.numCols + c];
        }
        // to is deeper; from is a column of to's table
        const int32_t tt = tableOf[to];
        const auto &t = tables[tt];
        const int32_t c = colIndexIn(tt, from);
        KASSERT(c >= 0);
        return t.inArr()[static_cast<size_t>(rowOf[to]) * t.numCols + c];
    }

    // Region-scheme structures (see FHLRegions.h / FHLStore.h).
    fhl::Values ladder;

    // Frees the hierarchical tables and their index structures (hub mode keeps only labels).
    // Mutable table access for the super-region rebuild (values change in place; the
    // layout and index maps stay valid).
    std::vector<HubTableBuilder::TableData> &mutableTables() {
        return tables;
    }

    void releaseTables() {
        tables.clear();
        tables.shrink_to_fit();
        tableOf.clear();
        tableOf.shrink_to_fit();
        rowOf.clear();
        rowOf.shrink_to_fit();
        colOf.clear();
        colOf.shrink_to_fit();
        bandOfLevel.clear();
        bandMask.clear();
        bandTableAtPrefix.clear();
    }

    uint64_t sizeInBytes() const {
        return sizeof(FHLData) + baseData.sizeInBytes() + sizeTablesInBytes() +
               ladder.sizeInBytes();
    }

    uint64_t sizeTablesInBytes() const {
        uint64_t size = tableOf.size() * sizeof(int32_t) + rowOf.size() * sizeof(int32_t);
        for (const auto &table : tables)
            size += table.sizeInBytes();
        for (const auto &cols : colOf)
            size += cols.size() * sizeof(int32_t);
        size += bandOfLevel.size() * sizeof(int32_t) + bandMask.size() * sizeof(uint64_t);
        for (const auto &prefixTables : bandTableAtPrefix)
            size += prefixTables.size() * sizeof(int32_t);
        return size;
    }

private:
    void buildIndexMaps() {
        tableOf.assign(numHubs, -1);
        rowOf.assign(numHubs, -1);
        for (int32_t ti = 0; ti < static_cast<int32_t>(tables.size()); ++ti) {
            const auto &table = tables[ti];
            for (int32_t i = 0; i < table.numRows; ++i) {
                tableOf[table.rowIds[i]] = ti;
                rowOf[table.rowIds[i]] = i;
            }
        }
        // Flat per-table column maps cost numTables x numHubs ints; with many small
        // tables (small spans) fall back to binary search / merges over the sorted column ids.
        const uint64_t flatBytes =
                tables.size() * static_cast<uint64_t>(numHubs) * sizeof(int32_t);
        flatCols = flatBytes <= FLAT_COL_MAP_BUDGET;
        colOf.clear();
        if (flatCols) {
            colOf.assign(tables.size(), {});
            for (int32_t ti = 0; ti < static_cast<int32_t>(tables.size()); ++ti) {
                const auto &table = tables[ti];
                auto &cols = colOf[ti];
                cols.assign(numHubs, -1);
                for (int32_t j = 0; j < table.numCols; ++j)
                    cols[table.colIds[j]] = j;
            }
        }
    }

    // Groups tables by band (depthBase) and indexes each band's tables by subtree prefix, so
    // the table responsible for a given (level, side id) is found in O(1).
    void buildBandIndex() {
        bandOfLevel.clear();
        bandMask.clear();
        bandTableAtPrefix.clear();

        std::vector<int32_t> bandDepthBases;
        for (const auto &table : tables)
            if (std::find(bandDepthBases.begin(), bandDepthBases.end(), table.depthBase) ==
                bandDepthBases.end())
                bandDepthBases.push_back(table.depthBase);
        std::sort(bandDepthBases.begin(), bandDepthBases.end());

        bandMask.assign(bandDepthBases.size(), 0);
        bandTableAtPrefix.assign(bandDepthBases.size(), {});
        int32_t maxLevel = 0;
        for (size_t b = 0; b < bandDepthBases.size(); ++b) {
            int32_t prefixLevel = 0;
            int32_t span = 0;
            for (const auto &table : tables)
                if (table.depthBase == bandDepthBases[b]) {
                    prefixLevel = table.prefixLevel;
                    span = table.span;
                }
            bandMask[b] = prefixLevel > 0 ? (1ULL << prefixLevel) - 1 : 0;
            bandTableAtPrefix[b].assign(bandMask[b] + 1, -1);
            maxLevel = std::max(maxLevel, bandDepthBases[b] + span);
        }

        bandOfLevel.assign(maxLevel, -1);
        for (size_t b = 0; b < bandDepthBases.size(); ++b) {
            int32_t span = 0;
            for (const auto &table : tables)
                if (table.depthBase == bandDepthBases[b])
                    span = table.span;
            for (int32_t l = bandDepthBases[b];
                 l < std::min(bandDepthBases[b] + span, maxLevel); ++l)
                bandOfLevel[l] = static_cast<int32_t>(b);
        }

        for (int32_t ti = 0; ti < static_cast<int32_t>(tables.size()); ++ti) {
            const auto &table = tables[ti];
            const auto b = std::find(bandDepthBases.begin(), bandDepthBases.end(),
                                     table.depthBase) - bandDepthBases.begin();
            bandTableAtPrefix[b][table.prefix & bandMask[b]] = ti;
        }
    }

    static constexpr uint64_t FLAT_COL_MAP_BUDGET = 256 * 1024 * 1024;

    HubAccessData baseData;
    int32_t numHubs;

    std::vector<HubTableBuilder::TableData> tables;
    std::vector<int32_t> tableOf; // hub id -> table owning its row
    std::vector<int32_t> rowOf;   // hub id -> row within its table
    bool flatCols = true;
    std::vector<std::vector<int32_t>> colOf; // per table: hub id -> column, or -1

    std::vector<int32_t> bandOfLevel;               // separator level -> band index
    std::vector<uint64_t> bandMask;                 // band -> mask for the subtree prefix
    std::vector<std::vector<int32_t>> bandTableAtPrefix; // band -> subtree prefix -> table
};
