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

// METRIC-DEPENDENT values for every scheme, in one place because a query reads across them:
// seeds, the pruned rectangle CSR with its skip bounds, per-hub rows, the boundary-point
// cliques, the
// in-region rows, and the materialized top labels. Customization refills exactly these.
namespace fhl {

// Move-only int32 array without value-initialization (avoids serial gigabyte memsets;
// every slot is written by the region-parallel builders before it is read).
struct RawArray {
    std::unique_ptr<int32_t[]> p;
    int64_t n = 0;
    void allocateForOverwrite(const int64_t m) {
        p = std::make_unique_for_overwrite<int32_t[]>(m);
        n = m;
    }
    void assignFill(const int64_t m, const int32_t v) {
        allocateForOverwrite(m);
        std::fill(p.get(), p.get() + m, v);
    }
    void clear() {
        p.reset();
        n = 0;
    }
    bool empty() const { return n == 0; }
    int64_t size() const { return n; }
    int32_t *data() { return p.get(); }
    const int32_t *data() const { return p.get(); }
};

// ------------------------------------------------------------------ metric-dependent part
struct Values {
    int32_t mode = 0; // 0 off, 1 flat, 2 overlay
    bool symmetric = false;

    // Seeds: exact distances to the own leaf boundary (aligned to the region's B list).
    std::vector<int64_t> seedPos; // by rank; empty rows for hubs
    RawArray seedFwd;
    RawArray seedBwd; // empty on symmetric metrics

    // Leaf cliques: pairwise boundary distances (global in flat mode, restricted in
    // overlay mode), aligned to the leaf B lists. Serve the local pair queries and the
    // flat CSR pruning.
    struct CliqueLevel {
        std::vector<int64_t> off; // by region
        std::vector<int32_t> M;   // d(x -> y)
    };
    CliqueLevel leafClique;

    // Flat mode: per-column CSR of the pruned rectangles plus column minima, hub
    // value rows, and the optional per-vertex top labels.
    struct FlatCsr {
        // A survivor's row index is a position within its OWN region's boundary point list, not a
        // global hub id, so one byte covers it -- and the narrow array fits more of the
        // inner scan's operands per cache line. Widen the alias and MAX_ROWS together.
        using RowIdx = uint8_t;
        static constexpr int32_t MAX_ROWS = 255;
        std::vector<int64_t> ptr; // size = total flat columns + 1
        std::vector<RowIdx> row;  // surviving row index within the region's B list
        std::vector<int32_t> val;
        // Relocatable mode (super-region rebuilds): explicit per-column end offsets, so
        // one region's segment can be rewritten (or moved to the tail) without shifting
        // every later column. Empty = the plain monotone layout, where end == ptr + 1.
        std::vector<int64_t> colEnd;
        std::vector<int64_t> segStart; // per leaf region: segment start and capacity
        std::vector<int64_t> segCap;
        bool active() const { return !ptr.empty(); }
        // One pointer select per query side; the inner loop stays branch-free.
        const int64_t *endArray() const {
            return colEnd.empty() ? ptr.data() + 1 : colEnd.data();
        }
        uint64_t sizeInBytes() const {
            return (ptr.size() + colEnd.size() + segStart.size() + segCap.size()) * 8 +
                   row.size() * sizeof(RowIdx) + val.size() * 4;
        }
    };
    FlatCsr csrOut;
    FlatCsr csrIn; // empty on symmetric metrics
    std::vector<int32_t> colMinOut;
    std::vector<int32_t> colMinIn;
    // Block minima over colMin (blocks of COL_BLOCK columns within a region). The far scan
    // skips almost every column one by one; a block bound skips a whole block at once.
    // Admissible: the block min is a lower bound for every column in it.
    static constexpr int32_t COL_BLOCK = 32;
    std::vector<int64_t> colBlockFirst; // per leaf region
    std::vector<int32_t> colBlockMinOut, colBlockMinIn;
    std::vector<int32_t> tRowFwd; // hub value rows, aligned to the column runs
    std::vector<int32_t> tRowBwd;
    int32_t topLabLevels = 0;
    std::vector<int64_t> topLabPos; // by rank
    RawArray topLabFwd;
    RawArray topLabBwd;

    // Overlay mode: restricted meet cliques and hub seeds.
    // In-region rows: d_R(v -> col) and d_R(col -> v) for the region's top-K columns.
    int32_t locLevels = 0;
    std::vector<int64_t> locPos; // by rank
    RawArray locFwd;
    RawArray locBwd; // empty on symmetric metrics

    std::vector<CliqueLevel> meetClique;
    std::vector<int64_t> tSeedPos;
    std::vector<int32_t> tSeedFwd;
    std::vector<int32_t> tSeedBwd;

    bool active() const { return mode != 0; }

    uint64_t sizeInBytes() const {
        uint64_t s = sizeof(Values) + seedPos.size() * 8 +
                     (seedFwd.size() + seedBwd.size()) * 4;
        s += leafClique.off.size() * 8 + leafClique.M.size() * 4;
        s += csrOut.sizeInBytes() + csrIn.sizeInBytes();
        s += (colMinOut.size() + colMinIn.size() + tRowFwd.size() + tRowBwd.size()) * 4;
        s += colBlockFirst.size() * 8 +
             (colBlockMinOut.size() + colBlockMinIn.size()) * 4;
        s += topLabPos.size() * 8 + (topLabFwd.size() + topLabBwd.size()) * 4;
        for (const auto &cl : meetClique)
            s += cl.off.size() * 8 + cl.M.size() * 4;
        s += tSeedPos.size() * 8 + (tSeedFwd.size() + tSeedBwd.size()) * 4;
        s += locPos.size() * 8 + (locFwd.size() + locBwd.size()) * 4;
        return s;
    }
};


} // namespace fhl
