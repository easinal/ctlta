#pragma once

// FHL's own hierarchy over the separator decomposition. Split from CTNR's because
// the two want different things from it: classic TNR needs a uniform top-k-levels
// hub-node set to fill a T x T table, while FHL needs the SIZE CUT (leaf regions of a
// target size, an optional inner cut) and a dense hub numbering for its per-hub rows.

#include <kassert/kassert.hpp>
#include <stack>

#include "Algorithms/CCH/CCH.h"
#include "DataStructures/Partitioning/SeparatorDecompositionWalk.h"
#include "DataStructures/Partitioning/SeparatorTree.h"
#include "Algorithms/CCH/CCHMetric.h"

#include "Tools/Constants.h"
#include "Tools/Timer.h"
#include "Algorithms/FHL/core/FHLConstants.h"

class HubHierarchy {

    using SideId = uint64_t;
    using HubId = fhl::HubId;
    using Level = fhl::Level;

    struct Identity {
        int operator[](const int x) const {
            return x;
        }
    };

public:

    HubHierarchy() = default;

    // Uniform mode: every separator vertex in the top `newHubLevelThreshold` levels is a hub.
    // This is the rule overlay cuts on; FHL itself prefers preprocessSizeCut below, because a
    // fixed depth gives wildly uneven region sizes and the seed count per vertex follows the
    // region's boundary.
    template<typename InputGraphT, typename PermuteVertexIdT = Identity>
    void preprocess(const InputGraphT &inputGraph, const int newHubLevelThreshold,
                    const SeparatorDecomposition &sepDecomp,
                    const PermuteVertexIdT &permuteVertexId = {}) {
        std::cout << "Depth of sepDecomp is " << sepdecomp::depth(sepDecomp) << std::endl;
        resetFor(sepDecomp, inputGraph.numVertices(), /*forSizeCut=*/false, 0);
        hubLevelThreshold = newHubLevelThreshold;
        walkDecomposition(sepDecomp, permuteVertexId,
                          [&](const int, const Level depth) { return depth >= hubLevelThreshold; },
                          [](const int) { return false; });
        finishHubs();
        KASSERT(std::all_of(packedSideIds.begin(), packedSideIds.end(),
                            [](const uint64_t &id) { return id != static_cast<uint64_t>(-1); }));
        std::cout << "FHL: Selected " << hubRanks.size() << " hubs from top "
                  << (int) hubLevelThreshold << " levels" << std::endl;
        std::cout << "Total number of vertices: " << packedSideIds.size() << std::endl;
    }

    // Size-cut mode: leaf regions are maximal subtrees holding at most targetRegionSize
    // vertices (a childless subtree always qualifies); every separator vertex strictly ABOVE a
    // leaf region is a hub. targetInnerSize > 0 makes a second, finer cut inside each leaf
    // region -- the vertices left above that cut are the region's own hub set.
    template<typename InputGraphT, typename PermuteVertexIdT = Identity>
    void preprocessSizeCut(const InputGraphT &inputGraph, const int32_t targetRegionSize,
                           const int32_t targetInnerSize,
                           const SeparatorDecomposition &sd,
                           const PermuteVertexIdT &permuteVertexId = {}) {
        resetFor(sd, inputGraph.numVertices(), /*forSizeCut=*/true, targetInnerSize);

        // Subtree vertex counts (own separator + descendants), bottom-up.
        std::vector<int64_t> subSize(sd.tree.size(), 0);
        for (size_t nd = 0; nd < sd.tree.size(); ++nd)
            subSize[nd] = sd.lastSeparatorVertex(nd) - sd.firstSeparatorVertex(nd);
        sepdecomp::forEachNodeInDfsOrder(
                sd, [](const int, const int) {},
                [&](const int child, const int parent) { subSize[parent] += subSize[child]; });

        const auto smallEnough = [&](const int node, const int32_t limit) {
            return subSize[node] > 0 && (subSize[node] <= limit || sd.leftChild(node) == 0);
        };
        walkDecomposition(sd, permuteVertexId,
                          [&](const int node, const Level) { return smallEnough(node, targetRegionSize); },
                          [&](const int node) { return smallEnough(node, targetInnerSize); });
        finishHubs();
        if (maxCutDepth > 63) // packed side ids hold one bit per level
            throw std::invalid_argument("size-cut: cut depth exceeds 63 levels; "
                                        "raise the target region size");
        // Not a valid level test in this mode; kept as the level-stride bound.
        hubLevelThreshold = static_cast<Level>(std::min(126, maxCutDepth + 1));
        std::cout << "FHL size-cut: " << numRegionsSizeCut << " leaf regions (target size "
                  << targetRegionSize << "), " << hubRanks.size()
                  << " hubs, max cut depth " << maxCutDepth << std::endl;
    }

    bool isSizeCut() const { return sizeCut; }
    int32_t numLeafRegions() const { return numRegionsSizeCut; }
    int32_t getMaxCutDepth() const { return maxCutDepth; }
    const std::vector<int32_t> &getLeafRegionRaw() const { return leafRegionRaw; }
    // Sub-region of a vertex inside its leaf region; -1 marks the region's local hub set
    // (the vertices left above the inner size cut).
    bool hasInnerCut() const { return innerCutUsed; }
    const std::vector<int32_t> &getInnerRegionRawVec() const { return innerRegionRaw; }

    size_t numVertices() const {
        return packedSideIds.size();
    }

    Level getHubLevelThreshold() const {
        return hubLevelThreshold;
    }

    // Given the rank of a vertex in the CCH-order, returns its level in the separator hierarchy.
    inline Level getVertexLevel(const int32_t &v) const {
        KASSERT(v >= 0 && v < vertexLevel.size());
        return vertexLevel[v];
    }

    // Given the ranks of two vertices in the CCH-order, returns the level of their lowest common ancestor
    // in the separator hierarchy.
    Level getLevelOfLowestCommonAncestor(const int32_t &s, const int32_t &t) const {

        const Level minInputLevel = std::min(getVertexLevel(s), getVertexLevel(t));

        // XOR packed side IDs to find out lowest common level in separator hierarchy.
        const Level l = static_cast<Level>(lowestOneBit(packedSideIds[s] ^ packedSideIds[t]));

        if (l >= 0)
            return std::min(l, minInputLevel);

        // If packed side IDs of s and t are exactly the same, the branch of s subsumes the branch of t or vice
        // versa. In this case, the lowest common ancestor is the lower one of the two vertices.
        return minInputLevel;
    }

    int32_t numHubs() const {
        return static_cast<int32_t>(hubRanks.size());
    }

    // Rank -> dense hub id with no n-sized map: a bitmap of the hub ranks plus a
    // popcount prefix per 512 ranks yields the RANK-order position, and selectToId turns
    // that into the (level desc, rank asc) id the rest of the code sorts by. The map used
    // to be one int32 per vertex, almost all of them -1.
    static constexpr size_t RANK_BLOCK_WORDS = 8; // 512 ranks per block
    void buildHubIndex(const int32_t numVerts) {
        const size_t words = (static_cast<size_t>(numVerts) + 63) / 64;
        hubBits.assign(words, 0);
        std::vector<std::pair<int32_t, HubId>> byRank(hubRanks.size());
        for (size_t i = 0; i < hubRanks.size(); ++i)
            byRank[i] = {hubRanks[i], static_cast<HubId>(i)};
        std::sort(byRank.begin(), byRank.end());
        selectToId.resize(byRank.size());
        for (size_t p = 0; p < byRank.size(); ++p) {
            selectToId[p] = byRank[p].second;
            hubBits[static_cast<size_t>(byRank[p].first) >> 6] |=
                    1ULL << (byRank[p].first & 63);
        }
        const size_t blocks = (words + RANK_BLOCK_WORDS - 1) / RANK_BLOCK_WORDS;
        blockRank.assign(blocks + 1, 0);
        int32_t running = 0;
        for (size_t b = 0; b < blocks; ++b) {
            blockRank[b] = running;
            for (size_t w = b * RANK_BLOCK_WORDS;
                 w < std::min(words, (b + 1) * RANK_BLOCK_WORDS); ++w)
                running += __builtin_popcountll(hubBits[w]);
        }
        blockRank[blocks] = running;
    }

    // Dense id, or -1 when v is not a hub.
    int32_t hubIdOrInvalid(const int32_t &v) const {
        return isHub(v) ? static_cast<int32_t>(hubIdOfRank(v)) : -1;
    }

    bool isHub(const int32_t &v) const {
        KASSERT(v >= 0 && v < vertexLevel.size());
        return (hubBits[static_cast<size_t>(v) >> 6] >> (v & 63)) & 1ULL;
    }

    int32_t rankOfHub(const HubId &index) const {
        KASSERT(index >= 0 && index < hubRanks.size());
        return hubRanks[index];
    }

    HubId hubIdOfRank(const int32_t &v) const {
        KASSERT(isHub(v));
        const size_t w = static_cast<size_t>(v) >> 6;
        const size_t b = w / RANK_BLOCK_WORDS;
        int32_t p = blockRank[b];
        for (size_t i = b * RANK_BLOCK_WORDS; i < w; ++i)
            p += __builtin_popcountll(hubBits[i]);
        p += __builtin_popcountll(hubBits[w] & ((1ULL << (v & 63)) - 1));
        return selectToId[p];
    }

    uint64_t sizeInBytes() const {
        return sizeof(HubHierarchy) +
               vertexLevel.size() * sizeof(decltype(vertexLevel)::value_type) +
               packedSideIds.size() * sizeof(decltype(packedSideIds)::value_type) +
               hubRanks.capacity() * sizeof(decltype(hubRanks)::value_type) +
               hubBits.size() * 8 + blockRank.size() * 4 +
               selectToId.size() * sizeof(decltype(selectToId)::value_type);
    }

    uint64_t getPackedSideId(const int32_t &v) const {
        KASSERT(v >= 0 && v < packedSideIds.size());
        return packedSideIds[v];
    }

private:
    void resetFor(const SeparatorDecomposition &sd, const size_t n, const bool forSizeCut,
                  const int32_t targetInnerSize) {
        if (!sepdecomp::hasStrictDissectionStructure(sd))
            throw std::invalid_argument("HubHierarchy requires strict dissection "
                                        "structure of separator decomposition.");
        packedSideIds.assign(n, static_cast<uint64_t>(-1));
        vertexLevel.assign(n, -1);
        hubRanks.clear();
        sizeCut = forSizeCut;
        innerCutUsed = targetInnerSize > 0;
        maxCutDepth = 0;
        numRegionsSizeCut = 0;
        if (forSizeCut) {
            leafRegionRaw.assign(n, -1);
            innerRegionRaw.assign(n, -1);
        }
    }

    // ONE recursive walk for both cut rules. It assigns every vertex its level and side-id
    // prefix; `openRegion` decides where the hub zone ends (everything strictly above an
    // opened region is a hub) and `openInner` does the same one scale down, inside a region.
    // Everything that flows DOWN the tree is a parameter, so there are no shadow stacks to
    // keep in sync. Recursion depth is the decomposition depth, bounded by the 63-level
    // guard in the size-cut caller.
    template<typename PermuteVertexIdT, typename OpenRegionT, typename OpenInnerT>
    void walkDecomposition(const SeparatorDecomposition &sd, const PermuteVertexIdT &permute,
                           OpenRegionT openRegion, OpenInnerT openInner) {
        int32_t numRegions = 0, numInner = 0;
        const auto walk = [&](auto &&self, const int node, const Level depth,
                              const SideId prefix, int32_t reg, int32_t inner) -> void {
            if (reg < 0 && openRegion(node, depth)) {
                reg = numRegions++;
                maxCutDepth = std::max<int32_t>(maxCutDepth, depth);
            }
            if (reg >= 0 && inner < 0 && openInner(node))
                inner = numInner++;
            for (auto vSd = sd.lastSeparatorVertex(node) - 1;
                 vSd >= sd.firstSeparatorVertex(node); --vSd) {
                const auto v = permute[vSd];
                vertexLevel[v] = depth;
                packedSideIds[v] = prefix;
                if (reg < 0) {
                    hubRanks.push_back(v);
                } else if (sizeCut) {
                    leafRegionRaw[v] = reg;
                    innerRegionRaw[v] = inner; // -1 = the region's own local hubs
                }
            }
            const int left = sd.leftChild(node);
            if (left == 0)
                return;
            self(self, left, depth + 1, prefix, reg, inner);
            const int right = sd.rightSibling(left);
            if (right != 0)
                self(self, right, depth + 1, prefix | (SideId{1} << depth), reg, inner);
        };
        walk(walk, 0, 0, 0, -1, -1);
        if (sizeCut)
            numRegionsSizeCut = numRegions;
    }

    // Hubs get their dense ids in (level desc, rank asc) order; the rank -> id index follows.
    void finishHubs() {
        if (hubRanks.size() > std::numeric_limits<HubId>::max())
            throw std::invalid_argument(
                    "Number of hubs (" + std::to_string(hubRanks.size()) +
                    ") exceeds HubId range; raise the cut size or lower the level threshold.");
        std::sort(hubRanks.begin(), hubRanks.end(), [&](const int32_t a, const int32_t b) {
            return vertexLevel[a] > vertexLevel[b] ||
                   (vertexLevel[a] == vertexLevel[b] && a < b);
        });
        buildHubIndex(static_cast<int32_t>(vertexLevel.size()));
    }


    std::vector<SideId> packedSideIds; // store which side each vertex is on in each level of separator hierarchy
    std::vector<Level> vertexLevel; // Map each vertex to its level in the SD

    // The small subset of vertices in the hubLevelThreshold highest levels make up the hubs.
    // Every hub gets an internal index in [0, numHubs-1] for distance table lookup.
    // We identify hubs by this index.
    // To map a hub's CCH rank r to its dense id i, use hubIdOfRank(r).
    // To map a hub index i to the CCH-rank r of the associated vertex, use r = hubRanks[i].
    Level hubLevelThreshold;
    bool sizeCut = false;             // size-cut mode: hub-ness is a per-vertex property
    std::vector<int32_t> leafRegionRaw; // size-cut only: vertex -> leaf region (-1 for hubs)
    std::vector<int32_t> innerRegionRaw; // vertex -> sub-region inside its leaf region
    bool innerCutUsed = false;
    int32_t numRegionsSizeCut = 0;
    int32_t maxCutDepth = 0;
    std::vector<int32_t> hubRanks; // The hubs' CCH ranks, ordered by (level desc, rank asc)
    std::vector<uint64_t> hubBits;      // one bit per rank: is it a hub?
    std::vector<int32_t> blockRank;         // set bits before each 512-rank block
    std::vector<HubId> selectToId;  // rank-order position -> sorted hub id
};
