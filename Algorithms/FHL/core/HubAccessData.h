#pragma once

// FHL's per-vertex hub access lists: for every vertex, which hubs its shortest paths can
// first reach and how far they are. This is SCAFFOLDING, not a query structure --
// customization folds it into per-region seeds and then releases it (see
// FHLMetric::releaseAccessNodes), which is why there is no distance table in
// here at all. CTNR's CTNRData keeps the same arrays alive AS its query structure and
// carries the T x T table alongside; forked so that neither side constrains the other.

#include <cstdint>
#include <utility>
#include <vector>

#include "Algorithms/FHL/core/FHLConstants.h"

namespace fhl { struct SeedRectangle; struct Overlay; }
#include "DataStructures/Labels/BasicLabelSet.h"
#include "DataStructures/Labels/ParentInfo.h"
#include "DataStructures/Labels/SimdLabelSet.h"
#include "DataStructures/Utilities/IteratorRange.h"
#include "Tools/Simd/AlignedVector.h"

class HubAccessData {

    // Access-hub distance runs are padded to this width so the label-set scans stay aligned.
    static constexpr int LOGK = 3; // K = 8
    static constexpr int K = 1 << LOGK;
    static constexpr bool USE_SIMD = (LOGK > 0);
    static_assert(LOGK != 1, "LOGK=1 (K=2) not supported for SIMD label sets.");

    using LabelSet = std::conditional_t<USE_SIMD,
            SimdLabelSet<LOGK, ParentInfo::NO_PARENT_INFO>,
            BasicLabelSet<0, ParentInfo::NO_PARENT_INFO>>;
    using DistanceLabel = typename LabelSet::DistanceLabel;
    using LabelMask = typename LabelSet::LabelMask;

    template<typename T>
    using DistanceVector = AlignedVector<T>;

public:

    explicit HubAccessData(const int numVertices) : numVertices(numVertices) {}

    ConstantVectorRange<int32_t, DistanceVector> getForwardDistances(const int v) const {
        return {forwardDistances.begin() + pos[rankToIdx(v)],
                forwardDistances.begin() + pos[rankToIdx(v) + 1]};
    }

    ConstantVectorRange<int32_t, DistanceVector> getBackwardDistances(const int v) const {
        // Symmetric metrics store no backward distances; forward doubles for both directions.
        const auto &arr = backwardDistances.empty() ? forwardDistances : backwardDistances;
        return {arr.begin() + pos[rankToIdx(v)], arr.begin() + pos[rankToIdx(v) + 1]};
    }

    uint64_t sizeInBytes() const {
        return sizeof(HubAccessData) + sizeAccessHubsInBytes();
    }

    uint64_t sizeAccessHubsInBytes() const {
        return pos.size() * sizeof(decltype(pos)::value_type) +
               accessHubs.size() * sizeof(decltype(accessHubs)::value_type) +
               forwardDistances.size() * sizeof(decltype(forwardDistances)::value_type) +
               backwardDistances.size() * sizeof(decltype(backwardDistances)::value_type);
    }

    std::pair<double, double> computeAverageNumberOfNonInftyAccessHubs() const {
        if (pos.empty()) // scaffolding already released
            return {0.0, 0.0};
        int64_t fSum = 0, bSum = 0, count = 0;
        for (int idx = 0; idx < numVertices; ++idx) {
            int localF = 0, localB = 0;
            for (auto i = pos[idx]; i < pos[idx + 1]; ++i) {
                if (forwardDistances[i] != FHL_INFTY)
                    ++localF;
                if ((backwardDistances.empty() ? forwardDistances : backwardDistances)[i] !=
                    FHL_INFTY)
                    ++localB;
            }
            fSum += localF;
            bSum += localB;
            ++count;
        }
        return {static_cast<double>(fSum) / count, static_cast<double>(bSum) / count};
    }

private:

    // Internal vertex index for CCH rank r. Ranks are inverted so the top-down access-hub
    // construction writes sequentially.
    inline int rankToIdx(const int r) const {
        return numVertices - 1 - r;
    }

    friend class HubAccessPreprocessor;

    template<typename>
    friend class FHLMetric;

    friend struct fhl::SeedRectangle;
    friend struct fhl::Overlay;

    int numVertices;

    // Vertex v's access hubs live in accessHubs[pos[v] .. pos[v+1]-1].
    std::vector<int32_t> pos;
    std::vector<fhl::HubId> accessHubs;

    DistanceVector<int32_t> forwardDistances;
    DistanceVector<int32_t> backwardDistances;
};

// A CCH edge (u,v) where u is not a hub and v is -- hence an access hub of u.
struct HubAccessEdge {
    int32_t edge = INVALID_EDGE;     // edge ID in the CCH up-graph
    int32_t position = INVALID_INDEX; // index of v's distance entry for u in HubAccessData
};
