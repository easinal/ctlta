#pragma once

#include <vector>
#include <unordered_map>
#include <cstdint>

#include "DataStructures/Utilities/IteratorRange.h"
#include "DataStructures/Partitioning/SeparatorDecomposition.h"
#include "Algorithms/CTL/BalancedTopologyCentricTreeHierarchy.h"
#include "CTNRConstants.h"
#include "DataStructures/Labels/SimdLabelSet.h"
#include "DataStructures/Labels/BasicLabelSet.h"

class CTNRData {
    static constexpr int LOGK = 3; // K=8
    static constexpr int K = 1 << LOGK;
    static_assert(LOGK != 1, "LOGK=1 (K=2) not supported for SIMD label sets.");
    static constexpr bool USE_SIMD = (LOGK > 0);

    using LabelSet = std::conditional_t<USE_SIMD, SimdLabelSet<LOGK, ParentInfo::NO_PARENT_INFO>, BasicLabelSet<0, ParentInfo::NO_PARENT_INFO>>;
    using DistanceLabel = typename LabelSet::DistanceLabel;
    using LabelMask = typename LabelSet::LabelMask;

    template<typename T>
    using DistanceVector = AlignedVector<T>;

public:

    explicit CTNRData(const int numTransitNodes, const int numVertices) :
    numTransitNodes(numTransitNodes),
    numVertices(numVertices),
    distanceTable(numTransitNodes * numTransitNodes, CTNR_INFTY) {}

    // Input: Internal transit node index of access nodes.
    DEBUG_NOINLINE
    int getDistanceBetweenTransitNodes(int32_t indexS, int32_t indexT) const {
        return distanceTable[indexS * numTransitNodes + indexT];
    }

    ConstantVectorRange<int32_t> getAccessNodes(const int v) const {
        return {accessNodes.begin() + pos[rankToIdx(v)], accessNodes.begin() + pos[rankToIdx(v) + 1]};
    }

    ConstantVectorRange<int32_t, DistanceVector> getForwardDistances(const int v) const {
        return {forwardDistances.begin() + pos[rankToIdx(v)], forwardDistances.begin() + pos[rankToIdx(v) + 1]};
    }

    ConstantVectorRange<int32_t, DistanceVector> getBackwardDistances(const int v) const {
        return {backwardDistances.begin() + pos[rankToIdx(v)], backwardDistances.begin() + pos[rankToIdx(v) + 1]};
    }

    // Memory usage calculation including node levels
    uint64_t sizeInBytes() const {
        uint64_t size = sizeof(CTNRData);

        size += pos.size() * sizeof(int32_t);
        size += accessNodes.size() * sizeof(int32_t);
        size += forwardDistances.size() * sizeof(int32_t);
        size += backwardDistances.size() * sizeof(int32_t);
        size += distanceTable.size() * sizeof(int32_t);

        return size;
    }

    uint64_t sizeDistanceTableInBytes() const {
        return distanceTable.size() * sizeof(int32_t);
    }

    uint64_t sizeAccessNodesInBytes() const {
        return pos.size() * sizeof(int32_t) +
               accessNodes.size() * sizeof(int32_t) +
               forwardDistances.size() * sizeof(int32_t) +
               backwardDistances.size() * sizeof(int32_t);
    }

private:

    // Get internal vertex index for CCH rank r.
    // Invert ranks for sequential writing order during top-down access node construction.
    inline int rankToIdx(const int r) const {
        return numVertices - 1 - r ;
    }

    void resetDistanceTable() {
        distanceTable.assign(numTransitNodes * numTransitNodes, CTNR_INFTY);
    }

    // Input: Internal transit node index of access nodes.
    void setDistanceBetweenTransitNodes(int32_t indexS, int32_t indexT, int32_t distance) {
        distanceTable[indexS * numTransitNodes + indexT] = distance;
    }

    friend class CTNRMetric;
    friend class CTNRPreprocessor;

    int numTransitNodes;
    int numVertices;

    // Range of access nodes for vertex v is stored in accessNodes[pos[v]..pos[v+1]-1]
    std::vector<int32_t> pos;
    std::vector<int32_t> accessNodes;

    DistanceVector<int32_t> forwardDistances;
    DistanceVector<int32_t> backwardDistances;

    // distanceTable[i * numTransitNodes + j] = distance from transit node i to transit node j, where i and j are
    // internal transit node indices.
    std::vector<int32_t> distanceTable;


};


// Contains information on a CCH-edge (u,v) where u is not a transit node and v is a transit node (and as such an access node of u).
struct AccessNodeEdge {
    int32_t edge = INVALID_EDGE; // edge ID in CCH upGraph
    int32_t position = INVALID_INDEX; // index of distance entry for access node v of u in CTNRData
};