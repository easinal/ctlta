#pragma once

#include <vector>
#include <unordered_map>
#include <cstdint>

#include "DataStructures/Utilities/IteratorRange.h"
#include "DataStructures/Partitioning/SeparatorDecomposition.h"
#include "Algorithms/CTL/BalancedTopologyCentricTreeHierarchy.h"

class CTNRData {

public:

    struct AccessNode {
        int32_t nodeIndex = INVALID_INDEX; // internal transit node index
        int32_t distance = INFTY; // distance to/from access node
    };

    explicit CTNRData(const int numTransitNodes, const int numVertices) :
    numTransitNodes(numTransitNodes),
    numVertices(numVertices),
    distanceTable(numTransitNodes * numTransitNodes, INFTY) {}

    // Input: Internal transit node index of access nodes.
    int getDistanceBetweenTransitNodes(int32_t indexS, int32_t indexT) const {
        return distanceTable[indexS * numTransitNodes + indexT];
    }

    ConstantVectorRange<AccessNode> getForwardAccessNodes(const int v) const {
        return {forwardAccess.begin() + forwardPos[rankToIdx(v)], forwardAccess.begin() + forwardPos[rankToIdx(v) + 1]};
    }

    ConstantVectorRange<AccessNode> getBackwardAccessNodes(const int v) const {
        return {backwardAccess.begin() + backwardPos[rankToIdx(v)], backwardAccess.begin() + backwardPos[rankToIdx(v) + 1]};
    }

    // Memory usage calculation including node levels
    uint64_t sizeInBytes() const {
        uint64_t size = sizeof(CTNRData);

        size += forwardPos.size() * sizeof(int32_t);
        size += forwardAccess.size() * sizeof(AccessNode);
        size += backwardPos.size() * sizeof(int32_t);
        size += backwardAccess.size() * sizeof(AccessNode);
        size += distanceTable.size() * sizeof(int32_t);

        return size;
    }

private:

    // Get internal vertex index for CCH rank r.
    // Invert ranks for sequential writing order during top-down access node construction.
    inline int rankToIdx(const int r) const {
        return numVertices - 1 - r ;
    }

    void resetDistanceTable() {
        distanceTable.assign(numTransitNodes * numTransitNodes, INFTY);
    }

    // Input: Internal transit node index of access nodes.
    void setDistanceBetweenTransitNodes(int32_t indexS, int32_t indexT, int32_t distance) {
        distanceTable[indexS * numTransitNodes + indexT] = distance;
    }

    friend class CTNRMetric;
    friend class TransitDistanceTableBuilder;

    int numTransitNodes;
    int numVertices;

    // Range of forward access nodes for vertex v is stored in forwardAccess[forwardPos[v]..forwardPos[v+1]-1]
    std::vector<int32_t> forwardPos;
    std::vector<AccessNode> forwardAccess;

    // Range of backward access nodes for vertex v is stored in backwardAccess[backwardPos[v]..backwardPos[v+1]-1]
    std::vector<int32_t> backwardPos;
    std::vector<AccessNode> backwardAccess;

    // distanceTable[i * numTransitNodes + j] = distance from transit node i to transit node j, where i and j are
    // internal transit node indices.
    std::vector<int32_t> distanceTable;


};