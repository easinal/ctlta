#pragma once

#include <vector>
#include <unordered_map>
#include <cstdint>

#include "DataStructures/Utilities/IteratorRange.h"
#include "DataStructures/Partitioning/SeparatorDecomposition.h"
#include "Algorithms/CTL/BalancedTopologyCentricTreeHierarchy.h"

class CTNRData {

public:

    explicit CTNRData(const int numTransitNodes, const int numVertices) :
    numTransitNodes(numTransitNodes),
    numVertices(numVertices),
    distanceTable(numTransitNodes * numTransitNodes, INFTY) {}

    // Input: Internal transit node index of access nodes.
    int getDistanceBetweenTransitNodes(int32_t indexS, int32_t indexT) const {
        return distanceTable[indexS * numTransitNodes + indexT];
    }

    ConstantVectorRange<int32_t> getForwardAccessNodes(const int v) const {
        return {forwardNodes.begin() + forwardPos[rankToIdx(v)], forwardNodes.begin() + forwardPos[rankToIdx(v) + 1]};
    }

    ConstantVectorRange<int32_t> getForwardDistances(const int v) const {
        return {forwardDistances.begin() + forwardPos[rankToIdx(v)], forwardDistances.begin() + forwardPos[rankToIdx(v) + 1]};
    }

    ConstantVectorRange<int32_t> getBackwardAccessNodes(const int v) const {
        return {backwardNodes.begin() + backwardPos[rankToIdx(v)], backwardNodes.begin() + backwardPos[rankToIdx(v) + 1]};
    }

    ConstantVectorRange<int32_t> getBackwardDistances(const int v) const {
        return {backwardDistances.begin() + backwardPos[rankToIdx(v)], backwardDistances.begin() + backwardPos[rankToIdx(v) + 1]};
    }

    // Memory usage calculation including node levels
    uint64_t sizeInBytes() const {
        uint64_t size = sizeof(CTNRData);

        size += forwardPos.size() * sizeof(int32_t);
        size += forwardNodes.size() * sizeof(int32_t);
        size += forwardDistances.size() * sizeof(int32_t);
        size += backwardPos.size() * sizeof(int32_t);
        size += backwardNodes.size() * sizeof(int32_t);
        size += backwardDistances.size() * sizeof(int32_t);
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

    int numTransitNodes;
    int numVertices;

    // Range of forward access nodes for vertex v is stored in forwardAccess[forwardPos[v]..forwardPos[v+1]-1]
    std::vector<int32_t> forwardPos;
    std::vector<int32_t> forwardNodes;
    std::vector<int32_t> forwardDistances;

    // Range of backward access nodes for vertex v is stored in backwardAccess[backwardPos[v]..backwardPos[v+1]-1]
    std::vector<int32_t> backwardPos;
    std::vector<int32_t> backwardNodes;
    std::vector<int32_t> backwardDistances;

    // distanceTable[i * numTransitNodes + j] = distance from transit node i to transit node j, where i and j are
    // internal transit node indices.
    std::vector<int32_t> distanceTable;


};