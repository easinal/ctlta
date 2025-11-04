#pragma once

#include <vector>
#include <unordered_map>
#include <cstdint>

#include "DataStructures/Partitioning/SeparatorDecomposition.h"
#include "Algorithms/CTL/BalancedTopologyCentricTreeHierarchy.h"

class CTNRData {

public:

    CTNRData() = default;

    // Input: Internal transit node index of access nodes.
    int getDistanceBetweenTransitNodes(int32_t accessS, int32_t accessT) const {
        return distanceTable[accessS][accessT];
    }

    const std::vector<std::vector<int32_t>> &getForwardAccessNodes() const { return forwardAccessNodes; }

    const std::vector<std::vector<int32_t>> &getForwardAccessDistances() const { return forwardAccessDistances; }

    const std::vector<std::vector<int32_t>> &getBackwardAccessNodes() const { return backwardAccessNodes; }

    const std::vector<std::vector<int32_t>> &getBackwardAccessDistances() const { return backwardAccessDistances; }

    const std::vector<std::vector<int32_t>> &getDistanceTable() const { return distanceTable; }

    // Memory usage calculation including node levels
    uint64_t sizeInBytes() const {
        uint64_t size = sizeof(CTNRData);

        size += forwardAccessNodes.size() * sizeof(std::vector<int32_t>);
        size += forwardAccessDistances.size() * sizeof(std::vector<int32_t>);
        size += backwardAccessNodes.size() * sizeof(std::vector<int32_t>);
        size += backwardAccessDistances.size() * sizeof(std::vector<int32_t>);
        size += distanceTable.size() * sizeof(std::vector<int32_t>);

        for (const auto &vec: forwardAccessNodes) {
            size += vec.capacity() * sizeof(int32_t);
        }
        for (const auto &vec: forwardAccessDistances) {
            size += vec.capacity() * sizeof(int32_t);
        }
        for (const auto &vec: backwardAccessNodes) {
            size += vec.capacity() * sizeof(int32_t);
        }
        for (const auto &vec: backwardAccessDistances) {
            size += vec.capacity() * sizeof(int32_t);
        }

        for (const auto &vec: distanceTable) {
            size += vec.capacity() * sizeof(int32_t);
        }

        return size;
    }

private:

    friend class CTNRMetric;

    // Access Nodes (indexed by rank IDs)

    // Given a rank r, forwardAccessNodes[r] contains the internal transit node indices of the forward access nodes
    // of the vertex with rank r.
    std::vector<std::vector<int32_t>> forwardAccessNodes;
    // Given a rank r, forwardAccessDistances[r] contains the distances to each forward access node of the vertex
    // with rank r.
    std::vector<std::vector<int32_t>> forwardAccessDistances;

    // Given a rank r, backwardAccessNodes[r] contains the internal transit node indices of the backward access nodes
    // of the vertex with rank r.
    std::vector<std::vector<int32_t>> backwardAccessNodes;
    // Given a rank r, backwardAccessDistances[r] contains the distances to each backward access node of the vertex
    // with rank r.
    std::vector<std::vector<int32_t>> backwardAccessDistances;

    // distanceTable[i][j] = distance from transit node i to transit node j, where i and j are internal transit node
    // indices.
    std::vector<std::vector<int32_t>> distanceTable;


};