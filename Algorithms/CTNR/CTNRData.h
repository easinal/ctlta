#pragma once

#include <vector>
#include <unordered_map>
#include <cstdint>

#include "DataStructures/Partitioning/SeparatorDecomposition.h"
#include "Algorithms/CTL/BalancedTopologyCentricTreeHierarchy.h"

class CTNRData {

public:

    CTNRData(const SeparatorDecomposition &sepDecomp,
             const BalancedTopologyCentricTreeHierarchy &hierarchy,
             int transitNodeThreshold)
            : sepDecomp(sepDecomp),
              hierarchy(hierarchy),
              transitNodeThreshold(transitNodeThreshold) {}

    void init() {

        selectTransitNodes();

        forwardAccessNodes.resize(hierarchy.numVertices());
        forwardAccessDistances.resize(hierarchy.numVertices());
        backwardAccessNodes.resize(hierarchy.numVertices());
        backwardAccessDistances.resize(hierarchy.numVertices());

    }

    // Input: Internal transit node index of access nodes.
    int getDistanceBetweenTransitNodes(int32_t accessS, int32_t accessT) const {
        return distanceTable[accessS][accessT];
    }

    int getTransitNodeThreshold() const { return transitNodeThreshold; }

    const std::vector<int32_t> &getTransitNodes() const { return transitNodes; }

    const std::vector<std::vector<int32_t>> &getForwardAccessNodes() const { return forwardAccessNodes; }

    const std::vector<std::vector<int32_t>> &getForwardAccessDistances() const { return forwardAccessDistances; }

    const std::vector<std::vector<int32_t>> &getBackwardAccessNodes() const { return backwardAccessNodes; }

    const std::vector<std::vector<int32_t>> &getBackwardAccessDistances() const { return backwardAccessDistances; }

    const std::vector<std::vector<int32_t>> &getDistanceTable() const { return distanceTable; }
    const std::unordered_map<int32_t, int32_t> &
    gettransitNodeToDistanceTableIndex() const { return transitNodeToDistanceTableIndex; }

    // TODO: Replace this with calls to hierarchy wherever it occurs
    int getVertexLevel(const int32_t &v) const {
        return separatorNodeToLevel.at(v);
    }

    // Memory usage calculation including node levels
    uint64_t sizeInBytes() const {
        uint64_t size = sizeof(CTNRData);
        size += transitNodes.capacity() * sizeof(int32_t);

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

    void selectTransitNodes() {
        separatorNodeToLevel.resize(sepDecomp.tree.size(), -1);
        std::function<void(int, int)> collectTransitNodes =
                [&](int node, int level) {
                    if (level <= transitNodeThreshold) {
                        separatorNodeToLevel[node] = level;
                        for (int v = sepDecomp.firstSeparatorVertex(node);
                             v < sepDecomp.lastSeparatorVertex(node); ++v) {
                            transitNodes.push_back(v);
                            transitVertexToLevel[v] = level;
                        }
                    } else {
                        // return;
                        separatorNodeToLevel[node] = level;
                    }

                    int child = sepDecomp.leftChild(node);
                    while (child != 0) {
                        collectTransitNodes(child, level + 1);
                        child = sepDecomp.rightSibling(child);
                    }
                };

        collectTransitNodes(0, 0);
        auto compareByLevel = [&](int32_t a, int32_t b) {
            return transitVertexToLevel[a] < transitVertexToLevel[b];
        };
        std::sort(transitNodes.begin(), transitNodes.end(), compareByLevel);
        for (int i = 0; i < transitNodes.size(); ++i) {
            // std::cout<<"transitNodes["<<i<<"]: "<<transitNodes[i]<<std::endl;
            transitNodeToDistanceTableIndex[transitNodes[i]] = i;
        }
        std::cout << "CTNR: Selected " << transitNodes.size() << " transit nodes from top " << transitNodeThreshold
                  << " levels" << std::endl;
        std::cout << "Total number of vertices: " << forwardAccessNodes.size() << std::endl;
    }

    friend class CTNRMetric;

    const SeparatorDecomposition &sepDecomp;
    const BalancedTopologyCentricTreeHierarchy &hierarchy;

    std::vector<int32_t> separatorNodeToLevel;  // separator node ID -> level

    // Transit Node related
    int transitNodeThreshold;
    std::vector<int32_t> transitNodes;
    std::unordered_map<int32_t, int32_t> transitVertexToLevel;  // vertex ID -> level
    std::unordered_map<int32_t, int32_t> transitNodeToDistanceTableIndex; // key: TN rank id

    // Access Nodes (indexed by rank IDs)
    std::vector<std::vector<int32_t>> forwardAccessNodes;  // forwardAccessNodes[rank(v)] = TN ranks
    std::vector<std::vector<int32_t>> forwardAccessDistances;  // distances
    // Backward Access (a -> v), indexed by rank IDs
    std::vector<std::vector<int32_t>> backwardAccessNodes;      // backwardAccessNodes[rank(v)] = TN ranks
    std::vector<std::vector<int32_t>> backwardAccessDistances;  // distances

    // Distance table
    std::vector<std::vector<int32_t>> distanceTable;  // distanceTable[i][j] = distance from node i to node j


};