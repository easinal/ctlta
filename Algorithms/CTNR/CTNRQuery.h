#pragma once

#include "Algorithms/CTNR/CTNRMetric.h"
#include "Algorithms/CCH/EliminationTreeQuery.h"
#include "Algorithms/Dijkstra/BiDijkstra.h"
#include "DataStructures/Labels/BasicLabelSet.h"
#include "DataStructures/Labels/ParentInfo.h"
#include "Tools/Constants.h"
#include <vector>
#include <cstdint>

template<typename InputGraphT>
class CTNRQuery {
public:
    using InputGraph = InputGraphT;
    using Metric = CTNRMetric<InputGraphT>;

    // Constructor
    CTNRQuery(const Metric& metric)
        : metric(metric),
          cchQuery(metric.getCCH().getUpwardGraph(), 
                   metric.getCCH().getUpwardGraph(),
                   nullptr, nullptr), // weights will be set during queries
          transitNodeDistances(),
          lastSource(INVALID_VERTEX),
          lastTarget(INVALID_VERTEX),
          lastDistance(INFTY) {
        
        // Initialize transit node distance matrix
        const auto& transitNodes = metric.getTransitNodes();
        const int numTransitNodes = transitNodes.size();
        transitNodeDistances.resize(numTransitNodes * numTransitNodes, INFTY);
    }

    // Main query method
    int32_t run(int32_t source, int32_t target) {
        lastSource = source;
        lastTarget = target;
        lastDistance = INFTY;

        // Check if both source and target are transit nodes
        const auto& transitNodes = metric.getTransitNodes();
        int sourceTransitIdx = findTransitNodeIndex(source, transitNodes);
        int targetTransitIdx = findTransitNodeIndex(target, transitNodes);

        if (sourceTransitIdx != -1 && targetTransitIdx != -1) {
            // Both are transit nodes - use precomputed distance
            const int numTransitNodes = transitNodes.size();
            lastDistance = transitNodeDistances[sourceTransitIdx * numTransitNodes + targetTransitIdx];
        } else if (sourceTransitIdx != -1 || targetTransitIdx != -1) {
            // One is a transit node - use hybrid approach
            lastDistance = queryWithTransitNode(source, target, sourceTransitIdx, targetTransitIdx);
        } else {
            // Neither is a transit node - use CCH query
            lastDistance = cchQuery.run(source, target);
        }

        return lastDistance;
    }

    // Get the distance from the last query
    int32_t getDistance() const {
        return lastDistance;
    }

    // Get the last source vertex
    int32_t getLastSource() const {
        return lastSource;
    }

    // Get the last target vertex
    int32_t getLastTarget() const {
        return lastTarget;
    }

    // Memory usage calculation
    uint64_t sizeInBytes() const {
        uint64_t size = sizeof(CTNRQuery);
        size += cchQuery.sizeInBytes();
        size += transitNodeDistances.capacity() * sizeof(int32_t);
        return size;
    }

private:
    const Metric& metric;
    EliminationTreeQuery<typename CCH::UpGraph, BasicLabelSet<0, ParentInfo::NO_PARENT_INFO>> cchQuery;
    std::vector<int32_t> transitNodeDistances; // Precomputed distances between transit nodes
    int32_t lastSource;
    int32_t lastTarget;
    int32_t lastDistance;

    // Find the index of a vertex in the transit nodes list
    int findTransitNodeIndex(int32_t vertex, const std::vector<int32_t>& transitNodes) const {
        auto it = std::find(transitNodes.begin(), transitNodes.end(), vertex);
        return (it != transitNodes.end()) ? (it - transitNodes.begin()) : -1;
    }

    // Query when one vertex is a transit node
    int32_t queryWithTransitNode(int32_t source, int32_t target, int sourceTransitIdx, int targetTransitIdx) {
        const auto& transitNodes = metric.getTransitNodes();
        int32_t minDist = INFTY;

        if (sourceTransitIdx != -1) {
            // Source is a transit node - find distance to target via other transit nodes
            for (size_t i = 0; i < transitNodes.size(); ++i) {
                int32_t transitNode = transitNodes[i];
                int32_t distToTransit = transitNodeDistances[sourceTransitIdx * transitNodes.size() + i];
                int32_t distFromTransit = cchQuery.run(transitNode, target);
                
                if (distToTransit != INFTY && distFromTransit != INFTY) {
                    minDist = std::min(minDist, distToTransit + distFromTransit);
                }
            }
            
            // Also try direct path from source to target
            int32_t directDist = cchQuery.run(source, target);
            minDist = std::min(minDist, directDist);
        } else {
            // Target is a transit node - find distance from source via other transit nodes
            for (size_t i = 0; i < transitNodes.size(); ++i) {
                int32_t transitNode = transitNodes[i];
                int32_t distToTransit = cchQuery.run(source, transitNode);
                int32_t distFromTransit = transitNodeDistances[i * transitNodes.size() + targetTransitIdx];
                
                if (distToTransit != INFTY && distFromTransit != INFTY) {
                    minDist = std::min(minDist, distToTransit + distFromTransit);
                }
            }
            
            // Also try direct path from source to target
            int32_t directDist = cchQuery.run(source, target);
            minDist = std::min(minDist, directDist);
        }

        return minDist;
    }
};