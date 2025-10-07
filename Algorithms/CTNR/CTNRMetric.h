#pragma once

#include "Algorithms/CCH/CCHMetric.h"
#include "Algorithms/CH/CH.h"
#include "Algorithms/CTL/BalancedTopologyCentricTreeHierarchy.h"
#include "DataStructures/Partitioning/SeparatorDecomposition.h"
#include "Tools/Constants.h"
#include <vector>
#include <cstdint>

template<typename InputGraphT>
class CTNRMetric {
public:
    using InputGraph = InputGraphT;

    // Constructor
    CTNRMetric(const SeparatorDecomposition& sepDecomp, int transitNodeThreshold)
        : sepDecomp(sepDecomp), 
          transitNodeThreshold(transitNodeThreshold),
          hierarchy(sepDecomp),
          cch(),
          cchMetric(),
          transitNodes() {
    }

    // Preprocessing phase - builds the hierarchy and identifies transit nodes
    void preprocess(const InputGraph& inputGraph) {
        // Build CCH from the input graph
        cch = CCH(sepDecomp, inputGraph);
        
        // Identify transit nodes based on threshold
        identifyTransitNodes(inputGraph);
        
        // Build the hierarchy for CTNR
        buildHierarchy();
    }

    // Customization phase - customizes weights
    void customize(const int32_t* inputWeights) {
        cchMetric = CCHMetric(cch, inputWeights);
        cchMetric.customize();
    }

    // Getters
    const BalancedTopologyCentricTreeHierarchy& getHierarchy() const { 
        return hierarchy; 
    }
    
    const CCH& getCCH() const { 
        return cch; 
    }
    
    const std::vector<int32_t>& getTransitNodes() const { 
        return transitNodes; 
    }

    // Memory usage calculation
    uint64_t sizeInBytes() const {
        uint64_t size = sizeof(CTNRMetric);
        size += hierarchy.sizeInBytes();
        size += cch.sizeInBytes();
        size += cchMetric.sizeInBytes();
        size += transitNodes.capacity() * sizeof(int32_t);
        return size;
    }

private:
    const SeparatorDecomposition& sepDecomp;
    int transitNodeThreshold;
    BalancedTopologyCentricTreeHierarchy hierarchy;
    CCH cch;
    CCHMetric cchMetric;
    std::vector<int32_t> transitNodes;

    void identifyTransitNodes(const InputGraph& inputGraph) {
        // Implementation to identify transit nodes based on separator decomposition
        // and threshold criteria
        transitNodes.clear();
        
        // For each separator in the decomposition, select nodes that exceed threshold
        for (int v = 0; v < inputGraph.numVertices(); ++v) {
            // Check if vertex is a potential transit node based on criteria
            if (isTransitNode(v, inputGraph)) {
                transitNodes.push_back(v);
            }
        }
    }

    bool isTransitNode(int vertex, const InputGraph& inputGraph) const {
        // Simplified transit node identification logic
        // In practice, this would involve more sophisticated criteria
        return inputGraph.outDegree(vertex) >= transitNodeThreshold;
    }

    void buildHierarchy() {
        // Build the balanced topology-centric tree hierarchy for CTNR
        // This would be customized for CTNR's specific needs
    }
};