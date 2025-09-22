#pragma once

#include <unordered_set>
#include <iostream>
#include "DataStructures/Partitioning/SeparatorDecomposition.h"

// Selects transit nodes from the top k levels of separator decomposition.
// Transit nodes are vertices that appear in the upper levels of the separator tree
// and serve as important routing hubs for non-local queries.
class TransitNodeSelector {
public:
    using TransitNodeSet = std::unordered_set<int>;

    // Default constructor.
    TransitNodeSelector() = default;

    // Constructs a transit node selector for the given separator decomposition.
    TransitNodeSelector(const SeparatorDecomposition& sepDecomp, int numLevels)
        : decomp(sepDecomp), numLevels(numLevels) {}

    // Selects transit nodes from the top k levels of the separator decomposition.
    TransitNodeSet selectTransitNodes() const {
        TransitNodeSet transitNodes;
        
        std::cout << "Selecting transit nodes from top " << numLevels << " levels..." << std::endl;
        
        // Select vertices from the top numLevels levels of the separator tree
        int actualLevels = std::min(numLevels, static_cast<int>(decomp.tree.size()));
        
        for (int level = 0; level < actualLevels; ++level) {
            int firstVertex = decomp.firstSeparatorVertex(level);
            int lastVertex = decomp.lastSeparatorVertex(level);
            
            for (int i = firstVertex; i < lastVertex; ++i) {
                if (i >= 0 && i < static_cast<int>(decomp.order.size())) {
                    int vertex = decomp.order[i];
                    transitNodes.insert(vertex);
                }
            }
        }
        
        std::cout << "Selected " << transitNodes.size() << " transit nodes." << std::endl;
        return transitNodes;
    }

    // Returns the number of levels to consider for transit node selection.
    int getNumLevels() const {
        return numLevels;
    }

    // Returns the separator decomposition.
    const SeparatorDecomposition& getDecomp() const {
        return decomp;
    }

private:
    SeparatorDecomposition decomp;
    int numLevels = 0;
};