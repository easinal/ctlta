#pragma once

#include <unordered_set>
#include "DataStructures/Partitioning/SeparatorDecomposition.h"

// Selects transit nodes from the top k levels of separator decomposition
class TransitNodeSelector {
public:
    using TransitNodeSet = std::unordered_set<int>;

    // Default constructor
    TransitNodeSelector() = default;

    // Constructs selector with separator decomposition and number of levels
    TransitNodeSelector(const SeparatorDecomposition& sepDecomp, int numLevels)
        : decomp(sepDecomp), numLevels(numLevels) {}

    // Selects transit nodes from top k levels
    TransitNodeSet selectTransitNodes() const {
        TransitNodeSet transitNodes;
        
        // Select vertices from the top numLevels levels of the separator decomposition
        for (int level = 0; level < std::min(numLevels, static_cast<int>(decomp.tree.size())); ++level) {
            int firstVertex = decomp.firstSeparatorVertex(level);
            int lastVertex = decomp.lastSeparatorVertex(level);
            
            for (int i = firstVertex; i < lastVertex; ++i) {
                if (i < decomp.order.size()) {
                    transitNodes.insert(decomp.order[i]);
                }
            }
        }
        
        return transitNodes;
    }

    // Returns the number of levels used
    int getNumLevels() const {
        return numLevels;
    }

private:
    SeparatorDecomposition decomp;
    int numLevels = 0;
};