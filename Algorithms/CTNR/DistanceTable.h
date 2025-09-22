#pragma once

#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <fstream>
#include <limits>
#include "Tools/Constants.h"

// Stores distances between all pairs of transit nodes for CTNR queries.
// The table is filled during customization phase and used during query phase.
class DistanceTable {
public:
    using TransitNodeSet = std::unordered_set<int>;

    // Default constructor.
    DistanceTable() = default;

    // Constructs a distance table for the given set of transit nodes.
    explicit DistanceTable(const TransitNodeSet& transitNodes)
        : transitNodes(transitNodes) {
        // Initialize distance table with infinity values
        for (int u : transitNodes) {
            for (int v : transitNodes) {
                if (u == v) {
                    distances[u][v] = 0;
                } else {
                    distances[u][v] = INFTY;
                }
            }
        }
    }

    // Gets the distance between two transit nodes.
    int getDistance(int source, int target) const {
        auto it1 = distances.find(source);
        if (it1 != distances.end()) {
            auto it2 = it1->second.find(target);
            if (it2 != it1->second.end()) {
                return it2->second;
            }
        }
        return INFTY;
    }

    // Sets the distance between two transit nodes.
    void setDistance(int source, int target, int distance) {
        distances[source][target] = distance;
    }

    // Returns true if the table contains the given transit node.
    bool hasTransitNode(int node) const {
        return transitNodes.find(node) != transitNodes.end();
    }

    // Returns the set of transit nodes.
    const TransitNodeSet& getTransitNodes() const {
        return transitNodes;
    }

    // Returns the number of transit nodes.
    size_t numTransitNodes() const {
        return transitNodes.size();
    }

    // Returns true if the distance table is fully initialized.
    bool isComplete() const {
        for (int u : transitNodes) {
            for (int v : transitNodes) {
                if (u != v && getDistance(u, v) == INFTY) {
                    return false;
                }
            }
        }
        return true;
    }

    // Returns the size of the data structure in bytes.
    size_t sizeInBytes() const {
        size_t size = sizeof(*this);
        size += transitNodes.size() * sizeof(int);
        size += distances.size() * sizeof(std::unordered_map<int, int>);
        for (const auto& pair : distances) {
            size += pair.second.size() * (sizeof(int) + sizeof(int));
        }
        return size;
    }

    // Writes the distance table to a binary file.
    void writeTo(std::ofstream& out) const {
        // Write number of transit nodes
        int numTransitNodes = transitNodes.size();
        out.write(reinterpret_cast<const char*>(&numTransitNodes), sizeof(numTransitNodes));
        
        // Write transit node IDs
        for (int node : transitNodes) {
            out.write(reinterpret_cast<const char*>(&node), sizeof(node));
        }
        
        // Write distances in consistent order
        for (int u : transitNodes) {
            for (int v : transitNodes) {
                int distance = getDistance(u, v);
                out.write(reinterpret_cast<const char*>(&distance), sizeof(distance));
            }
        }
    }

    // Reads the distance table from a binary file.
    void readFrom(std::ifstream& in) {
        // Read number of transit nodes
        int numTransitNodes;
        in.read(reinterpret_cast<char*>(&numTransitNodes), sizeof(numTransitNodes));
        
        // Clear existing data
        transitNodes.clear();
        distances.clear();
        
        // Read transit node IDs
        std::vector<int> nodeList;
        for (int i = 0; i < numTransitNodes; ++i) {
            int node;
            in.read(reinterpret_cast<char*>(&node), sizeof(node));
            transitNodes.insert(node);
            nodeList.push_back(node);
        }
        
        // Read distances in the same order they were written
        for (int i = 0; i < numTransitNodes; ++i) {
            int u = nodeList[i];
            for (int j = 0; j < numTransitNodes; ++j) {
                int v = nodeList[j];
                int distance;
                in.read(reinterpret_cast<char*>(&distance), sizeof(distance));
                distances[u][v] = distance;
            }
        }
    }

private:
    TransitNodeSet transitNodes;
    std::unordered_map<int, std::unordered_map<int, int>> distances;
};