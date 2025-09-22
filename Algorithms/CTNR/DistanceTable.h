#pragma once

#include <unordered_set>
#include <unordered_map>
#include <fstream>

// Stores distances between all pairs of transit nodes
class DistanceTable {
public:
    using TransitNodeSet = std::unordered_set<int>;

    // Default constructor
    DistanceTable() = default;

    // Constructs distance table for the given transit nodes
    explicit DistanceTable(const TransitNodeSet& transitNodes)
        : transitNodes(transitNodes) {
        // Initialize distance table (will be filled during customization)
        for (int u : transitNodes) {
            for (int v : transitNodes) {
                if (u == v) {
                    distances[u][v] = 0;
                } else {
                    distances[u][v] = std::numeric_limits<int>::max();
                }
            }
        }
    }

    // Gets the distance between two transit nodes
    int getDistance(int source, int target) const {
        auto it1 = distances.find(source);
        if (it1 != distances.end()) {
            auto it2 = it1->second.find(target);
            if (it2 != it1->second.end()) {
                return it2->second;
            }
        }
        return std::numeric_limits<int>::max();
    }

    // Sets the distance between two transit nodes
    void setDistance(int source, int target, int distance) {
        distances[source][target] = distance;
    }

    // Returns the size of the data structure in bytes
    size_t sizeInBytes() const {
        size_t size = sizeof(*this);
        size += transitNodes.size() * sizeof(int);
        size += distances.size() * sizeof(std::unordered_map<int, int>);
        for (const auto& pair : distances) {
            size += pair.second.size() * (sizeof(int) + sizeof(int));
        }
        return size;
    }

    // Writes the distance table to a binary file
    void writeTo(std::ofstream& out) const {
        int numTransitNodes = transitNodes.size();
        out.write(reinterpret_cast<const char*>(&numTransitNodes), sizeof(numTransitNodes));
        
        for (int node : transitNodes) {
            out.write(reinterpret_cast<const char*>(&node), sizeof(node));
        }
        
        for (int u : transitNodes) {
            for (int v : transitNodes) {
                int distance = getDistance(u, v);
                out.write(reinterpret_cast<const char*>(&distance), sizeof(distance));
            }
        }
    }

    // Reads the distance table from a binary file
    void readFrom(std::ifstream& in) {
        int numTransitNodes;
        in.read(reinterpret_cast<char*>(&numTransitNodes), sizeof(numTransitNodes));
        
        transitNodes.clear();
        distances.clear();
        
        for (int i = 0; i < numTransitNodes; ++i) {
            int node;
            in.read(reinterpret_cast<char*>(&node), sizeof(node));
            transitNodes.insert(node);
        }
        
        for (int u : transitNodes) {
            for (int v : transitNodes) {
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