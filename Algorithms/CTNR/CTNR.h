#pragma once

#include <cassert>
#include <cstdint>
#include <vector>
#include <unordered_set>
#include <unordered_map>

#include "DataStructures/Graph/Graph.h"
#include "DataStructures/Partitioning/SeparatorDecomposition.h"
#include "Algorithms/CTNR/TransitNodeSelector.h"
#include "Algorithms/CTNR/AccessNodeComputer.h"
#include "Algorithms/CTNR/LocalityFilter.h"
#include "Algorithms/CTNR/DistanceTable.h"
#include "Tools/Constants.h"

// A customizable transit node routing algorithm based on separator decomposition.
// Uses top k levels of separator decomposition as transit nodes and computes
// access nodes systematically for efficient shortest path queries.
class CTNR {
public:
    using TransitNodeSet = std::unordered_set<int>;
    using AccessNodeMap = std::unordered_map<int, std::vector<int>>; // vertex -> access nodes
    using DistanceMap = std::unordered_map<int, int>; // access node -> distance

    // Constructs an empty CTNR.
    CTNR() = default;

    // Constructs a CTNR from the specified binary file.
    explicit CTNR(std::ifstream &in) {
        readFrom(in);
    }

    // Builds the metric-independent CTNR for the specified graph and separator decomposition.
    template<typename InputGraphT>
    void preprocess(const InputGraphT &inputGraph, const SeparatorDecomposition &sepDecomp, 
                   int numTransitLevels = 5) {
        assert(inputGraph.numVertices() == sepDecomp.order.size());
        
        decomp = sepDecomp;
        numVertices = inputGraph.numVertices();
        numEdges = inputGraph.numEdges();
        
        // Select transit nodes from top k levels of separator decomposition
        transitNodeSelector = TransitNodeSelector(sepDecomp, numTransitLevels);
        transitNodes = transitNodeSelector.selectTransitNodes();
        
        // Initialize locality filter with bit strings for LCA computation
        localityFilter = LocalityFilter(sepDecomp, transitNodes);
        
        // Compute access nodes for all vertices systematically
        AccessNodeComputer<InputGraphT> accessNodeComputer(inputGraph, sepDecomp, transitNodes);
        accessNodes = accessNodeComputer.computeAccessNodes();
        
        // Initialize distance table (will be filled during customization)
        distanceTable = DistanceTable(transitNodes);
    }

    // Returns the number of vertices in the graph.
    int numVerticesInGraph() const {
        return numVertices;
    }

    // Returns the number of edges in the graph.
    int numEdgesInGraph() const {
        return numEdges;
    }

    // Returns the number of transit nodes.
    int numTransitNodes() const {
        return transitNodes.size();
    }

    // Returns the set of transit nodes.
    const TransitNodeSet& getTransitNodes() const {
        return transitNodes;
    }

    // Returns access nodes for a given vertex.
    const std::vector<int>& getAccessNodes(int vertex) const {
        auto it = accessNodes.find(vertex);
        if (it != accessNodes.end()) {
            return it->second;
        }
        static const std::vector<int> empty;
        return empty;
    }

    // Returns the distance table.
    const DistanceTable& getDistanceTable() const {
        return distanceTable;
    }

    // Returns the distance table (non-const version for customization).
    DistanceTable& getDistanceTable() {
        return distanceTable;
    }

    // Returns the locality filter.
    const LocalityFilter& getLocalityFilter() const {
        return localityFilter;
    }

    // Checks if a query between two vertices is local.
    bool isLocalQuery(int source, int target) const {
        return localityFilter.isLocal(source, target);
    }

    // Returns the size of the data structure in bytes.
    size_t sizeInBytes() const {
        size_t size = sizeof(*this);
        size += transitNodes.size() * sizeof(int);
        size += accessNodes.size() * (sizeof(int) + sizeof(std::vector<int>));
        for (const auto& pair : accessNodes) {
            size += pair.second.size() * sizeof(int);
        }
        size += distanceTable.sizeInBytes();
        size += localityFilter.sizeInBytes();
        return size;
    }

    // Writes the CTNR to a binary file.
    void writeTo(std::ofstream &out) const {
        // Write basic information
        out.write(reinterpret_cast<const char*>(&numVertices), sizeof(numVertices));
        out.write(reinterpret_cast<const char*>(&numEdges), sizeof(numEdges));
        
        // Write separator decomposition
        decomp.writeTo(out);
        
        // Write transit nodes
        int numTransitNodes = transitNodes.size();
        out.write(reinterpret_cast<const char*>(&numTransitNodes), sizeof(numTransitNodes));
        for (int node : transitNodes) {
            out.write(reinterpret_cast<const char*>(&node), sizeof(node));
        }
        
        // Write access nodes
        int numAccessNodeEntries = accessNodes.size();
        out.write(reinterpret_cast<const char*>(&numAccessNodeEntries), sizeof(numAccessNodeEntries));
        for (const auto& pair : accessNodes) {
            int vertex = pair.first;
            int numAccessNodes = pair.second.size();
            out.write(reinterpret_cast<const char*>(&vertex), sizeof(vertex));
            out.write(reinterpret_cast<const char*>(&numAccessNodes), sizeof(numAccessNodes));
            for (int accessNode : pair.second) {
                out.write(reinterpret_cast<const char*>(&accessNode), sizeof(accessNode));
            }
        }
        
        // Write distance table
        distanceTable.writeTo(out);
        
        // Write locality filter
        localityFilter.writeTo(out);
    }

    // Reads the CTNR from a binary file.
    void readFrom(std::ifstream &in) {
        // Read basic information
        in.read(reinterpret_cast<char*>(&numVertices), sizeof(numVertices));
        in.read(reinterpret_cast<char*>(&numEdges), sizeof(numEdges));
        
        // Read separator decomposition
        decomp.readFrom(in);
        
        // Read transit nodes
        int numTransitNodes;
        in.read(reinterpret_cast<char*>(&numTransitNodes), sizeof(numTransitNodes));
        transitNodes.clear();
        for (int i = 0; i < numTransitNodes; ++i) {
            int node;
            in.read(reinterpret_cast<char*>(&node), sizeof(node));
            transitNodes.insert(node);
        }
        
        // Read access nodes
        int numAccessNodeEntries;
        in.read(reinterpret_cast<char*>(&numAccessNodeEntries), sizeof(numAccessNodeEntries));
        accessNodes.clear();
        for (int i = 0; i < numAccessNodeEntries; ++i) {
            int vertex, numAccessNodes;
            in.read(reinterpret_cast<char*>(&vertex), sizeof(vertex));
            in.read(reinterpret_cast<char*>(&numAccessNodes), sizeof(numAccessNodes));
            std::vector<int> accessNodeList;
            for (int j = 0; j < numAccessNodes; ++j) {
                int accessNode;
                in.read(reinterpret_cast<char*>(&accessNode), sizeof(accessNode));
                accessNodeList.push_back(accessNode);
            }
            accessNodes[vertex] = accessNodeList;
        }
        
        // Read distance table
        distanceTable.readFrom(in);
        
        // Read locality filter
        localityFilter.readFrom(in);
    }

private:
    SeparatorDecomposition decomp;
    int numVertices = 0;
    int numEdges = 0;
    
    TransitNodeSet transitNodes;
    AccessNodeMap accessNodes;
    DistanceTable distanceTable;
    LocalityFilter localityFilter;
    TransitNodeSelector transitNodeSelector;
    // Note: AccessNodeComputer is a template class, so we can't store it as a member
    // We'll create it when needed in the preprocess method
};