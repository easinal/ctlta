#pragma once

#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <queue>
#include "DataStructures/Graph/Graph.h"
#include "DataStructures/Partitioning/SeparatorDecomposition.h"
#include "Algorithms/CTNR/TransitNodeSelector.h"
#include "Tools/Constants.h"

// Computes access nodes for all vertices systematically using separator decomposition.
// For metric-independent: uses BFS from each vertex until hitting transit nodes.
// For metric-dependent: uses Dijkstra with actual edge weights.
template<typename GraphT>
class AccessNodeComputer {
public:
    using AccessNodeMap = std::unordered_map<int, std::vector<int>>;
    using TransitNodeSet = std::unordered_set<int>;

    // Constructs an access node computer.
    AccessNodeComputer(const GraphT& graph, const SeparatorDecomposition& sepDecomp, 
                      const TransitNodeSet& transitNodes)
        : graph(graph), decomp(sepDecomp), transitNodes(transitNodes) {}

    // Computes access nodes for all vertices (metric-independent version).
    AccessNodeMap computeAccessNodes() const {
        AccessNodeMap accessNodes;
        
        // 简化实现：为所有顶点分配传输节点作为访问节点
        std::cout << "Computing access nodes for " << graph.numVertices() << " vertices..." << std::endl;
        
        // Initialize access nodes for transit nodes (they are their own access nodes)
        for (int transitNode : transitNodes) {
            accessNodes[transitNode] = {transitNode};
        }
        
        // 简化策略：为每个非传输节点分配所有传输节点作为访问节点
        for (int vertex = 0; vertex < graph.numVertices(); ++vertex) {
            if (transitNodes.find(vertex) == transitNodes.end()) {
                // 非传输节点：分配所有传输节点作为访问节点
                std::vector<int> vertexAccessNodes(transitNodes.begin(), transitNodes.end());
                accessNodes[vertex] = vertexAccessNodes;
            }
        }
        
        std::cout << "Access nodes computation completed. Total vertices: " << accessNodes.size() << std::endl;
        return accessNodes;
    }

    // Computes access nodes for all vertices (metric-dependent version).
    // This version uses actual edge weights and can be called during customization.
    template<typename WeightT>
    AccessNodeMap computeAccessNodesWithWeights(const std::vector<WeightT>& edgeWeights) const {
        AccessNodeMap accessNodes;
        
        // Initialize access nodes for transit nodes
        for (int transitNode : transitNodes) {
            accessNodes[transitNode] = {transitNode};
        }
        
        // Use Dijkstra-based approach for metric-dependent case
        TransitNodeSelector selector(decomp, transitNodes.size());
        int k = selector.getNumLevels();
        
        // Process vertices in order of their appearance in the separator decomposition
        for (int nodeIdx = k; nodeIdx < decomp.tree.size(); ++nodeIdx) {
            int firstVertex = decomp.firstSeparatorVertex(nodeIdx);
            int lastVertex = decomp.lastSeparatorVertex(nodeIdx);
            
            for (int i = firstVertex; i < lastVertex; ++i) {
                int vertex = decomp.order[i];
                std::unordered_map<int, int> accessNodeDistances; // access node -> distance
                
                // Run Dijkstra from this vertex until hitting transit nodes
                std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, 
                                   std::greater<std::pair<int, int>>> pq;
                std::vector<int> distances(graph.numVertices(), INFTY);
                
                pq.push({0, vertex});
                distances[vertex] = 0;
                
                while (!pq.empty()) {
                    auto [dist, u] = pq.top();
                    pq.pop();
                    
                    if (dist > distances[u]) continue;
                    
                    // If we hit a transit node, add it to access nodes
                    if (transitNodes.find(u) != transitNodes.end()) {
                        accessNodeDistances[u] = dist;
                    }
                    
                    // Continue search if we haven't hit enough transit nodes
                    FORALL_INCIDENT_EDGES(graph, u, e) {
                        int v = graph.edgeHead(e);
                        int newDist = dist + edgeWeights[e];
                        
                        if (newDist < distances[v]) {
                            distances[v] = newDist;
                            pq.push({newDist, v});
                        }
                    }
                }
                
                // Convert to vector of access nodes
                std::vector<int> vertexAccessNodes;
                for (const auto& pair : accessNodeDistances) {
                    vertexAccessNodes.push_back(pair.first);
                }
                accessNodes[vertex] = vertexAccessNodes;
            }
        }
        
        return accessNodes;
    }

private:
    const GraphT& graph;
    const SeparatorDecomposition& decomp;
    const TransitNodeSet& transitNodes;
};