#pragma once

#include "Algorithms/CTL/BalancedTopologyCentricTreeHierarchy.h"
#include "Algorithms/CCH/CCH.h"
#include "Algorithms/CCH/CCHMetric.h"
#include "Algorithms/CH/CH.h"
#include "DataStructures/Partitioning/SeparatorDecomposition.h"
#include "Tools/Constants.h"
#include <vector>
#include <unordered_map>
#include <map>
#include <memory>

template<typename InputGraphT>
class CTNRMetric {
public:
    using InputGraph = InputGraphT;
    
    // Constructor
    CTNRMetric(const SeparatorDecomposition& sepDecomp, int transitNodeThreshold)
        : sepDecomp(sepDecomp), hierarchy(), cch(), cchMetric(nullptr), 
          transitNodeThreshold(transitNodeThreshold) {
    }

    // Preprocessing phase
    void preprocess(const InputGraph& inputGraph) {
        cch.preprocess(inputGraph, sepDecomp);
        hierarchy.preprocess(inputGraph, sepDecomp);
        
        selectTransitNodes();
        
        forwardAccessNodes.resize(inputGraph.numVertices());
        forwardAccessDistances.resize(inputGraph.numVertices());
        backwardAccessNodes.resize(inputGraph.numVertices());
        backwardAccessDistances.resize(inputGraph.numVertices());
        
    }

    // Customization phase
    void customize(const int32_t* inputWeights) {
        
        cchMetric = std::make_unique<CCHMetric>(cch, inputWeights);
        cchMetric->customize();
        minCH = cchMetric->buildMinimumWeightedCH();
        computeAccessNodes();
        computeDistanceTable();
        pruneAccessNodesByDominance();
    }

    // Getters
    const BalancedTopologyCentricTreeHierarchy& getHierarchy() const { return hierarchy; }
    const CCH& getCCH() const { return cch; }
    const std::vector<int32_t>& getTransitNodes() const { return transitNodes; }
    const std::vector<std::vector<int32_t>>& getForwardAccessNodes() const { return forwardAccessNodes; }
    const std::vector<std::vector<int32_t>>& getForwardAccessDistances() const { return forwardAccessDistances; }
    const std::vector<std::vector<int32_t>>& getBackwardAccessNodes() const { return backwardAccessNodes; }
    const std::vector<std::vector<int32_t>>& getBackwardAccessDistances() const { return backwardAccessDistances; }
    const std::vector<std::vector<int32_t>>& getDistanceTable() const { return distanceTable; }
    const std::unordered_map<int32_t, int32_t>& gettransitNodeToDistanceTableIndex() const { return transitNodeToDistanceTableIndex; }
    const CH& getMinCH() const { return minCH; }
    int getTransitNodeThreshold() const { return transitNodeThreshold; }
    // Memory usage calculation including node levels
    uint64_t sizeInBytes() const {
        uint64_t size = sizeof(CTNRMetric);
        size += hierarchy.sizeInBytes();
        size += cch.sizeInBytes();
        if (cchMetric) size += cchMetric->sizeInBytes();
        size += minCH.sizeInBytes();
        size += transitNodes.capacity() * sizeof(int32_t);
        
        for (const auto& vec : forwardAccessNodes) {
            size += vec.capacity() * sizeof(int32_t);
        }
        for (const auto& vec : forwardAccessDistances) {
            size += vec.capacity() * sizeof(int32_t);
        }
        for (const auto& vec : backwardAccessNodes) {
            size += vec.capacity() * sizeof(int32_t);
        }
        for (const auto& vec : backwardAccessDistances) {
            size += vec.capacity() * sizeof(int32_t);
        }
        
        for (const auto& vec : distanceTable) {
            size += vec.capacity() * sizeof(int32_t);
        }
        
        return size;
    }
    std::unordered_map<int32_t, int32_t> transitNodesIDToLevel;  // separator node ID -> level
    std::unordered_map<int32_t, int32_t> transitVertexToLevel;  // vertex ID -> level

private:
    // Core data structures
    SeparatorDecomposition sepDecomp;
    BalancedTopologyCentricTreeHierarchy hierarchy;
    CCH cch;
    std::unique_ptr<CCHMetric> cchMetric;


    // Transit Node related
    int transitNodeThreshold;
    std::vector<int32_t> transitNodes;
    std::unordered_map<int32_t, int32_t> transitNodeToDistanceTableIndex; // key: TN rank id

    // Access Nodes (indexed by rank IDs)
    std::vector<std::vector<int32_t>> forwardAccessNodes;  // forwardAccessNodes[rank(v)] = TN ranks
    std::vector<std::vector<int32_t>> forwardAccessDistances;  // distances
    // Backward Access (a -> v), indexed by rank IDs
    std::vector<std::vector<int32_t>> backwardAccessNodes;      // backwardAccessNodes[rank(v)] = TN ranks
    std::vector<std::vector<int32_t>> backwardAccessDistances;  // distances

    // Distance table
    std::vector<std::vector<int32_t>> distanceTable;  // distanceTable[i][j] = distance from node i to node j
    CH minCH;


    // Helper methods
    void selectTransitNodes();
    void computeAccessNodes();
    void computeDistanceTable();
    void pruneAccessNodesByDominance();
    int32_t getTransitNodeDistance(int32_t accessS, int32_t accessT) const;
    
};

// Template implementation
#include "Algorithms/CCH/EliminationTreeQuery.h"
#include "Algorithms/CH/CHQuery.h"
#include "DataStructures/Labels/BasicLabelSet.h"
#include "DataStructures/Labels/ParentInfo.h"
#include <algorithm>
#include <iostream>

template<typename InputGraphT>
void CTNRMetric<InputGraphT>::selectTransitNodes() {
    std::function<void(int, int)> collectTransitNodes = 
        [&](int node, int level) {
            if (level <= transitNodeThreshold) {
                transitNodesIDToLevel[node] = level;
                for (int v = sepDecomp.firstSeparatorVertex(node);
                     v < sepDecomp.lastSeparatorVertex(node); ++v) {
                        transitNodes.push_back(v);
                        transitVertexToLevel[v] = level; 
                }
            } else {
                // return;
                transitNodesIDToLevel[node] = level;
                for (int v = sepDecomp.firstSeparatorVertex(node);
                     v < sepDecomp.lastSeparatorVertex(node); ++v) {
                        transitVertexToLevel[v] = level;  
                }
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
    std::sort(transitNodes.begin(), transitNodes.end(),compareByLevel);
    for(int i = 0; i < transitNodes.size(); ++i) {
        transitNodeToDistanceTableIndex[transitNodes[i]] = i;
    }
    std::cout << "CTNR: Selected " << transitNodes.size() << " transit nodes from top " << transitNodeThreshold << " levels" << std::endl;
}

template<typename InputGraphT>
void CTNRMetric<InputGraphT>::computeAccessNodes() {

    auto compareByLevel = [&](int32_t a, int32_t b) {
        return transitVertexToLevel[a] < transitVertexToLevel[b];
    };
    cch.forEachVertexTopDown([&](int32_t rv) {
        std::unordered_map<int, int> fMin;
        std::unordered_map<int, int> bMin;

        if(transitVertexToLevel.find(rv) != transitVertexToLevel.end()) {
            fMin[rv] = 0;
            bMin[rv] = 0;
        }

        FORALL_INCIDENT_EDGES(cch.getUpwardGraph(), rv, e) {
            const int neighbor = cch.getUpwardGraph().edgeHead(e);
            const int wUp = cchMetric->upwardWeights()[e];
            const int wDown = cchMetric->downwardWeights()[e];

            const auto& fa = forwardAccessNodes[neighbor];
            const auto& fd = forwardAccessDistances[neighbor];
            for (size_t i = 0; i < fa.size(); ++i) {
                const int rTN = fa[i];
                const int dist = fd[i] + wUp;
                auto it = fMin.find(rTN);
                if (it == fMin.end() || dist < it->second) fMin[rTN] = dist;
            }

            const auto& ba = backwardAccessNodes[neighbor];
            const auto& bd = backwardAccessDistances[neighbor];
            for (size_t i = 0; i < ba.size(); ++i) {
                const int rTN = ba[i];
                const int dist = bd[i] + wDown;
                auto it = bMin.find(rTN);
                if (it == bMin.end() || dist < it->second) bMin[rTN] = dist;
            }
        }
        forwardAccessNodes[rv].clear();
        forwardAccessDistances[rv].clear();
        backwardAccessNodes[rv].clear();
        backwardAccessDistances[rv].clear();

        forwardAccessNodes[rv].reserve(fMin.size());
        forwardAccessDistances[rv].reserve(fMin.size());
        for (const auto& kv : fMin) { forwardAccessNodes[rv].push_back(kv.first); }

        backwardAccessNodes[rv].reserve(bMin.size());
        backwardAccessDistances[rv].reserve(bMin.size());
        for (const auto& kv : bMin) { backwardAccessNodes[rv].push_back(kv.first); }
        sort(forwardAccessNodes[rv].begin(), forwardAccessNodes[rv].end(), compareByLevel);
        sort(backwardAccessNodes[rv].begin(), backwardAccessNodes[rv].end(), compareByLevel);

        for (auto &node: forwardAccessNodes[rv]) {
            forwardAccessDistances[rv].push_back(fMin[node]);
        }
        for (auto &node: backwardAccessNodes[rv]) {
            backwardAccessDistances[rv].push_back(bMin[node]);
        }
    });
}

template<typename InputGraphT>
void CTNRMetric<InputGraphT>::computeDistanceTable() {
    const int n = (int)transitNodes.size();
    distanceTable.assign(n, std::vector<int32_t>(n, INFTY));
    using LabelSet = BasicLabelSet<0, ParentInfo::NO_PARENT_INFO>;
    #pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        CHQuery<LabelSet> chq(minCH);
        for (int j = 0; j < n; ++j) {
            if (i == j) { distanceTable[i][j] = 0; continue; }
            chq.run(transitNodes[i], transitNodes[j]);
            distanceTable[i][j] = chq.getDistance();
            // std::cout<<"distanceTable["<<i<<"]["<<j<<"]: "<<distanceTable[i][j]<<std::endl;
        }
    }
}

// template<typename InputGraphT>
// void CTNRMetric<InputGraphT>::pruneAccessNodesByDominance(auto &accessNodes, auto &accessDistances) {
    
//     for (int32_t v = 0; v < forwardAccessNodes.size(); ++v) {
//         auto& an = accessNodes[v];
//         auto& ad = accessDistances[v];
        
//         if (an.size() <= 1) continue;
//         std::vector<int> keep(an.size(), 1);
        
//         for (size_t i = 0; i < an.size(); ++i) {
//             if (!keep[i]) continue;
//             for (size_t j = 0; j < an.size(); ++j) {
//                 if (i == j || !keep[j]) continue;
//                 int32_t dij = this->getTransitNodeDistance(an[i], an[j]);
//                 if (dij == INFTY) continue;
//                 if (ad[i] + dij <= ad[j]) keep[j] = 0;
//             }
//         }
//         size_t w = 0;
//         for (size_t i = 0; i < an.size(); ++i) if (keep[i]) {
//             an[w] = an[i];
//             ad[w] = ad[i];
//             ++w;
//         }
//         an.resize(w);
//         ad.resize(w);
//     }
// }

template<typename InputGraphT>
void CTNRMetric<InputGraphT>::pruneAccessNodesByDominance() {
    auto pruneOne = [&](std::vector<int32_t>& nodes, std::vector<int32_t>& dists) {
        if (nodes.size() <= 1) return;
        std::vector<int> keep(nodes.size(), 1);
        for (size_t i = 0; i < nodes.size(); ++i) {
            if (!keep[i]) continue;
            for (size_t j = 0; j < nodes.size(); ++j) {
                if (i == j || !keep[j]) continue;
                int32_t dij = this->getTransitNodeDistance(nodes[i], nodes[j]);
                if (dij == INFTY) continue;
                if (dists[i] + dij <= dists[j]) keep[j] = 0;
            }
        }
        size_t w = 0;
        for (size_t i = 0; i < nodes.size(); ++i) if (keep[i]) {
            nodes[w] = nodes[i];
            dists[w] = dists[i];
            ++w;
        }
        nodes.resize(w);
        dists.resize(w);
    };

    for (int32_t v = 0; v < (int32_t)forwardAccessNodes.size(); ++v) {
        pruneOne(forwardAccessNodes[v], forwardAccessDistances[v]);
    }
    for (int32_t v = 0; v < (int32_t)backwardAccessNodes.size(); ++v) {
        pruneOne(backwardAccessNodes[v], backwardAccessDistances[v]);
    }
}

template<typename InputGraphT>
int32_t CTNRMetric<InputGraphT>::getTransitNodeDistance(int32_t accessS, int32_t accessT) const {
    auto itS = transitNodeToDistanceTableIndex.find(accessS);
    auto itT = transitNodeToDistanceTableIndex.find(accessT);
    if (itS == transitNodeToDistanceTableIndex.end() || itT == transitNodeToDistanceTableIndex.end()) {
        return INFTY;
    }
    return distanceTable[itS->second][itT->second];
}