#pragma once

#include "Algorithms/CTL/BalancedTopologyCentricTreeHierarchy.h"
#include "Algorithms/CCH/CCH.h"
#include "Algorithms/CCH/CCHMetric.h"
#include "Algorithms/CH/CH.h"
#include "DataStructures/Partitioning/SeparatorDecomposition.h"
#include "Tools/Constants.h"
#include "Algorithms/CTNR/CTNRData.h"
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <memory>

class CTNRMetric {
public:
    
    // Constructor
    CTNRMetric(const CCH& cch, const int32_t *const inputWeights)
        : cch(cch), cchMetric(cch, inputWeights) {}

//    // Preprocessing phase
//    void preprocess(const InputGraph& inputGraph) {
//        cch.preprocess(inputGraph, sepDecomp);
//        hierarchy.preprocess(inputGraph, sepDecomp);
//
//        selectTransitNodes();
//
//        forwardAccessNodes.resize(inputGraph.numVertices());
//        forwardAccessDistances.resize(inputGraph.numVertices());
//        backwardAccessNodes.resize(inputGraph.numVertices());
//        backwardAccessDistances.resize(inputGraph.numVertices());
//
//    }

    // Customization phase
    void customize(CTNRData& data) {
        int64_t dummy1, dummy2, dummy3, dummy4;
        customizeWithMeasurements(data, dummy1, dummy2, dummy3, dummy4);
    }

    // Sets measurement parameters to times for each step in microseconds.
    void customizeWithMeasurements(CTNRData& data, int64_t& cchCustomizationTime, int64_t& accessNodeComputationTime,
                                  int64_t& distanceTableComputationTime, int64_t& accessNodePruningTime) {
        KASSERT(data.forwardAccessNodes.size() == cch.getUpwardGraph().numVertices(), "Data not initialized.");
        Timer timer;
        minCH = cchMetric.buildMinimumWeightedCH();
        cchCustomizationTime = timer.elapsed<std::chrono::microseconds>();
        timer.restart();
        computeAccessNodes(data);
        accessNodeComputationTime = timer.elapsed<std::chrono::microseconds>();
        timer.restart();
        computeDistanceTable(data);
        distanceTableComputationTime = timer.elapsed<std::chrono::microseconds>();
        timer.restart();
        pruneAccessNodesByDominance(data);
        accessNodePruningTime = timer.elapsed<std::chrono::microseconds>();
    }

    // Getters
    const CH& getMinCH() const { return minCH; }
    // Memory usage calculation including node levels
    uint64_t sizeInBytes() const {
        uint64_t size = sizeof(CTNRMetric);
        size += cchMetric.sizeInBytes();
        size += minCH.sizeInBytes();
        
        return size;
    }

private:
    // Core data structures
    const CCH& cch;
    CCHMetric cchMetric;
    CH minCH;

    // Helper methods
    void computeAccessNodes(CTNRData& data);
    void computeDistanceTable(CTNRData &data);
    void pruneAccessNodesByDominance(CTNRData &data);
//    int32_t getTransitNodeDistance(int32_t accessS, int32_t accessT) const;
    
};

// Template implementation
#include "Algorithms/CCH/EliminationTreeQuery.h"
#include "Algorithms/CH/CHQuery.h"
#include "DataStructures/Labels/BasicLabelSet.h"
#include "DataStructures/Labels/ParentInfo.h"
#include <algorithm>
#include <iostream>



void CTNRMetric::computeAccessNodes(CTNRData& data) {

    auto compareByLevel = [&](int32_t a, int32_t b) {
        return data.transitVertexToLevel[a] < data.transitVertexToLevel[b];
    };
    cch.forEachVertexTopDown([&](int32_t rv) {
        if(data.transitVertexToLevel.find(rv) != data.transitVertexToLevel.end()) {
            data.forwardAccessNodes[rv].push_back(rv);
            data.forwardAccessDistances[rv].push_back(0);
            data.backwardAccessNodes[rv].push_back(rv);
            data.backwardAccessDistances[rv].push_back(0);
        }else{
            std::unordered_map<int, int> fMin;
            std::unordered_map<int, int> bMin;

            FORALL_INCIDENT_EDGES(cch.getUpwardGraph(), rv, e) {
                const int neighbor = cch.getUpwardGraph().edgeHead(e);
                const int wUp = cchMetric.upwardWeights()[e];
                const int wDown = cchMetric.downwardWeights()[e];

                if(wUp != INFTY) {
                    const auto& fa = data.forwardAccessNodes[neighbor];
                    const auto& fd = data.forwardAccessDistances[neighbor];
                    for (size_t i = 0; i < fa.size(); ++i) {
                        const int rTN = fa[i];
                            const int dist = fd[i] + wUp;
                            auto it = fMin.find(rTN);
                            if (it == fMin.end() || dist < it->second) fMin[rTN] = dist;
                    }
                }

                if(wDown != INFTY) {
                    const auto& ba = data.backwardAccessNodes[neighbor];
                    const auto& bd = data.backwardAccessDistances[neighbor];
                    for (size_t i = 0; i < ba.size(); ++i) {
                        const int rTN = ba[i];
                            const int dist = bd[i] + wDown;
                            auto it = bMin.find(rTN);
                            if (it == bMin.end() || dist < it->second) bMin[rTN] = dist;   
                    }
                }
            }

            data.forwardAccessNodes[rv].clear();
            data.forwardAccessDistances[rv].clear();
            data.backwardAccessNodes[rv].clear();
            data.backwardAccessDistances[rv].clear();

            data.forwardAccessNodes[rv].reserve(fMin.size());
            data.forwardAccessDistances[rv].reserve(fMin.size());
            data.backwardAccessNodes[rv].reserve(bMin.size());
            data.backwardAccessDistances[rv].reserve(bMin.size());
            for (const auto& kv : fMin) { data.forwardAccessNodes[rv].push_back(kv.first); }
            for (const auto& kv : bMin) { data.backwardAccessNodes[rv].push_back(kv.first); }
            sort(data.forwardAccessNodes[rv].begin(), data.forwardAccessNodes[rv].end(), compareByLevel);
            sort(data.backwardAccessNodes[rv].begin(), data.backwardAccessNodes[rv].end(), compareByLevel);

            for (auto &node: data.forwardAccessNodes[rv]) {
                data.forwardAccessDistances[rv].push_back(fMin[node]);
            }
            for (auto &node: data.backwardAccessNodes[rv]) {
                data.backwardAccessDistances[rv].push_back(bMin[node]);
            }
            // if(rv%100000 == 0) {
            //     std::cout<<"Forward Last Level of "<<rv<<": "<<forwardLastLevel<<std::endl;
            //     std::cout<<"Backward Last Level of "<<rv<<": "<<backwardLastLevel<<std::endl;
            //     std::cout<<"Forward Access Nodes of "<<rv<<": "<<forwardAccessNodes[rv].size()<<std::endl;
            //     for(auto &node: forwardAccessNodes[rv]) {
            //         std::cout<<node<<" "<<"level: "<<transitVertexToLevel[node]<<" distance: "<<fMin[node]<<std::endl;
            //     }
            //     std::cout<<std::endl;   
            //     std::cout<<"Backward Access Nodes of "<<rv<<": "<<backwardAccessNodes[rv].size()<<std::endl;
            //     for(auto &node: backwardAccessNodes[rv]) {
            //         std::cout<<node<<" "<<"level: "<<transitVertexToLevel[node]<<" distance: "<<bMin[node]<< std::endl;
            //     }
            //     std::cout<<std::endl;
            // }
        }
    });
}
//TODO: use PHAST to accelerate distance table computation
void CTNRMetric::computeDistanceTable(CTNRData &data) {
    const int n = (int)data.transitNodes.size();
    data.distanceTable.assign(n, std::vector<int32_t>(n, INFTY));
    using LabelSet = BasicLabelSet<0, ParentInfo::NO_PARENT_INFO>;
#pragma omp parallel
    {
        EliminationTreeQuery<LabelSet> chq(minCH, cch.getEliminationTree());
#pragma omp for
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i == j) {
                    data.distanceTable[i][j] = 0;
                    continue;
                }
                chq.run(data.transitNodes[i], data.transitNodes[j]);
                data.distanceTable[i][j] = chq.getDistance();
                // std::cout<<"distanceTable["<<i<<"]["<<j<<"]: "<<distanceTable[i][j]<<std::endl;
            }
        }
    }
}

void CTNRMetric::pruneAccessNodesByDominance(CTNRData &data) {
    auto pruneOne = [&](std::vector<int32_t>& nodes, std::vector<int32_t>& dists, bool isForward) {
        for(int i = 0; i < nodes.size(); ++i) {
            nodes[i] = data.transitNodeToDistanceTableIndex[nodes[i]];
        }
        if (nodes.size() <= 1) return;
        std::vector<bool> keep(nodes.size(), true);
        for (size_t i = 0; i < nodes.size(); ++i) {
            if (!keep[i]) continue;
            for (size_t j = 0; j < nodes.size(); ++j) {
                if (i == j || !keep[j]) continue;
                int32_t dij;
                dij = data.getDistanceBetweenTransitNodes(nodes[i], nodes[j]);
                if (dij == INFTY) continue;
                if (isForward) {
                    if (dists[i] + dij <= dists[j]) keep[j] = false;
                }else{
                    if (dists[j] + dij <= dists[i]) keep[i] = false;
                }
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

    for (int32_t v = 0; v < (int32_t)data.forwardAccessNodes.size(); ++v) {
        pruneOne(data.forwardAccessNodes[v], data.forwardAccessDistances[v], true);
    }
    for (int32_t v = 0; v < (int32_t)data.backwardAccessNodes.size(); ++v) {
        pruneOne(data.backwardAccessNodes[v], data.backwardAccessDistances[v], false);
    }
}

//template<typename InputGraphT>
//int32_t CTNRMetric<InputGraphT>::getTransitNodeDistance(int32_t accessS, int32_t accessT) const {
//    return distanceTable[accessS][accessT];
//    // auto itS = transitNodeToDistanceTableIndex.find(accessS);
//    // auto itT = transitNodeToDistanceTableIndex.find(accessT);
//    // if (itS == transitNodeToDistanceTableIndex.end() || itT == transitNodeToDistanceTableIndex.end()) {
//    //     return INFTY;
//    // }
//    // return distanceTable[itS->second][itT->second];
//}
