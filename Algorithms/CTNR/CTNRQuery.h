#pragma once

#include "Algorithms/CTNR/CTNRMetric.h"
#include "Algorithms/CCH/EliminationTreeQuery.h"
#include "DataStructures/Labels/BasicLabelSet.h"
#include "DataStructures/Labels/ParentInfo.h"
#include "Tools/Constants.h"
#include <memory>
#include <algorithm>
#include <climits>

template<typename InputGraphT>
class CTNRQuery {
public:
    using InputGraph = InputGraphT;
    using Metric = CTNRMetric<InputGraphT>;
    using LabelSet = BasicLabelSet<0, ParentInfo::NO_PARENT_INFO>;

    // Constructor
    CTNRQuery(const Metric& metric)
        : metric(metric),
          ETquery(metric.getMinCH(), metric.getCCH().getEliminationTree()) {
        forwardAccessNodes = &metric.getForwardAccessNodes();
        forwardAccessDistances = &metric.getForwardAccessDistances();
        backwardAccessNodes = &metric.getBackwardAccessNodes();
        backwardAccessDistances = &metric.getBackwardAccessDistances();
        transitNodeToDistanceTableIndex = &metric.gettransitNodeToDistanceTableIndex();
        distanceTable = &metric.getDistanceTable();
    }

    // Main query method (s, t are rank IDs)
    int32_t run(int32_t s, int32_t t) {
        const uint32_t lca = metric.getHierarchy().getLowestCommonHub(s, t);
        int lcaDepth = metric.separatorNodeToLevel[lca];
        if( lcaDepth > metric.getTransitNodeThreshold()) {
            lastModeIsLocal = true;
            std::cout<<"lca: "<<lca<<", depth: "<<lcaDepth<<", is not a transit node"<<std::endl;
            return localQuery(s, t);
        }else{
            lastModeIsLocal = false;
            std::cout<<"lca: "<<lca<<", depth: "<<lcaDepth<<", is a transit node"<<std::endl;
            return transitNodeQuery(s, t, lcaDepth);
        }
    }
    int32_t getDistance() const { return lastDistance; }
    const char* getLastMode() const { return lastModeIsLocal ? "local" : "transit"; }

private:
    const Metric& metric;
    int32_t lastDistance = INFTY;
    bool lastModeIsLocal = true;
    
    const std::vector<std::vector<int32_t>>* forwardAccessNodes;
    const std::vector<std::vector<int32_t>>* forwardAccessDistances;
    const std::vector<std::vector<int32_t>>* backwardAccessNodes;
    const std::vector<std::vector<int32_t>>* backwardAccessDistances;
    const std::unordered_map<int32_t, int32_t>* transitNodeToDistanceTableIndex;
    const std::vector<std::vector<int32_t>>* distanceTable;
    EliminationTreeQuery<LabelSet> ETquery;
    // Local query using elimination tree
    int32_t localQuery(int32_t s, int32_t t) {
        //TODO: use distance bound from transit node query
        ETquery.run(s, t);
        lastDistance = ETquery.getDistance();
        return lastDistance;
    }

    // Transit node query using three-hop approach
    int32_t transitNodeQuery(int32_t s, int32_t t, int lcaNodeLevel) {
        int32_t minDist = INFTY;
        
        // Access arrays are indexed by rank IDs
        const auto &aS = (*forwardAccessNodes)[s]; 
        const auto &dS = (*forwardAccessDistances)[s];
        const auto &aT = (*backwardAccessNodes)[t];
        const auto &dT = (*backwardAccessDistances)[t];
        size_t candS = aS.size();
        size_t candT = aT.size();
        size_t evaluated = 0;
        std::cout<<"lcaNodeLevel: "<<lcaNodeLevel<<std::endl;

        // const auto byLevel = [&](int level, int v){
        //     auto it = metric.transitVertexToLevel.find(v);
        //     const int vLevel = (it == metric.transitVertexToLevel.end()) ? INT_MAX : it->second;
        //     return level < vLevel;
        // };
        // const size_t sBound = std::upper_bound(aS.begin(), aS.end(), lcaNodeLevel, byLevel) - aS.begin();
        // const size_t tBound = std::upper_bound(aT.begin(), aT.end(), lcaNodeLevel, byLevel) - aT.begin();
        const size_t sBound = aS.size();
        const size_t tBound = aT.size();

        for (int i = 0; i<sBound; ++i) {
            if (dS[i] >= minDist) continue;
            auto itS = transitNodeToDistanceTableIndex->find(aS[i]);
            if (itS == transitNodeToDistanceTableIndex->end()) continue;
            const int idxS = itS->second;


            for (int j = 0; j<tBound; ++j) {
                if (dT[j] >= minDist) continue;            
                auto itT = transitNodeToDistanceTableIndex->find(aT[j]);
                if (itT == transitNodeToDistanceTableIndex->end()) continue;
                const int idxT = itT->second;

                const int32_t mid = (*distanceTable)[idxS][idxT];
                if (mid >= minDist) { 
                    continue;
                }
                const int32_t total = dS[i] + mid + dT[j];
                ++evaluated;
                if (total < minDist) minDist = total;
            }
        }
        std::cout << "candS: " << candS << ", candT: " << candT << ", validS: " << sBound << ", validT: " << tBound << ", evaluated: " << evaluated << std::endl;
        lastDistance = minDist;
        return minDist;
    }
};

