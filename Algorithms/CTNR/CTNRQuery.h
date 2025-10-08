#pragma once

#include "Algorithms/CTNR/CTNRMetric.h"
#include "Algorithms/CCH/EliminationTreeQuery.h"
#include "DataStructures/Labels/BasicLabelSet.h"
#include "DataStructures/Labels/ParentInfo.h"
#include "Tools/Constants.h"
#include <memory>

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
        const int lcaDepth = metric.getHierarchy().getVertexDepth(lca);
        std::cout<<"lca: "<<lca<<", lcaDepth: "<<lcaDepth<<", s_depth: "<<metric.getHierarchy().getVertexDepth(s)<<", t_depth: "<<metric.getHierarchy().getVertexDepth(t)<<std::endl;
        const bool local = (lcaDepth > metric.getTransitNodeThreshold());
        if (!local) {
            lastModeIsLocal = false;
            return transitNodeQuery(s, t, lcaDepth);
        } else {
            lastModeIsLocal = true;
            return localQuery(s, t);
        }
    }
    int32_t getDistance() const { return lastDistance; }
    const char* getLastMode() const { return lastModeIsLocal ? "local" : "transit"; }

    // Memory usage calculation
    uint64_t sizeInBytes() const {
        uint64_t size = sizeof(CTNRQuery);
        size += ETquery.sizeInBytes();
        return size;
    }

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

    // Local query using CH elimination tree (s, t are rank IDs)
    int32_t localQuery(int32_t s, int32_t t) {
        // Use CHQuery on minCH; CH uses original vertex IDs
        ETquery.run(s, t);
        lastDistance = ETquery.getDistance();
        return lastDistance;
    }

    // Transit node query using three-hop approach
    int32_t transitNodeQuery(int32_t s, int32_t t, uint32_t lca) {
        int32_t minDist = INFTY;
        
        // Access arrays are indexed by rank IDs
        const auto &aS = (*forwardAccessNodes)[s];   // contents: original TN ids
        const auto &dS = (*forwardAccessDistances)[s];
        const auto &aT = (*backwardAccessNodes)[t];  // contents: original TN ids
        const auto &dT = (*backwardAccessDistances)[t];
        std::cout<<"s: "<<s<<" depth: "<<metric.getHierarchy().getVertexDepth(s)<<", t: "<<t<<" depth: "<<metric.getHierarchy().getVertexDepth(t)<<std::endl;
        // const int lcaDepth = metric.getHierarchy().getNumHubs(lca);
        std::cout<<"lca: "<<lca<<std::endl;
        // Diagnostics counters
        size_t candS = aS.size();
        size_t candT = aT.size();
        size_t prunedSByLCA = 0;
        size_t prunedTByLCA = 0;
        size_t missingIdx = 0;
        size_t infMid = 0;
        size_t evaluated = 0;
        
        for (size_t i = 0; i < aS.size(); ++i) {
            std::cout<<"aS[i]: "<<aS[i]<<", depth: "<<metric.getHierarchy().getNumHubs(aS[i])<<std::endl;
            if (metric.getHierarchy().getNumHubs(aS[i]) > lca) {
                ++prunedSByLCA;
                continue;
            }
            if (dS[i] == INFTY) continue;
            auto itS = transitNodeToDistanceTableIndex->find(aS[i]);
            if (itS == transitNodeToDistanceTableIndex->end()) { ++missingIdx; continue; }
            const int idxS = itS->second;


            for (size_t j = 0; j < aT.size(); ++j) {
                std::cout<<"aT[j]: "<<aT[j]<<", depth: "<<metric.getHierarchy().getNumHubs(aT[j])<<std::endl;
                if (metric.getHierarchy().getNumHubs(aT[j]) > lca) {
                    ++prunedTByLCA; continue;
                }
                if (dT[j] == INFTY) continue;
                auto itT = transitNodeToDistanceTableIndex->find(aT[j]);
                if (itT == transitNodeToDistanceTableIndex->end()) { ++missingIdx; continue; }
                const int idxT = itT->second;

                const int32_t mid = (*distanceTable)[idxS][idxT];
                if (mid == INFTY) { ++infMid; continue; }
                const int32_t total = dS[i] + mid + dT[j];
                ++evaluated;
                if (total < minDist) minDist = total;
            }
        }
        std::cout << "candS: " << candS << ", candT: " << candT << ", prunedSByLCA: " << prunedSByLCA << ", prunedTByLCA: " << prunedTByLCA << ", missingIdx: " << missingIdx << ", infMid: " << infMid << ", evaluated: " << evaluated << std::endl;
        (void)candS; (void)candT; (void)prunedSByLCA; (void)prunedTByLCA; (void)missingIdx; (void)infMid; (void)evaluated;

        lastDistance = minDist;
        return minDist;
    }
};