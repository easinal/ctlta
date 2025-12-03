#pragma once

#include "Algorithms/CTNR/CTNRMetric.h"
#include "Algorithms/CCH/EliminationTreeQuery.h"
#include "DataStructures/Labels/BasicLabelSet.h"
#include "DataStructures/Labels/ParentInfo.h"
#include "Tools/Constants.h"
#include "CTNRConstants.h"
#include <memory>
#include <algorithm>
#include <climits>

template<typename InputGraphT>
class CTNRQuery {
public:
    using LabelSet = BasicLabelSet<0, ParentInfo::NO_PARENT_INFO>;

    // Constructor
    CTNRQuery(const TransitNodeHierarchy &hierarchy, CTNRData &data,
              const std::vector<int>& localEliminationTree,
              const CCH::UpGraph &cchGraph,
              int const * const cchUpWeights,
              int const * const cchDownWeights)
            : hierarchy(hierarchy),
              data(data),
              localQuery(cchGraph, cchUpWeights, cchDownWeights, localEliminationTree) {}

    // Main query method (s, t are rank IDs)
    int32_t run(int32_t s, int32_t t) {
        const ctnr::Level lcaLevel = hierarchy.getLevelOfLowestCommonAncestor(s, t);
        int32_t dist = runTransitNodeQuery(s, t);
        if (lcaLevel >= hierarchy.getTransitNodeThreshold()) {
            lastModeIsLocal = true;
            dist = std::min(dist, runLocalQuery(s, t));
        } else {
            lastModeIsLocal = false;
        }

        lastDistance = dist;
        return dist;
    }

    int32_t getDistance() const { return lastDistance; }

    const char *getLastMode() const { return lastModeIsLocal ? "local" : "transit"; }

    uint64_t sizeInBytes() const {
        uint64_t size = sizeof(CTNRQuery);
        size += localQuery.sizeInBytes();
        return size;
    }

private:

    // Local query using elimination tree
    int32_t runLocalQuery(int32_t s, int32_t t) {
        localQuery.run(s, t);
        return localQuery.getDistance();
    }

    // Transit node query using three-hop approach
    int32_t runTransitNodeQuery(const int32_t s, const int32_t t) {
        int32_t minDist = CTNR_INFTY;


        // Access arrays are indexed by rank IDs
        const auto &accessNodesS = data.getAccessNodes(s);
        const auto &accessDistancesS = data.getForwardDistances(s);
        const auto &accessNodesT = data.getAccessNodes(t);
        const auto &accessDistancesT = data.getBackwardDistances(t);

        const bool sIsTransit = hierarchy.isTransitNode(s);
        const bool tIsTransit = hierarchy.isTransitNode(t);
        if (sIsTransit && tIsTransit) {
            const int32_t dist = data.getDistanceBetweenTransitNodes(
                    hierarchy.getTransitNodeIndexOfRank(s),
                    hierarchy.getTransitNodeIndexOfRank(t));
            return dist;
        }
        if (sIsTransit) {
            const ctnr::TransitNodeId nodeS = hierarchy.getTransitNodeIndexOfRank(s);
            const auto numAccessT = accessNodesT.size();
            for (auto j = 0; j < numAccessT; ++j) {
                const ctnr::TransitNodeId nodeT = accessNodesT[j];
                const int32_t distT = accessDistancesT[j];
                const int32_t mid = data.getDistanceBetweenTransitNodes(nodeS, nodeT);
                const int32_t total = mid + distT;
                if (total < minDist)
                    minDist = total;
            }
            return minDist;
        }
        if (tIsTransit) {
            const ctnr::TransitNodeId nodeT = hierarchy.getTransitNodeIndexOfRank(t);
            const auto numAccessS = accessNodesS.size();
            for (auto i = 0; i < numAccessS; ++i) {
                const ctnr::TransitNodeId nodeS = accessNodesS[i];
                const int32_t distS = accessDistancesS[i];
                const int32_t mid = data.getDistanceBetweenTransitNodes(nodeS, nodeT);
                const int32_t total = distS + mid;
                if (total < minDist)
                    minDist = total;
            }
            return minDist;
        }

        const auto numAccessS = accessNodesS.size();
        const auto numAccessT = accessNodesT.size();
        for (auto i = 0; i < numAccessS; ++i) {
            const ctnr::TransitNodeId nodeS = accessNodesS[i];
            const int32_t distS = accessDistancesS[i];
            if (distS >= minDist) continue;

            for (auto j = 0; j < numAccessT; ++j) {
                const ctnr::TransitNodeId nodeT = accessNodesT[j];
                const int32_t distT = accessDistancesT[j];
//                if (distT >= CTNR_INFTY) continue;
//                if (distS + distT >= minDist) continue;
                const int32_t mid = data.getDistanceBetweenTransitNodes(nodeS, nodeT);
                const int32_t total = distS + mid + distT;
                if (total < minDist)
                    minDist = total;
            }
        }
        return minDist;
    }


    const TransitNodeHierarchy &hierarchy;
    const CTNRData &data;
    int32_t lastDistance = CTNR_INFTY;
    bool lastModeIsLocal = true;

    EliminationTreeQuery<LabelSet, false> localQuery;
};

