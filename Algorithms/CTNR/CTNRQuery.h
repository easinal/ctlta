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
    using LabelSet = BasicLabelSet<0, ParentInfo::NO_PARENT_INFO>;

    // Constructor
    CTNRQuery(const TransitNodeHierarchy &hierarchy, CTNRData &data, const CCH &cch,
              const CH &minimumWeightedCH)
            : hierarchy(hierarchy),
              data(data),
              ETquery(minimumWeightedCH, cch.getEliminationTree()) {}

    // Main query method (s, t are rank IDs)
    int32_t run(int32_t s, int32_t t) {
        const auto lcaLevel = hierarchy.getLevelOfLowestCommonAncestor(s, t);
        if (lcaLevel >= hierarchy.getTransitNodeThreshold()) {
            lastModeIsLocal = true;
            return localQuery(s, t);
        } else {
            lastModeIsLocal = false;
            const int32_t tnDist = transitNodeQuery(s, t);
//            KASSERT(tnDist == localQuery(s,t));
            return tnDist;
        }
    }

    int32_t getDistance() const { return lastDistance; }

    const char *getLastMode() const { return lastModeIsLocal ? "local" : "transit"; }

private:
    const TransitNodeHierarchy &hierarchy;
    const CTNRData &data;
    int32_t lastDistance = INFTY;
    bool lastModeIsLocal = true;

    EliminationTreeQuery<LabelSet> ETquery;

    // Local query using elimination tree
    int32_t localQuery(int32_t s, int32_t t) {
        //TODO: use distance bound from transit node query
        ETquery.run(s, t);
        lastDistance = ETquery.getDistance();
        return lastDistance;
    }

    // Transit node query using three-hop approach
    int32_t transitNodeQuery(const int32_t s, const int32_t t) {
        int32_t minDist = INFTY;

        // Access arrays are indexed by rank IDs
        const auto &accessNodesS = data.getForwardAccessNodes(s);
        const auto &accessNodesT = data.getBackwardAccessNodes(t);

        for (const auto& as : accessNodesS) {
            if (as.distance >= minDist) continue;

            for (const auto& at : accessNodesT) {
                const int32_t mid = data.getDistanceBetweenTransitNodes(as.nodeIndex, at.nodeIndex);
                const int32_t total = as.distance + mid + at.distance;
                if (total < minDist)
                    minDist = total;
            }
        }
        lastDistance = minDist;
        return minDist;
    }
};

