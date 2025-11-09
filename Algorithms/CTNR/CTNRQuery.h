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
    CTNRQuery(const TransitNodeHierarchy &hierarchy, CTNRData &data,
              const std::vector<int>& localEliminationTree,
              const CH &localMinimumWeightedCH)
            : hierarchy(hierarchy),
              data(data),
              localQuery(localMinimumWeightedCH, localEliminationTree) {}

    // Main query method (s, t are rank IDs)
    int32_t run(int32_t s, int32_t t) {
        const auto lcaLevel = hierarchy.getLevelOfLowestCommonAncestor(s, t);
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

private:
    const TransitNodeHierarchy &hierarchy;
    const CTNRData &data;
    int32_t lastDistance = INFTY;
    bool lastModeIsLocal = true;

    EliminationTreeQuery<LabelSet> localQuery;

    // Local query using elimination tree
    int32_t runLocalQuery(int32_t s, int32_t t) {
        localQuery.run(s, t);
        return localQuery.getDistance();
    }

    // Transit node query using three-hop approach
    int32_t runTransitNodeQuery(const int32_t s, const int32_t t) {
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
        return minDist;
    }
};

