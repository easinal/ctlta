#pragma once

#include "Algorithms/CTL/TruncatedTreeLabelling.h"
#include "Algorithms/CTL/BalancedTopologyCentricTreeHierarchy.h"
#include "Algorithms/CCH/CCHMetric.h"
#include <vectorclass.h>
#include "Tools/Constants.h"

#ifndef CTL_USE_PERFECT_CUSTOMIZATION
#define CTL_USE_PERFECT_CUSTOMIZATION false
#endif

// An individual metric for truncated tree labelling. Uses a CCH to build customized hub labelling.
template<typename LabellingT, typename LabelSetT, bool USE_PERFECT_CUSTOMIZATION = CTL_USE_PERFECT_CUSTOMIZATION>
class CTLMetric {

public:

    using SearchGraph = std::conditional_t<USE_PERFECT_CUSTOMIZATION, CH::SearchGraph, CCH::UpGraph>;

    // Constructs an individual metric incorporating the specified input weights in the specified
    // BalancedTopologyCentricTreeHierarchy on the basis of the specified CCH.
    CTLMetric(const BalancedTopologyCentricTreeHierarchy &hierarchy, const CCH &cch, const int32_t *const inputWeights)
            : hierarchy(hierarchy), cch(cch), cchMetric(cch, inputWeights), minimumWeightedCH() {}

    void buildCustomizedCTL(LabellingT &ctl) {
        if constexpr (USE_PERFECT_CUSTOMIZATION)
            minimumWeightedCH = cchMetric.buildMinimumWeightedCH();
        else
            cchMetric.customize();
        customizeLabelling(ctl);
    }

    // Returns the upward graph of the metric, which is either the CCH graph or the upward graph after perfect
    // customization.
    const SearchGraph &upwardGraph() const {
        if constexpr (USE_PERFECT_CUSTOMIZATION)
            return minimumWeightedCH.upwardGraph();
        else
            return cch.getUpwardGraph();
    }

    // Returns the downward graph of the metric, which is either the CCH graph or the downward graph after perfect
    // customization.
    const SearchGraph &downwardGraph() const {
        if constexpr (USE_PERFECT_CUSTOMIZATION)
            return minimumWeightedCH.downwardGraph();
        else
            return cch.getUpwardGraph();
    }

    // Returns current weights on edges of upward graph.
    int const *upwardWeights() const {
        if constexpr (USE_PERFECT_CUSTOMIZATION)
            return &minimumWeightedCH.upwardGraph().template get<CH::Weight>(0);
        else
            return cchMetric.upWeights.data();
    }

    // Returns current weights on edges of downward graph.
    int const *downwardWeights() const {
        if constexpr (USE_PERFECT_CUSTOMIZATION)
            return &minimumWeightedCH.downwardGraph().template get<CH::Weight>(0);
        else
            return cchMetric.downWeights.data();
    }

    uint64_t sizeInBytes() const {
        return sizeof(*this) + cchMetric.sizeInBytes() + minimumWeightedCH.sizeInBytes();
    }

private:

    using LabelSet = LabelSetT;
    static constexpr int K = LabelSet::K;
    using Batch = typename LabelSet::DistanceLabel;
    using BatchMask = typename LabelSet::LabelMask;

    void customizeLabelling(LabellingT &ctl) {
        ctl.reset();

        const auto &upGraph = upwardGraph();
        const auto &downGraph = downwardGraph(); // reverse downward graph
        const auto upWeights = upwardWeights();
        const auto downWeights = downwardWeights();

        cch.forEachVertexTopDown([&](const int u) {
            // Do not build labels for truncated vertices.
            if (hierarchy.isVertexTruncated(u))
                return;

            BatchMask improved;
            Batch dU, dV, eU;
            const auto numHubsU = hierarchy.getNumHubs(u);

            // Customize upward label of u using upper neighbors
            auto uUpLabel = ctl.upLabel(u);
            uUpLabel.initializeLastHubDist(); // distance to self
            if constexpr (LabelSet::KEEP_PARENT_EDGES)
                uUpLabel.initializeLastHubPathEdge(); // edge to self
            int32_t* startUUp = uUpLabel.startDists();
            int32_t* startEdgesUUp = uUpLabel.startEdges();
            FORALL_INCIDENT_EDGES(upGraph, u, e) {
                const auto v = upGraph.edgeHead(e);
                const auto upWeight = upWeights[e];
                const auto numHubsV = hierarchy.getNumHubs(v);
                KASSERT(numHubsV < numHubsU);
                KASSERT(numHubsV == hierarchy.getLowestCommonHub(u, v));
                const auto vUpLabel = ctl.cUpLabel(v);
                int const * const startVUp = vUpLabel.startDists();
                const uint32_t numBlocks = numHubsV / K + (numHubsV % K != 0);
                const Batch weightAsBatch(upWeight);
                const Batch eAsBatch(e);
                for (uint32_t b = 0; b < numBlocks; ++b) {
                    dV.load(startVUp + b * K);
                    dU.load(startUUp + b * K);
                    const auto dNew = weightAsBatch + dV;
                    improved = dNew < dU;
                    dU = select(improved, dNew, dU);
                    dU.store(startUUp + b * K);
                    if constexpr (LabelSet::KEEP_PARENT_EDGES) {
                        eU.load(startEdgesUUp + b * K);
                        eU = select(improved, eAsBatch, eU);
                        eU.store(startEdgesUUp + b * K);
                    }
                }
            }

            // Customize (reverse) downward label of u using upper neighbors
            auto uDownLabel = ctl.downLabel(u);
            uDownLabel.initializeLastHubDist(); // distance to self
            if constexpr (LabelSet::KEEP_PARENT_EDGES)
                uDownLabel.initializeLastHubPathEdge(); // edge to self
            int32_t* startUDown = uDownLabel.startDists();
            int32_t* startEdgesUDown = uDownLabel.startEdges();
            FORALL_INCIDENT_EDGES(downGraph, u, e) {
                const auto v = downGraph.edgeHead(e);
                const auto downWeight = downWeights[e];
                const auto numHubsV = hierarchy.getNumHubs(v);
                KASSERT(numHubsV < numHubsU);
                KASSERT(numHubsV == hierarchy.getLowestCommonHub(u, v));
                const auto vDownLabel = ctl.cDownLabel(v);
                int32_t const * const startVDown = vDownLabel.startDists();
                const uint32_t numBlocks = numHubsV / K + (numHubsV % K != 0);
                const Batch weightAsBatch(downWeight);
                const Batch eAsBatch(e);
                for (uint32_t b = 0; b < numBlocks; ++b) {
                    dV.load(startVDown + b * K);
                    dU.load(startUDown + b * K);
                    const auto dNew = weightAsBatch + dV;
                    improved = dNew < dU;
                    dU = select(improved, dNew, dU);
                    dU.store(startUDown + b * K);
                    if constexpr (LabelSet::KEEP_PARENT_EDGES) {
                        eU.load(startEdgesUDown + b * K);
                        eU = select(improved, eAsBatch, eU);
                        eU.store(startEdgesUDown + b * K);
                    }
                }
            }
        });
    }

    const BalancedTopologyCentricTreeHierarchy &hierarchy;
    const CCH &cch;
    CCHMetric cchMetric;

    CH minimumWeightedCH;
};

