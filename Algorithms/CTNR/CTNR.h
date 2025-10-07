#pragma once

#include "Algorithms/CTNR/CTNRMetric.h"
#include "Algorithms/CTNR/CTNRQuery.h"
#include "DataStructures/Partitioning/SeparatorDecomposition.h"
#include <memory>

template<typename InputGraphT>
class CTNR {
public:
    using InputGraph = InputGraphT;
    using Metric = CTNRMetric<InputGraphT>;
    using Query = CTNRQuery<InputGraphT>;

    // Constructor
    CTNR(const SeparatorDecomposition& sepDecomp, int transitNodeThreshold)
        : metric(sepDecomp, transitNodeThreshold), queryEngine(metric) {
    }

    // Preprocessing phase
    void preprocess(const InputGraph& inputGraph) {
        metric.preprocess(inputGraph);
    }

    // Customization phase
    void customize(const int32_t* inputWeights) {
        metric.customize(inputWeights);
    }

    // Query phase
    int32_t query(int32_t s, int32_t t) {
        return queryEngine.run(s, t);
    }

    // Getters
    const Metric& getMetric() const { return metric; }
    const Query& getQuery() const { return queryEngine; }
    const BalancedTopologyCentricTreeHierarchy& getHierarchy() const { return metric.getHierarchy(); }
    const CCH& getCCH() const { return metric.getCCH(); }
    const std::vector<int32_t>& getTransitNodes() const { return metric.getTransitNodes(); }
    
    // Memory usage
    uint64_t sizeInBytes() const {
        return metric.sizeInBytes() + queryEngine.sizeInBytes() + sizeof(CTNR<InputGraphT>);
    }

private:
    Metric metric;
    Query queryEngine;
};