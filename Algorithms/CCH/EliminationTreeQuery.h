#pragma once

#include <array>
#include <cassert>
#include <cstdint>
#include <utility>
#include <vector>

#include "Algorithms/CH/CH.h"
#include "Algorithms/CCH/UpwardEliminationTreeSearch.h"
#include "Algorithms/Dijkstra/Dijkstra.h"
#include "Tools/Constants.h"

// An implementation of an elimination tree query, which computes shortest paths in a customizable
// contraction hierarchy without using priority queues. Depending on the used label set, it keeps
// parent vertices and/or edges, and computes multiple shortest paths simultaneously, optionally
// using SSE or AVX instructions.
//
// While the run() methods interleave the forward and reverse elimination tree search, the methods
// pinForwardSearch() and runReverseSearch() allow us to perform the forward and reverse search one
// after another. This is useful when performing multiple queries from the same source (or sources)
// in succession. In that case, it suffices to perform the forward search once and use its distance
// labels for multiple reverse searches.
template<typename LabelSetT, bool UseCH = true>
class EliminationTreeQuery {
private:
    using DistanceLabel = typename LabelSetT::DistanceLabel;
    using ParentLabel = typename LabelSetT::ParentLabel;

    // The pruning criterion for an elimination tree query that computes k shortest paths
    // simultaneously. We can prune the search at v if d_i(v) >= mu_i for all i = 1, ..., k.
    struct PruningCriterion {
        // Constructs a pruning criterion for an elimination tree query.
        PruningCriterion(const DistanceLabel &tentativeDistances) noexcept
                : tentativeDistances(&tentativeDistances) {}

        // Returns true if the search can be pruned at v.
        template<typename DistanceLabelContT>
        bool operator()(const int, const DistanceLabel &distToV, const DistanceLabelContT &) const {
            return allSet(distToV >= *tentativeDistances);
        }

        const DistanceLabel *tentativeDistances; // One tentative distance per simultaneous search.
    };

    static constexpr int K = LabelSetT::K; // The number of simultaneous shortest-path computations.

    static constexpr bool DO_NOT_USE_FAST_QUERY =
#ifdef NO_FAST_ELIMINATION_TREE_QUERY
            true;
#else
            false;
#endif

public:

    // Constructs an elimination tree query instance for a given CH (useful for using minimum weighted CH obtained
    // from perfect customization of CCH).
    EliminationTreeQuery(const CH &ch, const std::vector<int32_t> &eliminTree) requires (UseCH && DO_NOT_USE_FAST_QUERY)
            : forwardSearch(ch.upwardGraph(), eliminTree),
              reverseSearch(ch.downwardGraph(), eliminTree) {
        assert(ch.upwardGraph().numVertices() == eliminTree.size());
    }

    // Constructs an elimination tree query instance for a given CH (useful for using minimum weighted CH obtained
    // from perfect customization of CCH).
    EliminationTreeQuery(const CH &ch, const std::vector<int32_t> &eliminTree) requires (UseCH && !DO_NOT_USE_FAST_QUERY)
            : forwardSearch(ch.upwardGraph(), eliminTree, {tentativeDistances}),
              reverseSearch(ch.downwardGraph(), eliminTree, {tentativeDistances}) {
        assert(ch.upwardGraph().numVertices() == eliminTree.size());
    }

    // Constructs an elimination tree query instance for a CCH after basic customization, using the whole CCH graph with
    // metric-independent shortcuts.
    EliminationTreeQuery(const CCH::UpGraph& graph,
                         int const * const upWeights,
                         int const * const downWeights,
                         const std::vector<int32_t> &eliminTree) requires (!UseCH && DO_NOT_USE_FAST_QUERY)
            : forwardSearch(graph, upWeights, eliminTree),
              reverseSearch(graph, downWeights, eliminTree) {}

    // Constructs an elimination tree query instance for a CCH after basic customization, using the whole CCH graph with
    // metric-independent shortcuts.
    EliminationTreeQuery(const CCH::UpGraph& graph,
                         int const * const upWeights,
                         int const * const downWeights,
                         const std::vector<int32_t> &eliminTree) requires (!UseCH && !DO_NOT_USE_FAST_QUERY)
            : forwardSearch(graph, upWeights, eliminTree, {tentativeDistances}),
              reverseSearch(graph, downWeights, eliminTree, {tentativeDistances}) {}

    // Move constructor.
    EliminationTreeQuery(EliminationTreeQuery &&other) noexcept
            : forwardSearch(std::move(other.forwardSearch)),
              reverseSearch(std::move(other.reverseSearch)) {
        if constexpr (!DO_NOT_USE_FAST_QUERY) {
            forwardSearch.pruneSearch = {tentativeDistances};
            reverseSearch.pruneSearch = {tentativeDistances};
        }
    }

    // Runs an elimination tree query from s to t.
    void run(const int s, const int t) {
        std::array<int, K> sources;
        std::array<int, K> targets;
        sources.fill(s);
        targets.fill(t);
        run(sources, targets);
    }

    // Runs an elimination tree query that computes multiple shortest paths simultaneously.
    void run(const std::array<int, K> &sources, const std::array<int, K> &targets) {
        forwardSearch.init(sources);
        reverseSearch.init(targets);
        tentativeDistances = INFTY;
        while (forwardSearch.nextVertices.minKey() != INVALID_VERTEX &&
               reverseSearch.nextVertices.minKey() != INVALID_VERTEX) {
            if (forwardSearch.nextVertices.minKey() <= reverseSearch.nextVertices.minKey()) {
                updateTentativeDistances(forwardSearch.nextVertices.minKey());
                forwardSearch.distanceLabels[forwardSearch.settleNextVertex()] = INFTY;
            } else {
                reverseSearch.distanceLabels[reverseSearch.settleNextVertex()] = INFTY;
            }
        }

        // Reset distance labels of rest of branches in case searches did not meet.
        // This is only needed if the query is run on subgraphs of the CH and elimination tree.
        // Otherwise, only the label of the root in the reverse search would remain to be reset.
        //        reverseSearch.distanceLabels[reverseSearch.searchGraph.numVertices() - 1] = INFTY;
        int v = forwardSearch.nextVertices.minKey();
        while (v != INVALID_VERTEX) {
            forwardSearch.distanceLabels[v] = INFTY;
            v = forwardSearch.nextVertex();
        }
        v = reverseSearch.nextVertices.minKey();
        while (v != INVALID_VERTEX) {
            reverseSearch.distanceLabels[v] = INFTY;
            v = reverseSearch.nextVertex();
        }
    }

    // Runs a forward search from s and pins (stores) its distance labels.
    void pinForwardSearch(const int s) {
        pinForwardSearch(&s, &s + 1);
    }

    // Runs a forward search from multiple sources and pins (stores) its distance labels.
    template<typename IteratorT>
    void pinForwardSearch(const IteratorT firstSource, const IteratorT lastSource) {
        tentativeDistances = INFTY;
        forwardSearch.run(firstSource, lastSource);
    }

    // Runs a reverse search from t, which considers the pinned forward labels.
    void runReverseSearch(const int t) {
        runReverseSearch(&t, &t + 1);
    }

    // Runs a reverse search from multiple targets, which considers the pinned forward labels.
    template<typename IteratorT>
    void runReverseSearch(const IteratorT firstTarget, const IteratorT lastTarget) {
        reverseSearch.init(firstTarget, lastTarget);
        tentativeDistances = INFTY;
        while (reverseSearch.nextVertices.minKey() != INVALID_VERTEX) {
            updateTentativeDistances(reverseSearch.nextVertices.minKey());
            reverseSearch.distanceLabels[reverseSearch.settleNextVertex()] = INFTY;
        }
    }

    // Returns the length of the i-th shortest path.
    int getDistance(const int i = 0) {
        return tentativeDistances[i];
    }

    // Returns the edges in the upward graph on the up segment of the up-down path (in reverse order).
    const std::vector<int32_t> &getUpEdgePath(const int i = 0) {
        assert(tentativeDistances[i] != INFTY);
        return forwardSearch.getReverseEdgePath(meetingVertices.vertex(i), i);
    }

    // Returns the edges in the downward graph on the down segment of the up-down path.
    const std::vector<int32_t> &getDownEdgePath(const int i = 0) {
        assert(tentativeDistances[i] != INFTY);
        return reverseSearch.getReverseEdgePath(meetingVertices.vertex(i), i);
    }

    uint64_t sizeInBytes() const {
        return forwardSearch.sizeInBytes()
               + reverseSearch.sizeInBytes()
               + sizeof(tentativeDistances)
               + sizeof(meetingVertices);
    }

private:
    // Checks whether the path via v improves the tentative distance for any search.
    void updateTentativeDistances(const int v) {
        const auto distances = forwardSearch.distanceLabels[v] + reverseSearch.distanceLabels[v];
        meetingVertices.setVertex(v, distances < tentativeDistances);
        tentativeDistances.min(distances);
    }

    template<typename Container>
    static bool checkAllInfty(const Container &cont, const size_t size) {
        for (int i = 0; i < size; ++i)
            if (!allSet(cont[i] == INFTY)) {
                KASSERT(false);
                return false;
            }
        return true;
    }

    using SearchGraphT = std::conditional_t<UseCH, CH::SearchGraph, CCH::UpGraph>;

    using UpwardSearch = std::conditional_t<DO_NOT_USE_FAST_QUERY,
            UpwardEliminationTreeSearch<LabelSetT, elimintree::PruningCriterion, SearchGraphT>,
            UpwardEliminationTreeSearch<LabelSetT, PruningCriterion, SearchGraphT>
    >;

    UpwardSearch forwardSearch;       // The forward search from the source(s).
    UpwardSearch reverseSearch;       // The reverse search from the target(s).

    DistanceLabel tentativeDistances; // One tentative distance per simultaneous search.
    ParentLabel meetingVertices;      // One meeting vertex per simultaneous search.
};
