#pragma once

#include "Algorithms/CTL/TruncatedTreeLabelling.h"
#include "Algorithms/Dijkstra/DagShortestPaths.h"

template<typename SearchGraphT, typename LabellingT, typename LabelSetT>
class CTLQuery {

    using LabelSet = LabelSetT;
    static constexpr uint32_t K = LabelSet::K;
    using BatchMask = LabelSet::LabelMask;
    using Batch = typename LabelSet::DistanceLabel;

    static constexpr bool NoTruncatedVertices = BalancedTopologyCentricTreeHierarchy::NoTruncatedVertices;

    // Temporary label being constructed for truncated vertices.
    struct TemporaryLabel {

        TemporaryLabel() = default;

        void init(const size_t numHubs) {
            _numHubs = numHubs;
            dists.assign(numHubs, INFTY);
            if constexpr (LabelSet::KEEP_PARENT_EDGES)
                accessVertices.assign(numHubs, INVALID_VERTEX);
        }

        const int32_t& dist(const uint32_t &hubIdx) const {
            KASSERT(hubIdx < numHubs());
            return dists[hubIdx];
        }

        int32_t& dist(const uint32_t &hubIdx) {
            KASSERT(hubIdx < _numHubs);
            return dists[hubIdx];
        }

        const int &accessVertex(const uint32_t &hubIdx) const requires LabelSet::KEEP_PARENT_EDGES {
            KASSERT(hubIdx < _numHubs);
            return accessVertices[hubIdx];
        }

        int &accessVertex(const uint32_t &hubIdx) requires LabelSet::KEEP_PARENT_EDGES {
            KASSERT(hubIdx < _numHubs);
            return accessVertices[hubIdx];
        }

        uint32_t numHubs() const {
            return _numHubs;
        }

        uint64_t sizeInBytes() const {
            return sizeof(TemporaryLabel) + dists.capacity() * sizeof(Batch) + accessVertices.capacity() * sizeof(Batch);
        }

    private:
        uint32_t _numHubs;
        AlignedVector<int32_t> dists;
        AlignedVector<int32_t> accessVertices;
    };

    template<bool UP>
    struct PruneSearchAtUntruncatedVertices {

        PruneSearchAtUntruncatedVertices(TemporaryLabel &temporaryLabel,
                                         std::vector<int32_t> &truncatedSearchSpace,
                                         const BalancedTopologyCentricTreeHierarchy &hierarchy,
                                         const LabellingT &ctl)
                : temporaryLabel(temporaryLabel), truncatedSearchSpace(truncatedSearchSpace), hierarchy(hierarchy),
                  ctl(ctl) {}

        template<typename DistanceLabelT, typename DistanceLabelContainerT>
        bool operator()(const int v, DistanceLabelT &distToV, DistanceLabelContainerT &) {

            // If v is truncated, store it in truncated search space, and do not prune.
            if (hierarchy.isVertexTruncated(v)) {
                truncatedSearchSpace.push_back(v);
                return false;
            }

            // TODO: manually vectorize (?)
            // Update temporary label of source / target with label of v.
            const auto labelOfV = UP ? ctl.upLabel(v) : ctl.downLabel(v);
            const auto endHub = std::min(temporaryLabel.numHubs(), labelOfV.numHubs);
            for (int i = 0; i < endHub; ++i) {
                const auto distViaV = distToV[0] + labelOfV.dist(i);
                auto& distTemp = temporaryLabel.dist(i);
                if (distViaV < distTemp) {
                    distTemp = distViaV;
                    if constexpr (LabelSet::KEEP_PARENT_EDGES)
                        temporaryLabel.accessVertex(i) = v;
                }
            }

            // Prune at untruncated vertices.
            return true;
        }

        TemporaryLabel &temporaryLabel;
        std::vector<int32_t> &truncatedSearchSpace;
        const BalancedTopologyCentricTreeHierarchy &hierarchy;
        const LabellingT &ctl;
    };

    struct MockSearch {
        MockSearch(const SearchGraphT &, const int32_t *, const PruneSearchAtUntruncatedVertices<true> &) {}

        uint64_t sizeInBytes() const {
            return 0; // No memory used in mock search.
        }
    };

    using TruncatedVertexUpwardSearch = std::conditional_t<NoTruncatedVertices, MockSearch, DagShortestPaths<SearchGraphT, BasicLabelSet<0, ParentInfo::FULL_PARENT_INFO>, PruneSearchAtUntruncatedVertices<true>>>;
    using TruncatedVertexDownwardSearch = std::conditional_t<NoTruncatedVertices, MockSearch, DagShortestPaths<SearchGraphT, BasicLabelSet<0, ParentInfo::FULL_PARENT_INFO>, PruneSearchAtUntruncatedVertices<false>>>;



public:

    CTLQuery(const BalancedTopologyCentricTreeHierarchy &hierarchy,
             const SearchGraphT &upGraph,
             const SearchGraphT &downGraph,
             int const *const upWeights,
             int const *const downWeights,
             const LabellingT &ctl)
            : hierarchy(hierarchy), upGraph(upGraph), downGraph(downGraph), ctl(ctl),
              buildUpLabelSearch(upGraph, upWeights, {tempUpLabel, upTruncatedSearchSpace, hierarchy, ctl}),
              buildDownLabelSearch(downGraph, downWeights, {tempDownLabel, downTruncatedSearchSpace, hierarchy, ctl}) {}

    // Expects ranks in the underlying separator decomposition order as inputs.
    void run(const int32_t s, const int32_t t) {
        lastDistance = INFTY;
        lastMeetingHubIdx = INVALID_INDEX;
        minCCHDistMeetingVertex = INVALID_VERTEX;
        lastS = s;
        lastT = t;
        const auto lch = hierarchy.getLowestCommonHub(s, t);

        if constexpr (NoTruncatedVertices) {
            const auto sUpLabel = ctl.cUpLabel(s);
            const auto tDownLabel = ctl.cDownLabel(t);
            computeMinDistanceInLabels(sUpLabel, tDownLabel, lch);
        } else {
            if (!hierarchy.isVertexTruncated(s) && !hierarchy.isVertexTruncated(t)) {
                const auto sUpLabel = ctl.cUpLabel(s);
                const auto tDownLabel = ctl.cDownLabel(t);
                computeMinDistanceInLabels(sUpLabel, tDownLabel, lch);
            } else if (!hierarchy.isVertexTruncated(s)) {
                const auto sUpLabel = ctl.cUpLabel(s);
                buildTempDownLabel(t, lch);
                computeMinDistanceInLabels(sUpLabel, tempDownLabel, lch);
            } else if (!hierarchy.isVertexTruncated(t)) {
                buildTempUpLabel(s, lch);
                const auto tDownLabel = ctl.cDownLabel(t);
                computeMinDistanceInLabels(tempUpLabel, tDownLabel, lch);
            } else {
                buildTempUpLabel(s, lch);
                buildTempDownLabel(t, lch);
                computeMinDistanceInLabels(tempUpLabel, tempDownLabel, lch);

                // If both s and t are truncated, they may be within the same truncated separator subtree, which we can
                // check by comparing their number of hubs with the lch.
                // In this case, their shortest distance may not use a vertex that is high enough in the hierarchy
                // to be non-truncated. In this case, the shortest path is found using just the distances found during
                // the two elimination tree searches in the CCH. We only have to consider the truncated vertices in
                // the search space of the higher ranked vertex between s and t.
                if (hierarchy.getNumHubs(s) == lch && hierarchy.getNumHubs(t) == lch) {
                    const auto &searchSpace = s > t ? upTruncatedSearchSpace : downTruncatedSearchSpace;
                    for (const auto &v: searchSpace) {
                        const auto cchDist = buildUpLabelSearch.getDistance(v) + buildDownLabelSearch.getDistance(v);
                        if (cchDist < lastDistance) {
                            lastDistance = cchDist;
                            minCCHDistMeetingVertex = v;
                        }
                    }
                }
            }
        }

    }

    // Returns the distance of the last query.
    int32_t getDistance() const {
        return lastDistance;
    }

    int32_t getLastMeetingHubIdx() const {
        return lastMeetingHubIdx;
    }

    // Returns the CCH edges in the upward graph on the up segment of the up-down path (in reverse order to conform to
    // default orientation in graph-traversal-based searches).
    template<bool hasPathEdges = LabelSet::KEEP_PARENT_EDGES, std::enable_if_t<hasPathEdges, bool> = true>
    const std::vector<int32_t> &getUpEdgePath() {
        lastUpPath.clear();

        if constexpr (!NoTruncatedVertices) {
            // If the best distance was found using not labels but the CCH directly, return the path to the meeting
            // vertex in the CCH:
            if (minCCHDistMeetingVertex != INVALID_VERTEX) {
                lastUpPath = buildUpLabelSearch.getReverseEdgePath(minCCHDistMeetingVertex);
                return lastUpPath;
            }
        }

        int32_t v;
        if constexpr (NoTruncatedVertices) {
            v = lastS;
        } else {
            if (hierarchy.isVertexTruncated(lastS)) {
                // If s was truncated, get path to access vertex, i.e., the first non-truncated vertex used on the up path,
                // using parent pointers in topo search.
                const auto accVertex = tempUpLabel.accessVertex(lastMeetingHubIdx);
                lastUpPath = buildUpLabelSearch.getReverseEdgePath(accVertex);
                std::reverse(lastUpPath.begin(), lastUpPath.end());
                v = accVertex;
            } else {
                // If s was not truncated, the access vertex is simply s itself.
                v = lastS;
            }
        }


        // Enumerate path from access vertex to meeting hub using path edges stored in labels.
        int32_t e;
        while ((e = ctl.upPathEdge(v, lastMeetingHubIdx)) != INVALID_EDGE) {
            KASSERT(e >= 0 && e < upGraph.numEdges());
            lastUpPath.push_back(e);
            v = upGraph.edgeHead(e);
        }
        std::reverse(lastUpPath.begin(), lastUpPath.end());
        return lastUpPath;
    }

    // Returns the CCH edges in the upward graph on the up segment of the up-down path. Edges are returned in
    // arbitrary order. To retrieve edges in order of path, use getUpEdgePath() instead.
    template<bool hasPathEdges = LabelSet::KEEP_PARENT_EDGES, std::enable_if_t<hasPathEdges, bool> = true>
    const std::vector<int32_t> &getEdgesOnUpPathUnordered() {
        lastUpPath.clear();

        if constexpr (!NoTruncatedVertices) {
            // If the best distance was found using not labels but the CCH directly, return the path to the meeting
            // vertex in the CCH:
            if (minCCHDistMeetingVertex != INVALID_VERTEX) {
                lastUpPath = buildUpLabelSearch.getReverseEdgePath(minCCHDistMeetingVertex);
                return lastUpPath;
            }
        }

        int32_t v;
        if constexpr (NoTruncatedVertices) {
            v = lastS;
        } else {
            if (hierarchy.isVertexTruncated(lastS)) {
                // If s was truncated, get path to access vertex, i.e., the first non-truncated vertex used on the up path,
                // using parent pointers in topo search.
                const auto accVertex = tempUpLabel.accessVertex(lastMeetingHubIdx);
                lastUpPath = buildUpLabelSearch.getReverseEdgePath(accVertex);
                v = accVertex;
            } else {
                // If s was not truncated, the access vertex is simply s itself.
                v = lastS;
            }
        }

        // Enumerate path from access vertex to meeting hub using path edges stored in labels.
        int32_t e;
        while ((e = ctl.upPathEdge(v, lastMeetingHubIdx)) != INVALID_EDGE) {
            KASSERT(e >= 0 && e < upGraph.numEdges());
            lastUpPath.push_back(e);
            v = upGraph.edgeHead(e);
        }
        return lastUpPath;
    }

    // Returns the CCH edges in the upward graph on the down segment of the up-down path.
    template<bool hasPathEdges = LabelSet::KEEP_PARENT_EDGES, std::enable_if_t<hasPathEdges, bool> = true>
    const std::vector<int32_t> &getDownEdgePath() {
        lastDownPath.clear();

        if constexpr (!NoTruncatedVertices) {
            // If the best distance was found using not labels but the CCH directly, return the path to the meeting
            // vertex in the CCH:
            if (minCCHDistMeetingVertex != INVALID_VERTEX) {
                lastDownPath = buildDownLabelSearch.getReverseEdgePath(minCCHDistMeetingVertex);
                return lastDownPath;
            }
        }

        int32_t v;
        if constexpr (NoTruncatedVertices) {
            v = lastT;
        } else {
            if (hierarchy.isVertexTruncated(lastT)) {
                // If t was truncated, get path to access vertex, i.e., the first non-truncated vertex used on the down path,
                // using parent pointers in topo search.
                const auto accVertex = tempDownLabel.accessVertex(lastMeetingHubIdx);
                lastDownPath = buildDownLabelSearch.getReverseEdgePath(accVertex);
                std::reverse(lastDownPath.begin(), lastDownPath.end());
                v = accVertex;
            } else {
                // If t was not truncated, the access vertex is simply t itself.
                v = lastT;
            }
        }

        // Enumerate path from access vertex to meeting hub using path edges stored in labels.
        int32_t e;
        while ((e = ctl.downPathEdge(v, lastMeetingHubIdx)) != INVALID_EDGE) {
            KASSERT(e >= 0 && e < downGraph.numEdges());
            lastDownPath.push_back(e);
            v = downGraph.edgeHead(e);
        }
        std::reverse(lastDownPath.begin(), lastDownPath.end());
        return lastDownPath;
    }

    // Returns the CCH edges in the upward graph on the down segment of the up-down path.  Edges are returned in
    // arbitrary order. To retrieve edges in order of path, use getDownEdgePath() instead.
    template<bool hasPathEdges = LabelSet::KEEP_PARENT_EDGES, std::enable_if_t<hasPathEdges, bool> = true>
    const std::vector<int32_t> &getEdgesOnDownPathUnordered() {
        lastDownPath.clear();

        if constexpr (!NoTruncatedVertices) {
            // If the best distance was found using not labels but the CCH directly, return the path to the meeting
            // vertex in the CCH:
            if (minCCHDistMeetingVertex != INVALID_VERTEX) {
                lastDownPath = buildDownLabelSearch.getReverseEdgePath(minCCHDistMeetingVertex);
                return lastDownPath;
            }
        }

        int32_t v;
        if constexpr (NoTruncatedVertices) {
            v = lastT;
        } else {
            if (hierarchy.isVertexTruncated(lastT)) {
                // If t was truncated, get path to access vertex, i.e., the first non-truncated vertex used on the down path,
                // using parent pointers in topo search.
                const auto accVertex = tempDownLabel.accessVertex(lastMeetingHubIdx);
                lastDownPath = buildDownLabelSearch.getReverseEdgePath(accVertex);
                v = accVertex;
            } else {
                // If t was not truncated, the access vertex is simply t itself.
                v = lastT;
            }
        }


        // Enumerate path from access vertex to meeting hub using path edges stored in labels.
        int32_t e;
        while ((e = ctl.downPathEdge(v, lastMeetingHubIdx)) != INVALID_EDGE) {
            KASSERT(e >= 0 && e < downGraph.numEdges());
            lastDownPath.push_back(e);
            v = downGraph.edgeHead(e);
        }
        return lastDownPath;
    }

    uint64_t sizeInBytes() const {
        uint64_t size = sizeof(CTLQuery);
        size += sizeof(lastS) + sizeof(lastT) + sizeof(lastDistance) + sizeof(lastMeetingHubIdx) + sizeof(minCCHDistMeetingVertex);
        size += lastUpPath.capacity() * sizeof(int32_t) + lastDownPath.capacity() * sizeof(int32_t);
        size += tempUpLabel.sizeInBytes();
        size += tempDownLabel.sizeInBytes();
        size += buildUpLabelSearch.sizeInBytes();
        size += buildDownLabelSearch.sizeInBytes();
        size += upTruncatedSearchSpace.capacity() * sizeof(int32_t);
        size += downTruncatedSearchSpace.capacity() * sizeof(int32_t);
        return size;
    }

private:

    // Populates tempUpLabel with up label for truncated vertex v by running topological upwards search from v.
    // Whenever search runs into a non-truncated vertex w, the label of w is used to update tempUpLabel and the
    // search is pruned.
    // Populates distances only for hubs 0..lch (inclusive).
    void buildTempUpLabel(const int32_t v, const uint32_t lch) {
        KASSERT(lch <= hierarchy.getNumHubs(v));
        tempUpLabel.init(lch);
        upTruncatedSearchSpace.clear();
        buildUpLabelSearch.run(v);
    }

    // Populates tempDownLabel with down label for truncated vertex v by running topological reverse-downwards search
    // from v. Whenever search runs into a non-truncated vertex w, the label of w is used to update tempDownLabel and
    // the search is pruned.
    // Populates distances only for hubs 0..lch (inclusive).
    void buildTempDownLabel(const int32_t v, const uint32_t lch) {
        KASSERT(lch <= hierarchy.getNumHubs(v));
        tempDownLabel.init(lch);
        downTruncatedSearchSpace.clear();
        buildDownLabelSearch.run(v);
    }


    template<typename UpLabel, typename DownLabel>
    inline void computeMinDistanceInLabels(const UpLabel &up, const DownLabel &down, const uint32_t lowestCommonHub) {

        // Compute minimum distance by combining labels
        for (int i = 0; i < lowestCommonHub; ++i) {
            const auto newDist = up.dist(i) + down.dist(i);
            if (newDist < lastDistance) {
                lastDistance = newDist;
                lastMeetingHubIdx = i;
            }
        }
    }

    const BalancedTopologyCentricTreeHierarchy &hierarchy;
    const SearchGraphT &upGraph;
    const SearchGraphT &downGraph;
    const LabellingT &ctl;

    int32_t lastS;
    int32_t lastT;
    int32_t lastDistance;
    uint32_t lastMeetingHubIdx;
    int32_t minCCHDistMeetingVertex; // if shortest path does not use labels but CCH distances, this is set to meeting vertex in CCH
    std::vector<int32_t> lastUpPath;
    std::vector<int32_t> lastDownPath;

    TemporaryLabel tempUpLabel;
    TemporaryLabel tempDownLabel;
    TruncatedVertexUpwardSearch buildUpLabelSearch;
    TruncatedVertexDownwardSearch buildDownLabelSearch;
    std::vector<int32_t> upTruncatedSearchSpace; // All truncated vertices that the last up search visited
    std::vector<int32_t> downTruncatedSearchSpace; // All truncated vertices that the last down search visited

};

