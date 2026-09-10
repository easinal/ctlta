#pragma once

#include <kassert/kassert.hpp>
#include <stack>

#include "Algorithms/CCH/CCH.h"
#include "DataStructures/Partitioning/SeparatorDecompositionWalk.h"
#include "DataStructures/Partitioning/SeparatorTree.h"
#include "Algorithms/CCH/CCHMetric.h"

#include "Tools/Constants.h"

class TransitNodeHierarchy {

    using SideId = uint64_t;
    using TransitNodeId = ctnr::TransitNodeId;
    using Level = ctnr::Level;

    struct Identity {
        int operator[](const int x) const {
            return x;
        }
    };

public:

    TransitNodeHierarchy() = default;

    // Builds the transit node hierarchy for the specified graph and separator decomposition.
    // If specified, applies a permutation to the vertex IDs in the separator decomposition before processing them.
    // This is useful if the CCH graph uses a different vertex ordering.
    template<typename InputGraphT, typename PermuteVertexIdT = Identity>
    void preprocess(const InputGraphT &inputGraph, const int newTransitNodeThreshold,
                    const SeparatorDecomposition &sepDecomp, const PermuteVertexIdT& permuteVertexId = {}) {

        const auto sdDepth = sepdecomp::depth(sepDecomp);
        std::cout << "Depth of sepDecomp is " << sdDepth << std::endl;

        if (!sepdecomp::hasStrictDissectionStructure(sepDecomp))
            throw std::invalid_argument("TransitNodeHierarchy requires strict dissection "
                                        "structure of separator decomposition.");

        // Build labels and numCommonHubsComputer
        packedSideIds.clear();
        packedSideIds.resize(inputGraph.numVertices(), static_cast<uint64_t>(-1));
        vertexLevel.clear();
        vertexLevel.resize(inputGraph.numVertices(), -1);

        transitNodeThreshold = newTransitNodeThreshold;
        transitNodes.clear();
        transitNodeIndexOfRank.assign(inputGraph.numVertices(), -1);

        computeVertexLocationInSepDecomp(sepDecomp, permuteVertexId);

        KASSERT(std::all_of(packedSideIds.begin(), packedSideIds.end(),
                            [](const uint64_t &id) { return id != static_cast<uint64_t>(-1); }));

//        // Todo: remove debug
//        // Count number of vertices in each level and print
//        std::vector<size_t> levelCounts(sdDepth, 0);
//        for (const auto &level: vertexLevel) {
//            KASSERT(level < sdDepth);
//            ++levelCounts[level];
//        }
//        for (size_t level = 0; level < levelCounts.size(); ++level) {
//            std::cout << "SepDecomp Level " << level << ": " << levelCounts[level] << " vertices" << std::endl;
//        }
    }

    size_t numVertices() const {
        return packedSideIds.size();
    }

    Level getTransitNodeThreshold() const {
        return transitNodeThreshold;
    }

    // Given the rank of a vertex in the CCH-order, returns its level in the separator hierarchy.
    inline Level getVertexLevel(const int32_t &v) const {
        KASSERT(v >= 0 && v < vertexLevel.size());
        return vertexLevel[v];
    }

    // Given the ranks of two vertices in the CCH-order, returns the level of their lowest common ancestor
    // in the separator hierarchy.
    Level getLevelOfLowestCommonAncestor(const int32_t &s, const int32_t &t) const {

        const Level minInputLevel = std::min(getVertexLevel(s), getVertexLevel(t));

        // XOR packed side IDs to find out lowest common level in separator hierarchy.
        const Level l = static_cast<Level>(lowestOneBit(packedSideIds[s] ^ packedSideIds[t]));

        if (l >= 0)
            return std::min(l, minInputLevel);

        // If packed side IDs of s and t are exactly the same, the branch of s subsumes the branch of t or vice
        // versa. In this case, the lowest common ancestor is the lower one of the two vertices.
        return minInputLevel;
    }

    int32_t numTransitNodes() const {
        return static_cast<int32_t>(transitNodes.size());
    }

    bool isTransitNode(const int32_t &v) const {
        KASSERT(v >= 0 && v < vertexLevel.size());
        return vertexLevel[v] < transitNodeThreshold;
    }

    int32_t getRankOfTransitNodeIndex(const TransitNodeId &index) const {
        KASSERT(index >= 0 && index < transitNodes.size());
        return transitNodes[index];
    }

    TransitNodeId getTransitNodeIndexOfRank(const int32_t &v) const {
        KASSERT(isTransitNode(v));
        KASSERT(transitNodeIndexOfRank[v] != -1);
        return transitNodeIndexOfRank[v];
    }

    uint64_t sizeInBytes() const {
        return sizeof(TransitNodeHierarchy) +
               vertexLevel.size() * sizeof(decltype(vertexLevel)::value_type) +
               packedSideIds.size() * sizeof(decltype(packedSideIds)::value_type) +
               transitNodes.capacity() * sizeof(decltype(transitNodes)::value_type) +
               transitNodeIndexOfRank.size() * sizeof(decltype(transitNodeIndexOfRank)::value_type);
    }

    // TODO: remove debug getter for transit nodes
    const std::vector<int32_t> &getTransitNodes() const {
        return transitNodes;
    }

    const std::vector<TransitNodeId> &getTransitNodeIndexOfRankVector() const {
        return transitNodeIndexOfRank;
    }

private:
    // Finds depth, side bitvector, and truncation flag of each vertex.
    template<typename PermuteVertexIdT>
    void computeVertexLocationInSepDecomp(const SeparatorDecomposition &sd, const PermuteVertexIdT& permuteVertexId) {


        std::stack<bool> doneWithLeftChild;
        doneWithLeftChild.push(false);
        Level depth = 1;
        SideId packedSideId = 0;

        // Set location info for root node separator vertices
        for (auto vSd = sd.lastSeparatorVertex(0) - 1; vSd >= sd.firstSeparatorVertex(0); --vSd) {
            const auto v = permuteVertexId[vSd];
            vertexLevel[v] = 0;
            packedSideIds[v] = packedSideId;
            if (0 < transitNodeThreshold)
                transitNodes.push_back(v);
        }

        const auto recurse = [&](const int /*parent*/, const int child) {
            KASSERT(doneWithLeftChild.size() == depth);

            // The child will get a new side ID with a new bit stating which of the two children it is.
            // Set next bit in packedSideId to 0 for recursion to left child and to 1 for recursion to right child.
            setBit(packedSideId, depth - 1, doneWithLeftChild.top());

            // Set location info for separator vertices at child.
            for (auto vSd = sd.lastSeparatorVertex(child) - 1; vSd >= sd.firstSeparatorVertex(child); --vSd) {
            const auto v = permuteVertexId[vSd];
                vertexLevel[v] = depth;
                packedSideIds[v] = packedSideId;
                if (depth < transitNodeThreshold)
                    transitNodes.push_back(v);
            }

            // Memorize that we are not done with left child of child
            ++depth;
            doneWithLeftChild.push(false);

        };

        const auto backtrack = [&](const int /*child*/, const int /*parent*/) {
            KASSERT(doneWithLeftChild.size() == depth);

            // Remember that (at least) one child of parent is done. This way, after left child of parent is done,
            // we know that next recursion from parent is right child.
            doneWithLeftChild.pop();
            doneWithLeftChild.top() = true;
            --depth;

            // Packed ID of child has one bit more than that of parent. Reset this bit since we are going back to parent.
            setBit(packedSideId, depth - 1, false);
        };

        sepdecomp::forEachNodeInDfsOrder(sd, recurse, backtrack);

        // Re-order transit nodes by decreasing level and increasing rank, and build mapping from CCH rank to index within transit nodes
        auto compareByLevelAndRank = [&](int32_t a, int32_t b) {
            return vertexLevel[a] > vertexLevel[b] || (vertexLevel[a] == vertexLevel[b] && a < b);
        };
        if (transitNodes.size() > std::numeric_limits<TransitNodeId>::max())
            throw std::invalid_argument(
                    "Number of transit nodes (" + std::to_string(transitNodes.size()) +
                    ") exceeds TransitNodeId range; lower the transit node level threshold.");
        std::sort(transitNodes.begin(), transitNodes.end(), compareByLevelAndRank);
        for (int i = 0; i < transitNodes.size(); ++i) {
            transitNodeIndexOfRank[transitNodes[i]] = i;
        }

        std::cout << "CTNR: Selected " << transitNodes.size() << " transit nodes from top " << transitNodeThreshold
                  << " levels" << std::endl;
        std::cout << "Total number of vertices: " << packedSideIds.size() << std::endl;
    }


    std::vector<SideId> packedSideIds; // store which side each vertex is on in each level of separator hierarchy
    std::vector<Level> vertexLevel; // Map each vertex to its level in the SD

    // The small subset of vertices in the transitNodeThreshold highest levels make up the transit nodes.
    // Every transit node gets an internal index in [0, numTransitNodes-1] for distance table lookup.
    // We identify transit nodes by this index.
    // To map the CCH-rank of a transit node r  to its index i, use i = transitNodeToDistanceTableIndex[r].
    // To map a transit node index i to the CCH-rank r of the associated vertex, use r = transitNodes[i].
    Level transitNodeThreshold;
    std::vector<int32_t> transitNodes; // List of vertices in the top transitNodeThreshold levels which make up transit nodes
    std::vector<TransitNodeId > transitNodeIndexOfRank; // maps CCH rank to index in transitNodes
};

