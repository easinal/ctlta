#pragma once

#include <kassert/kassert.hpp>
#include <stack>

#include "Algorithms/CCH/CCH.h"
#include "DataStructures/Partitioning/SeparatorTree.h"
#include "Algorithms/CCH/CCHMetric.h"

#include "Tools/Constants.h"

class TransitNodeHierarchy {

    using PackedSideId = uint64_t;

public:

    using SdLevel = uint8_t; // type for level in separator decomposition

    TransitNodeHierarchy() = default;

    // Builds the metric-independent transit node hierarchy for the specified graph and separator decomposition.
    template<typename InputGraphT>
    void preprocess(const InputGraphT &inputGraph, const SeparatorDecomposition &sepDecomp,
                    const int newTransitNodeThreshold) {

        const auto sdDepth = computeSepDecompDepth(sepDecomp);
        std::cout << "Depth of sepDecomp is " << sdDepth << std::endl;

        if (!hasStrictDissectionStructure(sepDecomp))
            throw std::invalid_argument("TransitNodeHierarchy requires strict dissection "
                                        "structure of separator decomposition.");

        // Build labels and numCommonHubsComputer
        packedSideIds.clear();
        packedSideIds.resize(inputGraph.numVertices(), static_cast<PackedSideId>(-1));
        sdLevel.clear();
        sdLevel.resize(inputGraph.numVertices(), static_cast<SdLevel>(-1));

        transitNodeThreshold = newTransitNodeThreshold;
        transitNodes.clear();
        transitNodeToDistanceTableIndex.clear();

        computeVertexLocationInSepDecomp(sepDecomp);

        KASSERT(std::all_of(packedSideIds.begin(), packedSideIds.end(),
                            [](const PackedSideId &id) { return id != static_cast<PackedSideId>(-1); }));

        // Todo: remove debug
        // Count number of vertices in each level and print
        std::vector<size_t> levelCounts(sdDepth, 0);
        for (const auto &level: sdLevel) {
            KASSERT(level < sdDepth);
            ++levelCounts[level];
        }
        for (size_t level = 0; level < levelCounts.size(); ++level) {
            std::cout << "SepDecomp Level " << level << ": " << levelCounts[level] << " vertices" << std::endl;
        }
    }

    size_t numVertices() const {
        return packedSideIds.size();
    }

    int getTransitNodeThreshold() const {
        return transitNodeThreshold;
    }

    // Given the ranks of two vertices in the CCH-order, returns the level of their lowest common ancestor
    // in the separator decompisition.
    int32_t getSdLevelOfLowestCommonAncestor(const int32_t &s, const int32_t &t) const {

        const int minInputSdLevel = static_cast<int>(std::min(sdLevel[s], sdLevel[t]));

        // XOR packed side IDs to find out lowest common level in separator hierarchy.
        const int l = lowestOneBit(packedSideIds[s] ^ packedSideIds[t]);

        if (l >= 0)
            return std::min(l, minInputSdLevel);

        // If packed side IDs of s and t are exactly the same, the branch of s subsumes the branch of t or vice
        // versa. In this case, the lowest common ancestor is the lower one of the two vertices.
        return minInputSdLevel;
    }

    int32_t numTransitNodes() const {
        return static_cast<int32_t>(transitNodes.size());
    }

    bool isTransitNode(const int32_t &v) const {
        KASSERT(v >= 0 && v < sdLevel.size());
        return sdLevel[v] < transitNodeThreshold;
    }

    int32_t getRankOfTransitNodeIndex(const int32_t &index) const {
        KASSERT(index >= 0 && index < transitNodes.size());
        return transitNodes[index];
    }

    int32_t getTransitNodeIndexOfRank(const int32_t &v) const {
        KASSERT(isTransitNode(v));
        return transitNodeToDistanceTableIndex.at(v);
    }

    uint64_t sizeInBytes() const {
        return sizeof(TransitNodeHierarchy) +
               sdLevel.size() * sizeof(decltype(sdLevel)::value_type) +
               packedSideIds.size() * sizeof(decltype(packedSideIds)::value_type) +
               transitNodes.capacity() * sizeof(decltype(transitNodes)::value_type) +
               transitNodeToDistanceTableIndex.size() * (sizeof(int32_t) + sizeof(int32_t));
    }

    // TODO: remove debug getter for transit nodes
    const std::vector<int32_t> &getTransitNodes() const {
        return transitNodes;
    }

private:
    // Returns true if every separator node in decomposition has at most two children, or false otherwise.
    static bool hasStrictDissectionStructure(const SeparatorDecomposition &sd) {
        for (const auto &n: sd.tree)
            if (n.rightSibling != 0 && sd.tree[n.rightSibling].rightSibling != 0)
                return false;
        return true;
    }


    template<typename RecurseCallbackT,
            typename BacktrackCallbackT>
    static void forEachSepDecompNodeInDfsOrder(const SeparatorDecomposition &sd,
                                               RecurseCallbackT recurse,
                                               BacktrackCallbackT backtrack) {
        std::stack<uint32_t> sdNodesStack;
        sdNodesStack.push(0);
        bool returnedFromChildren = false;
        while (true) {
            const auto node = sdNodesStack.top();

            if (!returnedFromChildren && sd.leftChild(node) != 0) {
                recurse(node, sd.leftChild(node));
                sdNodesStack.push(sd.leftChild(node));
                continue;
            }

            // Done with this node. If there are siblings continue with siblings, otherwise return to parent.
            sdNodesStack.pop();
            if (sdNodesStack.empty())
                break; // Finished when stack becomes empty

            backtrack(node, sdNodesStack.top());
            if (sd.rightSibling(node) != 0) {
                recurse(sdNodesStack.top(), sd.rightSibling(node));
                sdNodesStack.push(sd.rightSibling(node));
                returnedFromChildren = false;
            } else {
                returnedFromChildren = true;
            }
        }
    }

    static size_t computeSepDecompDepth(const SeparatorDecomposition &sd) {
        uint32_t maxDepth = 0;
        uint32_t curDepth = 1;
        forEachSepDecompNodeInDfsOrder(sd,
                                       [&](const int, const int) {
                                           ++curDepth;
                                           maxDepth = std::max(maxDepth, curDepth);
                                       },
                                       [&](const int, const int) {
                                           --curDepth;
                                       });
        return maxDepth;
    }

    // Finds depth and side bitvector of each vertex.
    void computeVertexLocationInSepDecomp(const SeparatorDecomposition &sd) {


        std::stack<bool> doneWithLeftChild;
        doneWithLeftChild.push(false);
        SdLevel depth = 1;
        PackedSideId packedSideId = 0;

        int numSubtreesSeparatedByTransitNodes = 0;

        // Set location info for root node separator vertices
        for (auto v = sd.lastSeparatorVertex(0) - 1; v >= sd.firstSeparatorVertex(0); --v) {
            sdLevel[v] = 0;
            packedSideIds[v] = packedSideId;
            if (0 < transitNodeThreshold)
                transitNodes.push_back(v);
        }

        const auto recurse = [&](const int /*parent*/, const int child) {
            KASSERT(doneWithLeftChild.size() == depth);

            if (depth == transitNodeThreshold) {
                ++numSubtreesSeparatedByTransitNodes;
            }

            // The child will get a new side ID with a new bit stating which of the two children it is.
            // Set next bit in packedSideId to 0 for recursion to left child and to 1 for recursion to right child.
            setBit(packedSideId, depth - 1, doneWithLeftChild.top());

            // Set location info for separator vertices at child.
            for (auto v = sd.lastSeparatorVertex(child) - 1; v >= sd.firstSeparatorVertex(child); --v) {
                sdLevel[v] = depth;
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

        forEachSepDecompNodeInDfsOrder(sd, recurse, backtrack);

        // Re-order transit nodes by decreasing SD-level and increasing rank, and build mapping from CCH rank to index within transit nodes
        auto compareBySdLevelAndRank = [&](int32_t a, int32_t b) {
            return sdLevel[a] > sdLevel[b] || (sdLevel[a] == sdLevel[b] && a < b);
        };
        std::sort(transitNodes.begin(), transitNodes.end(), compareBySdLevelAndRank);
        for (int i = 0; i < transitNodes.size(); ++i) {
            transitNodeToDistanceTableIndex[transitNodes[i]] = i;
        }

        std::cout << "CTNR: Selected " << transitNodes.size() << " transit nodes from top " << static_cast<uint32_t>(transitNodeThreshold)
                  << " SD levels" << std::endl;
        std::cout << "CTNR: Number of subtrees separated by transit nodes: "
                  << numSubtreesSeparatedByTransitNodes << std::endl;
        std::cout << "Total number of vertices: " << packedSideIds.size() << std::endl;
    }


    std::vector<PackedSideId> packedSideIds; // Store which side each vertex is on in each level of separator hierarchy
    std::vector<SdLevel> sdLevel; // Map rank in separator decomposition of each vertex to its level in the SD

    // The small subset of vertices in the transitNodeThreshold highest SD-levels make up the transit nodes.
    // Every transit node gets an internal index in [0, numTransitNodes-1] for distance table lookup.
    // We identify transit nodes by this index.
    // To map the CCH-rank of a transit node r  to its index i, use i = transitNodeToDistanceTableIndex[r].
    // To map a transit node index i to the CCH-rank r of the associated vertex, use r = transitNodes[i].
    SdLevel transitNodeThreshold = 0;
    std::vector<int32_t> transitNodes; // List of vertices in the top transitNodeThreshold SD-levels which make up transit nodes
    std::unordered_map<int32_t, int32_t> transitNodeToDistanceTableIndex; // maps CCH rank to index in transitNodes
};

