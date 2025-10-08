#pragma once

#include <kassert/kassert.hpp>
#include <stack>

#include "Algorithms/CCH/CCH.h"
#include "DataStructures/Partitioning/SeparatorTree.h"
#include "Algorithms/CCH/CCHMetric.h"

#include "Tools/Constants.h"

class BalancedTopologyCentricTreeHierarchy {

    static constexpr uint32_t MaxTruncatedSubtreeSize = CTL_THETA;

public:

    static constexpr bool NoTruncatedVertices = MaxTruncatedSubtreeSize == 0;

    BalancedTopologyCentricTreeHierarchy() = default;

    // Builds the metric-independent CCH for the specified graph and separator decomposition.
    template<typename InputGraphT>
    void preprocess(const InputGraphT &inputGraph, const SeparatorDecomposition &sepDecomp) {

        const auto sdDepth = computeSepDecompDepth(sepDecomp);
        std::cout << "Depth of sepDecomp is " << sdDepth << std::endl;

        if (!hasStrictDissectionStructure(sepDecomp))
            throw std::invalid_argument("BalancedTopologyCentricTreeHierarchy requires strict dissection "
                                        "structure of separator decomposition.");

        uint32_t maxUntruncatedDepth = 1;
        const auto sdNodesToTruncateAt = getSepDecompNodesToTruncateAt(sepDecomp, MaxTruncatedSubtreeSize,
                                                                       maxUntruncatedDepth);
        std::cout << "Max depth of untruncated node in sepDecomp is " << maxUntruncatedDepth << std::endl;

        // Build labels and numCommonHubsComputer
        packedSideIds.clear();
        packedSideIds.resize(inputGraph.numVertices(), static_cast<uint64_t>(-1));
        truncateVertex.resize(inputGraph.numVertices(), false);

        std::vector<uint32_t> depthOfVertex(inputGraph.numVertices(), 0);
        computeVertexLocationInSepDecomp(sepDecomp, sdNodesToTruncateAt, depthOfVertex);

        firstSepSizeSum.clear();
        firstSepSizeSum.resize(inputGraph.numVertices() + 1);
        sepSizeSumsOnBranch.clear();
        computeSepSizeSumsOnBranches(sepDecomp, sdNodesToTruncateAt, depthOfVertex);

        KASSERT(std::all_of(packedSideIds.begin(), packedSideIds.end(),
                            [](const uint64_t &id) { return id != static_cast<uint64_t>(-1); }));
    }

    const uint32_t &getNumHubs(const int32_t& v) const {
        return sepSizeSumsOnBranch[firstSepSizeSum[v + 1] - 1];
    }

    size_t numVertices() const {
        return packedSideIds.size();
    }

    // Expects ranks in the CCH-order as inputs.
    uint32_t getLowestCommonHub(const int32_t &s, const int32_t &t) const {

        // XOR packed side IDs to find out lowest common level in separator hierarchy.
        const int l = lowestOneBit(packedSideIds[s] ^ packedSideIds[t]);

        // If packed side IDs of s and t are exactly the same, the branch of s subsumes the branch of t or vice
        // versa. In this case, the lowest common ancestor is the lower one of the two vertices.
        if (l < 0)
            return std::min(getNumHubs(s), getNumHubs(t));

        auto l2 = std::min(static_cast<uint32_t>(l), getVertexDepth(s));
        l2 = std::min(l2, getVertexDepth(t));
        KASSERT(l2 < firstSepSizeSum[s + 1] - firstSepSizeSum[s]);
        KASSERT(l2 < firstSepSizeSum[t + 1] - firstSepSizeSum[t]);
        return std::min(sepSizeSumsOnBranch[firstSepSizeSum[s] + l2], sepSizeSumsOnBranch[firstSepSizeSum[t] + l2]);
    }

    // Return whether vertex is truncated in which case it does not have a label.
    bool isVertexTruncated(const int v) const {
        return truncateVertex[v];
    }

    uint64_t sizeInBytes() const {
        return sizeof(BalancedTopologyCentricTreeHierarchy) +
        firstSepSizeSum.size() * sizeof(decltype(firstSepSizeSum)::value_type) +
               sepSizeSumsOnBranch.size() * sizeof(decltype(sepSizeSumsOnBranch)::value_type) +
               packedSideIds.size() * sizeof(decltype(packedSideIds)::value_type) +
               truncateVertex.numBlocks() * sizeof(BitVector::Block) ;
    }

    inline uint32_t getVertexDepth(const int32_t& v) const {
        KASSERT(v >= 0 && v < firstSepSizeSum.size() - 1);
        return firstSepSizeSum[v + 1] - firstSepSizeSum[v] - 1;
    }

private:
    // Returns true if every separator node in decomposition has at most two children, or false otherwise.
    static bool hasStrictDissectionStructure(const SeparatorDecomposition &sd) {
        for (const auto &n: sd.tree)
            if (n.rightSibling != 0 && sd.tree[n.rightSibling].rightSibling != 0)
                return false;
        return true;
    }


    template<typename RecurseCallbacKT,
            typename BacktrackCallbackT>
    static void forEachSepDecompNodeInDfsOrder(const SeparatorDecomposition &sd,
                                               RecurseCallbacKT recurse,
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

    // Traverses separator decomposition and computes nodes to truncate at s.t. the number of vertices in the subtree
    // rooted at each marked node is at most theta. Also sets maximum depth of an untruncated node.
    static BitVector
    getSepDecompNodesToTruncateAt(const SeparatorDecomposition &sd, const uint32_t theta,
                                  uint32_t &maxUntruncatedDepth) {
        BitVector doTruncate(sd.tree.size());
        std::stack<uint32_t> numVerticesOnBranch;
        numVerticesOnBranch.push(0);
        uint32_t curDepth = 1;
        maxUntruncatedDepth = 1;

        forEachSepDecompNodeInDfsOrder(
                sd,
                [&](const int, const int) {
                    numVerticesOnBranch.push(0);
                    ++curDepth;
                },
                [&](const int child, const int) {
                    // Backtrack from child to parent when child is done.
                    // Add number of vertices in separator at child to number of vertices
                    auto numVerticesChild = numVerticesOnBranch.top();
                    numVerticesOnBranch.pop();
                    numVerticesChild += sd.lastSeparatorVertex(child) - sd.firstSeparatorVertex(child);

                    // If number of vertices at child does not exceed theta, mark it to be treated as truncated.
                    if (numVerticesChild <= theta) {
                        for (int grandChild = sd.leftChild(child);
                             grandChild != 0; grandChild = sd.rightSibling(grandChild))
                            KASSERT(doTruncate[grandChild]);
                        doTruncate[child] = true;
                    } else {
                        maxUntruncatedDepth = std::max(maxUntruncatedDepth, curDepth);
                    }

                    --curDepth;

                    // Add number of vertices at child to number of vertices at parent
                    numVerticesOnBranch.top() += numVerticesChild;
                });

        // Potentially set doTruncate for root node
        KASSERT(numVerticesOnBranch.size() == 1);
        if (numVerticesOnBranch.top() <= theta) {
            KASSERT(maxUntruncatedDepth == 1);
            doTruncate[0] = true;
        }

        return doTruncate;
    }


    // Finds depth, side bitvector, and truncation flag of each vertex.
    void computeVertexLocationInSepDecomp(const SeparatorDecomposition &sd, const BitVector &truncateAtSdNode,
                                          std::vector<uint32_t>& depthOfVertex) {


        std::stack<bool> doneWithLeftChild;
        doneWithLeftChild.push(false);
        uint32_t depth = 1;
        uint64_t packedSideId = 0;

        // Set location info for root node separator vertices
        for (auto v = sd.lastSeparatorVertex(0) - 1; v >= sd.firstSeparatorVertex(0); --v) {
            depthOfVertex[v] = 0;
            packedSideIds[v] = packedSideId;
            truncateVertex[v] = truncateAtSdNode[0];
        }

        // If root is truncated, we are done.
        if (truncateAtSdNode[0])
            return;


        uint32_t lastNonTruncatedDepth = 1;

        const auto recurse = [&](const int parent, const int child) {
            KASSERT(!truncateAtSdNode[parent] || truncateAtSdNode[child]);
            KASSERT(doneWithLeftChild.size() == depth);

            // If the parent is not truncated, the child will get a new side ID with a new bit stating which
            // of the two children it is. If the parent is truncated, the child is also truncated and gets the same
            // side ID as the parent.
            if (!truncateAtSdNode[parent]) {
                KASSERT(depth == lastNonTruncatedDepth);
                // Set next bit in packedSideId to 0 for recursion to left child and to 1 for recursion to right child.
                setBit(packedSideId, depth - 1, doneWithLeftChild.top());
            }


            if (truncateAtSdNode[child]) {
                // If child is truncated, the vertices in the subtree rooted at child will have no labels.
                // Set vertex location info up to lowest non-truncated ancestor. Mark vertices as truncated.
                for (auto v = sd.lastSeparatorVertex(child) - 1; v >= sd.firstSeparatorVertex(child); --v) {
                    depthOfVertex[v] = lastNonTruncatedDepth;
                    packedSideIds[v] = packedSideId;
                    truncateVertex[v] = true;
                }
            } else {
                // If child is not truncated, set location info for separator vertices at child.
                for (auto v = sd.lastSeparatorVertex(child) - 1; v >= sd.firstSeparatorVertex(child); --v) {
                    depthOfVertex[v] = lastNonTruncatedDepth;
                    packedSideIds[v] = packedSideId;
                }

                ++lastNonTruncatedDepth;
            }

            // Memorize that we are not done with left child of child
            ++depth;
            doneWithLeftChild.push(false);

        };

        const auto backtrack = [&](const int child, const int parent) {
            KASSERT(doneWithLeftChild.size() == depth);

            // Remember that (at least) one child of parent is done. This way, after left child of parent is done,
            // we know that next recursion from parent is right child.
            doneWithLeftChild.pop();
            doneWithLeftChild.top() = true;
            --depth;

            // If parent was not truncated, packed ID of child has one bit more than that of parent. Reset this bit
            // since we are going back to parent.
            if (!truncateAtSdNode[parent]) {
                setBit(packedSideId, depth - 1, false);
            }

            // If child was not truncated, also reduce lastNonTruncatedDepth.
            if (!truncateAtSdNode[child])
                --lastNonTruncatedDepth;
        };

        forEachSepDecompNodeInDfsOrder(sd, recurse, backtrack);
    }

    // Populates firstSepSizeSum and sepSizeSumsOnBranches by traversing the separator decomposition.
    // Requires known number of hubs per vertex for offsets in 2D-array.
    void computeSepSizeSumsOnBranches(const SeparatorDecomposition &sd, const BitVector &truncateAtSdNode, const std::vector<uint32_t>& depthOfVertex) {

        // Set firstSepSizeSum to prefix sum over numHubs
        uint32_t offset = 0;
        for (uint32_t i = 0; i < depthOfVertex.size(); ++i) {
            firstSepSizeSum[i] = offset;
            offset += depthOfVertex[i] + !static_cast<bool>(truncateVertex[i]);
        }
        firstSepSizeSum[depthOfVertex.size()] = offset;

        // Initialize sepSizeSumsOnBranches
        sepSizeSumsOnBranch.resize(offset, std::numeric_limits<uint32_t>::max());

        // If root is truncated, assign 0 sep size sum to all vertices and return.
        if (truncateAtSdNode[0]) {
            for (int32_t v = 0; v < firstSepSizeSum.size() - 1; ++v) {
                KASSERT(firstSepSizeSum[v + 1] - firstSepSizeSum[v] == 1);
                sepSizeSumsOnBranch[firstSepSizeSum[v]] = 0;
            }
            return;
        }

        // Populate root node separator vertices
        uint32_t inSepIdx = 0;
        for (auto v = sd.lastSeparatorVertex(0) - 1; v >= sd.firstSeparatorVertex(0); --v) {
            sepSizeSumsOnBranch[firstSepSizeSum[v]] = inSepIdx + 1;
            ++inSepIdx;
        }

        std::vector<uint32_t> sepSizeSumsOnCurBranch;
        sepSizeSumsOnCurBranch.push_back(sd.lastSeparatorVertex(0) - sd.firstSeparatorVertex(0));


        const auto recurse = [&](const int parent, const int child) {
            KASSERT(!truncateAtSdNode[parent] || truncateAtSdNode[child]);


            if (truncateAtSdNode[child]) {
                // If child is truncated, the vertices in the subtree rooted at child will have no labels.
                // Set sepSizeSumsOnBranch up to last non-truncated ancestor.
                for (auto v = sd.lastSeparatorVertex(child) - 1; v >= sd.firstSeparatorVertex(child); --v) {
                    KASSERT(firstSepSizeSum[v + 1] - firstSepSizeSum[v] == sepSizeSumsOnCurBranch.size());
                    std::copy(sepSizeSumsOnCurBranch.begin(), sepSizeSumsOnCurBranch.end(), sepSizeSumsOnBranch.begin() + firstSepSizeSum[v]);
                }
            } else {


                // If child is not truncated, set number of hubs and side IDs for separator vertices at child.
                uint32_t inSepIdx = 0;
                for (auto v = sd.lastSeparatorVertex(child) - 1; v >= sd.firstSeparatorVertex(child); --v) {
                    KASSERT(firstSepSizeSum[v + 1] - firstSepSizeSum[v] == sepSizeSumsOnCurBranch.size() + 1);
                    // Populate up to parent
                    std::copy(sepSizeSumsOnCurBranch.begin(), sepSizeSumsOnCurBranch.end(), sepSizeSumsOnBranch.begin() + firstSepSizeSum[v]);
                    // Populate location within current separator
                    sepSizeSumsOnBranch[firstSepSizeSum[v + 1] - 1] = sepSizeSumsOnCurBranch.back() + inSepIdx + 1;
                    ++inSepIdx;
                }

                // Memorize separator size up to and including child
                sepSizeSumsOnCurBranch.push_back(sepSizeSumsOnCurBranch.back() + inSepIdx);
            }
        };

        const auto backtrack = [&](const int child, const int) {
            // If child was not truncated, it has an entry on sepSizeSumsOnCurBranch. Remove this entry.
            if (!truncateAtSdNode[child]) {
                sepSizeSumsOnCurBranch.pop_back();
            }
        };

        forEachSepDecompNodeInDfsOrder(sd, recurse, backtrack);

        for (const auto& s : sepSizeSumsOnBranch)
            KASSERT(s != std::numeric_limits<uint32_t>::max());
    }




    std::vector<uint64_t> packedSideIds; // store which side each vertex is on in each level of separator hierarchy
    BitVector truncateVertex; // If truncateVertex[v], vertex v is truncated and should not get a label

    // Let v be a vertex which is a separator vertex at level l. We store the sum of separator sizes, i.e., number of
    // hubs, at every level between 0 and l on the branch of the decomposition that leads to v. This allows us to
    // retrieve the number of common hubs between two vertices by finding their lowest common level in the
    // decomposition (using packedSideIds) and reading the number of hubs up to that level on the common branch.
    // The number at level l only includes all vertices in the separator that have higher (or equal) rank in the order.
    // Thus, level stores the total number of hubs for v.
    // The level of truncated vertices is considered to be the level of the lowest non-truncated ancestor.
    // We use an index array firstSepSizeSum s.t. the firstSepSizeSum[v + 1] - firstSepSizeSum[v] = l + 1 and the
    // sep size sums for vertex v are stored at sepSizeSumsOnBranch[firstSepSizeSum[v]..firstSepSizeSum[v+1]].
    std::vector<uint32_t> firstSepSizeSum;
    std::vector<uint32_t> sepSizeSumsOnBranch;
};

