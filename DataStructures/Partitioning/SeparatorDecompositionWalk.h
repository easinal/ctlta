#pragma once

#include <algorithm>
#include <cstdint>
#include <stack>

#include "DataStructures/Partitioning/SeparatorDecomposition.h"

// Walking a (strict) separator decomposition. Every hub/transit hierarchy built on top of a
// decomposition needs the same three things -- the structure check, one DFS order, and the
// resulting depth -- so they live here instead of once per hierarchy.
namespace sepdecomp {

// True when every node has at most two children, which is what a strict nested dissection
// produces and what the packed side ids (one bit per level) assume.
inline bool hasStrictDissectionStructure(const SeparatorDecomposition &sd) {
    for (const auto &n: sd.tree)
        if (n.rightSibling != 0 && sd.tree[n.rightSibling].rightSibling != 0)
            return false;
    return true;
}

// Depth-first walk. recurse(parent, child) fires when descending into child; backtrack(child,
// parent) when leaving child. Iterative on purpose: continental decompositions are ~35 deep
// but the recursion would carry the whole separator state.
template<typename RecurseCallbackT, typename BacktrackCallbackT>
void forEachNodeInDfsOrder(const SeparatorDecomposition &sd, RecurseCallbackT recurse,
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

        // Done with this node: continue with a sibling if there is one, else return to parent.
        sdNodesStack.pop();
        if (sdNodesStack.empty())
            break;

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

// Number of levels, root separator included.
inline size_t depth(const SeparatorDecomposition &sd) {
    size_t maxDepth = 0, curDepth = 1;
    forEachNodeInDfsOrder(sd,
                          [&](const int, const int) {
                              ++curDepth;
                              maxDepth = std::max(maxDepth, curDepth);
                          },
                          [&](const int, const int) { --curDepth; });
    return maxDepth;
}

} // namespace sepdecomp
