#pragma once

#include <cstdint>
#include <limits>

// Room to add three distances without overflow.
static constexpr const int32_t FHL_INFTY = std::numeric_limits<int>::max() / 3;

namespace fhl {
    // Dense id of a HUB -- a separator vertex above the region cut. Not a vertex id and not
    // a CCH rank: it indexes the per-hub rows and column runs, so it must stay dense in
    // [0, numHubs). int32: continental graphs hold hundreds of thousands of hubs, and -1
    // is the "not a hub" sentinel (see HubHierarchy::hubIdOrInvalid).
    using HubId = int32_t;
    // Depth in the separator decomposition: 0 is the root separator, larger is deeper.
    using Level = int8_t;
}
