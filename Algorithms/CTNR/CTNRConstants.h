#pragma once

#include <cstdint>
#include <limits>

// INFTY for CTNR needs to be robust for adding three distances without overflow.
static constexpr const int32_t CTNR_INFTY = std::numeric_limits<int>::max() / 3;

namespace ctnr {
    // int32, not uint16: USA holds 86178 transit nodes at thresh 14 (guard throws past the range).
    using TransitNodeId = int32_t;
    using Level = int8_t;
}

