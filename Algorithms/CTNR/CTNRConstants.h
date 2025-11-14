#pragma once

#include <cstdint>
#include <limits>

// INFTY for CTNR needs to be robust for adding three distances without overflow.
static constexpr const int32_t CTNR_INFTY = std::numeric_limits<int>::max() / 3;