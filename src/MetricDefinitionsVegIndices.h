#pragma once
#include "MetricDefinitionsSentinel.h"

#define VEG_INDICES(X) \
    X(NDVI, (p.B08 - p.B04) / (p.B08 + p.B04)) /* normalized difference vegetation index */ \
    X(NDMI, (p.B08 - p.B11) / (p.B08 + p.B11)) /*normalized difference moisture index */

namespace VI {
enum Index : std::size_t {
#define X(n, e) n,
    VEG_INDICES(X)
#undef X
    Count };

inline constexpr std::array<const char*, Count> names = {
#define X(n, e) #n,
    VEG_INDICES(X)
#undef X
};

inline constexpr std::array<double (*)(const Pixel&), Count> funcs = {
#define X(n, e) [](const Pixel& p) -> double { return e; },
    VEG_INDICES(X)
#undef X
};
}