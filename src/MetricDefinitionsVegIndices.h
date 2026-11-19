#pragma once
#include "MetricDefinitionsSentinel.h"
#include "RasterIO.h"

#define VEG_INDICES(X) \
    X(NDVI, (p.B08 - p.B04) / (p.B08 + p.B04)) /* Normalized Difference Vegetation Index */ \
    X(gNDVI,(p.B08 - p.B03) / (p.B08 + p.B03)) /* Green Normalized Difference Vegetation Index */ \
    X(NDMI, (p.B08 - p.B11) / (p.B08 + p.B11)) /* Normalized Difference Moisture Index */                  

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

inline double cleanValue(double v, double noData) {
    return (!std::isfinite(v) || v == noData) ? NAN : v;
}

inline Pixel makePixel(const RasterData& r10, const RasterData& r20, size_t i10, size_t i20)
{
    Pixel p;
#define X(name) p.name = cleanValue(r10.bands[sentinel::name##_10m][i10], r10.noDataValues[sentinel::name##_10m]);
    SENTINEL_BANDS_10M(X)
#undef X
#define X(name) p.name = cleanValue(r20.bands[sentinel::name##_20m][i20], r20.noDataValues[sentinel::name##_20m]);
    SENTINEL_BANDS_20M(X)
#undef X
    return p;
}