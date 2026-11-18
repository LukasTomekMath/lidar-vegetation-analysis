#pragma once
#include <array>
#include <cmath>

#define SENTINEL_BANDS_10M(X) X(B02) X(B03) X(B04) X(B08)
#define SENTINEL_BANDS_20M(X) X(B05) X(B06) X(B07) X(B8A) X(B11) X(B12) X(SCL)

namespace sentinel
{
enum Index10 : std::size_t
{
#define X(name) name##_10m,
    SENTINEL_BANDS_10M(X)
#undef X
    Count10
};

enum Index20 : std::size_t
{
#define X(name) name##_20m,
    SENTINEL_BANDS_20M(X)
#undef X
    Count20 // 7, SCL is the last band (index SCL_20m)
};

inline constexpr std::array<const char *, Count10> names10 = {
#define X(name) #name,
    SENTINEL_BANDS_10M(X)
#undef X
};

inline constexpr std::array<const char *, Count20> names20 = {
#define X(name) #name,
    SENTINEL_BANDS_20M(X)
#undef X
};
}

struct Pixel{
#define X(n) double n = NAN;
    SENTINEL_BANDS_10M(X)
    SENTINEL_BANDS_20M(X)
#undef X
};