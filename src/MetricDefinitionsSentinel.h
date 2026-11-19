#pragma once
#include <array>
#include <cmath>

#define SENTINEL_BANDS_10M(X) \
    X(B02) /* Blue       :   492.3 nm - 492.7 nm */ \
    X(B03) /* Green      :   558.9 nm - 559.8 nm*/ \
    X(B04) /* Red        :   664.6 nm - 664.9 nm */ \
    X(B08) /* NIR wide   :   832.8 nm - 832.9 nm */ 
#define SENTINEL_BANDS_20M(X) \
    X(B05) /* Red-edge 1 :   703.8 nm - 704.1 nm*/ \
    X(B06) /* Red-edge 2 :   739.1 nm - 740.5 nm*/ \
    X(B07) /* Red-edge 3 :   779.7 nm - 782.8 nm*/ \
    X(B8A) /* NIR narrow :   864.0 nm - 864.7 nm*/ \
    X(B11) /* SWIR 1     :   1610.4 nm - 1613.7 nm */ \
    X(B12) /* SWIR 2     :   2185.7 nm - 2204.4 nm*/ \
    X(SCL) /* Scene classification layer */ 

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