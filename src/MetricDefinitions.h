#pragma once
#include <array>
#include <cstddef>

// Zoznam metrik definovany cez X-macro, aby to bolo/a dalo sa menit na jendom mieste.
// rovnake poradie sa pouzije vo vypocte metrik, TIFF exportoch aj v stlpcoch CSV s feature vektormi.
// Pri zmene poradia sa meni format dat; existujuce tiff rastre musia mat rovnake poradie.
// z je normalizovana vyska po filtracii. BR maju v menovateli pocet vegetacnych bodov.
// m2, m3, m4 su centralne momenty s menovatelom n.
#define VEGETATION_METRIC_LIST(X) \
    X(Hmax) /* Maximalna normalizovana vyska vegetacie [m]. */ \
    X(Hmean) /* Priemerna normalizovana vyska vegetacie [m]. */ \
    X(Hmedian) /* 50. percentil vysok vegetacie [m]. */ \
    X(Hp25) /* 25. percentil vysok vegetacie [m]. */ \
    X(Hp75) /* 75. percentil vysok vegetacie [m]. */ \
    X(Hp95) /* 95. percentil vysok vegetacie [m]. */ \
    X(PPR) /* Podiel bodov zeme zo vsetkych bodov zeme+vegetacie. */ \
    X(DAM_z) /* Podiel vegetacnych bodov s vyskou >= Hmean. */ \
    X(BR_below_1) /* Podiel vegetacnych bodov s vyskou <= 1 m. */ \
    X(BR_1_2) /* Podiel vegetacnych bodov vo vrstve 1 az 2 m. */ \
    X(BR_2_3) /* Podiel vegetacnych bodov vo vrstve 2 az 3 m. */ \
    X(BR_above_3) /* Podiel vegetacnych bodov s vyskou >= 3 m. */ \
    X(BR_3_4) /* Podiel vegetacnych bodov vo vrstve 3 az 4 m. */ \
    X(BR_4_5) /* Podiel vegetacnych bodov vo vrstve 4 az 5 m. */ \
    X(BR_below_5) /* Podiel vegetacnych bodov s vyskou <= 5 m. */ \
    X(BR_5_20) /* Podiel vegetacnych bodov vo vrstve 5 az 20 m. */ \
    X(BR_above_20) /* Podiel vegetacnych bodov s vyskou >= 20 m. */ \
    X(Coeff_var_z) /* Variacny koeficient vysok vegetacie: Hstd/Hmean. */ \
    X(Hkurt) /* Spicatost vysok vegetacie: m4 / m2^2, bez odcitania 3 a korekcie vyberu. */ \
    X(Hskew) /* Sikmost vysok vegetacie: m3 / m2^(3/2), bez korekcie vyberu. */ \
    X(Hstd) /* Vyberova smerodajna odchylka vysok vegetacie [m]: sqrt(Hvar). */ \
    X(Hvar) /* Vyberovy rozptyl vysok vegetacie [m^2], menovatel n - 1. */ \
    X(Shannon) /* Shannonov index vysok vegetacie v 0.5 m vrstvach: -sum(p * ln(p)). */ \
    X(DTM) /* Priemer referencnych vysok normalizacie z 1 m buniek vo vystupnom pixeli [m]. */ \
    X(CRR) /* Canopy Relief Ratio: (Hmean - Hmin) / (Hmax - Hmin); pri rovnosti extremov 0. */ \
    X(VCI) /* Index vertikalnej komplexnosti: entropia 10 tried / ln(10); trieda = int(9*z/Hmax). */ \
    X(Rumple) /* Pomer 3D plochy trojuholnikov z Hmax a susedov k ich horizontalnemu priemetu. */

namespace metrics
{
// Z rovnakeho zoznamu vytvor cislicovane indexy a textove nazvy.
enum Index : std::size_t
{
#define VEGETATION_METRIC_INDEX(name) name,
    VEGETATION_METRIC_LIST(VEGETATION_METRIC_INDEX)
#undef VEGETATION_METRIC_INDEX
    Count
};

inline constexpr std::array<const char*, Count> names = {
#define VEGETATION_METRIC_NAME(name) #name,
    VEGETATION_METRIC_LIST(VEGETATION_METRIC_NAME)
#undef VEGETATION_METRIC_NAME
};
}

// Pomocne makra neunikaju do ostatnych suborov.
#undef VEGETATION_METRIC_LIST
