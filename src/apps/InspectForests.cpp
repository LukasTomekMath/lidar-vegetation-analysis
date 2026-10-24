#include "AppSupport.h"
#include "ForestManager.h"
#include "RasterIO.h"
#include "gdal_priv.h"
#include "ogrsf_frmts.h"
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <omp.h>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace fs = std::filesystem;

// Nacita lokality z KML a zapise ich atributy a hranice do pralesy.txt.
// Do tiles.txt zapise nazvy LAZ suborov pre jednotlive lokality a v konzole vypise
// pocty najdenych a chybajucich ciest v priecinku LAZ; obsah LAZ suborov nenacitava.
int runInspectForests(const AppOptions &options)
{

    std::string kmlDir = options.kmlDir;
    std::string lazDir = options.lazDir;

    if (!fs::exists(kmlDir))
    {
        std::cerr << "Directory does not exist: " << kmlDir << "\n";
        return -1;
    }

    GDALAllRegister();
    ForestManager forestManager(spatial::directoryLasCRS(lazDir));
    forestManager.loadDirectory(kmlDir);
    auto &forests = forestManager.getForests();

    // pralesy.txt obsahuje geometriu lokalit
    std::ofstream file(fs::path(options.outputDir) / "pralesy.txt");
    if (!file.is_open())
    {
        std::cerr << "Could not create pralesy.txt\n";
        return 1;
    }
    // tiles.txt obsahuje zoznam nazvov LAZ pre kazdu lokalitu
    std::ofstream tiles(fs::path(options.outputDir) / "tiles.txt");
    if (!tiles)
        throw std::runtime_error("Cannot create tiles.txt");
    int tile_counter = 0;

    for (size_t i = 0; i < forests.size(); ++i)
    {
        Forest &forest = forests[i];

        // Vypocet plochy, zaroven odstrani vnutorne hranice s plochou mensou nez 1 m2.
        forest.calculateForestArea();
        // Z opisaneho obdlznika (+padding) odvod nazvy suborov pre dlazdice 2 x 2 km.
        forest.findBoundingBox();
        forest.findTiles();

        // Zapis vsetky nazvy z findTiles(), vratane variantov pomenovania a duplicit.
        tiles << forest.getName() << ": \n";
        std::vector<std::string> &subory = forest.getTiles();
        tile_counter += subory.size();
        for (int k = 0; k < subory.size(); k++)
        {
            tiles << subory[k] << "\n";
        }
        tiles << "\n";

        auto &polygons = forest.getPolygons();

        // Atributy lokality a plocha prepoctena z m2 na hektare.
        file << "=== Forest #" << i << " ===\n";
        file << "Name: " << forest.getName() << "\n";
        file << "Segment: " << forest.getSegment() << "\n";
        file << "Kod prales: " << forest.getKodPrales() << "\n";
        file << std::fixed << std::setprecision(6);
        file << "Hectares: " << forest.getForestArea() / 10000.0 << "\n";
        file << "Polygons: " << polygons.size() << "\n\n";

        for (size_t j = 0; j < polygons.size(); ++j)
        {
            ForestPolygon &pg = polygons[j];

            auto &outerPoints = pg.outer.getPoints();

            file << "  Polygon #" << j << ":\n";
            file << "    Outer boundary points: " << outerPoints.size() << "\n";
            file << "    Outer perimeter: " << pg.outer.getPerimeter() << " m\n";

            // Pre kontrolu transformacie vypis povodne aj JTSK suradnice vrcholov.
            for (size_t k = 0; k < outerPoints.size(); ++k)
            {
                Point &p = outerPoints[k];

                file << std::fixed << std::setprecision(6);
                file << "      [" << k << "]\n";
                file << "         Original:  Lon=" << p.lon << ", Lat=" << p.lat << ", Alt=" << p.alt << "\n";

                file << std::fixed << std::setprecision(0);
                file << "         JTSK:      X=" << p.X << ", Y=" << p.Y << ", Z=" << p.Z << "\n";
            }

            // Rovnakym sposobom vypis zachovane vnutorne hranice (diery polygonu).
            if (!pg.inners.empty())
            {
                file << "    Inner boundaries: " << pg.inners.size() << "\n";

                for (size_t m = 0; m < pg.inners.size(); ++m)
                {
                    auto &innerPoints = pg.inners[m].getPoints();

                    file << "      Inner #" << m << " points = " << innerPoints.size() << "\n";

                    for (size_t k = 0; k < innerPoints.size(); ++k)
                    {
                        Point &p = innerPoints[k];

                        file << std::fixed << std::setprecision(6);
                        file << "        [" << k << "]\n";
                        file << "           Original:  Lon=" << p.lon << ", Lat=" << p.lat << ", Alt=" << p.alt << "\n";

                        file << std::fixed << std::setprecision(0);
                        file << "           JTSK:      X=" << p.X << ", Y=" << p.Y << ", Z=" << p.Z << "\n";
                    }
                }
            }

            file << "\n";
        }

        file << "\n";
    }

    file.close();
    tiles.close();
    std::cout << "Treba " << tile_counter << " suborov laz\n";

    // Over len existenciu ciest; pocty zahrnaju kazdy vyskyt nazvu, nie unikatne subory.
    int found = 0;
    int missing = 0;

    for (size_t i = 0; i < forests.size(); ++i)
    {
        Forest &forest = forests[i];

        auto &files = forest.getTiles();

        for (int j = 0; j < files.size(); j++)
        {
            fs::path fullPath = fs::path(lazDir) / files[j];
            if (fs::exists(fullPath))
            {
                found++;
            }
            else
            {
                // std::cout << forest.getName() << "\n";
                // std::cout << "[MISSING] " << files[j] << "\n";
                missing++;
            }
        }
    }
    std::cout << "\nSummary:\n";
    std::cout << "  Found:   " << found << "\n";
    std::cout << "  Missing: " << missing << "\n";

    return 0;
}
