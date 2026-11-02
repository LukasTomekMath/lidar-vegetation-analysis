#include "AppSupport.h"
#include "FeatureVector.h"
#include "MetricDefinitions.h"
#include "ForestManager.h"
#include "RasterIO.h"
#include "gdal_priv.h"
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

struct SatTilePair
{
    std::string path10;
    std::string path20;
};

static std::map<std::string, SatTilePair> collectSatPairs(const std::string &dir)
{
    std::map<std::string, SatTilePair> pairs;
    for (const auto &entry : fs::directory_iterator(dir))
    {
        if (!entry.is_regular_file())
            continue;

        auto ext = entry.path().extension().string();
        std::transform(ext.begin(), ext.end(), ext.begin(),
                       [](unsigned char c) { return std::tolower(c); });
        if (ext != ".tif" && ext != ".tiff")
            continue;

        const std::string stem = entry.path().stem().string();
        auto endsWith = [&](const std::string &suffix) {
            return stem.size() > suffix.size() &&
                   stem.compare(stem.size() - suffix.size(), suffix.size(), suffix) == 0;
        };

        if (endsWith("_10m"))
            pairs[stem.substr(0, stem.size() - 4)].path10 = entry.path().string();
        else if (endsWith("_20m"))
            pairs[stem.substr(0, stem.size() - 4)].path20 = entry.path().string();
    }
    return pairs;
}

int segmentFeaturesLaz(std::string kmlDir, std::string tifDirLaz, std::string output){
        std::cout << "Use case: Load tif files in a specified directory coresponding to kml files from a specified "
                 "directory and calculate feature vectors."
              << std::endl;

    if (!fs::exists(kmlDir))
    {
        std::cerr << "Directory does not exist: " << kmlDir << "\n";
        return -1;
    }

    // tu zacinaju features
    std::vector<std::string> files;

    // vytvorenie vektora so subormi
    for (const auto &entry : fs::directory_iterator(tifDirLaz))
    {
        if (!entry.is_regular_file())
            continue;

        auto path = entry.path();
        if (path.extension() == ".tif")
        {
            files.push_back(path.string());
        }
    }

    std::string csvFile = output;
    initializeFeatureVectorCSV(csvFile);

    // cyklus cez tiffy jednotlivych lokalit
    for (const auto &tifPath : files)
    {
        RasterData raster = readGeoTiff(tifPath);
        if (raster.nBands != static_cast<int>(metrics::Count))
            throw std::runtime_error("Expected " + std::to_string(metrics::Count) + " metric bands: " + tifPath);

        // Kazdy raster moze mat iny CRS; polygony aj masky vytvor v jeho sustave.
        OGRSpatialReference rasterCRS;
        if (raster.crsWkt.empty() || rasterCRS.importFromWkt(raster.crsWkt.c_str()) != OGRERR_NONE)
            throw std::runtime_error("Segment aggregation needs a known raster CRS: " + tifPath);
        ForestManager forestManager(rasterCRS);
        forestManager.loadDirectory(kmlDir);
        auto& forests = forestManager.getForests();
        for (auto& forest : forests)
        {
            forest.calculateForestArea();
            forest.findBoundingBox();
            forest.createMask();
        }

        // rozsah nacitaneho tiffu
        double tifXMin = raster.xMin;
        double tifXMax = raster.xMin + raster.nx * raster.pixelSize;
        double tifYMax = raster.yMax;
        double tifYMin = raster.yMax - raster.ny * raster.pixelSize;

        // auto stem = std::filesystem::path(tifPath).stem().string();

        for (auto &forest : forests)
        {
            bool matched = false; //(forest.getName() == stem);

            // najprv sa zistuje, ci bounding box forestu lezi v nacitanom tiffe

            // TODO: tu sa teoreticky moznomoze stat, ze ak pracujeme so segmentami, tak sa niektory
            // segment nachadza vo viacerych tiffoch (kedze tie su urobene cez bounding boxy lokalit)
            // potom sa preno vypocita feature vektor viackrat a v tabulke bude viackrat
            // to sa mi nepaci, vyzera, ze Maria mala aj snahu s tym nieco urobit skontrolovat to cez
            // premennu "stem", ale je to tu len zakomentovane a nedokoncene
            if (!matched)
            {
                matched = (forest.getRawMinX() >= tifXMin && forest.getRawMaxX() <= tifXMax &&
                           forest.getRawMinY() >= tifYMin && forest.getRawMaxY() <= tifYMax);
            }

            // Ak lezi vnutri tiffu, vypocita sa feature vektor
            // Z pixelov v maske sa vypocitaju statistiky jednotlivych pasiem a prida sa riadok do CSV.
            // NoData a nekonecne ci NaN hodnoty sa vynechavaju; bez platnych hodnot zostanu nuly.
            if (matched)
            {
                FeatureVector fv = computeFeatureVector(raster, forest);
                exportFeatureVectorCSV(fv, csvFile);
            }
        }
    }

    return 0;
}

int segmentFeaturesSat(const std::string shpDir, const std::string tifDirSat, const std::string outputSat){
    std::cout << "Use case: Load Sentinel-2 tif files from a specified directory, match them to shapefile "
                 "polygons from a specified directory and calculate feature vectors."
              << std::endl;

    if (!fs::exists(shpDir))
    {
        std::cerr << "Directory does not exist: " << shpDir << "\n";
        return -1;
    }
    if (!fs::exists(tifDirSat))
    {
        std::cerr << "Directory does not exist: " << tifDirSat << "\n";
        return -1;
    }

    // collect tif files
    std::vector<std::string> files;
    for (const auto &entry : fs::directory_iterator(tifDirSat))
    {
        if (!entry.is_regular_file())
            continue;

        auto ext = entry.path().extension().string();
        std::transform(ext.begin(), ext.end(), ext.begin(),
                       [](unsigned char c) { return std::tolower(c); });
        if (ext == ".tif" || ext == ".tiff")
            files.push_back(entry.path().string());
    }

    const std::string csvFile = outputSat;
    initializeFeatureVectorCSV(csvFile);

    for (const auto &[key, pair] : collectSatPairs(tifDirSat))
    {
        if (pair.path10.empty() || pair.path20.empty())
        {
            std::cerr << "Skipping " << key << ": missing 10m or 20m file\n";
            continue;
        }

        RasterData r10 = readGeoTiff(pair.path10);
        RasterData r20 = readGeoTiff(pair.path20);

        if (r10.nBands != static_cast<int>(sentinel::Count10))
            throw std::runtime_error("Expected " + std::to_string(sentinel::Count10) + " bands: " + pair.path10);
        if (r20.nBands != static_cast<int>(sentinel::Count20))
            throw std::runtime_error("Expected " + std::to_string(sentinel::Count20) + " bands: " + pair.path20);

        if (r10.crsWkt != r20.crsWkt)
            throw std::runtime_error("10m and 20m CRS differ for " + key);

        OGRSpatialReference rasterCRS;
        if (r10.crsWkt.empty() || rasterCRS.importFromWkt(r10.crsWkt.c_str()) != OGRERR_NONE)
            throw std::runtime_error("Unknown raster CRS: " + pair.path10);

        // extent from the 10m raster; the 20m one should cover the same footprint
        const double xMin = r10.xMin;
        const double xMax = r10.xMin + r10.nx * r10.pixelSize;
        const double yMax = r10.yMax;
        const double yMin = r10.yMax - r10.ny * r10.pixelSize;

        ForestManager forestManager(rasterCRS);
        forestManager.loadDirectory(shpDir);

        for (auto &forest : forestManager.getForests())
        {
            forest.calculateForestArea();
            forest.findBoundingBox();

            if (!(forest.getRawMinX() >= xMin && forest.getRawMaxX() <= xMax &&
                forest.getRawMinY() >= yMin && forest.getRawMaxY() <= yMax))
                continue;

            forest.createMask();
            //FeatureVector fv = computeFeatureVector(r10, r20, forest);
            //exportFeatureVectorCSV(fv, csvFile);
            std::cout << "Matched: " << forest.getName() << " in tile " << key << "\n";
        }
    }

    return 0;
}

int runSegmentFeatures(const AppOptions &options)
{
    // This aggregation loop is sequential.
    GDALAllRegister();

    std::string kmlDir = options.kmlDir;
    std::string tifDirLaz = options.tifDirLaz;
    std::string outputLaz = options.outputLaz;

    std::string shpDir = options.shpDir;
    std::string tifDirSat = options.tifDirSat;
    std::string outputSat = options.outputSat;

    if(!tifDirLaz.empty())
        segmentFeaturesLaz(kmlDir, tifDirLaz, outputLaz);
    if(!tifDirSat.empty())
        segmentFeaturesSat(shpDir, tifDirSat, outputSat);

    return 0;
}