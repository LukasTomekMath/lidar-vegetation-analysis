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

int segmentFeaturesSat(std::string shpDir, std::string tifDirSat, std::string outputSat){
    std::cout << "satellite segmenting";

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