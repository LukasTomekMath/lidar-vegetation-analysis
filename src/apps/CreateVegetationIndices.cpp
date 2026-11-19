#include "AppSupport.h"
#include "RasterIO.h"
#include "gdal_priv.h"
#include "MetricDefinitionsVegIndices.h"
#include <iostream>
#include <string>
#include <filesystem>

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

static bool writeIndexTiff(const std::string& path,
                    const std::vector<std::vector<float>>& idx,   // [vi::Count][nx*ny]
                    int nx, int ny,
                    const double geoTransform[6],                 // copy from source raster
                    const char* projWkt)
{
    GDALAllRegister();
    GDALDriver* drv = GetGDALDriverManager()->GetDriverByName("GTiff");
    char** opts = nullptr;
    opts = CSLSetNameValue(opts, "COMPRESS", "DEFLATE");
    opts = CSLSetNameValue(opts, "TILED", "YES");
    opts = CSLSetNameValue(opts, "BIGTIFF", "IF_SAFER");

    GDALDataset* ds = drv->Create(path.c_str(), nx, ny, (int)VI::Count, GDT_Float32, opts);
    CSLDestroy(opts);
    if(!ds) return false;

    ds->SetGeoTransform(const_cast<double*>(geoTransform));
    ds->SetProjection(projWkt);

    for(size_t i = 0; i < VI::Count; i++){
        GDALRasterBand* b = ds->GetRasterBand((int)i + 1);
        b->SetDescription(VI::names[i]);
        b->SetNoDataValue(-9999.0);
        if(b->RasterIO(GF_Write, 0, 0, nx, ny,
                        const_cast<float*>(idx[i].data()), nx, ny,
                        GDT_Float32, 0, 0) != CE_None){
            GDALClose(ds);
            return false;
        }
    }

    GDALClose(ds);
    return true;
}

static void printCurrentIndices(){
    std::cout << "Current Vegetation Indices:\n";
    for(size_t i = 0; i < VI::Count; ++i){
        std::cout << VI::names[i] << "\t";

        if(i % 4 == 0)
            std::cout << std::endl;
    }

    std::cout << std::endl << std::endl;
}

int createVegetationIndices(const std::string tifDir, const std::string outputDir){
    std::cout << "Use case: Load Sentinel-2 tif files from a specified directory, create "
                 "different vegetation indices and save them as tif files.\n";

    if(!fs::exists(tifDir)){
        std::cerr << "Directory does not exist: " << tifDir << std::endl;
        return -1;
    }

    printCurrentIndices();

    std::cout << "Reading tifs.\n";
    for(const auto&[key, pair] : collectSatPairs(tifDir)){
        if (pair.path10.empty() || pair.path20.empty()){
            std::cerr << "Skipping " << key << ": missing 10m or 20m file\n";
            continue;
        }

        RasterData r10 = readGeoTiff(pair.path10);
        RasterData r20 = readGeoTiff(pair.path20);

        if (r10.nBands != static_cast<int>(sentinel::Count10))
            throw std::runtime_error("Expected " + std::to_string(sentinel::Count10) + " bands: " + pair.path10);
        if (r20.nBands != static_cast<int>(sentinel::Count20))
            throw std::runtime_error("Expected " + std::to_string(sentinel::Count20) + " bands: " + pair.path20);
        if (r10.nx != 2 * r20.nx || r10.ny != 2 * r20.ny)
            throw std::runtime_error("10m/20m grids not 2:1 for " + key);

        const size_t n = (size_t)r10.nx * r10.ny;
        std::vector<std::vector<float>> out(VI::Count, std::vector<float>(n, -9999.f));

        for (int y = 0; y < r10.ny; y++){
            for (int x = 0; x < r10.nx; x++){
                    const size_t i10 = (size_t)y * r10.nx + x;
                    const size_t i20 = (size_t)(y / 2) * r20.nx + (x / 2);

                    Pixel p = makePixel(r10, r20, i10, i20);

                    for(size_t k = 0; k < VI::Count; k++){
                        double v = VI::funcs[k](p);
                        if(std::isfinite(v))
                            out[k][i10] = (float)v;
                    }
            }
        }

        const double gt[6] = { r10.xMin, r10.pixelSize, 0.0, r10.yMax, 0.0, -r10.pixelSize };
        const std::string outPath = (fs::path(outputDir) / (key + ".tif")).string();
        
        if(!writeIndexTiff(outPath, out, r10.nx, r10.ny, gt, r10.crsWkt.c_str()))
            throw std::runtime_error("Write failed: " + outPath);
        else
            std::cout << "Written: " << key << std::endl;
    }

    return 0;
}

int runCreateVegetationIndices(const AppOptions& options){
    GDALAllRegister();

    std::string tifDir = options.tifDirSat;
    std::string outputDir = options.outputDir;

    if(!tifDir.empty())
        createVegetationIndices(tifDir, outputDir);

    return 0;
}