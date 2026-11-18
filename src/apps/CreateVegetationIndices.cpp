#include "AppSupport.h"
#include "RasterIO.h"
#include "gdal_priv.h"
#include "MetricDefinitionsVegIndices.h"
#include <iostream>
#include <string>
#include <filesystem>

namespace fs = std::filesystem;

int createVegetationIndices(const std::string tifDir, const std::string outputDir){
    std::cout << "Use case: Load Sentinel-2 tif files from a specified directory, create "
                 "different vegetation indices and save them as tif files.\n";

    if(!fs::exists(tifDir)){
        std::cerr << "Directory does not exist: " << tifDir << std::endl;
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