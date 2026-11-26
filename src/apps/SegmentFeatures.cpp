#include "AppSupport.h"

#include "FeatureVector.h"
#include "FeatureVectorSat.h"
#include "FeatureVectorVI.h"

#include "MetricDefinitions.h"
#include "MetricDefinitionsSentinel.h"
#include "MetricDefinitionsVegIndices.h"

#include "ForestManager.h"
#include "RasterIO.h"
#include "gdal_priv.h"
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <optional>
#include <omp.h>
#include <atomic>

namespace fs = std::filesystem;

namespace {

struct SatTilePair
{
    std::string path10;
    std::string path20;
};

std::string normalizeKey(const std::string& s)
{
    std::string out;
    for (unsigned char c : s)
        if (c < 128 && std::isalnum(c))
            out += (char)std::tolower(c);
    return out;
}

// Returns cache key (EPSG code if known, else WKT). Fills `crs`.
std::string crsKeyFromWkt(const std::string &wkt, OGRSpatialReference &crs, const std::string &ctx)
{
    if (wkt.empty() || crs.importFromWkt(wkt.c_str()) != OGRERR_NONE)
        throw std::runtime_error("Unknown raster CRS: " + ctx);
    const char *code = crs.GetAuthorityCode(nullptr);
    return code ? std::string("EPSG:") + code : wkt;
}

// Header-only: no pixel read.
std::string readCrsWkt(const std::string &path)
{
    GDALDataset *ds = (GDALDataset *)GDALOpenEx(path.c_str(), GDAL_OF_RASTER | GDAL_OF_READONLY,
                                                nullptr, nullptr, nullptr);
    if (!ds) throw std::runtime_error("Cannot open: " + path);
    std::string wkt = ds->GetProjectionRef() ? ds->GetProjectionRef() : "";
    GDALClose(ds);
    return wkt;
}

struct ForestSet {
    std::unique_ptr<ForestManager> mgr;
    std::map<std::string, std::vector<Forest*>> groups;
};

ForestSet loadForestSet(const std::string &shpDir, const OGRSpatialReference &crs)
{
    ForestSet s;
    s.mgr = std::make_unique<ForestManager>(crs);
    s.mgr->loadDirectory(shpDir);
    // never load into mgr after this, Forest* would dangle
    for (auto &f : s.mgr->getForests())
    {
        f.calculateForestArea();
        f.findBoundingBox();
        s.groups[normalizeKey(f.getName() + "_" + f.getSegment())].push_back(&f);
    }
    return s;
}

std::map<std::string, SatTilePair> collectSatPairs(const std::string &dir)
{
    std::map<std::string, SatTilePair> pairs;
    for (const auto &entry : fs::directory_iterator(fs::u8path(dir)))
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

int segmentFeaturesSat(const std::string shpDir, const std::string tifDirSat, const std::string outputSat)
{
    std::cout << "Use case: Load Sentinel-2 tif files from a specified directory, match them to shapefile "
                 "polygons from a specified directory and calculate feature vectors."
              << std::endl;

    if (!fs::exists(shpDir))    { std::cerr << "Directory does not exist: " << shpDir << "\n";    return -1; }
    if (!fs::exists(tifDirSat)) { std::cerr << "Directory does not exist: " << tifDirSat << "\n"; return -1; }

    const std::string csvFile = getOutputName(tifDirSat, outputSat);
    initializeFeatureVectorSatCSV(csvFile);

    // serial setup: build work list + fill cache
    std::vector<std::pair<std::string, SatTilePair>> work;
    std::map<std::string, ForestSet> cache;

    for (const auto &[key, pair] : collectSatPairs(tifDirSat))
    {
        if (pair.path10.empty() || pair.path20.empty())
        {
            std::cerr << "Skipping " << key << ": missing 10m or 20m file\n";
            continue;
        }
        work.emplace_back(key, pair);

        OGRSpatialReference crs;
        const std::string crsKey = crsKeyFromWkt(readCrsWkt(pair.path10), crs, pair.path10);
        if (!cache.count(crsKey))
            cache.emplace(crsKey, loadForestSet(shpDir, crs));
    }
    // cache is read-only from here on

    std::exception_ptr firstError = nullptr;
    std::vector<std::optional<FeatureVectorSat>> results(work.size());
    std::atomic<size_t> done{0};
    const size_t total = work.size();

        // all cerr output goes through critical(io): one lock, no interleaving
    auto logLine = [&](const std::string &msg) {
#pragma omp critical(io)
        std::cerr << "\r\033[K" << msg << "\n";
    };
    auto progress = [&](const std::string &key) {
#pragma omp critical(io)
        std::cerr << "\r\033[K[" << ++done << "/" << total << "] " << key << std::flush;
    };

#pragma omp parallel for schedule(dynamic)
    for (long i = 0; i < work.size(); ++i)
    {
        try
        {
            const auto &[key, pair] = work[i];

            RasterData r10 = readGeoTiff(pair.path10);
            RasterData r20 = readGeoTiff(pair.path20);

            if (r10.nBands != static_cast<int>(sentinel::Count10))
                throw std::runtime_error("Expected " + std::to_string(sentinel::Count10) + " bands: " + pair.path10);
            if (r20.nBands != static_cast<int>(sentinel::Count20))
                throw std::runtime_error("Expected " + std::to_string(sentinel::Count20) + " bands: " + pair.path20);
            if (r10.crsWkt != r20.crsWkt)
                throw std::runtime_error("10m and 20m CRS differ for " + key);

            OGRSpatialReference crs;
            const std::string crsKey = crsKeyFromWkt(r10.crsWkt, crs, pair.path10);

            const auto cit = cache.find(crsKey);          // find read-only
            if (cit == cache.end())
                throw std::runtime_error("CRS not in cache: " + crsKey);

            const auto git = cit->second.groups.find(normalizeKey(key));
            if (git == cit->second.groups.end()){
                logLine("No segment for " + key);
                progress(key);
                continue;
            }
            const auto &parts = git->second;

            const double xMin = r10.xMin;
            const double xMax = r10.xMin + r10.nx * r10.pixelSize;
            const double yMax = r10.yMax;
            const double yMin = r10.yMax - r10.ny * r10.pixelSize;

            const bool inside = std::all_of(parts.begin(), parts.end(), [&](const Forest *f) {
                return f->getRawMinX() >= xMin && f->getRawMaxX() <= xMax &&
                       f->getRawMinY() >= yMin && f->getRawMaxY() <= yMax;
            });
            if (!inside){
                logLine("Outside raster: " + key);
                progress(key);
                continue;
            }

            for (Forest *f : parts) f->createMask();

            results[i] = computeFeatureVectorSat(r10, r20, parts);

#pragma omp critical(progress)
            {
                progress(key);
            }
        }
        catch (...)
        {
#pragma omp critical(err)
            if (!firstError) firstError = std::current_exception();
        }
    }

    std::cerr << "\n";
    
    if (firstError) std::rethrow_exception(firstError);

    for (const auto &r : results)
        if (r) exportFeatureVectorSatCSV(*r, csvFile);

    return 0;
}

int segmentFeaturesVI(const std::string shpDir, const std::string tifDirSat, const std::string outputSat)
{
    if (!fs::exists(fs::u8path(shpDir)))    { std::cerr << "Directory does not exist: " << shpDir << "\n";    return -1; }
    if (!fs::exists(fs::u8path(tifDirSat))) { std::cerr << "Directory does not exist: " << tifDirSat << "\n"; return -1; }

    const std::string csvFile = getOutputName(tifDirSat, outputSat);
    initializeFeatureVectorVICSV(csvFile);

    // serial setup: work list + CRS cache
    struct Job { std::string path, key; };
    std::vector<Job> work;
    std::map<std::string, ForestSet> cache;

    for (const auto &entry : fs::directory_iterator(fs::u8path(tifDirSat)))
    {
        if (!entry.is_regular_file() || entry.path().extension() != ".tif")
            continue;

        Job job{ entry.path().u8string(), normalizeKey(entry.path().stem().u8string()) };

        OGRSpatialReference crs;
        const std::string crsKey = crsKeyFromWkt(readCrsWkt(job.path), crs, job.path);
        if (!cache.count(crsKey))
            cache.emplace(crsKey, loadForestSet(shpDir, crs));

        work.push_back(std::move(job));
    }
    // directory_iterator order is unspecified: sort for deterministic CSV
    std::sort(work.begin(), work.end(), [](const Job &a, const Job &b) { return a.key < b.key; });
    // cache is read-only from here on

    std::vector<std::optional<FeatureVectorVI>> results(work.size());
    std::exception_ptr firstError = nullptr;
    std::atomic<size_t> done{0};
    const long long total = static_cast<long long>(work.size());

    // all cerr output goes through critical(io): one lock, no interleaving
    auto logLine = [&](const std::string &msg) {
#pragma omp critical(io)
        std::cerr << "\r\033[K" << msg << "\n";
    };
    auto progress = [&](const std::string &key) {
#pragma omp critical(io)
        std::cerr << "\r\033[K[" << ++done << "/" << total << "] " << key << std::flush;
    };

#pragma omp parallel for schedule(dynamic)
    for (long long i = 0; i < total; ++i)
    {
        const Job &job = work[i];
        try
        {
            RasterData raster = readGeoTiff(job.path);
            if (raster.nBands != static_cast<int>(VI::Count))
                throw std::runtime_error("Expected " + std::to_string(VI::Count) + " bands: " + job.path);

            OGRSpatialReference crs;
            const std::string crsKey = crsKeyFromWkt(raster.crsWkt, crs, job.path);

            const auto cit = cache.find(crsKey);
            if (cit == cache.end())
                throw std::runtime_error("CRS not in cache: " + crsKey);

            const auto git = cit->second.groups.find(job.key);
            if (git == cit->second.groups.end())
            {
                logLine("No segment for " + job.key);
                progress(job.key);
                continue;
            }
            const auto &parts = git->second;

            const double xMin = raster.xMin;
            const double xMax = raster.xMin + raster.nx * raster.pixelSize;
            const double yMax = raster.yMax;
            const double yMin = raster.yMax - raster.ny * raster.pixelSize;

            const bool inside = std::all_of(parts.begin(), parts.end(), [&](const Forest *f) {
                return f->getRawMinX() >= xMin && f->getRawMaxX() <= xMax &&
                       f->getRawMinY() >= yMin && f->getRawMaxY() <= yMax;
            });
            if (!inside)
            {
                logLine("Outside raster: " + job.key);
                progress(job.key);
                continue;
            }

            for (Forest *f : parts) f->createMask();

            results[i] = computeFeatureVectorVI(raster, parts);   // own slot, no race
            progress(job.key);
        }
        catch (...)
        {
#pragma omp critical(err)
            if (!firstError) firstError = std::current_exception();
            progress(job.key);   // count failed tifs too
        }
    }
    std::cerr << "\n";

    if (firstError) std::rethrow_exception(firstError);

    // serial write, tif order
    for (const auto &r : results)
        if (r) exportFeatureVectorVICSV(*r, csvFile);

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

    std::string mode = options.mode;

    if(!tifDirLaz.empty())
        segmentFeaturesLaz(kmlDir, tifDirLaz, outputLaz);
    if(!tifDirSat.empty()){
        if(mode == "vegetation_index")
            segmentFeaturesVI(shpDir, tifDirSat, outputSat);
        else
            segmentFeaturesSat(shpDir, tifDirSat, outputSat);
    }

    return 0;
}