#include "AppSupport.h"
#include "ForestManager.h"
#include "RasterIO.h"
#include "PointCloudMetrics.h"
#include <set>
#include <cctype>
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

// Vstupny subor a jeho velkost pre odhad priebehu spracovania.
struct FileInfo
{
    fs::directory_entry entry;
    std::uint64_t size;
};

// Spracuje cele dlazdice 2 x 2 km zo vstupneho priecinka na rastre metrik.
// Vystupne rozlisenie je 2 m; KML sa v tomto rezime nepouziva.
int runAllFiles(const AppOptions &options)
{
    // Predvolene pat vlakien; kazde spracuva vlastny DataHandler.
    int nprocs = options.threads ? options.threads : 5;
    int files_done = 0;
    std::atomic<int> failures{0};
    omp_set_num_threads(nprocs);

    int pixelsize = 2;

    std::cout << "Use case: Load all laz files in a specified directory and calculate tifs." << std::endl;
    std::cout << "Pixel size is set to: " << pixelsize << std::endl;
    std::cout << "Number of threads: " << nprocs << std::endl;

    // Prepis statistiky predosleho behu; kazdy uspesny vypocet prida jeden riadok.
    std::ofstream csv(fs::path(options.outputDir) / "processing_stats.csv", std::ios::trunc);
    if (!csv)
        throw std::runtime_error("Cannot create processing_stats.csv");
    csv << "Tile,File_MB,PointsInMesh,ReadTime_s,NormalizationTime_s,RedistributionTime_s,MetricsTime_s,ExportTime_s,"
           "DeallocationTime_s,OverallTime_s\n";
    csv.flush();

    std::string dataDirectoryPath = options.lazDir;
    std::cout << "Scanning for .laz files in '" << dataDirectoryPath << "'..." << std::endl;
    // vytvaranie vektora so subormi na paralelne spracovanie
    std::vector<FileInfo> files;
    for (const auto &file : fs::directory_iterator(dataDirectoryPath))
    {
        auto path = file.path();
        std::string filename = path.filename().string();

        // Vyber iba bezne .laz subory; COPC subory sa v tomto rezime vynechavaju.
        if (file.is_regular_file() && path.extension() == ".laz" && filename.find(".copc.laz") == std::string::npos)
        {
            validateTileName(path);
            auto sz = std::filesystem::file_size(file);
            files.push_back({file, sz});
        }
    }

    std::uint64_t totalBytes = 0;
    std::uint64_t processedBytes = 0;
    for (const auto &fe : files)
    {
        totalBytes += fe.size;
    }

    auto overall_start = std::chrono::high_resolution_clock::now();
    GDALAllRegister();

    // OpenMP prideluje pracu dynamicky po skupinach piatich poloziek.
    const int chunk = 5;

#pragma omp parallel for schedule(dynamic, chunk)
    for (int i = 0; i < static_cast<int>(files.size()); ++i)
    {
        auto file_start = std::chrono::high_resolution_clock::now(); // celkovy cas

        std::vector<std::string> fileNameParts;
        double upperLeftX = -1.0;
        double upperLeftY = -1.0;

        auto handler = std::make_unique<DataHandler>();
        handler->setOutputDirectory(options.outputDir);
        std::string lazFileName = files[i].entry.path().filename().string();
        handler->setAreaName(std::regex_replace(lazFileName, std::regex("(\\.laz)"), ""));
        splitString(handler->areaName(), fileNameParts);

        //		for (const auto& part : fileNameParts)
        //		{
        // #pragma omp critical
        //			std::cout << "file part: " << part << std::endl;
        //		}

        // Nazov tile_X_Y obsahuje lavy dolny roh; pre horny roh pridaj 2000 m k Y.
        upperLeftX = std::stod(fileNameParts[1]);
        upperLeftY = std::stod(fileNameParts[2]) + 2000.0;

        // upperLeftX = -450170.0;
        // upperLeftY = -1228390.0;

        std::ostringstream startMsg;
        startMsg << i + 1 << "/" << files.size() << " Processing file '" << lazFileName << "'\n";
#pragma omp critical
        {
            std::cout << startMsg.str() << std::flush;
        }

        // Nastav rozsah celej dlazdice a velkost vystupneho pixelu.
        handler->setupAreaInfo(files[i].entry.path().string(), upperLeftX, upperLeftY, upperLeftX + 2000,
                               upperLeftY - 2000, pixelsize);
        // handler->setupAreaInfo(files[i].entry.path().string(), upperLeftX, upperLeftY, -449300.0,
        // -1229310.0,pixelsize);

        // #pragma omp critical
        //		std::cout << "area name: " << handler->areaName() << std::endl;
        // Nacitanie, filtre, normalizaciu, redistribuciu, metriky a export vykona DataHandler.
        // Pri chybe vypoctu zaznamenaj zlyhanie a pokracuj dalsim suborom.
        try
        {
            if (!handler->performCalculation())
            {
                ++failures;
                continue;
            }
        }
        catch (const std::exception &error)
        {
            ++failures;
#pragma omp critical
            {
                std::cerr << "Calculation failed: " << error.what() << "\n";
            }
            continue;
        }

        // Uchovaj statistiky pred zrusenim handlera a zmeraj aj uvolnenie pamate.
        OutputData out = handler->getOutputData();

        auto dealloc_start = std::chrono::high_resolution_clock::now();
        handler.reset();
        auto dealloc_end = std::chrono::high_resolution_clock::now();

        auto file_end = std::chrono::high_resolution_clock::now();

        double deallocationTime =
            std::chrono::duration_cast<std::chrono::milliseconds>(dealloc_end - dealloc_start).count() / 1000.0;
        double fileTime = std::chrono::duration_cast<std::chrono::milliseconds>(file_end - file_start).count() / 1000.0;

        // Riadok CSV priprav lokalne, aby sa vypisy jednotlivych vlakien nemiesali.
        std::ostringstream local_csv;
        local_csv << lazFileName << ',' << (out.fileSizeBytes / (1024.0 * 1024.0)) << ',' << out.nPointsInMesh << ','
                  << out.readTime << ',' << out.normalizationTime << ',' << out.redistributionTime << ','
                  << out.metricsTime << ',' << out.exportTime << ',' << deallocationTime << ',' << fileTime << '\n';

        // Pocitadla zahrnaju iba uspesne vypocty; aktualizacia je chranena medzi vlaknami.
        int local_done;
        std::uint64_t local_processedBytes;
#pragma omp critical
        {
            files_done++;
            local_done = files_done;

            processedBytes += files[i].size;
            local_processedBytes = processedBytes;
        }

        auto overall_check = std::chrono::high_resolution_clock::now();
        double progressTime =
            std::chrono::duration_cast<std::chrono::milliseconds>(overall_check - overall_start).count() / 1000.0;
        // Odhad zostavajuceho casu vychadza z objemu vstupov, nie z poctu bodov.
        double bytesLeft = totalBytes - local_processedBytes;
        double timePerByte = local_processedBytes ? progressTime / local_processedBytes : 0.0;
        double estimatedTimeLeft = bytesLeft * timePerByte;

#pragma omp critical
        {
            csv << local_csv.str();
            csv.flush();
            std::cout << "\nDone file: " << lazFileName << "\n"
                      << "Progress: " << local_done << "/" << files.size() << " ("
                      << local_processedBytes / (1024.0 * 1024.0) << "MB / " << totalBytes / (1024.0 * 1024.0)
                      << "MB)\n"
                      << "Elapsed time: " << formatTime(progressTime) << "\n"
                      << "Estimated time left: " << formatTime(estimatedTimeLeft) << "\n\n"
                      << std::flush;
            ;
        }

        fileNameParts.clear();
    }
    // Nenulovy kod oznacuje aspon jeden neuspesny vypocet.
    return failures == 0 ? 0 : 1;
}

// Spracuje prieniky dlazdic s bounding boxami lokalit z KML (+padding).
// Kazdy prienik sa pocita samostatne; raster sa tu nemaskuje presnym polygonom.
int runForestFiles(const AppOptions &options)
{
    // paralelne sa spracuvaju lokality. Default 6 vlaked
    int nprocs = options.threads ? options.threads : 6;
    int files_done = 0;
    std::atomic<int> failures{0};
    omp_set_num_threads(nprocs);

    GDALAllRegister();

    std::cout << "Use case: Load laz files in a specified directory corespondin to kml files from a specified "
                 "directory and calculate tifs."
              << std::endl;
    std::string kmlDir = options.kmlDir;
    std::string lazDir = options.lazDir;

    if (!fs::exists(kmlDir))
    {
        std::cerr << "Directory does not exist: " << kmlDir << "\n";
        return -1;
    }
    if (!fs::exists(lazDir))
    {
        std::cerr << "Directory does not exist: " << lazDir << "\n";
        return -1;
    }

    // Prepis statistiky predosleho behu; kazdy uspesny vypocet prida 1 riadok.
    std::ofstream csv(fs::path(options.outputDir) / "processing_stats.csv", std::ios::trunc);
    if (!csv)
        throw std::runtime_error("Cannot create processing_stats.csv");
    csv << "Tile,File_MB,PointsInMesh,ReadTime_s,NormalizationTime_s,RedistributionTime_s,MetricsTime_s,ExportTime_s,"
           "DeallocationTime_s,OverallTime_s\n";

    ForestManager forestManager(spatial::directoryLasCRS(lazDir));
    forestManager.loadDirectory(kmlDir);
    auto &forests = forestManager.getForests();

    std::uint64_t totalBytes = 0;
    std::uint64_t processedBytes = 0;

    // Spocitaj existujuce vstupy pre odhad priebehu; opakovane nazvy sa zapocitaju znova.
    int num_files = 0;
    for (size_t i = 0; i < forests.size(); ++i)
    {
        Forest &forest = forests[i];

        forest.calculateForestArea();
        // bounding box s paddingom a kandidatske nazvy dlazdic pre lokalitu.
        forest.findBoundingBox();
        forest.findTiles();
        for (const auto &fe : forest.getTiles())
        {
            fs::path fullPath = fs::path(lazDir) / fe;
            if (fs::exists(fullPath))
            {
                validateTileName(fullPath);
                totalBytes += fs::file_size(fullPath);
                num_files++;
            }
        }
    }
    // lesy su nacitane, teraz by sa malo vsetko citat, cize cela laz sekcia kodu

    auto overall_start = std::chrono::high_resolution_clock::now();
    OGRRegisterAll();

    // OpenMP prideluje pracu dynamicky po skupinach piatich poloziek.
    const int chunk = 5;

    // Over jednoznacnost vystupnych nazvov este pred paralelnym vypoctom.
    std::set<std::string> mosaicNames;
    for (const auto& forest : forests)
    {
        const auto name = forest.getNameWithoutDiacritics();
        std::string key = name;
        std::transform(key.begin(), key.end(), key.begin(), [](unsigned char c) { return std::tolower(c); });
        if (name.empty() || name.find_first_of("/\\:*?\"<>|") != std::string::npos ||
            !mosaicNames.insert(key).second)
            throw std::runtime_error("Invalid or duplicate mosaic name: " + name);
    }
    const fs::path tilesDir = fs::path(options.outputDir) / "tiles";
    fs::create_directories(tilesDir);

    // Vlakna dostavaju lokality; subory jednej lokality sa spracuvaju postupne.
    // pre kazdy les musim najst tiles, ktore treba
#pragma omp parallel for schedule(dynamic, chunk)
    for (int i = 0; i < forests.size(); ++i)
    {
        Forest &forest = forests[i];
        auto &files = forest.getTiles(); // toto su tie subory, co chcem precitat pre dany les
        std::vector<std::string> generatedTifs;
        std::set<std::string> processedInputs;
        bool forestFailed = false;

        for (int j = 0; j < files.size(); j++) // prechadzanie cez jednotlive subory
        {
            fs::path fullPath = fs::path(lazDir) / files[j]; // cela cesta k suboru, tu chcem citat

            if (fs::exists(fullPath))
            {
                if (!processedInputs.insert(fullPath.string()).second)
                    continue;

                std::string lazFileName = files[j];
                // co sa ma spravit, ked sa konkretny subor najde - vypocita sa, co sa ma
                auto file_start = std::chrono::high_resolution_clock::now(); // celkovy cas

                std::vector<std::string> fileNameParts;
                double upperLeftX = -1.0;
                double upperLeftY = -1.0;
                double lowerRightX = -1.0;
                double lowerRightY = -1.0;

                auto handler = std::make_unique<DataHandler>();
                handler->setOutputDirectory(tilesDir.string());
                // std::string lazFileName = files[j];
                handler->setAreaName(std::regex_replace(lazFileName, std::regex("(\\.laz)"), ""));
                handler->setForestName(forest.getNameWithoutDiacritics());

                splitString(handler->areaName(), fileNameParts);

                // informacie o danej tile
                upperLeftX = std::stod(fileNameParts[1]);
                upperLeftY = std::stod(fileNameParts[2]) + 2000.0;
                lowerRightX = upperLeftX + 2000.0;
                lowerRightY = std::stod(fileNameParts[2]);

                // Vezmi prienik obdlznika lokality a aktualnej dlazdice.
                double clippedMinX = std::max(forest.getMinX(), upperLeftX);  // upper left X
                double clippedMaxX = std::min(forest.getMaxX(), lowerRightX); // lower right X
                double clippedMinY = std::max(forest.getMinY(), lowerRightY); // lower right Y
                double clippedMaxY = std::min(forest.getMaxY(), upperLeftY);  // upper left Y

                // Prazdny prienik preskoc bez vypoctu aj bez zaznamu do CSV.
                if (clippedMinX >= clippedMaxX || clippedMinY >= clippedMaxY)
                {
                    continue;
                }

                int pixelsize = forest.getPixelSize();

                std::ostringstream startMsg;
                startMsg << i + 1 << "/" << num_files << " Processing file '" << lazFileName << "' for "
                         << forest.getNameWithoutDiacritics() << "\n";
#pragma omp critical
                {
                    std::cout << startMsg.str();
                }

                // DataHandler dostane iba rozsah prieniku a rozlisenie z Forest (2 m).
                handler->setupAreaInfo(fullPath.string(), clippedMinX, clippedMaxY, clippedMaxX, clippedMinY,
                                       pixelsize);

                // Spusti rovnaky vypoctovy postup ako pri celych dlazdiciach.
                // Chyba tohto vypoctu nezastavi spracovanie dalsich vstupov.
                try
                {
                    if (!handler->performCalculation())
                    {
                        ++failures;
                        forestFailed = true;
                        continue;
                    }
                }
                catch (const std::exception &error)
                {
                    ++failures;
                    forestFailed = true;
#pragma omp critical
                    {
                        std::cerr << "Calculation failed: " << error.what() << "\n";
                    }
                    continue;
                }

                // Uchovaj statistiky pred zrusenim handlera a zmeraj uvolnenie pamate.
                generatedTifs.push_back((tilesDir / handler->outputFileName()).string());
                OutputData out = handler->getOutputData();

                auto dealloc_start = std::chrono::high_resolution_clock::now();
                handler.reset();
                auto dealloc_end = std::chrono::high_resolution_clock::now();

                auto file_end = std::chrono::high_resolution_clock::now();

                double deallocationTime =
                    std::chrono::duration_cast<std::chrono::milliseconds>(dealloc_end - dealloc_start).count() / 1000.0;
                double fileTime =
                    std::chrono::duration_cast<std::chrono::milliseconds>(file_end - file_start).count() / 1000.0;

                // Priprav riadok statistik pre tento vypocet lokality a dlazdice.
                std::ostringstream local_csv;
                local_csv << lazFileName << ',' << (out.fileSizeBytes / (1024.0 * 1024.0)) << ',' << out.nPointsInMesh
                          << ',' << out.readTime << ',' << out.normalizationTime << ',' << out.redistributionTime << ','
                          << out.metricsTime << ',' << out.exportTime << ',' << deallocationTime << ',' << fileTime
                          << '\n';

                // Zaznamenaj uspesny vypocet; ten isty LAZ moze patrit viacerym lokalitam.
                int local_done;
                std::uint64_t local_processedBytes;
#pragma omp critical
                {
                    files_done++;
                    local_done = files_done;

                    processedBytes += std::filesystem::file_size(fullPath);
                    local_processedBytes = processedBytes;
                }

                auto overall_check = std::chrono::high_resolution_clock::now();
                double progressTime =
                    std::chrono::duration_cast<std::chrono::milliseconds>(overall_check - overall_start).count() /
                    1000.0;
                // Odhad casu pouziva cele velkosti LAZ, aj ked sa pocita iba ich vyrez.
                double bytesLeft = totalBytes - local_processedBytes;
                double timePerByte = local_processedBytes ? progressTime / local_processedBytes : 0.0;
                double estimatedTimeLeft = bytesLeft * timePerByte;

#pragma omp critical
                {
                    csv << local_csv.str();
                    std::cout << "\nDone file: " << lazFileName << "\n"
                              << "Progress: " << local_done << "/" << num_files << " ("
                              << local_processedBytes / (1024.0 * 1024.0) << "MB / " << totalBytes / (1024.0 * 1024.0)
                              << "MB)\n"
                              << "Elapsed time: " << formatTime(progressTime) << "\n"
                              << "Estimated time left: " << formatTime(estimatedTimeLeft) << "\n\n";
                }

                fileNameParts.clear();
            }
            else
            {
                // Chybajuci kandidat sa iba vypise a nezvysuje pocet zlyhani.
                // co sa ma spravit, ked sa nenajde - nic
                std::cout << forest.getNameWithoutDiacritics() << "\n";
                std::cout << "[MISSING] " << files[j] << "\n";
            }
        }
        // Vsetky uspesne exportovane dlazdice lokality spoj do jedneho TIFF.
        if (!forestFailed)
        {
            try
            {
                const auto output = (fs::path(options.outputDir) /
                    (forest.getNameWithoutDiacritics() + ".tif")).string();
                mergeGeoTiffs(generatedTifs, output);
#pragma omp critical
                {
                    std::cout << "Merged " << generatedTifs.size() << " rasters: " << output << "\n";
                }
            }
            catch (const std::exception& error)
            {
                ++failures;
#pragma omp critical
                {
                    std::cerr << "Merging failed: " << error.what() << "\n";
                }
            }
        }
    }

    // Nenulovy kod oznacuje aspon jeden neuspesny vypocet.
    return failures == 0 ? 0 : 1;
}
