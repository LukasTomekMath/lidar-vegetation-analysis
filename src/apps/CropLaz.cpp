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
#include <sstream>
#include <stdexcept>
#include <vector>
#include <windows.h>

namespace fs = std::filesystem;

namespace
{
// Spusti las2las pre jeden vstupny subor a vyrez podla obdlznika XY.
int cropFile(const fs::path &executable, const fs::path &input, const fs::path &output, double minX, double minY,
             double maxX, double maxY)
{
    // Zostav argumenty -keep_xy s desatinnou bodkou a suradnicami na tri miesta.
    std::wostringstream command;
    command.imbue(std::locale::classic());
    const auto executablePath = executable.wstring();
    command << L'"' << executablePath << L"\" -i \"" << input.wstring() << L"\" -keep_xy " << std::fixed
            << std::setprecision(3) << minX << L' ' << minY << L' ' << maxX << L' ' << maxY << L" -o \""
            << output.wstring() << L'"';
    std::wstring arguments = command.str();
    STARTUPINFOW startup{};
    startup.cb = sizeof(startup);
    PROCESS_INFORMATION process{};
    // Spusti EXE priamo bez cmd.exe; cesty a nazvy sa nevyhodnocuju ako prikazy.
    if (!CreateProcessW(executablePath.c_str(), arguments.data(), nullptr, nullptr, FALSE, CREATE_NO_WINDOW, nullptr,
                        nullptr, &startup, &process))
        return -1;
    // Pockaj na dokoncenie vyrezu, zisti navratovy kod a uvolni handly procesu.
    const DWORD waitResult = WaitForSingleObject(process.hProcess, INFINITE);
    DWORD exitCode = 1;
    if (waitResult == WAIT_OBJECT_0)
        GetExitCodeProcess(process.hProcess, &exitCode);
    CloseHandle(process.hThread);
    CloseHandle(process.hProcess);
    return static_cast<int>(exitCode);
}
}

// Nacita lokality z KML suborov a pre kazdu lokalitu oreze prislusne LAZ tiles na prienik s bounding boxom (paddingovanym). Vyrezy zapise do suborov.
int runCropLaz(const AppOptions &options)
{
    // Cesta ku las2las voci pracovnemu adresaru
    const fs::path las2las = fs::absolute("libs/las2las64.exe");
    if (!fs::is_regular_file(las2las))
        throw std::runtime_error("Missing libs/las2las64.exe. Run crop_laz from the project root.");

    // Pocet vlakien pre paralelne spracovanie lokalit, default 6.
    int nprocs = options.threads ? options.threads : 6;
    omp_set_num_threads(nprocs);
    GDALAllRegister();

    std::cout << "Use case: Crop laz files according to the boundin boxes defined by forests in kml files."
              << std::endl;

    std::string kmlDir = options.kmlDir;
    std::string lazDir = options.lazDir;

    if (!fs::exists(kmlDir))
    {
        std::cerr << "Directory does not exist: " << kmlDir << "\n";
        return -1;
    }

    ForestManager forestManager(spatial::directoryLasCRS(lazDir));
    forestManager.loadDirectory(kmlDir);
    auto &forests = forestManager.getForests();

    // Pocet neuspesnych volani las2las zo vsetkych vlakien.
    std::atomic<int> failures{0};

// Lokality sa delia medzi vlakna priebezne; subory jednej lokality idu postupne.
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < forests.size(); i++)
    {
        Forest &forest = forests[i];

        // Vypocet plochy, zaroven odstrani vnutorne hranice s plochou mensou nez 1 m2.
        forest.calculateForestArea();
        // Z opisaneho obdlznika (+padding) odvod nazvy suborov pre dlazdice 2 x 2 km.
        forest.findBoundingBox();
        forest.findTiles();

        std::vector<std::string> &subory = forest.getTiles();

        // Prejdi vsetky nazvy z findTiles(), vratane variantov pomenovania a duplicit.
        for (int j = 0; j < subory.size(); j++)
        {
            fs::path fullPath = fs::path(lazDir) / subory[j];
            if (fs::exists(fullPath))
            {
#pragma omp critical
                {
                    std::cout << "Processing forest: " << forest.getNameWithoutDiacritics() << " and file " << subory[j] << std::endl;
                }

                fs::path outputDir = options.outputDir;

                // Zachovaj povodny nazov vstupu a pridaj nazov bez diakritiky.
                std::string originalName = fullPath.stem().string();
                const std::string forestName = forest.getNameWithoutDiacritics();
                fs::path outputPath = outputDir / fs::u8path(originalName + "_" + forestName + ".laz");

                // Usporiadaj hranice do poradia minX, minY, maxX, maxY pre las2las.
                double x_min = std::min(forest.getMinX(), forest.getMaxX());
                double x_max = std::max(forest.getMinX(), forest.getMaxX());
                double y_min = std::min(forest.getMinY(), forest.getMaxY());
                double y_max = std::max(forest.getMinY(), forest.getMaxY());

                // las2las prenasa CRS aj ostatne metadata; vystup uz neupravujeme.
                int result = cropFile(las2las, fullPath, outputPath, x_min, y_min, x_max, y_max);
                if (result != 0)
                    ++failures;

                if (result == 0)
                {
#pragma omp critical
                    {
                        std::cout << "Successfully saved to: " << outputPath.string() << "\n\n";
                    }
                }
                else
                {
#pragma omp critical
                    {
                        std::cerr << "Error: las2las64 returned code " << result << "\n\n";
                    }
                }
            }

            else
            {
                // Chybajuci kandidat sa iba vypise; nezvysuje pocet chyb orezania.
                std::cout << forest.getNameWithoutDiacritics() << "\n";
                std::cout << "[MISSING] " << subory[j] << "\n";
            }
        }
    }

    // Nenulovy kod znamena, ze aspon jedno spustenie las2las zlyhalo.
    return failures == 0 ? 0 : 1;
}
