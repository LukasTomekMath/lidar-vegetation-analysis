#include "AppSupport.h"
#include "gdal_priv.h"
#include "ogr_srs_api.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <windows.h>

namespace fs = std::filesystem;

void printUsage(const std::string &app)
{
    std::cout << "Usage: " << app;
    if (app == "lidar_metrics")
        std::cout << " --laz-dir DIR [--kml-dir DIR] --output-dir DIR [--threads N]\n"
                  << "Compute 2 m metric rasters for tiles or forest bounding boxes.\n";
    else if (app == "segment_features")
        std::cout << " --kml-dir DIR --tif-dir DIR --output FILE.csv\n"
                  << "Aggregate 27-band metric rasters inside segment masks.\n";
    else if (app == "crop_laz")
        std::cout << " --kml-dir DIR --laz-dir DIR --output-dir DIR [--threads N]\n"
                  << "Crop LAZ files to forest bounding boxes using libs/las2las64.exe. Run from the project root.\n";
    else
        std::cout << " --kml-dir DIR --laz-dir DIR --output-dir DIR\n"
                  << "Write pralesy.txt and tiles.txt and report file availability.\n";
    std::cout << "  -h, --help  Show this help. Quote paths containing spaces.\n";
}

AppOptions parseOptions(int argc, char **argv, const std::string &app)
{
    // help ma prednost aj pred kontrolou ostatnych argumentov.
    AppOptions result;
    for (int i = 1; i < argc; ++i)
    {
        if (std::string(argv[i]) == "--help" || std::string(argv[i]) == "-h")
        {
            result.help = true;
            return result;
        }
    }
    // Kazdy program prijima iba svoje podporovane prepinace.
    std::set<std::string> allowed;
    if (app == "lidar_metrics")
        allowed = {"--laz-dir", "--kml-dir", "--output-dir", "--threads"};
    else if (app == "segment_features")
        allowed = {"--kml-dir", "--tif-dir", "--output"};
    else if (app == "crop_laz")
        allowed = {"--kml-dir", "--laz-dir", "--output-dir", "--threads"};
    else
        allowed = {"--kml-dir", "--laz-dir", "--output-dir"};
    // Argumenty citaj ako dvojice --prepinac hodnota; duplicity odmietni.
    std::map<std::string, std::string> values;
    for (int i = 1; i < argc; ++i)
    {
        std::string key = argv[i];
        if (!allowed.count(key))
            throw std::runtime_error("Unknown option: " + key);
        if (values.count(key))
            throw std::runtime_error("Duplicate option: " + key);
        if (++i == argc || std::string(argv[i]).rfind("--", 0) == 0 || !*argv[i])
            throw std::runtime_error("Missing value for " + key);
        values[key] = argv[i];
    }
    // Povinne cesty vztiahni k pracovnemu adresaru a odstran z nich . a .. .
    // Ide o lexikalnu normalizaciu, nie o rozlisenie symbolickych odkazov.
    auto required = [&](const std::string &key) {
        if (!values.count(key))
            throw std::runtime_error("Required option: " + key);
        return fs::absolute(values.at(key)).lexically_normal().string();
    };
    // KML je volitelne iba pre lidar_metrics; jeho pritomnost vybera spracovanie lokalit.
    if (app != "segment_features")
        result.lazDir = required("--laz-dir");
    if (app != "lidar_metrics" || values.count("--kml-dir"))
        result.kmlDir = required("--kml-dir");
    if (app == "segment_features")
    {
        result.tifDir = required("--tif-dir");
        result.output = required("--output");
    }
    else
        result.outputDir = required("--output-dir");
    // Vyziadany pocet vlakien musi byt kladne cele cislo bez dalsich znakov.
    if (values.count("--threads"))
    {
        const auto &value = values.at("--threads");
        size_t end = 0;
        result.threads = std::stoi(value, &end);
        if (end != value.size() || result.threads < 1)
            throw std::runtime_error("--threads must be a positive integer");
    }
    // Vstupne priecinky musia existovat; chybajuce vystupne priecinky vytvor.
    for (const auto &dir : {result.lazDir, result.kmlDir, result.tifDir})
        if (!dir.empty() && !fs::is_directory(dir))
            throw std::runtime_error("Input directory does not exist: " + dir);
    if (!result.outputDir.empty())
    {
        // Odmietni zhodnu normalizovanu cestu pre vstupne LAZ a vystupy.
        if (!result.lazDir.empty() && fs::path(result.outputDir) == fs::path(result.lazDir))
            throw std::runtime_error("Output directory must differ from LAZ input directory");
        fs::create_directories(result.outputDir);
    }
    if (!result.output.empty())
        fs::create_directories(fs::path(result.output).parent_path());
    return result;
}

fs::path executableDirectory()
{
    // Windows API vracia cestu k EXE v UTF-16; odmietni aj prilis dlhu cestu.
    std::vector<wchar_t> path(32768);
    const DWORD length = GetModuleFileNameW(nullptr, path.data(), static_cast<DWORD>(path.size()));
    if (!length || length >= path.size())
        throw std::runtime_error("Cannot locate executable directory");
    return fs::path(std::wstring(path.data(), length)).parent_path();
}

void initializeGdal()
{
    // Pouzi proj.db dodanu s programom, aby transformacie nezaviseli od pracovneho adresara.
    const auto directory = executableDirectory();
    if (!fs::is_regular_file(directory / "proj.db"))
        throw std::runtime_error("proj.db is missing beside the executable");
    const std::string projDirectory = directory.u8string();
    const char *paths[] = {projDirectory.c_str(), nullptr};
    OSRSetPROJSearchPaths(paths);
    GDALAllRegister();
}

void splitString(std::string input, std::vector<std::string> &output)
{
    std::istringstream stream(input);
    std::string temp;
    while (std::getline(stream, temp, '_'))
        output.push_back(temp);
}

void validateTileName(const fs::path &path)
{
    // Z nazvu bez pripony ocakavaj aspon tri casti; dalsie casti su povolene.
    std::vector<std::string> parts;
    splitString(path.stem().string(), parts);
    if (parts.size() < 3)
        throw std::runtime_error("Expected tile_X_Y[...].laz: " + path.string());
    for (int i : {1, 2})
    {
        size_t end = 0;
        const double coordinate = std::stod(parts[i], &end);
        if (end != parts[i].size() || !std::isfinite(coordinate))
            throw std::runtime_error("Invalid tile coordinates: " + path.string());
    }
}

std::string formatTime(double seconds)
{
    int h = static_cast<int>(seconds) / 3600;
    int m = (static_cast<int>(seconds) % 3600) / 60;
    int s = static_cast<int>(seconds) % 60;
    std::ostringstream oss;
    if (h > 0)
        oss << h << "h ";
    if (m > 0 || h > 0)
        oss << m << "m ";
    oss << s << "s";
    return oss.str();
}

