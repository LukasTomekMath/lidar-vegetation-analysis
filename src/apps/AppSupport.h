#pragma once
#include <filesystem>
#include <string>
#include <vector>

// Spolocne argumenty programov; nepouzite cesty zostavaju prazdne.
struct AppOptions
{
    std::string lazDir;
    std::string shpDir;
    std::string kmlDir;
    std::string tifDirSat;
    std::string tifDirLaz;
    // Priecinok vystupov; segment_features namiesto neho pouziva output (cestu k CSV).
    std::string outputDir;
    std::string outputLaz;
    std::string outputSat;
    // Nula ponecha predvoleny pocet vlakien konkretneho rezimu.
    int threads = 0;
    bool help = false;
};

// Overi argumenty pre dany program, prevedie cesty na absolutne a vytvori vystupne priecinky.
// Pri chybe vyhodi vynimku; pri --help vrati iba nastavene help bez dalsich kontrol.
AppOptions parseOptions(int argc, char **argv, const std::string &app);
// Vypise syntax spustenia a kratky popis daneho programu do konzoly.
void printUsage(const std::string &app);
// Nastavi vyhladavanie PROJ databazy pri EXE a zaregistruje ovladace GDAL.
void initializeGdal();
// Rozdeli text podla podciarkovnika a prida casti do output; povodny obsah nemaze.
void splitString(std::string input, std::vector<std::string> &output);
// Cas "seconds" v sekundach prepocita na "hodiny minuty sekundy" avypise
std::string formatTime(double seconds);
// Overi konecne ciselne suradnice v druhej a tretej casti nazvu oddelenych podciarkovnikom.
void validateTileName(const std::filesystem::path &path);
// Vrati priecinok spusteneho EXE, nezavisle od aktualneho pracovneho adresara.
std::filesystem::path executableDirectory();

// Vstupne funkcie rezimov: metriky po dlazdiciach alebo lokalitach, agregacia,
// orezanie LAZ a kontrolny vypis lokalit.
int runAllFiles(const AppOptions &options);
int runForestFiles(const AppOptions &options);
int runSegmentFeatures(const AppOptions &options);
int runCropLaz(const AppOptions &options);
int runInspectForests(const AppOptions &options);
