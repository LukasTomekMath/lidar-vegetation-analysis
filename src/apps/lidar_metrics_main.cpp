#include "AppSupport.h"
#include <exception>
#include <iostream>

int main(int argc, char **argv)
{
    try
    {
        // Spracuj argumenty; help nevyzaduje inicializaciu GDAL/PROJ.
        const auto options = parseOptions(argc, argv, "lidar_metrics");
        if (options.help)
        {
            printUsage("lidar_metrics");
            return 0;
        }
        // Nastav PROJ podla umiestnenia EXE a zaregistruj ovladace GDAL.
        initializeGdal();
        // Bez KML pocitaj cele dlazdice, s KML iba prieniky s bounding boxami lokalit.
        return options.kmlDir.empty() ? runAllFiles(options) : runForestFiles(options);
    }
    // Chyby argumentov aj spracovania preved na chybovy navratovy kod.
    catch (const std::exception &error)
    {
        std::cerr << "Error: " << error.what() << "\n";
        printUsage("lidar_metrics");
        return 1;
    }
}
