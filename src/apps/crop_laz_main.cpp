#include "AppSupport.h"
#include <exception>
#include <iostream>

int main(int argc, char **argv)
{
    try
    {
        // Spracuj argumenty; help nevyzaduje inicializaciu GDAL/PROJ.
        const auto options = parseOptions(argc, argv, "crop_laz");
        if (options.help)
        {
            printUsage("crop_laz");
            return 0;
        }
        // Nastav PROJ podla umiestnenia EXE a zaregistruj ovladace GDAL.
        initializeGdal();
        return runCropLaz(options);
    }
    // Chyby argumentov aj spracovania preved na chybovy navratovy kod.
    catch (const std::exception &error)
    {
        std::cerr << "Error: " << error.what() << "\n";
        printUsage("crop_laz");
        return 1;
    }
}
