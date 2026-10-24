#include "AppSupport.h"
#include <exception>
#include <iostream>

int main(int argc, char **argv)
{
    try
    {
        const auto options = parseOptions(argc, argv, "segment_features");
        if (options.help)
        {
            printUsage("segment_features");
            return 0;
        }
        initializeGdal();
        return runSegmentFeatures(options);
    }
    catch (const std::exception &error)
    {
        std::cerr << "Error: " << error.what() << "\n";
        printUsage("segment_features");
        return 1;
    }
}
