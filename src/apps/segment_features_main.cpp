#include "AppSupport.h"
#include <cstring>
#include <exception>
#include <iostream>
#include <string>

namespace
{
bool hasFlag(int argc, char **argv, const char *flag)
{
    const std::size_t len = std::strlen(flag);
    for (int i = 1; i < argc; ++i)
    {
        // matches "--kml-dir <v>" and "--kml-dir=<v>"
        if (std::strncmp(argv[i], flag, len) == 0 &&
            (argv[i][len] == '\0' || argv[i][len] == '='))
            return true;
    }
    return false;
}
} // namespace

int main(int argc, char **argv)
{
    try
    {
        const bool wantLaz = hasFlag(argc, argv, "--kml-dir");
        const bool wantSat = hasFlag(argc, argv, "--shp-dir");
        const bool wantHelp = hasFlag(argc, argv, "--help") || hasFlag(argc, argv, "-h");

        if (wantHelp || (!wantLaz && !wantSat))
        {
            printUsage("segment_features_laz");
            printUsage("segment_features_sat");
            return wantHelp ? 0 : 1;
        }

        initializeGdal();

        int rc = 0;

        if (wantLaz)
        {
            const auto options = parseOptions(argc, argv, "segment_features_laz");
            rc = runSegmentFeatures(options);
        }
        if (wantSat)
        {
            const auto options = parseOptions(argc, argv, "segment_features_sat");
            const int satRc = runSegmentFeatures(options);
            if (rc == 0)
                rc = satRc;
        }

        return rc;
    }
    catch (const std::exception &error)
    {
        std::cerr << "Error: " << error.what() << "\n";
        printUsage("segment_features_laz");
        printUsage("segment_features_sat");
        return 1;
    }
}