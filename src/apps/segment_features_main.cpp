#include "AppSupport.h"
#include <cstring>
#include <exception>
#include <iostream>
#include <string>

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