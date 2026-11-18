#include "AppSupport.h"

#include <exception>
#include <iostream>

int main(int argc, char **argv)
{
    try{
        const bool wantHelp = hasFlag(argc, argv, "--help") || hasFlag(argc, argv, "-h");

        if (wantHelp)
        {
            printUsage("create_vegetation_indices");
            return wantHelp ? 0 : 1;
        }

        initializeGdal();

        int rc = 0;

        const auto options = parseOptions(argc, argv, "create_vegetation_indices");
        rc = runCreateVegetationIndices(options);

        return rc;
    }
    catch(const std::exception &error){
        std::cerr << "Error: " << error.what() << std::endl;
        printUsage("create_vegetation_indices");
        return 1;
    }
}