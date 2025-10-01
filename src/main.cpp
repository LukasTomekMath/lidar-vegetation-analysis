#define _USE_MATH_DEFINES

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <memory>
#include <chrono>
#include <iostream>
#include <filesystem>
#include <sstream>
#include <regex>

// GDAL include
#include "libs/gdal_include/gdal_priv.h"
#include "libs/gdal_include/gdal.h"

// LASTools include
#include "libs/LAStools_include/LASlib/lasreader.hpp"
#include "libs/LAStools_include/LASlib/laswriter.hpp"

// Data handler include
#include "PointCloudMetrics.h"

namespace fs = std::filesystem;

void splitString(std::string input, std::vector<std::string>& output)
{
	std::istringstream stream(input);
	std::string temp;
	while (std::getline(stream, temp, '_'))
	{
		output.push_back(temp);
	}
}

int main(int argc, char** argv)
{
	int nprocs = 6;
	omp_set_num_threads(nprocs);

	if (argc != 2)
	{
		std::cout << "Expected path to data directory or -h/--help argument." << std::endl;
		return -1;
	}

	std::string inputArg = std::string(argv[1]);

	if (inputArg == "-h" || inputArg == "--help")
	{
		std::cout << "Specify full path to the directory with PC data. For example:\n" << std::endl;
		std::cout << "    'process_PC_to_vegetation_metrics.exe <path to PC data directory>'" << std::endl;
		std::cout << std::endl;
		exit(0);
	}

	std::string dataDirectoryPath = inputArg;
	std::cout << "Scanning for .laz files in '" << dataDirectoryPath << "'..." << std::endl;
	std::vector<std::string> fileNameParts;
	double upperLeftX = -1.0;
	double upperLeftY = -1.0;
	for (const auto& file : fs::directory_iterator(dataDirectoryPath))
	{
		DataHandler* handler = new DataHandler;
		std::string lazFileName = file.path().filename().string();
		handler->setAreaName(std::regex_replace(lazFileName, std::regex("(\\.laz)"), ""));
		splitString(handler->areaName(), fileNameParts);

		for (const auto& part : fileNameParts)
		{
			std::cout << "file part: " << part << std::endl;
		}
		upperLeftX = std::stod(fileNameParts[1]);
		upperLeftY = std::stod(fileNameParts[2]) + 2000.0;
		std::cout << "Processing file '" << lazFileName << "', upper left coords [" << static_cast<int>(upperLeftX) << ", " << static_cast<int>(upperLeftY) << "]..." << std::endl;

		handler->setupAreaInfo(file.path().string(), upperLeftX, upperLeftY);
		std::cout << "area name: " << handler->areaName() << std::endl;
		handler->performCalculation();

		delete handler;
		fileNameParts.clear();
	}
}