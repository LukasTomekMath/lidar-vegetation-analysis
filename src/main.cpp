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

	std::ofstream stats_file("processing_stats.txt", std::ios::trunc);

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
	//vytvaranie vektora so subormi na paralelne spracovanie
	std::vector<fs::directory_entry> files;
	for (const auto& file : fs::directory_iterator(dataDirectoryPath))
	{
		files.push_back(file);
	}

#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < static_cast<int>(files.size()); ++i)
	{
		auto overall_start = std::chrono::high_resolution_clock::now(); //celkovy cas

		std::vector<std::string> fileNameParts;
		double upperLeftX = -1.0;
		double upperLeftY = -1.0;

		DataHandler* handler = new DataHandler;
		std::string lazFileName = files[i].path().filename().string();
		handler->setAreaName(std::regex_replace(lazFileName, std::regex("(\\.laz)"), ""));
		splitString(handler->areaName(), fileNameParts);

		for (const auto& part : fileNameParts)
		{
#pragma omp critical
			std::cout << "file part: " << part << std::endl;
		}

		upperLeftX = std::stod(fileNameParts[1]);
		upperLeftY = std::stod(fileNameParts[2]) + 2000.0;

#pragma omp critical
		std::cout << "Processing file '" << lazFileName << "', upper left coords [" << static_cast<int>(upperLeftX) << ", " << static_cast<int>(upperLeftY) << "]..." << std::endl;

		handler->setupAreaInfo(files[i].path().string(), upperLeftX, upperLeftY);

#pragma omp critical
		std::cout << "area name: " << handler->areaName() << std::endl;
		handler->performCalculation();

		ProcessingTimes times = handler->getProcessingTimes();

		auto dealloc_start = std::chrono::high_resolution_clock::now();
		delete handler;
		auto dealloc_end = std::chrono::high_resolution_clock::now();


		auto overall_end = std::chrono::high_resolution_clock::now();

		double deallocationTime = std::chrono::duration_cast<std::chrono::milliseconds>(dealloc_end - dealloc_start).count() / 1000.0;
		double overallTime = std::chrono::duration_cast<std::chrono::milliseconds>(overall_end - overall_start).count() / 1000.0;

#pragma omp critical
		{
			std::ofstream stats_file("processing_stats.txt", std::ios::app);
			stats_file << "Tile: " << lazFileName << "\n";
			stats_file << "Reading time: " << times.readTime << " s\n";
			stats_file << "Normalization time: " << times.normalizationTime << " s\n";
			stats_file << "Redistribution time: " << times.redistributionTime << " s\n";
			stats_file << "Metrics time: " << times.metricsTime << " s\n";
			stats_file << "Export time: " << times.exportTime << " s\n";
			stats_file << "Deallocation time: " << deallocationTime << " s\n";
			stats_file << "Overall time: " << overallTime << " s\n";
			stats_file << "-----------------------------\n";
		}

		fileNameParts.clear();
	}
}