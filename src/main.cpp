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

std::string formatTime(double seconds) {
	int h = static_cast<int>(seconds) / 3600;
	int m = (static_cast<int>(seconds) % 3600) / 60;
	int s = static_cast<int>(seconds) % 60;

	std::ostringstream oss;
	if (h > 0) oss << h << "h ";
	if (m > 0 || h > 0) oss << m << "m ";
	oss << s << "s";
	return oss.str();
}


int main(int argc, char** argv)
{
	int nprocs = 6;
	int files_done = 0;
	omp_set_num_threads(nprocs);

	std::ofstream csv("processing_stats.csv", std::ios::trunc);
	csv << "Tile,File_MB,PointsInMesh,ReadTime_s,NormalizationTime_s,RedistributionTime_s,MetricsTime_s,ExportTime_s,DeallocationTime_s,OverallTime_s\n";

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
	std::uint64_t totalBytes = 0;
	std::uint64_t processedBytes = 0;
	for (const auto& fe : files)
	{
		totalBytes+=std::filesystem::file_size(fe);
	}

	auto overall_start = std::chrono::high_resolution_clock::now();
#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < static_cast<int>(files.size()); ++i)
	{
		auto file_start = std::chrono::high_resolution_clock::now(); //celkovy cas

		std::vector<std::string> fileNameParts;
		double upperLeftX = -1.0;
		double upperLeftY = -1.0;

		DataHandler* handler = new DataHandler;
		std::string lazFileName = files[i].path().filename().string();
		handler->setAreaName(std::regex_replace(lazFileName, std::regex("(\\.laz)"), ""));
		splitString(handler->areaName(), fileNameParts);

//		for (const auto& part : fileNameParts)
//		{
//#pragma omp critical
//			std::cout << "file part: " << part << std::endl;
//		}

		upperLeftX = std::stod(fileNameParts[1]);
		upperLeftY = std::stod(fileNameParts[2]) + 2000.0;

#pragma omp critical
		std::cout << i+1 <<"/"<<files.size() << " Processing file '" << lazFileName << std::endl;

		handler->setupAreaInfo(files[i].path().string(), upperLeftX, upperLeftY);

//#pragma omp critical
//		std::cout << "area name: " << handler->areaName() << std::endl;
		handler->performCalculation();

#pragma omp atomic
		files_done++;

#pragma omp atomic
		processedBytes += std::filesystem::file_size(files[i]);

		OutputData out = handler->getOutputData();

		auto dealloc_start = std::chrono::high_resolution_clock::now();
		delete handler;
		auto dealloc_end = std::chrono::high_resolution_clock::now();


		auto file_end = std::chrono::high_resolution_clock::now();

		double deallocationTime = std::chrono::duration_cast<std::chrono::milliseconds>(dealloc_end - dealloc_start).count() / 1000.0;
		double fileTime = std::chrono::duration_cast<std::chrono::milliseconds>(file_end - file_start).count() / 1000.0;


#pragma omp critical
		csv << lazFileName << ','
			<< (out.fileSizeBytes / (1024.0 * 1024.0)) << ','
			<< out.nPointsInMesh << ','
			<< out.readTime << ','
			<< out.normalizationTime << ','
			<< out.redistributionTime << ','
			<< out.metricsTime << ','
			<< out.exportTime << ','
			<< deallocationTime << ','
			<< fileTime << '\n';

		auto overall_check = std::chrono::high_resolution_clock::now();
		double progressTime = std::chrono::duration_cast<std::chrono::milliseconds>(overall_check - overall_start).count() / 1000.0;

		double bytesLeft = totalBytes - processedBytes;
		double timePerByte = progressTime / processedBytes;
		double estimatedTimeLeft = bytesLeft * timePerByte;

		std::cout << std::endl;
		std::cout << "Done file: " << lazFileName << std::endl;
		std::cout << "Progress: "<<files_done<<"/"<<files.size()
			<< " ("<<processedBytes/(1024.0*1024.0) <<"MB / " 
			<< totalBytes / (1024.0 * 1024.0) << "MB)"
			<< std::endl ;
		std::cout << "Elapsed time: " << formatTime(progressTime) << std::endl;
		std::cout << "Estimated time left: " << formatTime(estimatedTimeLeft) << std::endl << std::endl;

		fileNameParts.clear();
	}
}