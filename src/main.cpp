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
#include "PolygonManager.h"

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

struct FileInfo { 
	fs::directory_entry entry; 
	std::uint64_t size; 
};

int main_all_files(int argc, char** argv)
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
	std::vector<FileInfo> files;
	for (const auto& file : fs::directory_iterator(dataDirectoryPath))
	{
		auto path = file.path();
		std::string filename = path.filename().string();

		if (path.extension() == ".laz" && filename.find(".copc.laz") == std::string::npos)
		{
			auto sz = std::filesystem::file_size(file);
			files.push_back({ file, sz });
		}
	}

	std::uint64_t totalBytes = 0;
	std::uint64_t processedBytes = 0;
	for (const auto& fe : files)
	{
		totalBytes += fe.size;
	}


	auto overall_start = std::chrono::high_resolution_clock::now();
	GDALAllRegister();

	const int chunk = 5;

#pragma omp parallel for schedule(dynamic, chunk)
	for (int i = 0; i < static_cast<int>(files.size()); ++i)
	{
		auto file_start = std::chrono::high_resolution_clock::now(); //celkovy cas

		std::vector<std::string> fileNameParts;
		double upperLeftX = -1.0;
		double upperLeftY = -1.0;

		DataHandler* handler = new DataHandler;
		std::string lazFileName = files[i].entry.path().filename().string();
		handler->setAreaName(std::regex_replace(lazFileName, std::regex("(\\.laz)"), ""));
		splitString(handler->areaName(), fileNameParts);

		//		for (const auto& part : fileNameParts)
		//		{
		//#pragma omp critical
		//			std::cout << "file part: " << part << std::endl;
		//		}

		upperLeftX = std::stod(fileNameParts[1]);
		upperLeftY = std::stod(fileNameParts[2]) + 2000.0;

		std::ostringstream startMsg;
		startMsg << i + 1 << "/" << files.size() << " Processing file '" << lazFileName << "'\n";
#pragma omp critical
		{
			std::cout << startMsg.str();
		}

		handler->setupAreaInfo(files[i].entry.path().string(), upperLeftX, upperLeftY);

		//#pragma omp critical
		//		std::cout << "area name: " << handler->areaName() << std::endl;
		handler->performCalculation();

		OutputData out = handler->getOutputData();

		auto dealloc_start = std::chrono::high_resolution_clock::now();
		delete handler;
		auto dealloc_end = std::chrono::high_resolution_clock::now();


		auto file_end = std::chrono::high_resolution_clock::now();

		double deallocationTime = std::chrono::duration_cast<std::chrono::milliseconds>(dealloc_end - dealloc_start).count() / 1000.0;
		double fileTime = std::chrono::duration_cast<std::chrono::milliseconds>(file_end - file_start).count() / 1000.0;

		std::ostringstream local_csv;
		local_csv << lazFileName << ','
			<< (out.fileSizeBytes / (1024.0 * 1024.0)) << ','
			<< out.nPointsInMesh << ','
			<< out.readTime << ','
			<< out.normalizationTime << ','
			<< out.redistributionTime << ','
			<< out.metricsTime << ','
			<< out.exportTime << ','
			<< deallocationTime << ','
			<< fileTime << '\n';

		int local_done;
		std::uint64_t local_processedBytes;
#pragma omp critical
		{
			files_done++;
			local_done = files_done;

			processedBytes += files[i].size;
			local_processedBytes = processedBytes;
		}

		auto overall_check = std::chrono::high_resolution_clock::now();
		double progressTime = std::chrono::duration_cast<std::chrono::milliseconds>(overall_check - overall_start).count() / 1000.0;
		double bytesLeft = totalBytes - local_processedBytes;
		double timePerByte = progressTime / local_processedBytes;
		double estimatedTimeLeft = bytesLeft * timePerByte;

#pragma omp critical
		{
			csv << local_csv.str();
			std::cout << "\nDone file: " << lazFileName << "\n"
				<< "Progress: " << local_done << "/" << files.size()
				<< " (" << local_processedBytes / (1024.0 * 1024.0) << "MB / "
				<< totalBytes / (1024.0 * 1024.0) << "MB)\n"
				<< "Elapsed time: " << formatTime(progressTime) << "\n"
				<< "Estimated time left: " << formatTime(estimatedTimeLeft) << "\n\n";
		}

		fileNameParts.clear();
	}
}

int main_curve(int argc, char** argv)
{
	if (argc != 3)
	{
		std::cout << "Usage: test_kml <KML directory>\n";
		return -1;
	}

	std::string kmlDir = argv[1];
	std::string lazDir = argv[2];

	if (!fs::exists(kmlDir))
	{
		std::cerr << "Directory does not exist: " << kmlDir << "\n";
		return -1;
	}

	ForestManager forestManager;
	int fileCount = 0;

	for (auto& entry : fs::directory_iterator(kmlDir))
	{
		if (entry.is_regular_file() && entry.path().extension() == ".kml")
		{
			std::string filename = entry.path().string();
			std::cout << "Loading KML: " << filename << std::endl;

			forestManager.loadFromKML(filename);
			fileCount++;
		}
	}

	auto& forests = forestManager.getForests();

	std::cout << "Finished reading " << fileCount << " KML files.\n";
	std::cout << "Total forests loaded: " << forests.size() << "\n\n";

	// Output file
	std::ofstream file("pralesy.txt");
	if (!file.is_open())
	{
		std::cerr << "Could not create pralesy.txt\n";
		return false;
	}
	std::ofstream tiles("tiles.txt");
	int tile_counter = 0;


	for (size_t i = 0; i < forests.size(); ++i)
	{
		Forest& forest = forests[i];

		//setup na konkretne subory
		forest.calculateForestArea();
		forest.findBoundingBox();
		forest.findTiles();

		//toto je len vypis do suboru ze tore subory treba laz
		tiles << forest.getName() << ": \n";
		std::vector<std::string>& subory = forest.getTiles();
		tile_counter += subory.size();
		for (int k = 0; k <subory.size(); k++)
		{
			tiles << subory[k]<<"\n";
		}
		tiles << "\n";

		auto& polygons = forest.getPolygons();

		file << "=== Forest #" << i << " ===\n";
		file << "Name: " << forest.getName() << "\n";
		file << std::fixed << std::setprecision(6);
		file << "Hectares: " << forest.getForestArea()/10000.0 << "\n";
		file << "Polygon groups: " << polygons.size() << "\n\n";

		for (size_t j = 0; j < polygons.size(); ++j)
		{
			PolygonGroup& pg = polygons[j];

			auto& outerPoints = pg.outer.getPoints();

			file << "  PolygonGroup #" << j << ":\n";
			file << "    Outer boundary points: " << outerPoints.size() << "\n";
			file << "    Outer perimeter: " << pg.outer.getPerimeter() << " m\n";

			for (size_t k = 0; k < outerPoints.size(); ++k)
			{
				Point& p = outerPoints[k];

				file << std::fixed << std::setprecision(6);
				file << "      [" << k << "]\n";
				file << "         Original:  Lon=" << p.lon
					<< ", Lat=" << p.lat
					<< ", Alt=" << p.alt << "\n";

				file << std::fixed << std::setprecision(0);
				file << "         JTSK:      X=" << p.X
					<< ", Y=" << p.Y
					<< ", Z=" << p.Z << "\n";
			}

			if (!pg.inners.empty())
			{
				file << "    Inner boundaries: " << pg.inners.size() << "\n";

				for (size_t m = 0; m < pg.inners.size(); ++m)
				{
					auto& innerPoints = pg.inners[m].getPoints();

					file << "      Inner #" << m
						<< " points = " << innerPoints.size() << "\n";

					for (size_t k = 0; k < innerPoints.size(); ++k)
					{
						Point& p = innerPoints[k];

						file << std::fixed << std::setprecision(6);
						file << "        [" << k << "]\n";
						file << "           Original:  Lon=" << p.lon
							<< ", Lat=" << p.lat
							<< ", Alt=" << p.alt << "\n";

						file << std::fixed << std::setprecision(0);
						file << "           JTSK:      X=" << p.X
							<< ", Y=" << p.Y
							<< ", Z=" << p.Z << "\n";
					}
				}
			}

			file << "\n";
		}

		file << "\n";
	}

	file.close();
	tiles.close();
	std::cout << "Treba " << tile_counter << " suborov laz\n";


	int found = 0;
	int missing = 0;

	for (size_t i = 0; i < forests.size(); ++i)
	{
		Forest& forest = forests[i];
		
		auto& files = forest.getTiles();

		for (int j = 0; j < files.size(); j++)
		{
			fs::path fullPath = fs::path(lazDir) / files[j];
			if (fs::exists(fullPath))
			{
				found++;
			}
			else
			{
				std::cout << forest.getName() << "\n";
				std::cout << "[MISSING] " << files[j] << "\n";
				missing++;
			}
		}

		
	}
	std::cout << "\nSummary:\n";
	std::cout << "  Found:   " << found << "\n";
	std::cout << "  Missing: " << missing << "\n";


}

int main(int argc, char** argv)
{
	//main_all_files(argc, argv);
	main_curve(argc, argv);
}