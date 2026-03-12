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
#include "libs/gdal_include/ogrsf_frmts.h"
#include "libs/gdal_include/cpl_conv.h"
#include "libs/gdal_include/ogr_srs_api.h"

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

//structura na precitany tif subor
struct RasterData
{
	int nx = 0;
	int ny = 0;
	int nBands = 0;
	std::vector<std::vector<double>> bands; // [band][y*nx + x]
	std::vector<double> noDataValues;
};

//len na zednodusenie prace s metrikami a co mam z nich ratat
struct MetricStats
{
	double mean;
	double std;
	double min;
	double max;
};

//samotny vysledny feature vector pre dany prales
struct FeatureVector
{
	std::string forestName;
	std::vector<MetricStats> metrics; 
};

void exportFeatureVectorCSV(const FeatureVector& fv, const std::string& csvPath)
{
	std::ofstream csv(csvPath,std::ios::app);
	if (!csv.is_open()) {
		std::cerr << "Cannot open CSV: " << csvPath << "\n";
		return;
	}

	csv << fv.forestName;
	for (auto& m : fv.metrics) {
		csv << "," << m.mean << "," << m.std << "," << m.min << "," << m.max;
	}
	csv << "\n";
	csv.flush();
}


//citanie tif subboru na spracovanie pomocou masky
RasterData readGeoTiff(const std::string& filename)
{
	RasterData out;

	GDALDataset* ds = (GDALDataset*)GDALOpen(filename.c_str(), GA_ReadOnly);
	if (!ds)
		throw std::runtime_error("Failed to open " + filename);

	out.nx = GDALGetRasterXSize(ds);
	out.ny = GDALGetRasterYSize(ds);
	out.nBands = GDALGetRasterCount(ds);

	out.bands.resize(out.nBands);
	out.noDataValues.resize(out.nBands);

	for (int b = 0; b < out.nBands; b++)
	{
		GDALRasterBandH band = GDALGetRasterBand(ds, b + 1);

		int hasNoData = 0;
		double noDataVal = GDALGetRasterNoDataValue(band, &hasNoData);
		out.noDataValues[b] = hasNoData ? noDataVal : std::numeric_limits<double>::quiet_NaN();

		out.bands[b].resize(out.nx * out.ny);

		GDALRasterIO(band, GF_Read,
			0, 0,
			out.nx, out.ny,
			out.bands[b].data(),
			out.nx, out.ny,
			GDT_Float64,
			0, 0);
	}

	GDALClose(ds);
	return out;
}

//testovacia funkcia
void applyMask(RasterData& r, const std::vector<std::vector<bool>>& mask, double noData)
{
	for (int y = 0; y < r.ny; y++)
	{
		for (int x = 0; x < r.nx; x++)
		{
			if (!mask[y][x])   // mimo lesa
			{
				int idx = y * r.nx + x;
				for (int b = 0; b < r.nBands; ++b)
				{
					r.bands[b][idx] = noData;
				}
			}
		}
	}
}
//testovacia funkcia
void writeMaskedTif(const std::string& outPath, const RasterData& r, GDALDataset* srcDS, double noData)
{
	GDALDriver* driver = GetGDALDriverManager()->GetDriverByName("GTiff");

	GDALDataset* out = driver->Create(outPath.c_str(), r.nx, r.ny, r.nBands, GDT_Float64, nullptr);

	// copy geotransform + projection
	double gt[6];
	srcDS->GetGeoTransform(gt);
	out->SetGeoTransform(gt);

	out->SetProjection(srcDS->GetProjectionRef());

	for (int b = 0; b < r.nBands; b++)
	{
		GDALRasterBand* band = out->GetRasterBand(b + 1);
		band->SetNoDataValue(noData);
		band->RasterIO(GF_Write, 0, 0, r.nx, r.ny,
			(void*)r.bands[b].data(), r.nx, r.ny,
			GDT_Float64, 0, 0);
	}

	GDALClose(out);
}
//testovacia funkcia
void debugMaskCheck(const std::string& inTif, const std::string& outTif, const std::vector<std::vector<bool>>& mask)
{
	GDALAllRegister();

	GDALDataset* src = (GDALDataset*)GDALOpen(inTif.c_str(), GA_ReadOnly);
	if (!src) throw std::runtime_error("Cannot open input TIF");

	RasterData r = readGeoTiff(inTif);

	double no_data = -9999.0;
	applyMask(r, mask, no_data);

	writeMaskedTif(outTif, r, src, no_data);

	GDALClose(src);
}

MetricStats computeStatsForBand(const std::vector<double>& band, const std::vector<std::vector<bool>>& mask, int nx, int ny,double noData)
{
	MetricStats s{};

	double sum = 0.0;
	double sumSq = 0.0;
	int count = 0;

	bool first = true;

	for (int y = 0; y < ny; y++)
	{
		for (int x = 0; x < nx; x++)
		{
			if (!mask[y][x]) continue; //mimo lesa

			int i = y * nx + x;
			double v = band[i];

			if (!std::isfinite(v) || v == noData|| std::isnan(v)) continue;

			if (first) {
				s.min = s.max = v;
				first = false;
			}

			if (v < s.min) s.min = v;
			if (v > s.max) s.max = v;

			sum += v;
			sumSq += v * v;
			count++;
		}
	}
	s.mean = sum / count;
	double var=0.0;

	for (int y = 0; y < ny; y++)
	{
		for (int x = 0; x < nx; x++)
		{
			if (!mask[y][x]) continue;

			int i = y * nx + x;
			double v = band[i];

			if (!std::isfinite(v) || v == noData || std::isnan(v)) continue;

			var += (v - s.mean) * (v - s.mean);
		}
	}
	var /= count;
	s.std = std::sqrt(var);

	return s;
}

FeatureVector computeFeatureVector(const RasterData& raster, const Forest& forest)
{
	FeatureVector fv;
	fv.forestName = forest.getName();
	fv.metrics.resize(raster.nBands);

	for (int i = 0; i < raster.nBands;i++)
	{
		fv.metrics[i] = computeStatsForBand(raster.bands[i],forest.getMask(),raster.nx,raster.ny,raster.noDataValues[i]);
	}
	return fv;
}


//povodne spracovanie vsetkych suborov (bez pralesov) + filtracia
int main_all_files(int argc, char** argv)
{
	int nprocs = 10;
	int files_done = 0;
	omp_set_num_threads(nprocs);

	int pixelsize = 10;

	std::cout << "Use case: Load all laz files in a specified directory and calculate tifs." << std::endl;
	std::cout << "Pixel size is set to: "<<pixelsize << std::endl;
	std::cout << "Number of threads: " << nprocs << std::endl;

	std::ofstream csv("processing_stats.csv", std::ios::trunc);
	csv << "Tile,File_MB,PointsInMesh,ReadTime_s,NormalizationTime_s,RedistributionTime_s,MetricsTime_s,ExportTime_s,DeallocationTime_s,OverallTime_s\n";
	csv.flush();

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
		startMsg << i + 1 << "/" << files.size() << " Processing file '" << lazFileName << "'\n" ;
#pragma omp critical
		{
			std::cout << startMsg.str() << std::flush;
		}

		handler->setupAreaInfo(files[i].entry.path().string(), upperLeftX, upperLeftY, upperLeftX+2000, upperLeftY-2000,pixelsize);

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
			csv.flush();
			std::cout << "\nDone file: " << lazFileName << "\n"
				<< "Progress: " << local_done << "/" << files.size()
				<< " (" << local_processedBytes / (1024.0 * 1024.0) << "MB / "
				<< totalBytes / (1024.0 * 1024.0) << "MB)\n"
				<< "Elapsed time: " << formatTime(progressTime) << "\n"
				<< "Estimated time left: " << formatTime(estimatedTimeLeft) << "\n\n" << std::flush;;
		}

		fileNameParts.clear();
	}
}

//testovanie nacitania kriviek a check ktore subory treba
int main_curve(int argc, char** argv)
{
	if (argc != 3)
	{
		std::cout << "Input is: (kml directory) (laz directory)\n";
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
	GDALAllRegister();

	for (auto& entry : fs::directory_iterator(kmlDir))
	{
		if (entry.is_regular_file() && entry.path().extension() == ".kml")
		{
			std::string filename = entry.path().string();
			std::cout << "Loading KML: " << filename << std::endl;

			//forestManager.loadFromKML(filename);
			forestManager.parseKML_GDAL(filename);
			fileCount++;
		}
	}

	auto& forests = forestManager.getForests();

	std::cout << "Finished reading " << fileCount << " KML files.\n";
	std::cout << "Total forests loaded: " << forests.size() << "\n\n";

	//output file
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

//spracovanie vsetkych suborov S krivkami pralesov + filtracia
int main_curve_files(int argc, char** argv)
{
	int nprocs = 10;
	int files_done = 0;
	omp_set_num_threads(nprocs);

	if (argc != 3)
	{
		std::cout << "Input is: (kml directory) (laz directory)\n";
		return -1;
	}

	std::cout << "Use case: Load laz files in a specified directory corespondin to kml files from a specified directory and calculate tifs." << std::endl;
	std::string kmlDir = argv[1];
	std::string lazDir = argv[2];


	if (!fs::exists(kmlDir))
	{
		std::cerr << "Directory does not exist: " << kmlDir << "\n";
		return -1;
	}
	if (!fs::exists(lazDir))
	{
		std::cerr << "Directory does not exist: " << lazDir << "\n";
		return -1;
	}


	std::ofstream csv("processing_stats.csv", std::ios::trunc);
	csv << "Tile,File_MB,PointsInMesh,ReadTime_s,NormalizationTime_s,RedistributionTime_s,MetricsTime_s,ExportTime_s,DeallocationTime_s,OverallTime_s\n";

	ForestManager forestManager;
	int fileCount = 0;

	for (auto& entry : fs::directory_iterator(kmlDir))
	{
		if (entry.is_regular_file() && entry.path().extension() == ".kml")
		{
			std::string filename = entry.path().string();
			std::cout << "Loading KML: " << filename << std::endl;

			//forestManager.loadFromKML(filename);
			forestManager.parseKML_GDAL(filename);
			fileCount++;
		}
	}
	auto& forests = forestManager.getForests();
	std::cout << "Finished reading " << fileCount << " KML files.\n";
	std::cout << "Total forests loaded: " << forests.size() << "\n\n";


	std::uint64_t totalBytes = 0;
	std::uint64_t processedBytes = 0;

	int num_files = 0;
	for (size_t i = 0; i < forests.size(); ++i)
	{
		Forest& forest = forests[i];

		//setup na konkretne subory
		forest.calculateForestArea();//tuto sa aj vymazavaju tie male inners lebo e to viazane na rozlohu
		forest.findBoundingBox();
		forest.findTiles();
		for (const auto& fe : forest.getTiles())
		{
			fs::path fullPath = fs::path(lazDir) / fe;
			if (fs::exists(fullPath)) {
				totalBytes += fs::file_size(fullPath);
				num_files++;
			}
		}
	}
	//lesy su nacitane, teraz by sa malo vsetko citat, cize cela laz sekcia kodu

	auto overall_start = std::chrono::high_resolution_clock::now();
	GDALAllRegister();
	OGRRegisterAll();

	const int chunk = 5;

	//pre kazdy les musim najst tiles ktore treba 
#pragma omp parallel for schedule(dynamic, chunk)
	for (int i = 0; i < forests.size(); ++i)
	{
		Forest& forest = forests[i];
		auto& files = forest.getTiles(); //toto su tie subory co chcem precitat pre dany les

		for (int j = 0; j < files.size(); j++)//prechadzanie cez jednotlive subory
		{
			fs::path fullPath = fs::path(lazDir) / files[j];//cela cest k suboru, tu chcem citat
			

			if (fs::exists(fullPath))
			{

				std::string lazFileName = files[j];
				std::string areaName = std::regex_replace(lazFileName, std::regex("(\\.laz)"), "");
				//std::string forestName = std::regex_replace(forest.getName(), std::regex(" "), "_");
				std::string forestName = forest.getName();

				std::string outputPath = "../../../filtrovane_metriky_5x5/" + areaName + "_" + forestName + ".tif";

				if (fs::exists(outputPath)) {
#pragma omp critical
					{
						files_done++;
						processedBytes += fs::file_size(fullPath);
						std::cout << "[SKIPPING] " << outputPath << " already exists.\n";
					}
					continue;
				}

				//co sa ma spravit ked sa konkretny subor najde- vypocita sa co sa ma
				auto file_start = std::chrono::high_resolution_clock::now(); //celkovy cas

				std::vector<std::string> fileNameParts; 
				double upperLeftX = -1.0; 
				double upperLeftY = -1.0; 
				double lowerRightX = -1.0;
				double lowerRightY = -1.0;

				DataHandler* handler = new DataHandler; 
				//std::string lazFileName = files[j];
				handler->setAreaName(std::regex_replace(lazFileName, std::regex("(\\.laz)"), "")); 
				handler->setForestName(forest.getName());

				splitString(handler->areaName(), fileNameParts); 

				//informacie o danej tile
				upperLeftX = std::stod(fileNameParts[1]); 
				upperLeftY = std::stod(fileNameParts[2]) + 2000.0; 
				lowerRightX = upperLeftX + 2000.0;
				lowerRightY = std::stod(fileNameParts[2]);

				double clippedMinX = std::max(forest.getMinX(), upperLeftX);	//upper left X
				double clippedMaxX = std::min(forest.getMaxX(), lowerRightX);	//lower right X
				double clippedMinY = std::max(forest.getMinY(), lowerRightY);	//lower right Y
				double clippedMaxY = std::min(forest.getMaxY(), upperLeftY);	//upper left Y


				if (clippedMinX >= clippedMaxX || clippedMinY >= clippedMaxY) {
					continue;
				}

				int pixelsize = forest.getPixelSize();


				std::ostringstream startMsg; 
				startMsg << i + 1 << "/" << num_files << " Processing file '" << lazFileName << "' for "<< forest.getName()<<"\n";
#pragma omp critical
				{
					std::cout << startMsg.str();
				}

				handler->setupAreaInfo(fullPath.string(), clippedMinX, clippedMaxY, clippedMaxX, clippedMinY, pixelsize);

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

					processedBytes += std::filesystem::file_size(fullPath);
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
						<< "Progress: " << local_done << "/" << num_files
						<< " (" << local_processedBytes / (1024.0 * 1024.0) << "MB / "
						<< totalBytes / (1024.0 * 1024.0) << "MB)\n"
						<< "Elapsed time: " << formatTime(progressTime) << "\n"
						<< "Estimated time left: " << formatTime(estimatedTimeLeft) << "\n\n";
				}

				fileNameParts.clear();

			}
			else
			{
				//co sa ma spravit ked sa nenajde- nic
				files_done++;
				std::cout << forest.getName() << "\n";
				std::cout << "[MISSING] " << files[j] << "\n";
			}
		}


	}

}

//nacitanie kriviek, maska, pocitanie feature vektorov
int main_features(int argc, char** argv)
{
	int nprocs = 10;
	int files_done = 0;
	omp_set_num_threads(nprocs);
	GDALAllRegister();

	if (argc != 3)
	{
		std::cout << "Input is: (kml directory) (tif directory)\n";
		return -1;
	}
	std::string kmlDir = argv[1];
	std::string tifDir = argv[2];
	std::string outDir = "../../../masky";
	std::cout << "Use case: Load tif files in a specified directory coresponding to kml files from a specified directory and calculate feature vectors." << std::endl;


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

			//forestManager.loadFromKML(filename);
			forestManager.parseKML_GDAL(filename);
			fileCount++;
		}
	}
	auto& forests = forestManager.getForests();
	std::cout << "Finished reading " << fileCount << " KML files.\n";
	std::cout << "Total forests loaded: " << forests.size() << "\n\n";

	std::cout << "Processin masks " << "\n\n";
	int num_files = 0;
	for (size_t i = 0; i < forests.size(); ++i)
	{
		Forest& forest = forests[i];

		//setup na konkretne subory
		forest.calculateForestArea();
		forest.findBoundingBox();
		forest.createMask();

		//std::string outFilename = outDir + "/mask_forest_" + std::regex_replace(forest.getName(), std::regex(" "), "_") + ".tif";
		std::string outFilename = outDir + "/mask_forest_" + forest.getName() + ".tif";
		forest.exportMaskToGeoTIFF(outFilename);

	}
	std::cout << "Masks done" << "\n\n";


	//tu zacinaju features
	std::vector<std::string> files;

	//vytvorenie vektora so subormi
	for (const auto& entry : fs::directory_iterator(tifDir))
	{
		if (!entry.is_regular_file()) continue;

		auto path = entry.path();
		if (path.extension() == ".tif")
		{
			files.push_back(path.string());
		}
	}
	

	std::string csvFile = "feature_vectors.csv";
	std::ofstream csv(csvFile, std::ios::trunc);

	//zapisanie headra
	csv << "ForestName";
	std::vector<std::string> metricNames = {
		"Hmax","Hmean","Hmedian","Hp25","Hp75","Hp95","PPR","DAM_z",
		"BR_below_1","BR_1_2","BR_2_3","BR_above_3","BR_3_4","BR_4_5",
		"BR_below_5","BR_5_20","BR_above_20","Coeff_var_z","Hkurt",
		"Hskew","Hstd","Hvar","Shannon"
	};
	for (auto& name : metricNames) {
		csv << "," << name << "_mean"
			<< "," << name << "_std"
			<< "," << name << "_min"
			<< "," << name << "_max";
	}
	csv << "\n";
	csv.flush();

	//iterovanie cez jednotlive subory tif
	for (const auto& tifPath : files)
	{
		RasterData raster = readGeoTiff(tifPath);

		auto stem = std::filesystem::path(tifPath).stem().string();

		Forest* forest = forestManager.getForestByFileName(stem);

		if (!forest) {
			std::cout << "No mask for forest: " << stem << "\n";
			continue;
		}/*
		else
		{
			std::cout << "Mask found for: " << stem << "\n";
			debugMaskCheck(tifPath,
				"../../../masky/TEST_" + stem + ".tif",
				forest->getMask());
		}*/

		//tu bude to pocitanie features
		FeatureVector fv = computeFeatureVector(raster, *forest);
		exportFeatureVectorCSV(fv, csvFile);
	}


}

//orezanie laz suborov ia na casti v bounding boxoch a ulozenie zvlast
int main_crop_laz(int argc, char** argv)
{
	int nprocs = 6;
	omp_set_num_threads(nprocs);

	if (argc != 3)
	{
		std::cout << "Input is: (kml directory) (laz directory)\n";
		return -1;
	}

	std::cout << "Use case: Crop laz files according to the boundin boxes defined by forests in kml files." << std::endl;


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

			//forestManager.loadFromKML(filename);
			forestManager.parseKML_GDAL(filename);
			fileCount++;
		}
	}

	auto& forests = forestManager.getForests();

	std::cout << "Finished reading " << fileCount << " KML files.\n";
	std::cout << "Total forests loaded: " << forests.size() << "\n\n";

	int tile_counter = 0;

#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < forests.size(); i++)
	{
		Forest& forest = forests[i];

		//setup na konkretne subory
		forest.calculateForestArea();
		forest.findBoundingBox();
		forest.findTiles();

		std::vector<std::string>& subory = forest.getTiles();
		tile_counter += subory.size();

		auto& polygons = forest.getPolygons();

		fs::path outputDir = "D:\\blahova\\pointclouds_forests";


		for (int j = 0; j < subory.size(); j++)
		{
			fs::path fullPath = fs::path(lazDir) / subory[j];
			if (fs::exists(fullPath))
			{
#pragma omp critical
				{
					std::cout << "Processing forest: " << forest.getName() << " and file "<<subory[j] << std::endl;
				}

				fs::path outputDir = "D:\\blahova\\pointclouds_forests";
				if (!fs::exists(outputDir)) {
					fs::create_directories(outputDir);
				}

				std::string originalName = fullPath.stem().string();
				std::string forestName = forest.getName();
				//fs::path outputPath = outputDir / (originalName + "_" + std::regex_replace(forest.getName(), std::regex(" "), "_") + ".laz");
				fs::path outputPath = outputDir / (originalName + "_" + forest.getName());

				std::string lasExe = "D:\\blahova\\git\\libs\\las2las64.exe";

				double x_min = std::min(forest.getMinX(), forest.getMaxX());
				double x_max = std::max(forest.getMinX(), forest.getMaxX());
				double y_min = std::min(forest.getMinY(), forest.getMaxY());
				double y_max = std::max(forest.getMinY(), forest.getMaxY());

				std::stringstream cmd;
				cmd << "\" \"" << lasExe << "\""
					<< " -i \"" << fullPath.string() << "\""
					<< " -keep_xy "
					<< std::fixed << std::setprecision(3)
					<< x_min << " " << y_min << " " << x_max << " " << y_max
					<< " -o \"" << outputPath.string() << "\" \"";

				std::string finalCmd = cmd.str();

				//std::cout << "DEBUG COMMAND: " << finalCmd << std::endl;

				int result = std::system(finalCmd.c_str());

				if (result == 0) 
				{
#pragma omp critical
					{
						std::cout << "Successfully saved to: " << outputPath.string() << "\n\n";
					}
				}
				else 
				{
#pragma omp critical
					{
						std::cerr << "Error: las2las64 returned code " << result << "\n\n";
					}
				}

			}

			else
			{
				std::cout << forest.getName() << "\n";
				std::cout << "[MISSING] " << subory[j] << "\n";
			}
		}
	}

}

int main(int argc, char** argv)
{
	const char* proj_paths[] = { "libs/dlls_to_copy", "./", nullptr };
	OSRSetPROJSearchPaths(proj_paths);

	main_all_files(argc, argv);
	
	//main_curve(argc, argv);

	//main_curve_files(argc, argv);

	//main_features(argc, argv);

	//main_crop_laz(argc, argv);
}