#pragma once

#ifdef min
#undef min
#endif

#ifdef max
#undef max
#endif

#include "../libs/gdal_include/gdal_priv.h"
#include "../libs/gdal_include/gdal.h"
#include "../libs/gdal_include/ogr_api.h"

#include "LASlib/lasreader.hpp"
#include "LASlib/laswriter.hpp"
// #include "laszip/laszip_api.h"

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <chrono>
#include <omp.h>
#include <vector>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <random>
#include <ctime>
#include <filesystem>

#include "stat_functions.h"

#define GROUND 2
#define VEG_LOW 3
#define VEG_MEDIUM 4
#define VEG_HIGH 5


struct OutputData {
	//casy spracovania jednotlivych ukonov
	double readTime = 0.0;
	double normalizationTime = 0.0;
	double redistributionTime = 0.0;
	double metricsTime = 0.0;
	double exportTime = 0.0;

	//info o subore do vystupu
	int nPointsInMesh = 0;
	std::uint64_t fileSizeBytes = 0;
};

struct AreaInfo
{
	// Grid size for normalization (should be 1x1m per pixel).
	int width_n = -1;
	int height_n = -1;
	// Grid size of the output raster (should be 10x10m per pixel).
	int width = -1;
	int height = -1;
	// Size of 1 pixel of the output raster (should be 10m).
	int desiredPixelSize = -1;
	// Coordinates of the upper left corner of the output raster.
	double xLeft = -1.0;
	double yLeft = -1.0;
	// LAZ file name from which to read the point cloud data.
	std::string lazFileName;
};

struct VegetationPoints
{
	std::vector<double> xCoords = {};
	std::vector<double> yCoords = {};
	std::vector<double> zCoords = {};
	std::vector<int> ptClass = {};
};

struct GroundPoints
{
	std::vector<double> xCoords = {};
	std::vector<double> yCoords = {};
	std::vector<double> zCoords = {};
};

struct PixelPoints
{
	VegetationPoints vegetationPts;
	GroundPoints groundPts;
};

class DataHandler
{
public:
	DataHandler();
	~DataHandler();


	void setupAreaInfo(const std::string &lazFileName, const double upperLeftX, const double upperLeftY, const double lowerRightX, const double lowerRightY, const int pixelSize);
	void reset();
	bool performCalculation();
	void setAreaName(const std::string &name) { m_areaName = name; }
	void setOutputDirectory(const std::string& directory) { m_outputDirectory = directory; }
	void setForestName(const std::string& name) { m_forestName = name; }//toto je na nazov lesa. v jednom tile moze byt viac lesov
	const std::string areaName() { return m_areaName;}
	std::string outputFileName() const;
	const OutputData& getOutputData() { return m_output; }

private:
	OGRSpatialReference m_sourceCRS; // CRS nacitane z aktualneho LAS/LAZ.
	OutputData m_output;
	AreaInfo m_areaInfo;
	std::vector<double> m_DTM_raw;	//na povodnom gride
	std::vector<double> m_DTM_final;	//na redistribuovanom
	std::vector<PixelPoints> m_meshPixels;
	std::vector<PixelPoints> m_meshPixelsRedistributed;
	std::vector<double *> m_metrics;
	std::string m_areaName = "";
	std::string m_outputDirectory = "raster_metrics";
	std::string m_forestName="";

	bool readLazFile();

	void filterPointsMAD();

	void filterPointsHp95();

	void filterPointsHp02();

	void cleanVegetationBelowGround();

	void normalizePoints();

	void redistributePoints();

	void computeMetrics();

	void exportMetrics();

	bool exportLAS(std::string fileName);

	bool exportLAZ(std::string fileName);
};
