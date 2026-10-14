#pragma once
#include <string>
#include <vector>

//structura na precitany tif subor
struct RasterData
{
	int nx = 0;
	int ny = 0;
	int nBands = 0;
	std::vector<std::vector<double>> bands; // [band][y*nx + x]
	std::vector<double> noDataValues;

	double xMin = 0.0;     //vlavo hore
	double yMax = 0.0;     //vlavo hore
	double pixelSize = 0.0;
};

RasterData readGeoTiff(const std::string& filename);
