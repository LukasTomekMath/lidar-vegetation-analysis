#pragma once
#include <string>
#include <vector>
#include "ogr_spatialref.h"

class LASheader;

// CRS sa cita zo vstupu a zapisuje do TIFF; bez oprav WKT a bez reprojekcie bodov.
namespace spatial
{
OGRSpatialReference readLasCRS(const LASheader& header, const std::string& source);
OGRSpatialReference readLasCRS(const std::string& path);
OGRSpatialReference directoryLasCRS(const std::string& directory);
OGRSpatialReference comparisonCRS(const OGRSpatialReference& input);
bool sameCRS(const OGRSpatialReference& a, const OGRSpatialReference& b);
std::string toWkt(const OGRSpatialReference& srs);
}

//structura na precitany tif subor
struct RasterData
{
	std::string crsWkt; // Povodne CRS vstupu, prazdne ak chyba.
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

// merge rastrov s rovnakym CRS, pasmami a rozlisenim, bez interpolacie.
void mergeGeoTiffs(const std::vector<std::string>& inputs, const std::string& output);

// Diagnosticka pomocna funkcia: ulozi kopiu vstupneho TIFF s pixelmi mimo masky nastavenymi na NoData (-9999).
// Maska musi zodpovedat rozmerom aj poradiu pixelov vstupneho rastra.
// outTif musi byt ina cesta nez inTif; vystup je urceny na vizualnu kontrolu v QGIS.
// Pri beznych vypoctoch sa tato funkcia nikde nevola.
void maskCheck(const std::string &inTif, const std::string &outTif, const std::vector<std::vector<bool>> &mask);
