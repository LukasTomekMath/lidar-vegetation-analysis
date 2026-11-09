#pragma once
#include <string>
#include <vector>

class Forest;
struct RasterData;

//len na zjednodusenie prace s metrikami a co mam z nich ratat
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
	std::string segment;      
	std::string kodPrales;
	std::string hodnotSeg;
	std::string poznamka;
	double hectares;
	std::vector<MetricStats> metrics; 
};

struct FeatureVectorSat
{
	std::string forestName;
	std::string segment;      
	std::string kodPrales;
	std::string hodnotSeg;
	std::string poznamka;
	double hectares;
	std::vector<MetricStats> metrics; 
};

MetricStats computeStatsForBand(const std::vector<double>& band, Forest& forest, const RasterData& raster, int bandIdx);
FeatureVector computeFeatureVector(const RasterData& raster, Forest& forest);

MetricStats computeStatsForBandSat(const std::vector<double>& band, Forest& forest, const RasterData& raster, int bandIdx, double maskRes);
FeatureVectorSat computeFeatureVector(const RasterData& r10, const RasterData& r20, Forest& forest);

// Vytvori alebo prepise CSV a zapise hlavicku; riadky sa potom pridavaju cez exportFeatureVectorCSV.
// Pri neuspesnom otvoreni suboru vyhodi vynimku.
void initializeFeatureVectorCSV(const std::string& csvPath);
void exportFeatureVectorCSV(const FeatureVector& fv, const std::string& csvPath);

void initializeFeatureVectorSatCSV(const std::string& csvPath);
void exportFeatureVectorCSV(const FeatureVectorSat& fv, const std::string& csvPath);
