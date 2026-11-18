#pragma once
#include <string>
#include <vector>

class Forest;
struct RasterData;

struct MetricStatsSat
{
	double mean;
	double std;
	double median;
};

struct FeatureVectorSat
{
	std::string forestName;
	std::string segment;      
	std::string kodPrales;
	std::string hodnotSeg;
	std::string poznamka;
	double hectares;
	std::vector<MetricStatsSat> metrics; 
};

MetricStatsSat computeStatsForBandSat(const std::vector<double>& band, const std::vector<Forest*>& parts, const RasterData& raster, int bandIdx, double maskRes);
FeatureVectorSat computeFeatureVectorSat(const RasterData& r10, const RasterData& r20, const std::vector<Forest*>& parts);

void initializeFeatureVectorSatCSV(const std::string& csvPath);
void exportFeatureVectorSatCSV(const FeatureVectorSat& fv, const std::string& csvPath);