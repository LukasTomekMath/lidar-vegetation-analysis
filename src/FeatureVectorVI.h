#pragma once
#include <string>
#include <vector>

class Forest;
struct RasterData;

struct MetricStatsVI
{
	double mean;
	double std;
	double median;
};

struct FeatureVectorVI
{
	std::string forestName;
	std::string segment;      
	std::string kodPrales;
	std::string hodnotSeg;
	std::string poznamka;
	double hectares;
	std::vector<MetricStatsVI> metrics; 
};

MetricStatsVI computeStatsForBandVI(const std::vector<double>& band, const std::vector<Forest*>& parts, const RasterData& raster, int bandIdx, double maskRes);
FeatureVectorVI computeFeatureVectorVI(RasterData& raster, const std::vector<Forest*>& parts);

void initializeFeatureVectorVICSV(const std::string& csvPath);
void exportFeatureVectorVICSV(const FeatureVectorVI& fv, const std::string& csvPath);