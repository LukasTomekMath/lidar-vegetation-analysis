#include "FeatureVector.h"
#include "MetricDefinitions.h"
#include "Forest.h"
#include "RasterIO.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>



MetricStats computeStatsForBand(const std::vector<double>& band, Forest& forest, const RasterData& raster, int bandIdx)
{
	MetricStats s{};

	double noData = raster.noDataValues[bandIdx];
	const auto& mask = forest.getMask();

	int maskH = mask.size();
	if (maskH == 0) return s;
	int maskW = mask[0].size();

	double fMinX = forest.getMinX();
	double fMaxY = forest.getMaxY();


	double sum = 0.0;
	int count = 0;
	std::vector<double> values;

	for (int my = 0; my < maskH; my++) 
	{
		for (int mx = 0; mx < maskW; mx++) 
		{
			if (!mask[my][mx]) continue; //mimo masky

			//prepocet na svetove suradnice, resp. metre
			double worldX = fMinX + mx * raster.pixelSize;
			double worldY = fMaxY - my * raster.pixelSize;

			//prepocet na suradnice v rastri, teda ktory pixel v tife to je
			int rx = std::round((worldX - raster.xMin) / raster.pixelSize);
			int ry = std::round((raster.yMax - worldY) / raster.pixelSize);

			//ci to nahodou nevyletelo von
			if (rx < 0 || rx >= raster.nx || ry < 0 || ry >= raster.ny) continue;

			double v = band[ry * raster.nx + rx];

			if (!std::isfinite(v) || v == noData) continue;

			if (count == 0) { s.min = s.max = v; }
			if (v < s.min) s.min = v;
			if (v > s.max) s.max = v;

			sum += v;
			values.push_back(v); 
			count++;
		}
	}

	if (count > 0) {
		s.mean = sum / count;
		double var = 0.0;
		for (double v : values) {
			var += (v - s.mean) * (v - s.mean);
		}
		s.std = std::sqrt(var / count);
	}

	return s;
}


FeatureVector computeFeatureVector(const RasterData& raster, Forest& forest)
{
	FeatureVector fv;
	fv.forestName = forest.getNameWithoutDiacritics();
	fv.segment = forest.getSegment();
	fv.kodPrales = forest.getKodPrales();
	fv.metrics.resize(raster.nBands);
	fv.hodnotSeg = forest.getHodnotSeg();
	fv.hectares = forest.getHectares();
	fv.poznamka = (forest.getPoznamka().empty() ? "none" : "\"" + forest.getPoznamka() + "\"");

	for (int i = 0; i < raster.nBands;i++)
	{
		fv.metrics[i] = computeStatsForBand(raster.bands[i],forest,raster,i);
	}
	return fv;
}


void initializeFeatureVectorCSV(const std::string& csvPath)
{
    std::ofstream csv(csvPath, std::ios::trunc);
    if (!csv)
        throw std::runtime_error("Cannot create CSV: " + csvPath);

    // zapisanie hlavicky
    csv << "ForestName,Segment,KodPrales,HodnotSeg,Hectares,Poznamka";
    for (const auto &name : metrics::names)
    {
        csv << "," << name << "_mean"
            << "," << name << "_std"
            << "," << name << "_min"
            << "," << name << "_max";
    }
    csv << "\n";
    csv.close();
}

void exportFeatureVectorCSV(const FeatureVector& fv, const std::string& csvPath)
{
	std::ofstream csv(csvPath,std::ios::app);
	if (!csv.is_open()) {
		std::cerr << "Cannot open CSV: " << csvPath << "\n";
		return;
	}

	csv << fv.forestName << "," << fv.segment << "," << fv.kodPrales<<","<<fv.hodnotSeg<<","<<fv.hectares<<","<<fv.poznamka;
	for (auto& m : fv.metrics) {
		csv << "," << m.mean << "," << m.std << "," << m.min << "," << m.max;
	}
	csv << "\n";
	csv.flush();
}
