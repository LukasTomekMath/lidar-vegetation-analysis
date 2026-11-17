#include "FeatureVector.h"
#include "MetricDefinitions.h"
#include "Forest.h"
#include "RasterIO.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include <map>
#include <unordered_set>
#include <algorithm>

//#define _DEBUG

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

MetricStats computeStatsForBandSat(const std::vector<double>& band, const std::vector<Forest*>& parts, const RasterData& raster, int bandIdx, double maskRes)
{
    MetricStats s{};
    const double noData = raster.noDataValues[bandIdx];

    std::vector<double> values;
    std::unordered_set<long long> seen;   // raster pixels already counted (overlapping parts)
    double sum = 0.0;

    for (Forest* fp : parts)
    {
        const auto& mask = fp->getMask();
        if (mask.empty()) continue;

        const int maskH = mask.size(), maskW = mask[0].size();
        const double x0 = fp->getMinX(), y0 = fp->getMaxY();

        // raster pixel range covering the mask bbox
        int c0 = (int)std::floor((x0 - raster.xMin) / raster.pixelSize);
        int c1 = (int)std::floor((x0 + maskW * maskRes - raster.xMin) / raster.pixelSize);
        int r0 = (int)std::floor((raster.yMax - y0) / raster.pixelSize);
        int r1 = (int)std::floor((raster.yMax - (y0 - maskH * maskRes)) / raster.pixelSize);

        for (int ry = std::max(r0, 0); ry <= std::min(r1, raster.ny - 1); ry++)
            for (int rx = std::max(c0, 0); rx <= std::min(c1, raster.nx - 1); rx++)
            {
                double wx = raster.xMin + (rx + 0.5) * raster.pixelSize;  // pixel center
                double wy = raster.yMax - (ry + 0.5) * raster.pixelSize;

                int mx = (int)std::floor((wx - x0) / maskRes);
                int my = (int)std::floor((y0 - wy) / maskRes);
                if (mx < 0 || mx >= maskW || my < 0 || my >= maskH) continue;
                if (!mask[my][mx]) continue;

                const long long idx = (long long)ry * raster.nx + rx;
                if (!seen.insert(idx).second) continue;   // already counted by another part

                double v = band[idx];
                if (!std::isfinite(v) || v == noData) continue;

                if (values.empty()) s.min = s.max = v;
                s.min = std::min(s.min, v);
                s.max = std::max(s.max, v);
                sum += v;
                values.push_back(v);
            }
    }

    if (!values.empty()) {
        s.mean = sum / values.size();
        double var = 0.0;
        for (double v : values) var += (v - s.mean) * (v - s.mean);
        s.std = std::sqrt(var / values.size());
    }
    return s;
}

FeatureVectorSat computeFeatureVector(const RasterData& r10, const RasterData& r20, const std::vector<Forest*>& parts)
{
    FeatureVectorSat fv;
    const Forest& first = *parts.front();

    fv.forestName = first.getNameWithoutDiacritics();
    fv.segment    = first.getSegment();
    fv.kodPrales  = first.getKodPrales();
    fv.hodnotSeg  = first.getHodnotSeg();
    fv.poznamka   = first.getPoznamka().empty()
                        ? "none"
                        : "\"" + first.getPoznamka() + "\"";

    fv.hectares = 0.0;
    for (Forest* f : parts) fv.hectares += f->getHectares();

    fv.metrics.reserve(r10.nBands + r20.nBands);

    // 10 m bands (e.g. B02, B03, B04, B08)
    for (int i = 0; i < r10.nBands; i++)
        fv.metrics.push_back(computeStatsForBandSat(r10.bands[i], parts, r10, i, 2.0));

    // 20 m bands (e.g. B05, B06, B07, B8A, B11, B12)
    for (int i = 0; i < r20.nBands; i++)
        fv.metrics.push_back(computeStatsForBandSat(r20.bands[i], parts, r20, i, 2.0));

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

void exportFeatureVectorCSV(const FeatureVectorSat& fv, const std::string& csvPath)
{
    std::ofstream csv(csvPath, std::ios::app);
    if (!csv.is_open()) {
        std::cerr << "Cannot open CSV: " << csvPath << "\n";
        return;
    }

    csv << fv.forestName << "," << fv.segment << "," << fv.kodPrales << ","
        << fv.hodnotSeg << "," << fv.hectares << "," << fv.poznamka;

    for (const auto& m : fv.metrics)
        csv << "," << m.mean << "," << m.std << "," << m.min << "," << m.max;

    csv << "\n";
}

void initializeFeatureVectorSatCSV(const std::string& csvPath)
{
    std::ofstream csv(csvPath, std::ios::trunc);
    if (!csv)
        throw std::runtime_error("Cannot create CSV: " + csvPath);

    // zapisanie hlavicky
	csv << "ForestName,Segment,KodPrales,HodnotSeg,Hectares,Poznamka";

    auto writeBandHeader = [&](const auto& names) {
		for (const auto& name : names)
			csv << "," << name << "_mean"
				<< "," << name << "_std"
				<< "," << name << "_min"
				<< "," << name << "_max";
	};

	writeBandHeader(sentinel::names10);
	writeBandHeader(sentinel::names20);

    csv << "\n";
    csv.close();
}