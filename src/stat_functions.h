#pragma once

#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>

#define NO_DATA_VALUE -999999999.0

namespace StatFunctions
{
	double mean(const std::vector<double>& data);

	double percentile(std::vector<double> data, double p);

	double kurtosis(const std::vector<double>& data, const double mean, bool flag = true);

	double skewness(const std::vector<double>& data, const double mean, bool flag = true);

	double variance(const std::vector<double>& data, const double mean);

	double shannonIndex(const std::vector<double>& data);

	double pointsBelow(const std::vector<double>& data, double z2);

	double pointsAbove(const std::vector<double>& data, double z1);

	double pointsBetween(const std::vector<double>& data, double z1, double z2);

}