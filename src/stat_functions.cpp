#include "stat_functions.h"

double StatFunctions::mean(const std::vector<double>& data)
{
    if (data.empty())
        return NO_DATA_VALUE;

    double sum = std::accumulate(data.begin(), data.end(), 0.0);
    double mean = sum / data.size();
    
    return mean;
}

double StatFunctions::percentile(std::vector<double> data, double p)
{
    if (data.empty())
        return NO_DATA_VALUE;

    std::sort(data.begin(), data.end()); // sort data in ascending order
    int n = data.size();
    p = p / 100;
    double k_p = 0.0;

    int j = static_cast<int>(n * p);

    if (std::abs(static_cast<double>(j) - n * p) < DBL_EPSILON)
        k_p = (data[j - 1] + data[j]) / 2.0;
    else
        k_p = data[j];

    return k_p;
}

double StatFunctions::variance(const std::vector<double>& data, const double mean)
{
    if (data.size() < 2)
        return NO_DATA_VALUE;

    double sum = 0.0;
    for (double val : data)
        sum += (val - mean) * (val - mean);

    return sum / (data.size() - 1);
}

double StatFunctions::kurtosis(const std::vector<double>& data, const double mean, bool flag)
{
    if (!flag && data.size() <= 3)
        return NO_DATA_VALUE;

    double sumUpper = 0.0;
    double sumLower = 0.0;
    double n = static_cast<double>(data.size());
    double temp = 0.0;

    for (double val : data)
    {
        temp = val - mean;
        sumUpper += temp * temp * temp * temp;
        sumLower += temp * temp;
    }

    double k1 = (sumUpper / n) / ((sumLower / n) * (sumLower / n));

    if (flag)
        return k1;
    else
    {
        double k0 = ((n - 1.0) / ((n - 2.0) * (n - 3.0))) * ((n + 1.0) * k1 - 3.0 * (n - 1.0)) + 3.0;
        return k0;
    }
}

double StatFunctions::skewness(const std::vector<double>& data, const double mean, bool flag)
{
    if (!flag && data.size() <= 2)
        return NO_DATA_VALUE;

    double sumUpper = 0.0;
    double sumLower = 0.0;
    double n = static_cast<double>(data.size());
    double temp = 0.0;

    for (double val : data)
    {
        temp = val - mean;
        sumUpper += temp * temp * temp;
        sumLower += temp * temp;
    }

    temp = std::sqrt(sumLower / n);
    double s1 = (sumUpper / n) / (temp * temp * temp);

    if (flag)
        return s1;
    else
    {
        double s0 = ((std::sqrt(n * (n - 1.0))) / (n - 2.0)) * s1;
        return s0;
    }
}

double StatFunctions::shannonIndex(const std::vector<double>& data)
{
    // skusit aj Simpson index
    if (data.empty())
        return NO_DATA_VALUE;

    double shannon = 0.0;
    double bound1 = 0.0;
    double bound2 = 0.0;
    double sum_p = 0.0;
    double p_i = 0.0;
    double n = static_cast<double>(data.size());
    int nPts = 0;
    double maxZ = *std::max_element(data.begin(), data.end());

    for (double i = 0; i < maxZ; i += 0.5)
    {
        bound1 = i;
        bound2 = i + 0.5;
    
        nPts = std::count_if(data.begin(), data.end(), [bound1, bound2](double pt) {return (pt >= bound1) && (pt < bound2); });

        p_i = static_cast<double>(nPts) / n;

        if (p_i > 0.0)
            shannon -= p_i * std::log(p_i);
    }

    return shannon;
}

double StatFunctions::pointsBelow(const std::vector<double>& data, double z2)
{
    if (data.size() == 1 && data[0] == 0.0)
        return 0.0;

    int pts = std::count_if(data.begin(), data.end(), 
        [z2](double pt) {return (std::abs(pt - z2) < DBL_EPSILON || pt < z2); });
    
    return static_cast<double>(pts);
}

double StatFunctions::pointsAbove(const std::vector<double>& data, double z1)
{
    int pts = std::count_if(data.begin(), data.end(),
        [z1](double pt) {return (std::abs(pt - z1) < DBL_EPSILON || pt > z1); });

    return static_cast<double>(pts);
}

double StatFunctions::pointsBetween(const std::vector<double>& data, double z1, double z2)
{
    int pts = std::count_if(data.begin(), data.end(),
        [z1,z2](double pt) {return (std::abs(pt - z1) < DBL_EPSILON || pt > z1) && (std::abs(pt - z2) < DBL_EPSILON || pt < z2); });

    return static_cast<double>(pts);
}
