#include "RasterIO.h"
#include <limits>
#include <stdexcept>

#include "libs/gdal_include/gdal_priv.h"

//citanie tif suboru na spracovanie pomocou masky
RasterData readGeoTiff(const std::string& filename)
{
	RasterData out;

	GDALDataset* ds = (GDALDataset*)GDALOpen(filename.c_str(), GA_ReadOnly);
	if (!ds)
		throw std::runtime_error("Failed to open " + filename);

	out.nx = GDALGetRasterXSize(ds);
	out.ny = GDALGetRasterYSize(ds);
	out.nBands = GDALGetRasterCount(ds);

	double adfGeoTransform[6];
	if (ds->GetGeoTransform(adfGeoTransform) == CE_None) {
		out.xMin = adfGeoTransform[0]; 
		out.yMax = adfGeoTransform[3]; 
		out.pixelSize = adfGeoTransform[1]; 
	}

	out.bands.resize(out.nBands);
	out.noDataValues.resize(out.nBands);


	for (int b = 0; b < out.nBands; b++)
	{
		GDALRasterBandH band = GDALGetRasterBand(ds, b + 1);

		int hasNoData = 0;
		double noDataVal = GDALGetRasterNoDataValue(band, &hasNoData);

		out.noDataValues[b] = hasNoData ? noDataVal : std::numeric_limits<double>::quiet_NaN();

		out.bands[b].resize(out.nx * out.ny);

		GDALRasterIO(band, GF_Read,
			0, 0,
			out.nx, out.ny,
			out.bands[b].data(),
			out.nx, out.ny,
			GDT_Float64,
			0, 0);
	}

	GDALClose(ds);
	return out;
}

