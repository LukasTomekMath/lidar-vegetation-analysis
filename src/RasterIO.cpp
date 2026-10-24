#include "RasterIO.h"
#include "LASlib/lasreader.hpp"
#include <filesystem>
#include <limits>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <memory>
#include <cstring>
#include "libs/gdal_include/gdal_utils.h"

#include "libs/gdal_include/gdal_priv.h"

//citanie tif suboru na spracovanie pomocou masky
RasterData readGeoTiff(const std::string& filename)
{
	RasterData out;

	GDALDataset* ds = (GDALDataset*)GDALOpen(filename.c_str(), GA_ReadOnly);
	if (!ds)
		throw std::runtime_error("Failed to open " + filename);

	std::unique_ptr<GDALDataset, decltype(&GDALClose)> guard(ds, GDALClose);
	if (const auto* srs = ds->GetSpatialRef()) out.crsWkt = spatial::toWkt(*srs);
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

	return out;
}

void mergeGeoTiffs(const std::vector<std::string>& inputs, const std::string& output)
{
    if (inputs.empty())
        throw std::runtime_error("No rasters to merge: " + output);
    using Dataset = std::unique_ptr<GDALDataset, decltype(&GDALClose)>;
    std::vector<Dataset> sources;
    std::vector<GDALDatasetH> handles;
    double reference[6]{};
    auto equal = [](double a, double b) { return std::abs(a - b) < 1e-8; };
    for (const auto& path : inputs)
    {
        if (path == output)
            throw std::runtime_error("Merge output is also an input: " + output);
        Dataset ds(static_cast<GDALDataset*>(GDALOpen(path.c_str(), GA_ReadOnly)), GDALClose);
        if (!ds)
            throw std::runtime_error("Cannot open merge input: " + path);
        double gt[6]{};
        if (ds->GetGeoTransform(gt) != CE_None || gt[1] <= 0 || gt[5] >= 0 || gt[2] != 0 || gt[4] != 0)
            throw std::runtime_error("Unsupported raster grid: " + path);
        if (sources.empty())
            std::copy(gt, gt + 6, reference);
        else
        {
            auto& first = sources.front();
            const auto* srs = ds->GetSpatialRef();
            const auto* firstSrs = first->GetSpatialRef();
            double dx = (gt[0] - reference[0]) / reference[1];
            double dy = (gt[3] - reference[3]) / reference[5];
            if (!equal(gt[1], reference[1]) || !equal(gt[5], reference[5]) ||
                !equal(dx, std::round(dx)) || !equal(dy, std::round(dy)) ||
                ds->GetRasterCount() != first->GetRasterCount() || (bool(srs) != bool(firstSrs)) || (srs && !spatial::sameCRS(*srs, *firstSrs)))
                throw std::runtime_error("Incompatible merge grid, CRS or bands: " + path);
            for (int b = 1; b <= ds->GetRasterCount(); ++b)
            {
                auto* a = first->GetRasterBand(b);
                auto* z = ds->GetRasterBand(b);
                int hasA, hasZ;
                double na = a->GetNoDataValue(&hasA), nz = z->GetNoDataValue(&hasZ);
                if (a->GetRasterDataType() != z->GetRasterDataType() ||
                    std::strcmp(a->GetDescription(), z->GetDescription()) != 0 || hasA != hasZ ||
                    (hasA && na != nz && !(std::isnan(na) && std::isnan(nz))))
                    throw std::runtime_error("Incompatible merge band metadata: " + path);
            }
        }
        handles.push_back(ds.get());
        sources.push_back(std::move(ds));
    }
    char strict[] = "-strict";
    // Ekvivalenciu CRS sme overili vyssie, vratane starsieho ESRI Krovak WKT.
    // BuildVRT nesmie opakovat prisne porovnanie tychto ekvivalentnych zapisov.
    char allowEquivalent[] = "-allow_projection_difference";
    char* vrtArgs[] = {strict, allowEquivalent, nullptr};
    std::unique_ptr<GDALBuildVRTOptions, decltype(&GDALBuildVRTOptionsFree)>
        vrtOptions(GDALBuildVRTOptionsNew(vrtArgs, nullptr), GDALBuildVRTOptionsFree);
    int usageError = 0;
    Dataset vrt(static_cast<GDALDataset*>(GDALBuildVRT("", static_cast<int>(handles.size()),
        handles.data(), nullptr, vrtOptions.get(), &usageError)), GDALClose);
    if (!vrt || usageError)
        throw std::runtime_error("Cannot build mosaic: " + output);
    // VRT nemusi preniest popisy pasiem; zachovaj nazvy metrik explicitne.
    for (int b = 1; b <= vrt->GetRasterCount(); ++b)
        vrt->GetRasterBand(b)->SetDescription(sources.front()->GetRasterBand(b)->GetDescription());
    char of[] = "-of", format[] = "GTiff", co[] = "-co", big[] = "BIGTIFF=IF_SAFER";
    char* args[] = {of, format, co, big, nullptr};
    std::unique_ptr<GDALTranslateOptions, decltype(&GDALTranslateOptionsFree)>
        options(GDALTranslateOptionsNew(args, nullptr), GDALTranslateOptionsFree);
    Dataset result(static_cast<GDALDataset*>(GDALTranslate(output.c_str(), vrt.get(), options.get(),
        &usageError)), GDALClose);
    if (!result || usageError || result->FlushCache() != CE_None)
        throw std::runtime_error("Cannot write mosaic: " + output);
}

// Interne pomocne diagnosticke funkcie, anonymny namespace
namespace
{
// Nastavi pixely mimo masky na NoData vo vsetkych pasmach rastra.
// Maska musi mat rozmery rastra a rovnake poradie riadkov a stlpcov.
void applyMask(RasterData &r, const std::vector<std::vector<bool>> &mask, double noData)
{
    for (int y = 0; y < r.ny; y++)
    {
        for (int x = 0; x < r.nx; x++)
        {
            if (!mask[y][x]) // mimo lesa
            {
                int idx = y * r.nx + x;
                for (int b = 0; b < r.nBands; ++b)
                {
                    r.bands[b][idx] = noData;
                }
            }
        }
    }
}
// Zapise diagnosticky GeoTIFF vo Float64 s georeferenciou zo zdrojoveho datasetu.
// Kazdemu pasmu nastavi zadanu hodnotu NoData.
void writeMaskedTif(const std::string &outPath, const RasterData &r, GDALDataset *srcDS, double noData)
{
    GDALDriver *driver = GetGDALDriverManager()->GetDriverByName("GTiff");

    if (!driver) throw std::runtime_error("GTiff driver is unavailable");
    GDALDataset *out = driver->Create(outPath.c_str(), r.nx, r.ny, r.nBands, GDT_Float64, nullptr);

    if (!out) throw std::runtime_error("Cannot create raster: " + outPath);
    std::unique_ptr<GDALDataset, decltype(&GDALClose)> guard(out, GDALClose);

    // Zachovaj polohu pixelov a suradnicovy system vstupneho TIFF.
    double gt[6];
    if (srcDS->GetGeoTransform(gt) != CE_None || out->SetGeoTransform(gt) != CE_None ||
        out->SetSpatialRef(srcDS->GetSpatialRef()) != CE_None)
        throw std::runtime_error("Cannot copy raster georeference: " + outPath);

    for (int b = 0; b < r.nBands; b++)
    {
        GDALRasterBand *band = out->GetRasterBand(b + 1);
        band->SetNoDataValue(noData);
        band->RasterIO(GF_Write, 0, 0, r.nx, r.ny, (void *)r.bands[b].data(), r.nx, r.ny, GDT_Float64, 0, 0);
    }

    if (out->FlushCache() != CE_None)
        throw std::runtime_error("Cannot flush raster: " + outPath);
}
} // koniec anonymneho namespace

// Vytvori maskovany TIFF na vizualnu kontrolu hranic a dier segmentu, napr. v QGIS.
void maskCheck(const std::string &inTif, const std::string &outTif, const std::vector<std::vector<bool>> &mask)
{
    GDALAllRegister();

    GDALDataset *src = (GDALDataset *)GDALOpen(inTif.c_str(), GA_ReadOnly);
    if (!src)
        throw std::runtime_error("Cannot open input TIF");

    std::unique_ptr<GDALDataset, decltype(&GDALClose)> guard(src, GDALClose);
    RasterData r = readGeoTiff(inTif);

    // Diagnosticky vystup pouziva pevnu hodnotu NoData
    double no_data = -9999.0;
    applyMask(r, mask, no_data);

    writeMaskedTif(outTif, r, src, no_data);


}


namespace spatial
{
// Starsi ESRI WKT zapisuje orientaciu Krovaka aj tromi nestandardnymi
// parametrami. Pri osiach East/North je tato konkretna trojica redundantna.
// Normalizujeme iba pracovnu kopiu; povodny CRS pre TIFF zostava zachovany.
OGRSpatialReference comparisonCRS(const OGRSpatialReference& input)
{
    OGRSpatialReference result(input);
    const char* projection = result.GetAttrValue("PROJECTION");
    OGRAxisOrientation x = OAO_Other, y = OAO_Other;
    result.GetAxis(nullptr, 0, &x);
    result.GetAxis(nullptr, 1, &y);
    if (!projection || !EQUAL(projection, "Krovak") || x != OAO_East || y != OAO_North ||
        result.GetProjParm("X_Scale", 0) != -1 || result.GetProjParm("Y_Scale", 0) != 1 ||
        result.GetProjParm("XY_Plane_Rotation", 0) != 90)
        return result;
    auto* node = result.GetAttrNode("PROJCS");
    if (!node) return result;
    for (int i = node->GetChildCount() - 1; i >= 0; --i)
    {
        const auto* child = node->GetChild(i);
        if (!EQUAL(child->GetValue(), "PARAMETER") || child->GetChildCount() < 2) continue;
        const auto* name = child->GetChild(0)->GetValue();
        if (EQUAL(name, "X_Scale") || EQUAL(name, "Y_Scale") || EQUAL(name, "XY_Plane_Rotation"))
            node->DestroyChild(i);
    }
    return result;
}

bool sameCRS(const OGRSpatialReference& a, const OGRSpatialReference& b)
{
    if (a.IsSame(&b)) return true;
    auto normalizedA = comparisonCRS(a), normalizedB = comparisonCRS(b);
    return normalizedA.IsSame(&normalizedB);
}

std::string toWkt(const OGRSpatialReference& srs)
{
    if (srs.IsEmpty()) return {};
    char* text = nullptr;
    const char* options[] = {"FORMAT=WKT2_2019", nullptr};
    if (srs.exportToWkt(&text, options) != OGRERR_NONE)
        throw std::runtime_error("Cannot serialize CRS");
    std::string result(text);
    CPLFree(text);
    return result;
}

OGRSpatialReference readLasCRS(const LASheader& header, const std::string& source)
{
    OGRSpatialReference srs;
    if (header.vlr_geo_ogc_wkt && *header.vlr_geo_ogc_wkt)
    {
        // LASlib spristupnuje WKT z VLR aj EVLR. Prenes ho bez uprav textu.
        if (srs.importFromWkt(header.vlr_geo_ogc_wkt) != OGRERR_NONE)
            throw std::runtime_error("Invalid LAS WKT: " + source);
    }
    else if (header.vlr_geo_keys)
    {
        if (header.vlr_geo_keys->number_of_keys && !header.vlr_geo_key_entries)
            throw std::runtime_error("Invalid LAS GeoKeys: " + source);
        // Starsie LAS mozu namiesto WKT obsahovat priamo EPSG kody.
        int projected = 0, geographic = 0, vertical = 0;
        for (unsigned i = 0; i < header.vlr_geo_keys->number_of_keys; ++i)
        {
            const auto& key = header.vlr_geo_key_entries[i];
            if (key.key_id != 3072 && key.key_id != 2048 && key.key_id != 4096) continue;
            if (key.tiff_tag_location || key.count != 1)
                throw std::runtime_error("LAS custom GeoKeys require WKT: " + source);
            if (key.key_id == 3072) projected = key.value_offset;
            if (key.key_id == 2048) geographic = key.value_offset;
            if (key.key_id == 4096) vertical = key.value_offset;
        }
        const int horizontal = projected ? projected : geographic;
        if (!horizontal || horizontal == 32767 || vertical == 32767)
            throw std::runtime_error("LAS custom GeoKeys require WKT: " + source);
        if (srs.importFromEPSG(horizontal) != OGRERR_NONE)
            throw std::runtime_error("Invalid LAS EPSG code: " + source);
        if (vertical)
        {
            OGRSpatialReference height, compound;
            if (height.importFromEPSG(vertical) != OGRERR_NONE ||
                compound.SetCompoundCS("LAS CRS", &srs, &height) != OGRERR_NONE)
                throw std::runtime_error("Invalid LAS vertical CRS: " + source);
            srs = compound;
        }
    }
    srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    return srs;
}

OGRSpatialReference readLasCRS(const std::string& path)
{
    LASreadOpener opener;
    opener.set_file_name(path.c_str());
    std::unique_ptr<LASreader> reader(opener.open());
    if (!reader) throw std::runtime_error("Cannot open LAS/LAZ: " + path);
    auto srs = readLasCRS(reader->header, path);
    reader->close();
    return srs;
}

// KML vyber dlaždic pracuje s jednou spolocnou sietou v priecinku.
OGRSpatialReference directoryLasCRS(const std::string& directory)
{
    OGRSpatialReference result;
    bool found = false;
    for (const auto& entry : std::filesystem::directory_iterator(directory))
    {
        if (!entry.is_regular_file() || entry.path().extension() != ".laz" ||
            entry.path().filename().string().find(".copc.laz") != std::string::npos) continue;
        auto srs = readLasCRS(entry.path().string());
        if (srs.IsEmpty()) throw std::runtime_error("KML selection needs a known LAS CRS: " + entry.path().string());
        if (found && !sameCRS(result, srs))
            throw std::runtime_error("KML selection found incompatible CRS in: " + entry.path().string());
        result = comparisonCRS(srs);
        found = true;
    }
    if (!found) throw std::runtime_error("No LAZ inputs for KML selection: " + directory);
    return result;
}
}
