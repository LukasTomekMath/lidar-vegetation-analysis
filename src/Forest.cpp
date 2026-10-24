#include "RasterIO.h"
#include <memory>
#include "Forest.h"
#include "TextUtils.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>

#include "libs/gdal_include/gdal_priv.h"
#include "libs/gdal_include/ogr_spatialref.h"

void Forest::setName(const std::string& n)
{
    const auto plainName = text_utils::withoutDiacritics(n);
    name = n;
    nameWithoutDiacritics = plainName;
}

void Forest::findBoundingBox()
{
    if (polygons.empty() || polygons[0].outer.getPoints().empty()) {
        minX = maxX = minY = maxY = 0;
        rawMinX = rawMaxX = rawMinY = rawMaxY = 0;
        return;
    }

    int pad = padding*pixelSize;

    //inicializacia na prvy bod, s tym budem najprv porovnavat
    minX = polygons[0].outer.getPoints()[0].X;
    maxX= polygons[0].outer.getPoints()[0].X;
    minY = polygons[0].outer.getPoints()[0].Y;
    maxY = polygons[0].outer.getPoints()[0].Y;

    for (auto& pg : polygons)
    {
        for (auto& p : pg.outer.getPoints()) {
            minX = std::min(minX, p.X);
            maxX = std::max(maxX, p.X);
            minY = std::min(minY, p.Y);
            maxY = std::max(maxY, p.Y);
        }
    }
    // Uloz presne hranice pred rozsirenim bounding boxu.
    rawMinX = minX;
    rawMinY = minY;
    rawMaxX = maxX;
    rawMaxY = maxY;

    // Zaokruhli smerom von na nasobky 10 m a pridaj pas so sirkou padding * pixelSize.
    minX = std::floor(minX / 10.0) * 10.0 -pad;
    minY = std::floor(minY / 10.0) * 10.0-pad;
    maxX = std::ceil(maxX / 10.0) * 10.0 +pad;
    maxY = std::ceil(maxY / 10.0) * 10.0+pad;
}

void Forest::findTiles()
{
    tiles.clear();
    // floor urci lavy dolny roh dlazdice aj pri zapornych suradniciach.
    int minTileX = static_cast<int>(std::floor(minX / 2000) * 2000);
    int maxTileX = static_cast<int>(std::floor(maxX / 2000) * 2000);
    int minTileY = static_cast<int>(std::floor(minY / 2000) * 2000);
    int maxTileY = static_cast<int>(std::floor(maxY / 2000) * 2000);

    for (int tx = minTileX; tx <= maxTileX; tx += 2000)
    {
        for (int ty = minTileY; ty <= maxTileY; ty += 2000)
        {
            // Pre kazdu dlazdicu pridaj tri varianty; samotny nazov lokality sa opakuje.
            tiles.push_back("tile_" + std::to_string(tx) + "_" + std::to_string(ty) + ".laz");
            tiles.push_back("tile_" + std::to_string(tx) + "_" + std::to_string(ty)  +"_" + nameWithoutDiacritics + ".laz");
            tiles.push_back(nameWithoutDiacritics + ".laz");
        }
    }
}

void Forest::calculateForestArea()
{
    int num = polygons.size();
    double sum_area = 0.0;

    for (int i = 0; i < num; i++)
    {
        ForestPolygon& pg = polygons[i];

        double sum_area_inner = 0.0;
        pg.outer.calculateCurveArea();
        if (!polygons[i].inners.empty())
        {
            int num_inners = pg.inners.size();
            for (int j = 0; j < num_inners; j++)
            {
                pg.inners[j].calculateCurveArea();
                double area = std::abs(pg.inners[j].getcurveArea());
                if (area < 1)
                {
                    // Po vymazani posun index spat, aby sa nepreskocila nasledujuca diera.
                    pg.inners.erase(pg.inners.begin() + j);
                    num_inners--;
                    j--;
                }
                else
                {
                    sum_area_inner += area;
                }
            }
        }
        sum_area += pg.outer.getcurveArea() - sum_area_inner;
    }
    forestArea = sum_area;
}

void Forest::rasteriseCurve(Curve& curve, std::vector<std::vector<bool>>& mask, double gridX0, double gridY0, int nx,
    int ny, double pixelSize, bool fillValue)
{
    // Scan-line algoritmus: priesecniky hran s riadkom zorad a vyplnaj po dvojiciach.

    auto& pts = curve.getPoints();
    int n = pts.size();

    for (int j = 0; j < ny; j++)
    {
        //y-ova suradnica, pre ktoru idem hladat priesecniky
        double y = gridY0 - (j + 0.5) * pixelSize;//pozicia stredu pixela v pociatku danej oblasti
        std::vector<double> xIntersections;

        for (int i = 0; i < n - 1; i++)
        {
            //hrana, pri ktorej prave zistujem, ci ma s mojim riadkom priesecnik
            Point& a = pts[i];
            Point& b = pts[i + 1];

            // Ostre nerovnosti vynechaju vodorovne hrany aj vrcholy presne na riadku.
            if ((a.Y > y && b.Y < y) || (a.Y < y && b.Y > y))
            {
                //upravenie parametrizacie usecky, hladam t, pre ktore je priesecnik, potom že ake je to x
                //x(t) = x1 + t (x2 - x1)
                //y(t) = y1 + t(y2 - y1)
                double t = (y - a.Y) / (b.Y - a.Y);
                double x = a.X + t * (b.X - a.X);
                xIntersections.push_back(x);
            }
        }

        if (xIntersections.size() < 2)
            continue;

        std::sort(xIntersections.begin(), xIntersections.end());

        for (int k = 0; k < xIntersections.size() - 1; k += 2)
        {
            //ak to ma na danom riadku viac priesecnikov, tak
            //x1----vnutri---x2---vonku---x3---vnutri---x4 atd., takze zaciatocny pre "dnu" je kazdy druhy (0,2,..)
            double x0 = xIntersections[k];
            double x1 = xIntersections[k + 1];

            //pozeram, ktory pixel v poradi je zaciatocny a koncovy pre "vnutri"
            //tie -0.5 zabezpecia, ze len pixel, ktoreho stred je vnutri, bude oznaceny ako vnutri
            int iStart = std::ceil((x0 - gridX0) / pixelSize - 0.5);
            int iEnd = std::floor((x1 - gridX0) / pixelSize - 0.5);

            // Orez rozsah stlpcov na rozmery masky.
            iStart = std::max(iStart, 0);
            iEnd = std::min(iEnd, nx - 1);


            for (int i = iStart; i <= iEnd; ++i)
            {
                mask[j][i] = fillValue;
            }
        }
    }
}

void Forest::createMask()
{
    int nx = std::ceil((maxX - minX) / pixelSize);
    int ny = std::ceil((maxY - minY) / pixelSize);

	//"zaciatok" gridu je v lavom hornom rohu, ale y ide dole, takze y0 je maxY
    double gridX0 = minX;
    double gridY0 = maxY;

    mask.assign(ny, std::vector<bool>(nx, false));

    for (int i = 0; i < polygons.size(); i++)
    {
        ForestPolygon& pg = polygons[i];

        // Najprv vypln vonkajsiu hranicu aktualneho polygonu.
        rasteriseCurve(pg.outer, mask, gridX0, gridY0, nx, ny, pixelSize, true);

        // Potom vymaz vsetky jeho diery.
        for (int inner = 0; inner < pg.inners.size(); inner++)
        {
            Curve& innerCurve = pg.inners[inner];
            rasteriseCurve(innerCurve, mask, gridX0, gridY0, nx, ny, pixelSize, false);
        }
    }
}

void Forest::exportMaskToGeoTIFF(const std::string& filename)
{    
    int ny = mask.size();
    int nx = mask.empty() ? 0 : mask[0].size();
    if (nx == 0 || ny == 0) return;

    GDALDriver* driver = GetGDALDriverManager()->GetDriverByName("GTiff");

    if (!driver)
    {
        std::cerr << "GTiff driver not available\n";
        return;
    }

    GDALDataset* ds = driver->Create(
        filename.c_str(),
        nx,
        ny,
        1,
        GDT_Byte,
        nullptr
    );

    if (!ds)
    {
        std::cerr << "Failed to create GeoTIFF\n";
        return;
    }


    std::unique_ptr<GDALDataset, decltype(&GDALClose)> datasetGuard(ds, GDALClose);

    // Sever je hore: zaciatok je v lavom hornom rohu a krok Y je zaporny.
    double geoTransform[6] = { 0 };
    geoTransform[0] = minX;       // X laveho horneho rohu
    geoTransform[1] = static_cast<double>(pixelSize);  // sirka pixela
    geoTransform[2] = 0.0;
    geoTransform[3] = maxY;       // Y laveho horneho rohu
    geoTransform[4] = 0.0;
    geoTransform[5] = static_cast<double>( - pixelSize); // riadky pribudaju smerom nadol

    if (ds->SetGeoTransform(geoTransform) != CE_None) {
        throw std::runtime_error("Cannot write mask grid: " + filename);
    }

    if (!crsWkt.empty() && ds->SetProjection(crsWkt.c_str()) != CE_None)
        throw std::runtime_error("Cannot write mask CRS: " + filename);

    GDALRasterBand* band = ds->GetRasterBand(1);
    std::vector<uint8_t> row(nx);

    // Boolovsku masku preved na bajty 0/1 a zapisuj po riadkoch.
    for (int y = 0; y < ny; ++y)
    {
        for (int x = 0; x < nx; ++x)
            row[x] = mask[y][x] ? 1 : 0;

        band->RasterIO(
            GF_Write,
            0, y,
            nx, 1,
            row.data(),
            nx, 1,
            GDT_Byte,
            0, 0
        );
    }

    if (ds->FlushCache() != CE_None)
        throw std::runtime_error("Cannot flush mask: " + filename);
}
