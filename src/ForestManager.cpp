#include "ForestManager.h"
#include "RasterIO.h"
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <regex>
#include <memory>

#include "libs/gdal_include/ogrsf_frmts.h"

void ForestManager::loadDirectory(const std::string& directory)
{
    int fileCount = 0;
    for (const auto& entry : std::filesystem::directory_iterator(directory))
    {
        std::filesystem::path ext = entry.path().extension();
        if (entry.is_regular_file() && (ext == ".kml" || ext == ".shp"))
        {
            const std::string filename = entry.path().string();
            std::cout << "Loading " << ext << " :" << filename << std::endl;

            if (!parseVector_GDAL(filename))
                throw std::runtime_error("Failed to load KML: " + filename);
            ++fileCount;
        }
    }
    std::cout << "Finished reading " << fileCount << " files.\n";
    std::cout << "Total forests loaded: " << forests.size() << "\n\n";
}

Forest* ForestManager::getForestByFileName(const std::string& fileName)
{
    std::string modifiedName = fileName;
    //std::replace(modifiedName.begin(), modifiedName.end(), '_', ' ');

    for (auto& forest : forests)
    {
        if (forest.getName() == modifiedName)
            return &forest;
    }
    return nullptr;
}

bool ForestManager::parseVector_GDAL(const std::string& filename)
{
    GDALDataset* ds = (GDALDataset*)GDALOpenEx(
        filename.c_str(),
        GDAL_OF_VECTOR,
        NULL, NULL, NULL);
    if (!ds)
    {
        std::cout << "Failed to open dataset: " << filename << "\n";
        return false;
    }
    std::unique_ptr<GDALDataset, decltype(&GDALClose)> datasetGuard(ds, GDALClose);

    OGRSpatialReference dst(targetCRS);
    dst.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    const std::string dstWkt = spatial::toWkt(dst);

    for (int i = 0; i < ds->GetLayerCount(); i++)
    {
        OGRLayer* layer = ds->GetLayer(i);
        layer->ResetReading();

        // CRS vrstvy: SHP ho ma v .prj, KML je vzdy WGS84.
        OGRSpatialReference src;
        if (const OGRSpatialReference* layerSrs = layer->GetSpatialRef())
        {
            src = *layerSrs;
        }
        else
        {
            std::cout << "Layer has no CRS (missing .prj?), assuming EPSG:4326\n";
            if (src.importFromEPSG(4326) != OGRERR_NONE)
                throw std::runtime_error("Cannot load EPSG:4326");
        }
        src.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

        OGRCoordinateTransformation* tr = OGRCreateCoordinateTransformation(&src, &dst);
        if (!tr)
            throw std::runtime_error("Cannot create coordinate transformation for " + filename);
        std::unique_ptr<OGRCoordinateTransformation, decltype(&OCTDestroyCoordinateTransformation)>
            transformGuard(tr, OCTDestroyCoordinateTransformation);

        OGRFeature* feature;
        while ((feature = layer->GetNextFeature()) != nullptr)
        {
            std::unique_ptr<OGRFeature, decltype(&OGRFeature::DestroyFeature)>
                featureGuard(feature, OGRFeature::DestroyFeature);

            Forest forest;
            forest.setCRS(dstWkt);

            if (feature->GetFieldIndex("Name") >= 0)
            {
                std::string name = feature->GetFieldAsString("Name");
                if (!name.empty())
                    forest.setName(std::regex_replace(name, std::regex(" "), "_"));
            }

            for (int f = 0; f < feature->GetFieldCount(); f++)
            {
                std::string fieldName = feature->GetFieldDefnRef(f)->GetNameRef();

                if (fieldName == "NAZOV")
                {
                    std::string name = feature->GetFieldAsString(f);
                    if (!name.empty())
                        forest.setName(std::regex_replace(name, std::regex(" "), "_"));
                }
                else if (fieldName == "HODNOT_SEG") forest.setHodnotSeg(feature->GetFieldAsString(f));
                else if (fieldName == "HECTARES")   forest.setHectares(feature->GetFieldAsDouble(f));
                else if (fieldName == "SEGMENT")    forest.setSegment(feature->GetFieldAsString(f));
                else if (fieldName == "KOD_PRALES") forest.setKodPrales(feature->GetFieldAsString(f));
                else if (fieldName == "POZNAMKA")   forest.setPoznamka(feature->GetFieldAsString(f));
            }

            OGRGeometry* geom = feature->GetGeometryRef();
            if (!geom) continue;

            processPolygon(geom, forest, tr);
            forests.push_back(forest);
        }
    }

    std::cout << "Loaded " << forests.size() << " forests (total so far)\n";
    return true;
}

void ForestManager::processPolygon(OGRGeometry* geom, Forest& forest, OGRCoordinateTransformation* tr)
{
    if (!geom)
        return;

    OGRwkbGeometryType type = wkbFlatten(geom->getGeometryType());

    //std::cout << "Geometry type: " << geom->getGeometryName() << std::endl;

    if (type == wkbPolygon)
    {
        OGRPolygon* poly = (OGRPolygon*)geom;

        ForestPolygon pg;

        OGRLinearRing* outer = poly->getExteriorRing();

        if (!outer)
        {
            std::cout << "Polygon has NO exterior ring\n";
            return;
        }


        Curve outerCurve;

        // Cyklus cez vsetky body outer krivky
        for (int i = 0; i < outer->getNumPoints(); i++)
        {
            Point p;

            p.lon = outer->getX(i);
            p.lat = outer->getY(i);
            p.alt = outer->getZ(i);

            double x = p.lon;
            double y = p.lat;
            double z = p.alt;

            // Transformacia lon, lat, alt na X,Y,Z
            if (!OCTTransform(tr, 1, &x, &y, &z))
                throw std::runtime_error("KML coordinate transformation failed");

            p.X = x;
            p.Y = y;
            p.Z = z;

            outerCurve.addPoint(p);
        }

        // Vypocet dlzky krivky
        outerCurve.calculateEdgeLengths();
        outerCurve.calculatePerimeter();

        pg.outer = outerCurve;

        int innerCount = poly->getNumInteriorRings();


        // Diery (inner krivky) sa spracuju rovnako; male diery sa odstranuju az pri vypocte plochy Forest.
        // cyklus cez vsetky inner krivky
        for (int r = 0; r < innerCount; r++)
        {
            OGRLinearRing* ring = poly->getInteriorRing(r);

            if (!ring)
                continue;

            Curve inner;

            // Cyklus cez vsetky body inner krivky
            for (int i = 0; i < ring->getNumPoints(); i++)
            {
                Point p;

                p.lon = ring->getX(i);
                p.lat = ring->getY(i);
                p.alt = ring->getZ(i);

                double x = p.lon;
                double y = p.lat;
                double z = p.alt;

                // Transformacia lon, lat, alt na X,Y,Z
                if (!OCTTransform(tr, 1, &x, &y, &z))
                throw std::runtime_error("KML coordinate transformation failed");

                p.X = x;
                p.Y = y;
                p.Z = z;

                inner.addPoint(p);
            }
            // Vypocet dlzky krivky
            inner.calculateEdgeLengths();
            inner.calculatePerimeter();

            pg.inners.push_back(inner);
        }

        forest.addPolygon(pg);
    }

    else if (type == wkbMultiPolygon)
    {

        // Vsetky casti viacnasobneho polygonu patria do tej istej lokality.
        OGRMultiPolygon* mp = (OGRMultiPolygon*)geom;

        for (int i = 0; i < mp->getNumGeometries(); i++)
        {
            processPolygon(mp->getGeometryRef(i), forest, tr);
        }
    }

    else if (type == wkbGeometryCollection)
    {
        // Rekurzia prejde aj vnorene kolekcie; nepolygonove casti nic nepridaju.
        OGRGeometryCollection* gc = (OGRGeometryCollection*)geom;

        for (int i = 0; i < gc->getNumGeometries(); i++)
        {
            processPolygon(gc->getGeometryRef(i), forest, tr);
        }
    }
}

