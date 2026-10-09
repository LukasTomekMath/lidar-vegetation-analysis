#include "ForestManager.h"
#include <iostream>
#include <regex>

#include "libs/gdal_include/ogrsf_frmts.h"

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

bool ForestManager::parseKML_GDAL(const std::string& filename)
{
    GDALDataset* ds = (GDALDataset*)GDALOpenEx(
        filename.c_str(),
        GDAL_OF_VECTOR,
        NULL, NULL, NULL);
    if (!ds)
    {
        std::cout << "Failed to open dataset\n";
        return false;
    }

    OGRSpatialReference src, dst;
    src.importFromEPSG(4326);
    dst.importFromEPSG(8353);

    src.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    dst.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    OGRCoordinateTransformation* tr = OGRCreateCoordinateTransformation(&src, &dst);

    for (int i = 0; i < ds->GetLayerCount(); i++)
    {
        OGRLayer* layer = ds->GetLayer(i);
        layer->ResetReading();
        OGRFeature* feature;
        while ((feature = layer->GetNextFeature()) != nullptr)
        {
            Forest forest;


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
                    if (!name.empty())  //ak je NAZOV, tak ten sa zapise
                        forest.setName(std::regex_replace(name, std::regex(" "), "_"));
                }
                else if (fieldName == "HODNOT_SEG")
                {
                    forest.setHodnotSeg(feature->GetFieldAsString(f));
                }
                else if (fieldName == "HECTARES")
                {
                    forest.setHectares(feature->GetFieldAsDouble(f));
                }
                else if (fieldName == "SEGMENT")
                {
                    forest.setSegment(feature->GetFieldAsString(f));
                }
                else if (fieldName == "KOD_PRALES")
                {
                    forest.setKodPrales(feature->GetFieldAsString(f));
                }
                else if (fieldName == "POZNAMKA")
                {
                    forest.setPoznamka(feature->GetFieldAsString(f));
                }
            }

            OGRGeometry* geom = feature->GetGeometryRef();
            if (!geom)
                continue;

            processPolygon(geom, forest, tr);
            forests.push_back(forest);
            OGRFeature::DestroyFeature(feature);
        }

        GDALClose(ds);
        std::cout << "Loaded " << forests.size() << " forests\n";
        return true;
    }
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

        for (int i = 0; i < outer->getNumPoints(); i++)
        {
            Point p;

            p.lon = outer->getX(i);
            p.lat = outer->getY(i);
            p.alt = outer->getZ(i);

            double x = p.lon;
            double y = p.lat;
            double z = p.alt;

            OCTTransform(tr, 1, &x, &y, &z);

            p.X = x;
            p.Y = y;
            p.Z = z;

            outerCurve.addPoint(p);
        }

        outerCurve.calculateEdgeLengths();
        outerCurve.calculatePerimeter();

        pg.outer = outerCurve;

        int innerCount = poly->getNumInteriorRings();


        for (int r = 0; r < innerCount; r++)
        {
            OGRLinearRing* ring = poly->getInteriorRing(r);

            if (!ring)
                continue;

            Curve inner;

            for (int i = 0; i < ring->getNumPoints(); i++)
            {
                Point p;

                p.lon = ring->getX(i);
                p.lat = ring->getY(i);
                p.alt = ring->getZ(i);

                double x = p.lon;
                double y = p.lat;
                double z = p.alt;

                OCTTransform(tr, 1, &x, &y, &z);

                p.X = x;
                p.Y = y;
                p.Z = z;

                inner.addPoint(p);
            }

            inner.calculateEdgeLengths();
            inner.calculatePerimeter();

            pg.inners.push_back(inner);
        }

        forest.addPolygon(pg);
    }

    else if (type == wkbMultiPolygon)
    {

        OGRMultiPolygon* mp = (OGRMultiPolygon*)geom;

        for (int i = 0; i < mp->getNumGeometries(); i++)
        {
            processPolygon(mp->getGeometryRef(i), forest, tr);
        }
    }

    else if (type == wkbGeometryCollection)
    {
        OGRGeometryCollection* gc = (OGRGeometryCollection*)geom;

        for (int i = 0; i < gc->getNumGeometries(); i++)
        {
            processPolygon(gc->getGeometryRef(i), forest, tr);
        }
    }
}

