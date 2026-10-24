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
        if (entry.is_regular_file() && entry.path().extension() == ".kml")
        {
            const std::string filename = entry.path().string();
            std::cout << "Loading KML: " << filename << std::endl;

            if (!parseKML_GDAL(filename))
                throw std::runtime_error("Failed to load KML: " + filename);
            ++fileCount;
        }
    }
    std::cout << "Finished reading " << fileCount << " KML files.\n";
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

    std::unique_ptr<GDALDataset, decltype(&GDALClose)> datasetGuard(ds, GDALClose);
    OGRSpatialReference src;
    OGRSpatialReference dst(targetCRS);
    if (src.importFromEPSG(4326) != OGRERR_NONE)
    {
        throw std::runtime_error("Cannot load EPSG:4326");
    }

    // Zachovaj poradie GIS osi: zemepisna dlzka/sirka a vystupne X/Y.
    src.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    dst.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    OGRCoordinateTransformation* tr = OGRCreateCoordinateTransformation(&src, &dst);
    if (!tr)
    {
        throw std::runtime_error("Cannot create KML coordinate transformation");
    }

    std::unique_ptr<OGRCoordinateTransformation, decltype(&OCTDestroyCoordinateTransformation)>
        transformGuard(tr, OCTDestroyCoordinateTransformation);

    for (int i = 0; i < ds->GetLayerCount(); i++)
    {
        OGRLayer* layer = ds->GetLayer(i);
        layer->ResetReading();
        OGRFeature* feature;
        while ((feature = layer->GetNextFeature()) != nullptr)
        {
            Forest forest;
            forest.setCRS(spatial::toWkt(dst));


            // Zakladny nazov z Name moze neskor prepisat atribut NAZOV.
            if (feature->GetFieldIndex("Name") >= 0)
            {
                std::string name = feature->GetFieldAsString("Name");
                if (!name.empty())
                    forest.setName(std::regex_replace(name, std::regex(" "), "_"));
            }

            // Prevezmi podporovane atributy; ostatne polia zostanu nevyuzite.
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

            // Prvok bez geometrie sa neprida; nepodporovana geometria moze dat prazdny Forest.
            OGRGeometry* geom = feature->GetGeometryRef();
            if (!geom)
            {
                OGRFeature::DestroyFeature(feature);
                continue;
            }

            processPolygon(geom, forest, tr);
            forests.push_back(forest);
            OGRFeature::DestroyFeature(feature);
        }

        // Navrat je vo vnutri cyklu vrstiev, preto sa dalsie vrstvy uz nenacitaju.
        std::cout << "Loaded " << forests.size() << " forests\n";
        return true;
    }
    return false;
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

