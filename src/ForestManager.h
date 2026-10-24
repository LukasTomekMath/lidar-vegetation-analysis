#pragma once
#include <string>
#include <vector>
#include <stdexcept>

#include "Forest.h"
#include "RasterIO.h"

class OGRGeometry;
class OGRCoordinateTransformation;

// Uchovava lokality nacitane z KML a zabezpecuje prevod ich geometrie do Forest.
class ForestManager {
private:
	std::vector<Forest> forests;
    OGRSpatialReference targetCRS;
    // Prida polygony a ich diery do forest
    // Uchova povodne suradnice vrcholov (lon, lat, alt) a pocita aj transformovane X,Y,Z
    // Vypocita dlzky hran a obvody, nie plochu, obdlznik ani masku.
    void processPolygon(OGRGeometry* geom,Forest& forest,OGRCoordinateTransformation* tr);

public:
    explicit ForestManager(const OGRSpatialReference& target)
        : targetCRS(spatial::comparisonCRS(target))
    {
        if (targetCRS.IsEmpty()) throw std::runtime_error("KML transformation needs a known target CRS");
        if (targetCRS.IsCompound() && targetCRS.StripVertical() != OGRERR_NONE)
            throw std::runtime_error("Cannot extract KML target horizontal CRS");
        targetCRS.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    }
    // Vrati upravitelny zoznam vsetkych doteraz nacitanych lokalit.
    std::vector<Forest>& getForests() { return forests; }
    // Nacita bezne .kml subory priamo z priecinka (bez rekurzie) a vypise suhrn.
    // Pred volanim inicializuj GDAL. 
    void loadDirectory(const std::string& directory);
    // Prida prvky prvej vrstvy KML do zoznamu.
    // Predpoklada vstup EPSG:4326 a prevadza ho do CRS zadaneho konstruktorom.
    bool parseKML_GDAL(const std::string& filename);
    // Vrati prvu lokalitu s presnou zhodou getName()
    // Argument musi byt uz samotny nazov (cesta ani pripona suboru sa neodstranuju)
    Forest* getForestByFileName(const std::string& name);
};


