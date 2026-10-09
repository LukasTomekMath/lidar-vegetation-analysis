#pragma once
#include <cstddef>
#include <string>
#include <vector>

#include "PolygonGeometry.h"

class Forest {
private:
    std::string name;
    std::vector<ForestPolygon> polygons; //Segment sa moze skladat z viacerych uzemi - polygonov. A kazdy polygon moze mat dieru/y
    double hectares;
    double forestArea;
    double minX, maxX, minY, maxY;
    double rawMinX, rawMaxX, rawMinY, rawMaxY;  //hodnoty bez paddingu

    std::string kod_prales="";
    std::string segment="";
    std::string hodnot_seg = "";
    std::string poznamka = "";

    std::vector<std::vector< bool >> mask;
    void rasteriseCurve( Curve& curve, std::vector<std::vector<bool>>& mask, double gridX0, double gridY0, int nx,
        int ny, double pixelSize, bool fillValue);//univerzalna funkcia aj na inner aj outer
    //maska, ktore pixely su vnutri a vonku
    //ratanie feature vektora podla masky - min, max, priemer, std

    int pixelSize = 2; //velkost pixela pre Forest. Toto bude treba pre univerzalnu pracu s LIDAR+Satelit datami urobit inak
    int padding = 2; //padding okolo bounding boxu (sirka padding pasu v pixeloch)

    std::vector<std::string> tiles;

public:
    Forest() : hectares(0.0), forestArea(0),
        minX(0), maxX(0), minY(0), maxY(0),
        rawMinX(0), rawMaxX(0), rawMinY(0), rawMaxY(0),
        kod_prales(""), segment(""),hodnot_seg(""),poznamka("") {}

    double getMinX() { return minX; }
    double getMinY() { return minY; }
    double getMaxX() { return maxX; }
    double getMaxY() { return maxY; }
    double getRawMinX() const { return rawMinX; }
    double getRawMinY() const { return rawMinY; }
    double getRawMaxX() const { return rawMaxX; }
    double getRawMaxY() const { return rawMaxY; }


    double getPixelSize() { return pixelSize; }
    
    const std::string& getName() const { return name; }
    const std::string& getSegment() const { return segment; }
    const std::string& getKodPrales() const { return kod_prales; }
    const std::string& getHodnotSeg() const { return hodnot_seg; }
    const std::string& getPoznamka() const { return poznamka; }

    void setName(const std::string& n) { name = n; }
    void setSegment(const std::string& s) { segment = s; }
    void setKodPrales(const std::string& k) { kod_prales = k; }
    void setHodnotSeg(const std::string& k) { hodnot_seg = k; }
    void setPoznamka(const std::string& k) { poznamka = k; }

    std::vector<ForestPolygon>& getPolygons()  { return polygons; }
    std::vector<std::string>& getTiles() { return tiles; }
    std::vector<std::vector< bool >>& getMask() { return mask; }
    const std::vector<std::vector<bool>>& getMask() const { return mask; }

    void addPolygon(const ForestPolygon& pg) { polygons.push_back(pg); }
    std::size_t getPolygonCount() { return polygons.size(); }

    double getHectares() { return hectares; }
    void setHectares(double h) { hectares = h; }

    void calculateForestArea();
    double getForestArea() { return forestArea; }

    void findBoundingBox();
    void findTiles();

    void createMask();

    void exportMaskToGeoTIFF(const std::string& filename);
};
