#pragma once
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cmath>

#include "libs/gdal_include/ogr_spatialref.h"
#include "libs/gdal_include/gdal.h"
#include "libs/gdal_include/gdal_priv.h"
#include "libs/gdal_include/gdal.h"

struct Point {
    double lon, lat, alt;// povodne suradnice
    double X, Y, Z; //transormovane suradnice
};


class Curve {
private:
	std::vector<Point> points;
	std::vector<double> edgeLengths;
	int numPoints;
	int length;
    double area;
	int curvePixels;//kolko pixelov zaerie krivka
    double perimeter;
    double curveArea;
public:
    Curve() : numPoints(0), length(0), area(0.0), curvePixels(0),perimeter(0) {}
    void addPoint(Point& p) { points.push_back(p); numPoints = points.size(); }


    std::vector<Point>& getPoints()  { return points; }

    void calculateEdgeLengths();
    void calculatePerimeter();
    void calculateCurveArea();

    double getPerimeter()  { return perimeter; }
    double getcurveArea() { return curveArea; }
};

struct PolygonGroup {
    //jedna outer boundary moze mat aj viac inner boundaries, zaroven potrebujem vediet 
    // ku ktorej outer skutocne patri, kedze jeden les moze mat viacej outer a  
    // teoreticky kazda outer moze mat viac inner
    Curve outer;
    std::vector<Curve> inners;
};

class Forest {
private:
    std::string name;
    std::vector<PolygonGroup> polygons; //les moze mat viac polygonov a kazdy polygon moze mat dieru/y
    double hectares;
    double forestArea;
    double minX, maxX, minY, maxY;

    std::vector<std::vector< bool >> mask;
    void rasteriseCurve( Curve& curve, std::vector<std::vector<bool>>& mask, double gridX0, double gridY0, int nx,
        int ny, double pixelSize, bool fillValue);//univerzalna funkcia aj na inner aj outer
    //maska ktore pixely su vnutri a vonku
    //ratanie feature vektora podla masky-min,max,priemer, std

    //padding specs ale aj definicia velkosti pixela do pocitania. teray je to tu, ale teoreickz to moze bz v maine v pripade ze by royne lesy mali mat rozne velkosti i guess
    int pixelSize = 5; 
    int padding = 2;

    std::vector<std::string> tiles;

public:
    Forest() : hectares(0.0), forestArea(0),minX(0), maxX(0), minY(0), maxY(0) {}

    double getMinX() { return minX; }
    double getMinY() { return minY; }
    double getMaxX() { return maxX; }
    double getMaxY() { return maxY; }
    double getPixelSize() { return pixelSize; }
    
    const std::string& getName() const { return name; }
    void setName(const std::string& n) { name = n; }

    std::vector<PolygonGroup>& getPolygons()  { return polygons; }
    std::vector<std::string>& getTiles() { return tiles; }
    std::vector<std::vector< bool >>& getMask() { return mask; }

    void addPolygon(const PolygonGroup& pg) { polygons.push_back(pg); }
    size_t getPolygonCount() { return polygons.size(); }

    double getHectares() { return hectares; }
    void setHectares(double h) { hectares = h; }

    void calculateForestArea();
    double getForestArea() { return forestArea; }

    void findBoundingBox();
    void findTiles();

    void createMask();

    void exportMaskToGeoTIFF(const std::string& filename);
};

class ForestManager {
private:
	std::vector<Forest> forests;
	bool parseKML(const std::string& content);

public:
    ForestManager() {};
    std::vector<Forest>& getForests() { return forests; }
	bool loadFromKML(const std::string& filename);
};


