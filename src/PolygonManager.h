#pragma once
#include <string>
#include <vector>

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
public:
    Curve() : numPoints(0), length(0), area(0.0), curvePixels(0) {}
    void addPoint(const Point& p) { points.push_back(p); numPoints = points.size(); }
    const std::vector<Point>& getPoints() const { return points; }
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
    double hectares = 0.0;
    double minX, maxX, minY, maxY;

public:
    Forest() : hectares(0.0), minX(0), maxX(0), minY(0), maxY(0) {}
    const std::string& getName() const { return name; }
    void setName(const std::string& n) { name = n; }
    void addPolygon(const PolygonGroup& pg) { polygons.push_back(pg); }
    const std::vector<PolygonGroup>& getPolygons() const { return polygons; }
    size_t getPolygonCount() { return polygons.size(); }
    void setHectares(double h) { hectares = h; }
    double getHectares() const { return hectares; }
};

class ForestManager {
private:
	std::vector<Forest> forests;
	bool parseKML(const std::string& content);

public:
    ForestManager() {};
    const std::vector<Forest>& getForests() const { return forests; }
	bool loadFromKML(const std::string& filename);
};


