#pragma once
#include <vector>

struct Point {
    double lon, lat, alt;// povodne suradnice
    double X, Y, Z; //transformovane suradnice
};


class Curve {
private:
	std::vector<Point> points;
	std::vector<double> edgeLengths;
	int numPoints;
	int length;
    double area;
	int curvePixels;//kolko pixelov zaberie krivka
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


struct ForestPolygon {
    //Segment/les sa moze skladat z viacerych uzemi. Kazdemu uzemiu zodpoveda tento ForestPolygon. 
    //Kazdy ForestPolygon ma jednu outer boundary a moze (ale nemusi) mat aj inner boundaries
    Curve outer; //outer boundary
    std::vector<Curve> inners; // inner boundary (moze byt z viacerych kriviek)
};
