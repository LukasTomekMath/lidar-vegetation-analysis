#pragma once
#include <vector>

// Bod uchovava povodne aj transformovane suradnice toho isteho vrcholu.
struct Point {
    double lon, lat, alt;// povodne suradnice
    double X, Y, Z; //transformovane suradnice
};


// Usporiadane vrcholy jednej hranice polygonu; vypocty pouzivaju iba suradnice X, Y.
// Pri uzavretej krivke musi byt prvy bod zopakovany aj na konci zoznamu.
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
    // Prida vrchol na koniec krivky; dlzky hran, obvod ani plochu neprepocita.
    void addPoint(Point& p) { points.push_back(p); numPoints = points.size(); }

    std::vector<Point>& getPoints()  { return points; }

    // Prepocita dlzky useciek medzi susednymi vrcholmi v rovine XY.
    // Uzatvaraciu hranu zahrnie iba pri zopakovani prveho bodu na konci.
    void calculateEdgeLengths();

    // Vypocita obvod z ulozenych dlzok hran; najprv treba zavolat calculateEdgeLengths().
    void calculatePerimeter();

    // Plocha uzavreta krivkou (v rovine XY).
    // Prvy bod musi byt zopakovany na konci; dlzky hran netreba vopred vypocitat.
    void calculateCurveArea();

    // Vrati ulozeny obvod bez prepocitania (pred prvym vypoctom je nulovy).
    double getPerimeter()  { return perimeter; }
    
    // Vrati ulozenu plochu; pred prvym citanim treba zavolat calculateCurveArea().
    double getcurveArea() { return curveArea; }
};


// Geometria jedneho polygonu; odcitanie dier od plochy riesi nadradena trieda Forest.
// Segment/les sa moze skladat z viacerych uzemi. Kazdemu uzemiu zodpoveda tento ForestPolygon. 
// Kazdy ForestPolygon ma jednu outer boundary a moze (ale nemusi) mat aj inner boundaries
struct ForestPolygon {
    
    Curve outer; //outer boundary
    std::vector<Curve> inners; // inner boundary (moze byt z viacerych kriviek)
};
