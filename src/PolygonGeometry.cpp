#include "PolygonGeometry.h"
#include <cmath>

void Curve::calculateEdgeLengths()
{
    edgeLengths.clear();
    if (points.size() < 2) return;

    double total = 0.0;
    for (size_t i = 1; i < points.size(); i++)
    {
        double dx = points[i].X - points[i - 1].X;
        double dy = points[i].Y - points[i - 1].Y;
        double len = std::sqrt(dx * dx + dy * dy);

        edgeLengths.push_back(len);
    }
}

void Curve::calculatePerimeter() 
{
    //dlzka krivky (v rovine xy)
    double sum = 0.0;
    for (int i=0;i<edgeLengths.size();i++)
        sum += edgeLengths[i];
    perimeter = sum;
}

void Curve::calculateCurveArea()
{
    //plocha uzavreta krivkou (v rovine xy)
    double sum = 0.0;
    for (int i = 0; i < numPoints-1; i++)
    {
        sum += points[i].X * points[i + 1].Y - points[i + 1].X * points[i].Y;
    }
    curveArea = 0.5 * fabs(sum);
}
