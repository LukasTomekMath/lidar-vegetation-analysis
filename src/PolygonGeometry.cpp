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
        // Euklidovska vzdialenost v rovine XY; vyska Z sa nezohladnuje.
        double len = std::sqrt(dx * dx + dy * dy);

        edgeLengths.push_back(len);
    }
}

void Curve::calculatePerimeter() 
{ 
    double sum = 0.0;
    for (int i=0;i<edgeLengths.size();i++)
        sum += edgeLengths[i];
    perimeter = sum;
}

void Curve::calculateCurveArea()
{
    // Vypocet plochy uzavretej rovinnou krivkou pomocou suradnic vrcholov
    double sum = 0.0;
    for (int i = 0; i < numPoints-1; i++)
    {
        sum += points[i].X * points[i + 1].Y - points[i + 1].X * points[i].Y;
    }
    // Absolutna hodnota odstrani zavislost znamienka od smeru poradia vrcholov.
    curveArea = 0.5 * fabs(sum);
}
