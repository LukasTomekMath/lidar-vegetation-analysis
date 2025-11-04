#define _USE_MATH_DEFINES
#include "PolygonManager.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cmath>



bool ForestManager::loadFromKML(const std::string& filename)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Could not open KML file: " << filename << std::endl;
        return false;
    }

    std::stringstream buffer;
    buffer << file.rdbuf();
    file.close();

    return parseKML(buffer.str());
}


bool ForestManager::parseKML(const std::string& content)
{
    size_t pos = 0;

    while (true)
    {
        size_t start = content.find("<Placemark", pos);
        if (start == std::string::npos) break;

        size_t end = content.find("</Placemark>", start);
        if (end == std::string::npos) break;    //ak sa  nenajde koniec placemarku, skonci sa parsovanie, asi neni dobry subor

        std::string placemark = content.substr(start, end - start);
        std::string name;

        size_t nameStart = placemark.find("<name>");
        size_t nameEnd = placemark.find("</name>");
        if (nameStart != std::string::npos)
        {
            nameStart += 6;
            name = placemark.substr(nameStart, nameEnd - nameStart);
        }

        Forest forest;
		forest.setName(name);

        size_t descStart = placemark.find("<description>");
        size_t descEnd = placemark.find("</description>");
        if (descStart != std::string::npos) 
        {
            descStart += 13;
            std::string description = placemark.substr(descStart, descEnd - descStart);

            //toto nefunguje
            size_t hectStart = description.find("HECTARES");
            if (hectStart != std::string::npos) 
            {
                size_t tdStart = description.find("<td>", hectStart);
                if (tdStart != std::string::npos) {
                    tdStart = description.find("<td>", tdStart + 4);
                    size_t tdEnd = description.find("</td>", tdStart);
                    if (tdStart != std::string::npos && tdEnd != std::string::npos) 
                    {
                        tdStart += 4;
                        std::string hectStr = description.substr(tdStart, tdEnd - tdStart);
                        std::replace(hectStr.begin(), hectStr.end(), ',', '.');
                        double hectares = std::stod(hectStr);
                        forest.setHectares(hectares);
                    }
                }
            }
        }

        size_t polyPos = 0;
        while (true)
        {
            size_t startPoly = placemark.find("<Polygon>", polyPos);
            if (startPoly == std::string::npos) break;

            size_t endPoly = placemark.find("</Polygon>", startPoly);
            if (endPoly == std::string::npos) break;

            std::string polyBlock = placemark.substr(startPoly, endPoly - startPoly);

            PolygonGroup pg;

            size_t outerStart = polyBlock.find("<outerBoundaryIs>");
            size_t outerEnd = polyBlock.find("</outerBoundaryIs>");
            if (outerStart != std::string::npos && outerEnd != std::string::npos)
            {
                size_t coordStart = polyBlock.find("<coordinates>", outerStart);
                size_t coordEnd = polyBlock.find("</coordinates>", coordStart);
                if (coordStart != std::string::npos && coordEnd != std::string::npos)
                {
                    coordStart += 13;
                    std::string coordText = polyBlock.substr(coordStart, coordEnd - coordStart);

                    std::stringstream ss(coordText);
                    std::vector<Point> points;
                    std::string token;
                    Curve outerCurve;
                    while (ss >> token)
                    {
                        std::replace(token.begin(), token.end(), ',', ' ');
                        std::stringstream pt(token);

                        Point p;
                        pt >> p.lon >> p.lat >> p.alt;
                        outerCurve.addPoint(p);
                    }
                    pg.outer = outerCurve;
                }
            }
            size_t innerPos = 0;
            while (true)
            {
                size_t innerStart = polyBlock.find("<innerBoundaryIs>", innerPos);
                if (innerStart == std::string::npos) break;

                size_t innerEnd = polyBlock.find("</innerBoundaryIs>", innerStart);
                if (innerEnd == std::string::npos) break;

                std::string innerBlock = polyBlock.substr(innerStart, innerEnd - innerStart);

                // one <innerBoundaryIs> can contain multiple <LinearRing>
                size_t ringPos = 0;
                while (true)
                {
                    size_t ringStart = innerBlock.find("<LinearRing>", ringPos);
                    if (ringStart == std::string::npos) break;
                    size_t ringEnd = innerBlock.find("</LinearRing>", ringStart);
                    if (ringEnd == std::string::npos) break;

                    std::string ringBlock = innerBlock.substr(ringStart, ringEnd - ringStart);

                    size_t coordStart = ringBlock.find("<coordinates>");
                    size_t coordEnd = ringBlock.find("</coordinates>", coordStart);
                    if (coordStart != std::string::npos && coordEnd != std::string::npos)
                    {
                        coordStart += 13;
                        std::string coordText = ringBlock.substr(coordStart, coordEnd - coordStart);

                        std::stringstream ss(coordText);
                        std::string token;
                        Curve innerCurve;

                        while (ss >> token)
                        {
                            std::replace(token.begin(), token.end(), ',', ' ');
                            std::stringstream pt(token);
                            Point p;
                            p.alt = 0;
                            pt >> p.lon >> p.lat >> p.alt;
                            innerCurve.addPoint(p);
                        }

                        pg.inners.push_back(innerCurve);
                    }
                    ringPos = ringEnd + 13;
                }
                innerPos = innerEnd + 17;
            }
            forest.addPolygon(pg);
            polyPos = endPoly + 9;
        }
        forests.push_back(forest);
        pos = end + 13; 
    }
    std::cout << "Loaded " << forests.size() << " polygons total.\n";
    return true;
}