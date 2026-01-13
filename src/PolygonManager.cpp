#define _USE_MATH_DEFINES
#include "PolygonManager.h"


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

    OGRSpatialReference src, dst;
    src.importFromEPSG(4326);
    dst.importFromEPSG(8353);

    src.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    dst.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

    OGRCoordinateTransformation* tr = OGRCreateCoordinateTransformation(&src, &dst);


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

                        double x = p.lon;
                        double y = p.lat;
                        double z = p.alt;

                        tr->Transform(1, &x, &y, &z);
                        p.X = x;
                        p.Y = y;
                        p.Z = z;

                        outerCurve.addPoint(p);
                    }
                    outerCurve.calculateEdgeLengths();
                    outerCurve.calculatePerimeter();
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

                            double x = p.lon;
                            double y = p.lat;
                            double z = p.alt;

                            tr->Transform(1, &x, &y, &z);
                            p.X = x;
                            p.Y = y;
                            p.Z = z;
                            innerCurve.addPoint(p);
                        }
                        innerCurve.calculateEdgeLengths();
                        innerCurve.calculatePerimeter();
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

void Curve::calculateEdgeLengths()
{
    edgeLengths.clear();
    if (points.size() < 2) return;

    double total = 0.0;
    for (size_t i = 1; i < points.size(); i++)
    {
        double dx = points[i].X - points[i - 1].X;
        double dy = points[i].Y - points[i - 1].Y;
        double dz = points[i].Z - points[i - 1].Z;
        double len = std::sqrt(dx * dx + dy * dy + dz * dz);

        edgeLengths.push_back(len);
    }
}

void Forest::findBoundingBox()
{
    if (polygons.empty() || polygons[0].outer.getPoints().empty()) {
        minX = maxX = minY = maxY = 0;
        return;
    }

    int pad = padding*pixelSize;//tu bude potom nejaký padding, este neviem ako, asi podla velkosti mriezky aeo taak
    //aleo mozno to bude member classy, este nvm

    //inicializacia na prvy bod s tym udem najprv porovnavat
    minX = polygons[0].outer.getPoints()[0].X;
    maxX= polygons[0].outer.getPoints()[0].X;
    minY = polygons[0].outer.getPoints()[0].Y;
    maxY = polygons[0].outer.getPoints()[0].Y;

    for (auto& pg : polygons)
    {
        for (auto& p : pg.outer.getPoints()) {
            minX = std::min(minX, p.X);
            maxX = std::max(maxX, p.X);
            minY = std::min(minY, p.Y);
            maxY = std::max(maxY, p.Y);
        }
    }
    minX = std::floor(minX / 10.0) * 10.0 -pad;
    minY = std::floor(minY / 10.0) * 10.0-pad;
    maxX = std::ceil(maxX / 10.0) * 10.0 +pad;
    maxY = std::ceil(maxY / 10.0) * 10.0+pad;
}

void Forest::findTiles()
{
    tiles.clear();
    int minTileX = static_cast<int>(std::floor(minX / 2000) * 2000);
    int maxTileX = static_cast<int>(std::floor(maxX / 2000) * 2000);
    int minTileY = static_cast<int>(std::floor(minY / 2000) * 2000);
    int maxTileY = static_cast<int>(std::floor(maxY / 2000) * 2000);

    for (int tx = minTileX; tx <= maxTileX; tx += 2000)
    {
        for (int ty = minTileY; ty <= maxTileY; ty += 2000)
        {
            tiles.push_back("tile_" + std::to_string(tx) + "_" + std::to_string(ty) + ".laz");
        }
    }
}

void Forest::calculateForestArea()
{
    int num = polygons.size();
    double sum_area = 0.0;

    for (int i = 0; i < num; i++)
    {
        PolygonGroup& pg = polygons[i];

        double sum_area_inner = 0.0;
        pg.outer.calculateCurveArea();
        if (!polygons[i].inners.empty())
        {
            int num_inners = pg.inners.size();
            for (int j = 0; j < num_inners; j++)
            {
                pg.inners[j].calculateCurveArea();
                double area = std::abs(pg.inners[j].getcurveArea());
                if (area < 1)
                {
                    pg.inners.erase(pg.inners.begin() + j);
                    num_inners--;
                    j--;
                }
                else
                {
                    sum_area_inner += area;
                }
            }
        }
        sum_area += pg.outer.getcurveArea() - sum_area_inner;
    }
    forestArea = sum_area;
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
    //toto je na plochu jednej krivky, zvlast bude funkcia na plochu lesa lebo 1 les moze byt viac kriviek
    double sum = 0.0;
    for (int i = 0; i < numPoints-1; i++)
    {
        sum += points[i].X * points[i + 1].Y - points[i + 1].X * points[i].Y;
    }
    curveArea = 0.5 * fabs(sum);
}



void Forest::rasteriseCurve(Curve& curve, std::vector<std::vector<bool>>& mask, double gridX0, double gridY0, int nx,
    int ny, double pixelSize, bool fillValue)
{
    //chcem to mat vseoebcne, lebo rovnko musim riesit outer aj inner boundary- pre outer vyplnam true (je vnutri), pre inner vyplnam false(je vonku)
    //je to ale vlastne iba scan line alg jak na PG
    //curve= krivka ktoru riesim, mask= kam zapisujem, gridX0 a gridY0= pociatocna suradnica (vlavo dole), nx a ny=pocet pixelov ktore riesim, fillvalue= ci to budem vyplnat true alebo false

    auto& pts = curve.getPoints();
    int n = pts.size();

    for (int j = 0; j < ny; j++)
    {
        //y-ova suradnica ktorej idem hladat priesecniky
        double y = gridY0 + (j + 0.5) * pixelSize;//pozicia stredu pixela v pociatku danej olasti
        std::vector<double> xIntersections;

        for (int i = 0; i < n - 1; i++)
        {
            //hrana ktoru prave idem riesit ci ma s mojim riadku priesecnik
            Point& a = pts[i];
            Point& b = pts[i + 1];

            if ((a.Y > y && b.Y < y) || (a.Y < y && b.Y > y))
            {
                //upravenie parametrizacie usecky, hladam t pre ktore je priesecnik, potom že ake je to x
                //x(t) = x1 + t (x2 - x1)
                //y(t) = y1 + t(y2 - y1)
                double t = (y - a.Y) / (b.Y - a.Y);
                double x = a.X + t * (b.X - a.X);
                xIntersections.push_back(x);
            }
        }

        if (xIntersections.size() < 2)
            continue;

        std::sort(xIntersections.begin(), xIntersections.end());

        for (int k = 0; k < xIntersections.size() - 1; k += 2)
        {
            //ak to ma na danom riadku viac priesecnikov, tak
            //x1----vnutri---x2---vonku---x3---vnutri---x4 atd takze zaciatocny pre "dnu" je kazdy druhy (0,2,..)
            double x0 = xIntersections[k];
            double x1 = xIntersections[k + 1];

            //pozeram ze ktory v poradi pixel je zaciatocny a koncovy pre "vnutri"
            //tie -0.5 zabecpecia ze len stred pixela ktory je vnutri bude dany ako vnutri
            int iStart = std::ceil((x0 - gridX0) / pixelSize - 0.5);
            int iEnd = std::floor((x1 - gridX0) / pixelSize - 0.5);

            iStart = std::max(iStart, 0);
            iEnd = std::min(iEnd, nx - 1);


            for (int i = iStart; i <= iEnd; ++i)
            {
                mask[j][i] = fillValue;
            }
        }
    }
}

void Forest::createMask()
{
    int nx = std::ceil((maxX - minX) / pixelSize);
    int ny = std::ceil((maxY - minY) / pixelSize);

    double gridX0 = minX;
    double gridY0 = minY;

    mask.assign(ny, std::vector<bool>(nx, false));


    for (int i = 0; i < polygons.size(); i++)
    {
        PolygonGroup& pg = polygons[i];

        //outer je vzdy len jedna krivka:
        rasteriseCurve(pg.outer, mask, gridX0, gridY0, nx, ny, pixelSize, true);

        //inner cisto teoreticky moze yt aj viac, aj ked myslim ze neni ak dam prec tie smeti male
        for (int inner = 0; inner < pg.inners.size(); inner++)
        {
            Curve& innerCurve = pg.inners[inner];
            rasteriseCurve(innerCurve, mask, gridX0, gridY0, nx, ny, pixelSize, false);
        }
    }
}

void Forest::exportMaskToGeoTIFF(const std::string& filename)
{
    int ny = mask.size();
    int nx = mask.empty() ? 0 : mask[0].size();
    if (nx == 0 || ny == 0) return;

    GDALAllRegister(); 

    GDALDriver* driver = GetGDALDriverManager()->GetDriverByName("GTiff");
    if (!driver)
    {
        std::cerr << "GTiff driver not available\n";
        return;
    }

    GDALDataset* ds = driver->Create(
        filename.c_str(),
        nx,
        ny,
        1,
        GDT_Byte,
        nullptr
    );

    if (!ds)
    {
        std::cerr << "Failed to create GeoTIFF\n";
        return;
    }

    // GeoTransform: row 0 = top
    double geoTransform[6] = { 0 };
    geoTransform[0] = minX;       // top-left X
    geoTransform[1] = pixelSize;  // pixel width
    geoTransform[2] = 0.0;
    geoTransform[3] = maxY;       // top-left Y
    geoTransform[4] = 0.0;
    geoTransform[5] = -pixelSize; // negative to flip Y
    ds->SetGeoTransform(geoTransform);

    // Set projection via EPSG
    OGRSpatialReference srs;
    srs.importFromEPSG(5514); // JTSK / Krovak
    char* wkt = nullptr;
    srs.exportToWkt(&wkt);
    ds->SetProjection(wkt);
    CPLFree(wkt);

    // Write mask, flipping vertically so row 0 = north
    GDALRasterBand* band = ds->GetRasterBand(1);
    std::vector<uint8_t> row(nx);

    for (int y = 0; y < ny; ++y)
    {
        int maskRow = ny - 1 - y; // flip: mask row 0 (south) → band row ny-1
        for (int x = 0; x < nx; ++x)
            row[x] = mask[maskRow][x] ? 1 : 0;

        band->RasterIO(
            GF_Write,
            0, y,
            nx, 1,
            row.data(),
            nx, 1,
            GDT_Byte,
            0, 0
        );
    }

    GDALClose(ds);
}