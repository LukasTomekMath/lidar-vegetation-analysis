#pragma once
#include <string>
#include <vector>

#include "Forest.h"

class OGRGeometry;
class OGRCoordinateTransformation;

class ForestManager {
private:
	std::vector<Forest> forests;
    void processPolygon(OGRGeometry* geom,Forest& forest,OGRCoordinateTransformation* tr);

public:
    ForestManager() {};
    std::vector<Forest>& getForests() { return forests; }
    bool parseKML_GDAL(const std::string& filename);
    Forest* getForestByFileName(const std::string& name);
};


