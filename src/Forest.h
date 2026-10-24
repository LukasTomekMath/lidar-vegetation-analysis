#pragma once
#include <cstddef>
#include <string>
#include <vector>

#include "PolygonGeometry.h"

// Lokalita alebo segment z KML: atributy, polygony a odvodene udaje pre spracovanie.
class Forest {
private:
    std::string crsWkt;
    std::string name;
    std::string nameWithoutDiacritics;
    std::vector<ForestPolygon> polygons; //Forest sa moze skladat z viacerych uzemi - polygonov. A kazdy polygon moze mat dieru/y
    double hectares;
    double forestArea;
    double minX, maxX, minY, maxY;
    double rawMinX, rawMaxX, rawMinY, rawMaxY;  //hodnoty bez paddingu

    std::string kod_prales="";
    std::string segment="";
    std::string hodnot_seg = "";
    std::string poznamka = "";

    // Maska v poradi [riadok][stlpec]; true oznacuje vnutro lokality.
    std::vector<std::vector< bool >> mask;
    
    // Vyplni pixely vnutri uzavretej krivky hodnotou fillValue v uz vytvorenej maske.
    // gridX0, gridY0 je lavy horny roh; nx, ny su rozmery masky, pixelSize je krok siete.
    // Krivka musi koncit kopiou prveho bodu; true prida vnutro, false vymaze dieru.
    void rasteriseCurve( Curve& curve, std::vector<std::vector<bool>>& mask, double gridX0, double gridY0, int nx,
        int ny, double pixelSize, bool fillValue);

    int pixelSize = 2; //velkost pixela pre Forest. Toto bude treba pre univerzalnu pracu s LIDAR+Satelit datami urobit inak
    int padding = 2; //padding okolo bounding boxu (sirka padding pasu v pixeloch)

    std::vector<std::string> tiles;

public:
    Forest() : hectares(0.0), forestArea(0),
        minX(0), maxX(0), minY(0), maxY(0),
        rawMinX(0), rawMaxX(0), rawMinY(0), rawMaxY(0),
        kod_prales(""), segment(""),hodnot_seg(""),poznamka("") {}

    void setCRS(const std::string& wkt) { crsWkt = wkt; }
    const std::string& getCRS() const { return crsWkt; }

    // Ulozene hranice po zaokruhleni a pridani paddingu; nastavi ich findBoundingBox().
    double getMinX() { return minX; }
    double getMinY() { return minY; }
    double getMaxX() { return maxX; }
    double getMaxY() { return maxY; }
    // Povodne extremy vrcholov bez zaokruhlenia a paddingu.
    double getRawMinX() const { return rawMinX; }
    double getRawMinY() const { return rawMinY; }
    double getRawMaxX() const { return rawMaxX; }
    double getRawMaxY() const { return rawMaxY; }


    double getPixelSize() { return pixelSize; }
    
    // Povodny nazov a jeho ulozena verzia bez diakritiky.
    const std::string& getName() const { return name; }
    const std::string& getNameWithoutDiacritics() const { return nameWithoutDiacritics; }
    const std::string& getSegment() const { return segment; }
    const std::string& getKodPrales() const { return kod_prales; }
    const std::string& getHodnotSeg() const { return hodnot_seg; }
    const std::string& getPoznamka() const { return poznamka; }

    // Nastavi nazov aj verziu bez diakritiky; pri chybe prevodu ponecha oba povodne.
    void setName(const std::string& n);
    void setSegment(const std::string& s) { segment = s; }
    void setKodPrales(const std::string& k) { kod_prales = k; }
    void setHodnotSeg(const std::string& k) { hodnot_seg = k; }
    void setPoznamka(const std::string& k) { poznamka = k; }

    // Pristup k ulozenym datam ich neprepocitava; po zmene geometrie obnov odvodene udaje.
    std::vector<ForestPolygon>& getPolygons()  { return polygons; }
    std::vector<std::string>& getTiles() { return tiles; }
    std::vector<std::vector< bool >>& getMask() { return mask; }
    const std::vector<std::vector<bool>>& getMask() const { return mask; }

    // Prida polygon bez prepocitania plochy, hranic, dlazdic a masky.
    void addPolygon(const ForestPolygon& pg) { polygons.push_back(pg); }
    std::size_t getPolygonCount() { return polygons.size(); }

    // Plocha z atributu KML v hektaroch, nezavisla od vypoctu z geometrie.
    double getHectares() { return hectares; }
    void setHectares(double h) { hectares = h; }

    // Vypocita plochu v m2 ako sucet ploch polygonov po odcitani dier.
    // Diery mensie nez 1 m2 pritom trvalo odstrani z geometrie.
    // Hranice musia mat prvy bod zopakovany na konci.
    void calculateForestArea();
    // Vrati ulozenu vypocitanu plochu v m2 (pred vypoctom nula).
    double getForestArea() { return forestArea; }

    // Ulozi extremy vonkajsich hranic aj obdlznik zaokruhleny von na 10 m s paddingom.
    // Ak chyba prvy polygon alebo jeho vrcholy, nastavi vsetky hranice na nulu.
    void findBoundingBox();
    // Z aktualneho obdlznika zostavi kandidatske nazvy LAZ pre dlazdice 2 x 2 km.
    // Najprv zavolaj findBoundingBox(); existencia suborov ani duplicity sa nekontroluju.
    void findTiles();

    // Vytvori masku nad obdlznikom z findBoundingBox() s krokom pixelSize.
    // Vonkajsie hranice vyplni, diery vymaze; krivky musia byt uzavrete opakovanim bodu.
    void createMask();

    // Zapise existujucu masku ako jednopasmovy GeoTIFF (0/1, CRS geometrie).
    // Pred volanim vytvor masku a zaregistruj GDAL ovladace; prazdnu masku neexportuje.
    void exportMaskToGeoTIFF(const std::string& filename);
};
