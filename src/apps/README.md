# Návod na používanie programov

Projekt spracúva klasifikované mračná bodov leteckého laserového skenovania
(LAZ), počíta rastrové metriky vegetácie a agreguje ich pre segmenty z KML.

## Prehľad programov

Po builde v konfigurácii **Release** vzniknú v `build/bin/Release` štyri EXE.

| Program | Účel | Hlavný výstup |
| --- | --- | --- |
| `crop_laz.exe` | Vyreže LAZ podľa opísaných obdĺžnikov lokalít z KML s presahom. | Výrezy `.laz` |
| `inspect_forests.exe` | Vypíše lokality z KML a preverí existenciu kandidátnych LAZ súborov. | `pralesy.txt`, `tiles.txt`, konzolový prehľad |
| `lidar_metrics.exe` | Vykoná filtráciu, výškovú normalizáciu a výpočet metrík na sieti 2 × 2 m. | 27-pásmové GeoTIFF rastre, `processing_stats.csv` |
| `segment_features.exe` | Agreguje hodnoty rastrov vnútri polygónov segmentov. | CSV s atribútmi a štatistikami metrík |

## Príprava

Potrebný je Windows, CMake a Visual Studio 2022 s C++ nástrojmi.
Projekt používa pribalené knižnice z `libs`; paralelizáciu zabezpečuje OpenMP,
MPI sa nepoužíva.

Všetky príkazy nižšie spúšťaj z **koreňa repozitára**, ktorý obsahuje
súbor `CMakeLists.txt` a priečinky `src` a `libs`.

Ak `build` už obsahuje konfiguráciu bez explicitného `-A`, použi pri opätovnej
konfigurácii `cmake -S . -B build` bez zmeny generátora/platformy.
CTest overuje nápovedu a odmietnutie chýbajúcich argumentov; nejde o overenie
správnosti vedeckých výpočtov.

Build skopíruje DLL a `proj.db` vedľa EXE. Pri prenose programov zachovaj aj
tieto súbory. `crop_laz` navyše vyžaduje `libs/las2las64.exe` voči pracovnému
adresáru, preto ho spúšťaj z koreňa projektu.

### Vstupy pre príklady

Príklady nižšie predpokladajú toto rozloženie. Cesty nahraď vlastnými;
priečinky `data` a ich obsah si pripravuje používateľ.

```text
lidar-vegetation-analysis/
  data/
    laz/                 vstupné LAZ dlaždice
    kml_lokality/        KML hranice lokalít pre výber oblasti výpočtu
    kml_segmenty/        KML segmenty pre agregáciu metrík
  results/               výstupy vytvorené programami
  libs/
  build/bin/Release/
```

KML lokalít určuje rozsah výpočtu, KML segmentov určuje plochy agregácie.
Zvoľ lokality tak, aby výsledné rastre pokryli celé požadované segmenty;
segmenty môžu siahať aj mimo hraníc pralesových lokalít.
Použi KML so štruktúrou atribútov očakávanou projektom (napr. `NAZOV`,
`KOD_PRALES`, `SEGMENT`), nie ľubovoľný súbor s polygónmi.

Programy čítajú súbory priamo zo zadaných priečinkov, nie rekurzívne.
Používaj prípony `.kml`, `.laz` a `.tif` malými písmenami. LAZ dlaždice
musia zodpovedať priestorovej sieti a pomenovaniu projektu; pre režim celých
dlaždíc sa očakáva tvar `tile_X_Y[...].laz` so súradnicami dlaždice 2 × 2 km.
Samotné premenovanie ľubovoľného LAZ súboru nezabezpečí správny výpočet.

## Odporúčaný postup

1. Priprav vstupné LAZ a KML a zostav programy.
2. Over lokality a dostupnosť LAZ pomocou `inspect_forests`.
3. Vypočítaj rastre pomocou `lidar_metrics` s KML lokalít.
4. Z finálnych mozaík exportuj CSV pomocou `segment_features` s KML segmentov.

`crop_laz` je voliteľný nástroj na prípravu výrezov, napríklad na vizuálnu
kontrolu. Pre uvedený výpočet metrík ho netreba spúšťať.
Vyber si jeden z nasledujúcich návodov podľa používaného terminálu.

## Použitie v CMD (Príkazový riadok)

Otvor CMD (`cmd.exe`) v koreni repozitára. Použi pripravené vstupy podľa
časti Príprava. Najprv zostav programy:

```bat
cmake -S . -B build -G "Visual Studio 17 2022" -A x64
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

Kontrola lokalít a dostupnosti LAZ:

```bat
build\bin\Release\inspect_forests.exe --kml-dir ".\data\kml_lokality" --laz-dir ".\data\laz" --output-dir ".\results\kontrola"
```

Výpočet metrík pre lokality z KML a vytvorenie mozaík:

```bat
build\bin\Release\lidar_metrics.exe --laz-dir ".\data\laz" --kml-dir ".\data\kml_lokality" --output-dir ".\results\metriky_lokality" --threads 6
```

Alternatíva: výpočet celých dlaždíc bez KML:

```bat
build\bin\Release\lidar_metrics.exe --laz-dir ".\data\laz" --output-dir ".\results\metriky_dlazdice" --threads 5
```

Výpočet features a uloženie do CSV:

```bat
build\bin\Release\segment_features.exe --kml-dir ".\data\kml_segmenty" --tif-dir ".\results\metriky_lokality" --output ".\results\segment_features.csv"
```

Voliteľné vytvorenie výrezov LAZ:

```bat
build\bin\Release\crop_laz.exe --kml-dir ".\data\kml_lokality" --laz-dir ".\data\laz" --output-dir ".\results\vyrezy_laz" --threads 6
```

Nápoveda a kontrola návratového kódu:

```bat
build\bin\Release\lidar_metrics.exe --help
echo %ERRORLEVEL%
```

Príkaz `echo %ERRORLEVEL%` spusti na samostatnom riadku hneď po programe,
ktorého výsledok chceš overiť. Nula znamená, že program neoznámil chybu;
nepotvrdzuje úplnosť dát.
Všetky uvedené príkazy na spustenie EXE sú na jednom riadku. Ak ich v CMD
rozdelíš na viac riadkov, použi na konci pokračujúceho riadka znak `^`
bez medzier za ním.

## Použitie v PowerShelli

Otvor PowerShell v koreni repozitára. Postup a vstupy sú rovnaké ako pri CMD.

Zostavenie a kontrola programov:

```powershell
cmake -S . -B build -G "Visual Studio 17 2022" -A x64
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

Kontrola lokalít a dostupnosti LAZ:

```powershell
.\build\bin\Release\inspect_forests.exe --kml-dir ".\data\kml_lokality" --laz-dir ".\data\laz" --output-dir ".\results\kontrola"
```

Výpočet metrík pre lokality z KML a vytvorenie mozaík:

```powershell
.\build\bin\Release\lidar_metrics.exe --laz-dir ".\data\laz" --kml-dir ".\data\kml_lokality" --output-dir ".\results\metriky_lokality" --threads 6
```

Alternatíva: výpočet celých dlaždíc bez KML:

```powershell
.\build\bin\Release\lidar_metrics.exe --laz-dir ".\data\laz" --output-dir ".\results\metriky_dlazdice" --threads 5
```

Agregácia finálnych mozaík pre segmenty do CSV:

```powershell
.\build\bin\Release\segment_features.exe --kml-dir ".\data\kml_segmenty" --tif-dir ".\results\metriky_lokality" --output ".\results\segment_features.csv"
```

Voliteľné vytvorenie výrezov LAZ:

```powershell
.\build\bin\Release\crop_laz.exe --kml-dir ".\data\kml_lokality" --laz-dir ".\data\laz" --output-dir ".\results\vyrezy_laz" --threads 6
```

Nápoveda a kontrola návratového kódu:

```powershell
.\build\bin\Release\lidar_metrics.exe --help
$LASTEXITCODE
```

Hodnotu `$LASTEXITCODE` skontroluj hneď po programe, ktorého výsledok chceš overiť.

## Výstupy a obmedzenia

### Kontrola lokalít: inspect_forests

Výstupný priečinok obsahuje:

- `pralesy.txt` – atribúty, plochy, hranice a pôvodné aj transformované súradnice lokalít;
- `tiles.txt` – kandidátne názvy LAZ súborov pre každú lokalitu.

Konzola vypíše počty nájdených a chýbajúcich ciest. Počty zahŕňajú varianty
názvov a opakované výskyty; nejde o počet unikátnych dlaždíc. Hlásenie
`Missing` preto samo osebe neznamená chýbajúce priestorové pokrytie.
Program číta hlavičky LAZ na zistenie spoločného CRS pre transformáciu KML.
Pri kontrole kandidátnych súborov overuje existenciu ciest; bodové dáta nečíta.

### Rastrové metriky: lidar_metrics

#### Pre lokality z KML

Program spracuje prieniky dlaždíc s rozšírenými opísanými obdĺžnikmi lokalít.
Výstup má rozlíšenie **2 × 2 m** a **27 pásiem**. Výpočet nie je orezaný
na presnú hranicu polygónu; maskovanie segmentu sa vykonáva pri výpočte feature vektora.

```text
results/metriky_lokality/
  tiles/                 čiastkové TIFF dlaždice
  Sur.tif                príklad finálnej mozaiky lokality Šúr
  processing_stats.csv   časy a údaje o spracovaní dlaždíc
```

Čiastkové rastre sa automaticky spoja cez GDALBuildVRT a GDALTranslate
bez interpolácie. Kontroluje sa zhoda mriežky, CRS, pásiem a NoData.
Názov finálneho TIFF je názov lokality bez diakritiky. Názvy lokalít musia
byť jedinečné aj po odstránení diakritiky a bez ohľadu na veľkosť písmen.

Pri chybe výpočtu dlaždice sa mozaika lokality nevytvorí. Chýbajúce kandidátne
LAZ sa preskakujú; nepokryté miesta v mozaike zostávajú NoData. Existencia
finálneho TIFF teda nezaručuje úplné dátové pokrytie.
Ak pre lokalitu nie je dostupná žiadna dlaždica, vytvorenie mozaiky skončí
chybou a program vráti nenulový návratový kód.

#### Pre celé dlaždice bez výberu lokalít

Bez `--kml-dir` sa spracujú celé dlaždice a TIFF sa zapíšu priamo do výstupného
priečinka. Automatická mozaika sa v tomto režime nevytvára. Súbory obsahujúce
v názve `.copc.laz` sa v tomto režime vynechávajú.

### Výpočet feature vektora: segment_features

Pre každú z 27 metrík sa podľa masky segmentu vypočíta priemer (`_mean`),
smerodajná odchýlka (`_std`), minimum (`_min`) a maximum (`_max`). CSV teda
obsahuje **108 štatistických stĺpcov** a atribúty `ForestName`, `Segment`,
`KodPrales`, `HodnotSeg`, `Hectares`, `Poznamka`. Externé biotopové údaje tento
program nepripája. V diplomovej práci je opísané použitie priemeru a smerodajnej
odchýlky (54 príznakov); aktuálny export obsahuje navyše minimá a maximá.

Dôležité pravidlá aktuálnej implementácie:

- Raster musí mať rozlíšenie 2 × 2 m a obsahovať 27 pásiem v poradí
  definovanom v `MetricDefinitions.h`, ktoré používa export `PointCloudMetrics`.
  Program kontroluje počet pásiem; správnosť ich významu a poradia musí
  zabezpečiť používateľ. Rozlíšenie sa neoveruje; maska segmentu používa
  pevný krok 2 m. Staršie exporty s inou zostavou metrík nie sú vhodné.
- Celý opísaný obdĺžnik segmentu bez presahu musí ležať v jednom TIFF.
  Program dlaždice nespája. Segment, ktorý nepokryje žiadny TIFF, sa nezapíše.
- Priradenie sa robí podľa priestorového rozsahu, nie podľa názvu súboru.
  Ak segment pokrýva viac TIFF súborov, môže sa v CSV objaviť viackrát.
  Po exporte skontroluj počet a jedinečnosť segmentov.
- Do `--tif-dir` zadaj priečinok s finálnymi mozaikami, nie podpriečinok `tiles`.
  Čítajú sa iba TIFF priamo v tomto priečinku.
- NoData a nekonečné či NaN hodnoty sa pri agregácii vynechávajú.
  Bez platných hodnôt zostanú štatistiky nulové, preto samotná nula v CSV
  nepreukazuje platné meranie. Smerodajná odchýlka používa delenie počtom
  platných hodnôt, nie počtom zmenšeným o jeden.

CSV je vstupom pre ďalšiu analýzu, napríklad v R. Tento program nevykonáva
PCA, LDA ani QDA.
`ForestName` sa zapisuje bez diakritiky; medzery v názvoch sa pri načítaní KML
nahrádzajú podčiarkovníkmi.

### Výrezy mračien: crop_laz

Program používa `libs/las2las64.exe` a vytvorí pre dostupné vstupné dlaždice
výrezy podľa opísaného obdĺžnika lokality s presahom. Neorezáva na presnú
hranicu polygónu a výrezy jednotlivých dlaždíc nezlučuje. Názov výstupu skladá
z názvu vstupného súboru a názvu lokality bez diakritiky.

Tento krok je užitočný napríklad na prezeranie mračien v externom programe.
**Nie je povinným krokom pred `lidar_metrics`**: uvedený návod počíta metriky
priamo z pôvodných LAZ dlaždíc.

## Prehľad argumentov

| Argument | Význam | Programy |
| --- | --- | --- |
| `--laz-dir DIR` | Existujúci priečinok so vstupnými LAZ | `crop_laz`, `inspect_forests`, `lidar_metrics` |
| `--kml-dir DIR` | Existujúci priečinok s KML | Všetky; voliteľný iba pri `lidar_metrics` |
| `--output-dir DIR` | Výstupný priečinok | `crop_laz`, `inspect_forests`, `lidar_metrics` |
| `--tif-dir DIR` | Existujúci priečinok s rastrovými metrikami | `segment_features` |
| `--output FILE.csv` | Cesta k výstupnému CSV | `segment_features` |
| `--threads N` | Kladný celý počet vlákien | `crop_laz`, `lidar_metrics` |
| `--help` alebo `-h` | Nápoveda bez výpočtu | Všetky |

Predvolený počet vlákien je 6 pre `crop_laz`, 6 pre `lidar_metrics` s KML a 5
bez KML. `segment_features` agreguje sekvenčne a `--threads` nepodporuje.

## Bežné problémy

- **Cesty:** relatívne cesty sa vyhodnocujú voči pracovnému adresáru.
  Cesty s medzerami daj do úvodzoviek. Výstupné priečinky vrátane rodičovského
  priečinka CSV sa vytvoria automaticky; vstupné priečinky musia existovať.
- **Opakované spustenie:** existujúce výstupy rovnakého názvu sa prepisujú.
  Pre nový výpočet použi samostatný čistý výstupný priečinok; staré TIFF-y sa
  automaticky neodstraňujú. Výstupný priečinok musí byť odlišný od LAZ vstupu.
- **Chýbajúce `proj.db` alebo DLL:** skontroluj obsah `build/bin/Release`
  a zopakuj build. Nestačí preniesť samotné EXE.
- **Chýbajúci `libs/las2las64.exe`:** spusti `crop_laz` z koreňa projektu
  a over existenciu pribaleného nástroja.
- **CSV bez očakávaných segmentov:** skontroluj rozsah rastrov, správny priečinok
  TIFF, 27 pásiem a pokrytie celého segmentu jedným rastrom.
- **Chyba spustenia:** skontroluj konzolové hlásenia. Návratový kód posledného
  programu zistíš v CMD cez `echo %ERRORLEVEL%`, v PowerShelli cez `$LASTEXITCODE`;
  nenulová hodnota signalizuje
  chybu, ktorú program oznámil. Nulový kód sám nepotvrdzuje úplnosť dát.

## Orientácia v zdrojovom kóde

Každý program má vlastný `*_main.cpp`. Workflowy sú v `CropLaz.cpp`,
`InspectForests.cpp`, `LidarMetrics.cpp` a `SegmentFeatures.cpp`.
`AppSupport` zdieľa spracovanie argumentov a inicializáciu GDAL/PROJ.
Výpočtové moduly v `src` tvoria statickú knižnicu `vegetation_core`.
Rozdelenie na EXE nemení vedeckú metodiku ani vzorce metrík.

Pri vyhľadávaní LAZ podľa KML sa názov lokality používa bez diakritiky: napr. `tile_X_Y_Sur.laz` pre Šúr. Kandidáty bez názvu lokality (`tile_X_Y.laz`) zostávajú podporované.

### Súradnicový systém

TIFF s metrikami preberá CRS priamo zo vstupného LAS/LAZ, bez nahradenia pevným
EPSG a bez reprojekcie bodov. Číta WKT z VLR/EVLR alebo GeoTIFF GeoKeys
(s kódom EPSG). Vlastné CRS vyžaduje WKT. Ak CRS chýba,
TIFF ho tiež nemá a program vypíše upozornenie. WKT sa číta priamo cez GDAL;
úvodzovky sa neopravujú. Nečitateľné metadáta skončia chybou.
Mozaiky a diagnostické TIFF zachovávajú CRS vstupov. Dlaždice s navzájom rozdielnymi
CRS sa nezlučujú do jednej mozaiky. Ekvivalentný starší ESRI Krovak WKT1
a novší WKT2 sa porovnávajú po normalizácii redundantných parametrov orientácie
v pracovnej kópii; samotný zhodný EPSG štítok nestačí.

`crop_laz` ponecháva prenos CRS na `las2las`, ktorý skopíruje metadáta pri
orezaní. Výstupný LAZ už dodatočne neupravuje. Nezadáva pevné EPSG ani reprojekciu.
Pri exporte TIFF sa CRS musí zapísať explicitne cez GDAL `SetSpatialRef`.

Pri práci s KML sa polygóny transformujú do horizontálneho CRS dát:
`lidar_metrics --kml-dir`, `crop_laz` a `inspect_forests` ho načítajú z hlavičiek
LAZ v priečinku. Tento výber dlaždíc vyžaduje jeden spoločný známy CRS v priečinku;
samotný dlaždicový export bez KML môže spracovať rôzne CRS oddelene.
`segment_features` transformuje polygóny podľa CRS každého TIFF.
Maska lesa má CRS svojej geometrie.

Vzorce metrík, súradnice bodov a číselné rozmery mriežky sa nemenia. Existujúci
výpočtový postup stále predpokladá metrické súradnice a dlaždice pomenované podľa
siete 2 × 2 km; zachovanie CRS samo osebe neprispôsobuje algoritmus stupňom či stopám.
Existujúce dátové TIFF sa automaticky neopravujú.

Regresné testy CRS vyžadujú CMake 3.22 alebo novší a možno ich spustiť
z priečinka `src`:

```powershell
cmake -S tests -B .build-crs-tests
cmake --build .build-crs-tests --config Release --target crs_regression crop_crs_regression -- /m:1 /nr:false
ctest --test-dir .build-crs-tests -C Release -R crs_regression --output-on-failure
```
