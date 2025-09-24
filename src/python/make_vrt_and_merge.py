import os
import subprocess
from osgeo import gdal

# Zozbieraj .tif súbory v aktuálnom priečinku
current_file_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(current_file_dir)
input_files = [f for f in os.listdir(current_file_dir) if f.endswith('.tif')]

# Zapíš zoznam do súboru pre gdalbuildvrt
with open('input_files.txt', 'w') as f:
    for file in input_files:
        f.write(f"{file}\n")

# Vytvor dočasný VRT
subprocess.run(['gdalbuildvrt', '-input_file_list', 'input_files.txt', 'merged.vrt'], check=True)

# Načítaj popisy bandov **len z prvého súboru**
band_descriptions = []
ds = gdal.Open(input_files[0])
for b in range(1, ds.RasterCount + 1):
    band_descriptions.append(ds.GetRasterBand(b).GetDescription())
ds = None

# Načítaj dočasný VRT a nastav popisy bandov
vrt_ds = gdal.Open('merged.vrt', gdal.GA_Update)
if vrt_ds is None:
    raise RuntimeError("Nepodarilo sa otvoriť merged.vrt")

for i, desc in enumerate(band_descriptions):
    band = vrt_ds.GetRasterBand(i+1)
    band.SetDescription(desc)

vrt_ds.FlushCache()
vrt_ds = None

# Preveď VRT na finálny TIFF s Float32 + kompresiou
subprocess.run([
    'gdal_translate',
    'merged.vrt',
    'merged.tif',
    '-ot', 'Float32',
    '-co', 'COMPRESS=ZSTD',
    '-co', 'PREDICTOR=2',
    '-co', 'BIGTIFF=YES',
    '-stats'  # pridá výpočet štatistík a uloženie do TIFF
], check=True)

# Vypočítaj štatistiky pre QGIS
#subprocess.run(['gdalinfo', '-stats', 'merged.tif'], check=True)
#subprocess.run(['gdal_edit.py', '-stats', 'merged.tif'], check=True)
