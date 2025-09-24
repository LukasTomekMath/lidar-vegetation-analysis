import subprocess
import os
from osgeo import gdal

# cesta k vstupnému TIFFu
current_file_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(current_file_dir)
input_tif = "merged.tif"

# otvorenie súboru a zistenie počtu pásiem
ds = gdal.Open(input_tif)
num_bands = ds.RasterCount

# pre každé pásmo vytvoriť samostatný COG
for i in range(1, num_bands + 1):
    band = ds.GetRasterBand(i)
    desc = band.GetDescription() or f"band{i}"

    output_tif = f"band{i:02d}_{desc}.tif"
    cmd = [
        "gdal_translate",
        "-b", str(i),
        "-of", "COG",        
        "-co", "COMPRESS=ZSTD",
        input_tif,
        output_tif,
        '-stats'  # pridá výpočet štatistík a uloženie do TIFF
    ]
    print(f"Vytváram {output_tif}...")
    subprocess.run(cmd, check=True)
