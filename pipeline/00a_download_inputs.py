# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

import os, sys, struct, urllib.request, zipfile
from common import *

CHECK = "--check" in sys.argv
def have(rel): return os.path.exists(os.path.join(ROOT, rel))
def fetch(url, dest):
    dest = os.path.join(ROOT, dest)
    if os.path.exists(dest): print("  present ", os.path.basename(dest)); return dest
    if CHECK: print("  MISSING ", os.path.basename(dest), "<-", url); return None
    os.makedirs(os.path.dirname(dest), exist_ok=True); print("  download", url)
    urllib.request.urlretrieve(url, dest); return dest
def unzip(path, into):
    if path and path.endswith(".zip"):
        with zipfile.ZipFile(path) as z: z.extractall(os.path.join(ROOT, into))

print("BGD_DATA_ROOT =", ROOT)
print("\n[1] geoBoundaries BGD ADM1")
p = fetch("https://github.com/wmgeolab/geoBoundaries/raw/main/releaseData/gbOpen/BGD/ADM1/geoBoundaries-BGD-ADM1-all.zip", "geoBoundaries-BGD-ADM1-all.zip")
if p and not have("geoBoundaries-BGD-ADM1-all/geoBoundaries-BGD-ADM1.shp"): unzip(p, "geoBoundaries-BGD-ADM1-all")

print("\n[2] GHS-BUILT-S R2023A epochs 2000, 2015, 2020")
for yr in (2000, 2015, 2020):
    name = f"GHS_BUILT_S_E{yr}_GLOBE_R2023A_54009_100_V1_0"
    if os.path.exists(os.path.join(GHSL_DIR, name + ".tif")): print("  present ", name + ".tif"); continue
    p = fetch(f"https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_BUILT_S_GLOBE_R2023A/{name[:-5]}/V1-0/{name}.zip", f"GHSL/{name}.zip")
    unzip(p, "GHSL")

print("\n[3] GLC_FCS30D tiles for 85-95E, 20-30N")
print("  Zenodo: https://doi.org/10.5281/zenodo.8239305  ->  files 'GLC_FCS30D_19852000maps_E85-E90.zip', '..._E90-E95.zip' (5-year maps)")
print("          and 'GLC_FCS30D_20002022maps_E85-E90.zip', '..._E90-E95.zip' (annual maps); unzip into GLC_FCS30D_tiles/")
print("  Then run 00c_mosaic_glc.py to build Data/Processed/GLC_Early_BGD_30m.tif and GLC_Late_BGD_30m.tif")
print("  present" if have("Data/Processed/GLC_Late_BGD_30m.tif") else "  MISSING", "Data/Processed/GLC_Late_BGD_30m.tif")

print("\n[4] ESA CCI / C3S land cover, 300 m, area subset 88-93E 21-27N")
missing = [yr for yr, f in ESA_FILES.items() if not have(os.path.join("ESACCI", f))]
if not missing: print("  present  all 8 netCDF subsets")
else:
    print("  MISSING years", missing)
    print("  Request with cdsapi (pip install cdsapi; ~/.cdsapirc with your key), dataset 'satellite-land-cover',")
    print("  variable 'all', version 'v2_0_7cds' for 1992-2015 and 'v2_1_1' for 2016-2020, area [27, 88, 21, 93], format zip;")
    print("  unzip the .nc into ESACCI/ with the exact names listed in pipeline/common.py (ESA_FILES).")
    if not CHECK:
        try:
            import cdsapi
            c = cdsapi.Client()
            for yr in missing:
                ver = "v2_0_7cds" if yr <= 2015 else "v2_1_1"
                out = os.path.join(ROOT, "ESACCI", f"lc_{yr}.zip"); os.makedirs(os.path.dirname(out), exist_ok=True)
                c.retrieve("satellite-land-cover", {"variable": "all", "year": str(yr), "version": ver, "area": [27, 88, 21, 93], "format": "zip"}, out)
                unzip(out, "ESACCI")
        except Exception as ex:
            print("  cdsapi not available or request failed:", ex)

print("\n[5] Protected areas, KBAs, IUCN ranges (free registration required; cannot be fetched automatically)")
print("  WDPA (January 2026 release): https://www.protectedplanet.net/en/thematic-areas/wdpa  -> export Bangladesh polygons to BGD_PA.shp")
print("  KBAs (March 2021 release):    https://www.keybiodiversityareas.org/kba-data/request  -> KBAsGlobal_2021_March_01_POL.shp")
print("  IUCN Red List spatial data:   https://www.iucnredlist.org/resources/spatial-data-download -> MAMMALS_TERRESTRIAL_ONLY, REPTILES_PART1/2, AMPHIBIANS_PART1/2")
for f in ["BGD_PA.shp", "MAMMALS_TERRESTRIAL_ONLY.shp", "REPTILES_PART1.shp", "REPTILES_PART2.shp", "AMPHIBIANS_PART1.shp", "AMPHIBIANS_PART2.shp",
          "KBAsGlobal_2021_March_01-20210718T193202Z-001/KBAsGlobal_2021_March_01/KBAsGlobal_2021_March_01_POL.shp"]:
    ok = have(f); msg = "present " if ok else "MISSING "
    if ok:
        p = os.path.join(ROOT, f)
        with open(p, "rb") as fh: declared = struct.unpack(">i", fh.read(28)[24:28])[0] * 2
        if declared != os.path.getsize(p): msg = "TRUNCATED"
    print(" ", msg, f)

print("\n[6] GAUD annual urban extent 1985-2015, 30 m")
print("  https://doi.org/10.6084/m9.figshare.11513178  -> tiles urban_grid_116.tif, urban_grid_117.tif and grids_world.shp into GAUD/")
for f in ["urban_grid_116.tif", "urban_grid_117.tif", "grids_world.shp"]: print("  present " if os.path.exists(os.path.join(GAUD_DIR, f)) else "  MISSING ", os.path.join(GAUD_DIR, f))
print("\nDONE")
