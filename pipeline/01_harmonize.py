# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

import os, sys, json, argparse
import numpy as np, rasterio, netCDF4, geopandas as gpd
from rasterio.transform import from_origin, Affine
from rasterio.features import rasterize
from common import *

ap = argparse.ArgumentParser(); ap.add_argument("--budget", type=float, default=float(os.environ.get("BGD_BUDGET", 1e12))); args = ap.parse_args()
B = Budget(args.budget)
PROG = os.path.join(HARM, "progress_01.json")
prog = json.load(open(PROG)) if os.path.exists(PROG) else {}
def save(): json.dump(prog, open(PROG, "w"))

ESA_LUT = make_lut(ESA_XWALK); GLC_LUT = make_lut(GLC_XWALK)
PROFILE = dict(driver="GTiff", dtype="uint8", count=1, nodata=0, compress="LZW", tiled=True,
               blockxsize=512, blockysize=512, BIGTIFF="IF_SAFER")

bnd = gpd.read_file(BOUNDARY_SHP).to_crs(4326)
bnd_geom = [bnd.geometry.union_all()]

def build_mask(path, transform, height, width):
    if os.path.exists(path): return
    m = rasterize(bnd_geom, out_shape=(height, width), transform=transform, fill=0, default_value=1, dtype="uint8", all_touched=False)
    with rasterio.open(path, "w", crs="EPSG:4326", transform=transform, height=height, width=width, **PROFILE) as d:
        d.write(m, 1)
    log(f"mask written {path}: {int(m.sum())} pixels inside")

ds = netCDF4.Dataset(os.path.join(ESA_DIR, ESA_FILES[1992]))
lat = ds["lat"][:]; lon = ds["lon"][:]
dlat = abs(float(lat[1] - lat[0])); dlon = abs(float(lon[1] - lon[0]))
esa_transform = from_origin(float(lon[0]) - dlon / 2, float(lat[0]) + dlat / 2, dlon, dlat)
esa_h, esa_w = len(lat), len(lon); ds.close()
build_mask(os.path.join(HARM, "mask_esa300.tif"), esa_transform, esa_h, esa_w)

early = GLC_EARLY
late = GLC_LATE
with rasterio.open(early) as r:
    glc_transform, glc_h, glc_w = r.transform, r.height, r.width
build_mask(os.path.join(HARM, "mask_glc30.tif"), glc_transform, glc_h, glc_w)

with rasterio.open(os.path.join(HARM, "mask_esa300.tif")) as r: esa_mask = r.read(1)
for yr in ESA_YEARS:
    out = os.path.join(HARM, f"ESA300_{yr}.tif")
    if os.path.exists(out): continue
    ds = netCDF4.Dataset(os.path.join(ESA_DIR, ESA_FILES[yr]))
    a = np.asarray(ds["lccs_class"][0, :, :]).astype(np.uint8); ds.close()
    h = ESA_LUT[a] * esa_mask
    with rasterio.open(out, "w", crs="EPSG:4326", transform=esa_transform, height=esa_h, width=esa_w, **PROFILE) as d:
        d.write(h, 1)
    log(f"ESA {yr} harmonized; classes present {np.unique(h[h>0]).tolist()}")

BLOCK = 1024
with rasterio.open(os.path.join(HARM, "mask_glc30.tif")) as r: mask_src = r
for yr in GLC_YEARS:
    out = os.path.join(HARM, f"GLC30_{yr}.tif")
    key = f"GLC30_{yr}"
    if prog.get(key) == "done": continue
    src_path, band = (early, [1985, 1990, 1995].index(yr) + 1) if yr < 2000 else (late, yr - 1999)
    if not os.path.exists(out):
        with rasterio.open(out, "w", crs="EPSG:4326", transform=glc_transform, height=glc_h, width=glc_w, **PROFILE) as d: pass
        prog[key] = 0; save()
    row0 = prog[key]
    with rasterio.open(src_path) as src, rasterio.open(out, "r+") as dst, rasterio.open(os.path.join(HARM, "mask_glc30.tif")) as mk:
        while row0 < glc_h:
            if not B.ok(25):
                log(f"budget exhausted at {key} row {row0}"); save(); print("CONTINUE"); sys.exit(0)
            nrows = min(BLOCK, glc_h - row0)
            win = rasterio.windows.Window(0, row0, glc_w, nrows)
            a = src.read(band, window=win)
            a[a > 255] = 0
            h = GLC_LUT[a.astype(np.uint8)] * mk.read(1, window=win)
            dst.write(h, 1, window=win)
            row0 += nrows; prog[key] = row0; save()
    prog[key] = "done"; save(); log(f"GLC {yr} harmonized (band {band} of {os.path.basename(src_path)})")

print("DONE")
