# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

import os, sys, glob, numpy as np, rasterio, geopandas as gpd
from rasterio.merge import merge
from rasterio.mask import mask
from common import *
from config import GLC_TILE_DIR

bnd = gpd.read_file(BOUNDARY_SHP).to_crs(4326); bb = bnd.total_bounds
tiles = [t for t in glob.glob(os.path.join(GLC_TILE_DIR, "**", "*.tif"), recursive=True)]
def covers(path):
    with rasterio.open(path) as s: b = s.bounds
    return not (b.right < bb[0] or b.left > bb[2] or b.top < bb[1] or b.bottom > bb[3])
series = {"early": ("19852000", GLC_EARLY), "late": ("20002022", GLC_LATE)}
for key, (tag, out) in series.items():
    sel = sorted(t for t in tiles if tag in os.path.basename(t) and covers(t))
    print(f"{key}: {len(sel)} tiles ->", out)
    for t in sel: print("   ", os.path.basename(t))
    if "--dry-run" in sys.argv or not sel: continue
    if os.path.exists(out): print("   exists, skipping"); continue
    os.makedirs(os.path.dirname(out), exist_ok=True)
    srcs = [rasterio.open(t) for t in sel]
    arr, tr = merge(srcs, bounds=tuple(bb), nodata=0)
    prof = srcs[0].profile; [s.close() for s in srcs]
    prof.update(height=arr.shape[1], width=arr.shape[2], transform=tr, count=arr.shape[0], nodata=0, compress="deflate", tiled=True, BIGTIFF="YES")
    with rasterio.open(out, "w", **prof) as d: d.write(arr)
    with rasterio.open(out) as s:
        clipped, tr2 = mask(s, bnd.geometry, crop=True, nodata=0); prof2 = s.profile
    prof2.update(height=clipped.shape[1], width=clipped.shape[2], transform=tr2)
    with rasterio.open(out, "w", **prof2) as d: d.write(clipped)
    log(f"GLC mosaic {key}: {clipped.shape} written to {out}")
print("DONE")
