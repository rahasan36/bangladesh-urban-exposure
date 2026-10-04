# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

import os, sys
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path[:0] = [REPO, os.path.join(REPO, "figures"), os.path.join(REPO, "sdm")]
from config import HARM, MR_TAB, TAB, OUT_FIG, SDM_RASTERS, SDM_WORK, SDM_DATA, REV, SUPP
os.chdir(SDM_WORK)
import os, numpy as np, pandas as pd, rasterio
from rasterio.warp import reproject, Resampling
with rasterio.open(f"{HARM}/LandArea300_GLC30.tif") as r: LAND = r.read(1); T, C, H, W = r.transform, r.crs, r.height, r.width
LAND = np.where(LAND < 0, 0, LAND)
hot = np.zeros((H, W), np.uint8)
with rasterio.open("hotspots_sdm_GLC30.tif") as s:
    reproject(s.read(1), hot, src_transform=s.transform, src_crs=s.crs, dst_transform=T, dst_crs=C, resampling=Resampling.nearest, src_nodata=255, dst_nodata=0)
hot = hot == 1; mask = LAND > 0
rd = lambda f: rasterio.open(os.path.join(HARM, f)).read(1)
R = np.clip(rd("GHSL_frac300_2020.tif"), 0, 1)
prods = {"GLC30_fraction": np.clip(rd("UrbanFrac300_GLC30_2020.tif"), 0, None).astype("float32"),
         "GLC30_plurality": (rd("GLC300maj_2020.tif") == 1).astype("float32"), "ESA300": (rd("ESA300_2020.tif") == 1).astype("float32")}
rows = []
for name, P in prods.items():
    for stratum, m in (("national", mask), ("hotspot", mask & hot), ("outside_hotspot", mask & ~hot)):
        p, r, l = P[m], R[m], LAND[m]; ref = float((r * l).sum())
        row = dict(product=name, stratum=stratum, land_km2=float(l.sum()), ref_km2=ref, prod_km2=float((p * l).sum()), capture_any=float((r * l)[p > 0].sum()) / ref, pearson_r=float(np.corrcoef(p, r)[0, 1]))
        for t in (0.05, 0.25):
            rb = r >= t; pb = (p >= t) if name == "GLC30_fraction" else (p > 0.5)
            tp = float((rb & pb).sum()); fp = float((~rb & pb).sum()); fn = float((rb & ~pb).sum())
            row[f"producers_t{t}"] = tp / (tp + fn); row[f"users_t{t}"] = tp / (tp + fp)
        rows.append(row)
d = pd.DataFrame(rows); d.to_csv(os.path.join(SDM_DATA, "validation_sdm_hotspots_2020.csv"), index=False); print(d.round(3).to_string(index=False))
