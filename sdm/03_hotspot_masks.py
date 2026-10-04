# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

import os, sys
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path[:0] = [REPO, os.path.join(REPO, "figures"), os.path.join(REPO, "sdm")]
from config import HARM, MR_TAB, TAB, OUT_FIG, SDM_RASTERS, SDM_WORK, SDM_DATA, REV, SUPP
os.chdir(SDM_WORK)
import numpy as np, pandas as pd, rasterio
from scipy.ndimage import uniform_filter
z = np.load("urban_on_sdm_grid.npz"); valid = np.load("valid.npy"); U = np.load("urban_sdm_layers.npz"); irr = z["irr"]; CA = z["cellkm"]; L = z["GLC30_2020_l"]
A = valid & (irr > -1e30) & (L / CA > 0.5)
def gistar(x, m):
    x = np.where(m, x, 0.0).astype(float); n = m.sum(); xb = x[m].mean(); s = np.sqrt((x[m] ** 2).mean() - xb ** 2)
    sw = uniform_filter(x, 3, mode="constant") * 9; wi = uniform_filter(m.astype(float), 3, mode="constant") * 9
    with np.errstate(all="ignore"): zz = (sw - xb * wi) / (s * np.sqrt((n * wi - wi ** 2) / (n - 1)))
    return np.where(m, zz, np.nan)
src = rasterio.open(SDM_RASTERS + "/irreplacecable_areas_up2.tif"); prof = src.profile; prof.update(dtype="uint8", nodata=255, count=1, compress="lzw")
rows = []
for p in ["GLC30", "GHSL100", "ESA300"]:
    g = U[f"{p}_2020"] - U[f"{p}_2000"]; h = gistar(g, A) > 2.58
    out = np.full(h.shape, 255, np.uint8); out[A] = 0; out[h] = 1
    with rasterio.open(f"hotspots_sdm_{p}.tif", "w", **prof) as d: d.write(out, 1)
    rows.append(dict(product=p, hotspot_km2=L[h].sum(), gain_km2=g[A].sum(), gain_in_hotspots_km2=g[h].sum(), share_in_hotspots=g[h].sum() / g[A].sum(),
                     land_share_hotspots=L[h].sum() / L[A].sum()))
t = pd.DataFrame(rows); t.to_csv(os.path.join(SDM_DATA, "hotspots_sdm_gain_capture.csv"), index=False); print(t.round(3).to_string(index=False))
