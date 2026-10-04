# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

import os, sys, numpy as np, pandas as pd, rasterio
from rasterio.warp import reproject, Resampling
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path[:0] = [REPO, os.path.join(REPO, "pipeline")]
from config import HARM, MR, MR_HARM, SDM_RASTERS, SDM_WORK, SDM_DATA
from common import row_areas_km2

ref = rasterio.open(os.path.join(SDM_RASTERS, "irreplacecable_areas_up2.tif")); T, C, SH = ref.transform, ref.crs, ref.shape
z = np.load(os.path.join(SDM_WORK, "urban_sdm_grid.npz")); L, CA = z["L"], z["CA"]
nval = np.load(os.path.join(SDM_WORK, "nval.npy")); irr = np.load(os.path.join(SDM_WORK, "irr.npy"))
meta = pd.read_csv(os.path.join(SDM_RASTERS, "targetRepresentation_taxonomicDetails_up3.csv"))
A = (nval == len(meta)) & (irr > -1e30) & (L / CA > 0.5)

def to_grid(a, src):
    o = np.zeros(SH, np.float64)
    reproject(a.astype(np.float64), o, src_transform=src.transform, src_crs=src.crs, dst_transform=T, dst_crs=C, resampling=Resampling.sum)
    return o
with rasterio.open(os.path.join(HARM, "LandArea300_GLC30.tif")) as r: land300 = np.where(r.read(1) < 0, 0, r.read(1)); land_src = r
def frac300(path):
    with rasterio.open(path) as r: return to_grid(np.clip(r.read(1), 0, None) * land300, r)
def glc_ladder(f, yr):
    with rasterio.open(os.path.join(MR, f"GLC_f{f}_{yr}_class.tif")) as r, rasterio.open(os.path.join(MR, f"GLC_f{f}_landfrac.tif")) as lf:
        px = row_areas_km2(r.transform, r.height, r.width)[:, None]
        return to_grid((r.read(1) == 1) * lf.read(1) * px, r)
def ghsl_ladder(f, yr):
    with rasterio.open(os.path.join(MR, f"GHSL_f{f}_{yr}_majcat20.tif")) as r, rasterio.open(os.path.join(MR, f"GHSL_f{f}_landfrac.tif")) as lf:
        return to_grid(r.read(1) * lf.read(1) * 0.01 * f * f, r)

import time
BUDGET = float(os.environ.get("BGD_BUDGET", 1e12)); t_start = time.time()
def out_of_time(margin=20): return time.time() - t_start > BUDGET - margin
CACHE = os.path.join(SDM_WORK, "urban_layers.npz")
U = {}
if os.path.exists(CACHE):
    c = np.load(CACHE); U = {tuple([k.rsplit("_", 1)[0], int(k.rsplit("_", 1)[1])]): c[k] for k in c.files}
LAYERS = [("GLC30", y, lambda y=y: frac300(os.path.join(HARM, f"UrbanFrac300_GLC30_{y}.tif"))) for y in (2000, 2004, 2008, 2012, 2016, 2020)] + \
         [("ESA300", y, lambda y=y: frac300(os.path.join(HARM, f"UrbanFrac300_ESA300_{y}.tif"))) for y in (2000, 2004, 2008, 2012, 2016, 2020)] + \
         [("GHSL100", y, lambda y=y: frac300(os.path.join(MR_HARM, f"UrbanFrac300_GHSL100_{y}.tif"))) for y in (2000, 2005, 2010, 2015, 2020)] + \
         [(n, y, lambda f=f, y=y: glc_ladder(f, y)) for y in (2000, 2020) for f, n in [(3, "GLCp90"), (10, "GLCp300"), (17, "GLCp510"), (33, "GLCp990")]] + \
         [(n, y, lambda f=f, y=y: ghsl_ladder(f, y)) for y in (2000, 2020) for f, n in [(3, "GHSLm300"), (5, "GHSLm500"), (10, "GHSLm1000")]]
for n, y, fn in LAYERS:
    if (n, y) in U: continue
    if out_of_time(40): np.savez(CACHE, **{f"{a}_{b}": v for (a, b), v in U.items()}); print("CONTINUE"); sys.exit(0)
    U[(n, y)] = fn()
np.savez(CACHE, **{f"{a}_{b}": v for (a, b), v in U.items()})
REPS = ["GLC30", "GLCp90", "GLCp300", "GLCp510", "GLCp990", "GHSL100", "GHSLm300", "GHSLm500", "GHSLm1000", "ESA300"]
TRAJ = [("GLC30", (2000, 2004, 2008, 2012, 2016, 2020)), ("ESA300", (2000, 2004, 2008, 2012, 2016, 2020)), ("GHSL100", (2000, 2005, 2010, 2015, 2020))]
PART = os.path.join(SDM_WORK, "sdm_species_partial.csv")
rows = pd.read_csv(PART).to_dict("records") if os.path.exists(PART) else []
seen = {r["species"] for r in rows} | set(open(PART + ".empty").read().split()) if os.path.exists(PART + ".empty") else {r["species"] for r in rows}
for sp in meta.species:
    if sp in seen: continue
    if out_of_time(): pd.DataFrame(rows).to_csv(PART, index=False); print("CONTINUE"); sys.exit(0)
    with rasterio.open(os.path.join(SDM_RASTERS, f"{sp}.tif")) as r: m = (r.read(1) == 1) & A
    h = L[m].sum()
    if h <= 0:
        open(PART + ".empty", "a").write(sp + "\n"); continue
    row = dict(species=sp, habitat_km2=h)
    for p in REPS:
        row[f"{p}_2000"] = 100 * U[(p, 2000)][m].sum() / h; row[f"{p}_2020"] = 100 * U[(p, 2020)][m].sum() / h; row[f"{p}_gain"] = row[f"{p}_2020"] - row[f"{p}_2000"]
    for p, ys in TRAJ:
        for y in ys: row[f"traj_{p}_{y}"] = 100 * U[(p, y)][m].sum() / h
    rows.append(row)
out = meta.merge(pd.DataFrame(rows), on="species")
out["thr"] = (out.threatStatus == "Threatened").astype(int); out["log_area"] = np.log10(out.habitat_km2)
os.makedirs(SDM_DATA, exist_ok=True); out.to_csv(os.path.join(SDM_DATA, "sdm_species.csv"), index=False)
print(f"taxa with habitat {len(out)} of {len(meta)}; domain {L[A].sum():,.0f} km2")
