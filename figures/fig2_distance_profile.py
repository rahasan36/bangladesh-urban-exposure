# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

import os, sys, numpy as np, pandas as pd, geopandas as gpd, rasterio
from rasterio.features import rasterize
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path[:0] = [REPO]
from config import OUT, HARM, MR_HARM, TAB_MS
UTM = 32646
pa = gpd.read_file(os.path.join(OUT, "BGD_PA_Jan2026.geojson")).to_crs(UTM); pa["geometry"] = pa.geometry.make_valid()
core = pa.geometry.union_all()
with rasterio.open(os.path.join(HARM, "LandArea300_GLC30.tif")) as r: land = np.where(r.read(1) < 0, 0, r.read(1)).astype(float); T, crs, sh = r.transform, r.crs, r.shape
U = {}
for p, d in [("GLC30", HARM), ("GHSL100", MR_HARM), ("ESA300", HARM)]:
    with rasterio.open(os.path.join(d, f"UrbanFrac300_{p}_2020.tif")) as r: U[p] = np.clip(r.read(1), 0, None).astype(float)
def mask(g): return rasterize([(gpd.GeoSeries([g], crs=UTM).to_crs(crs).iloc[0], 1)], out_shape=sh, transform=T, fill=0, dtype="uint8").astype(bool)
bands = [("core", None, None, core)]
for i in range(50):
    d0, d1 = i * 0.5, (i + 1) * 0.5
    g = core.buffer(d1 * 1000).difference(core.buffer(d0 * 1000) if d0 > 0 else core)
    bands.append((f"{d0:g}-{d1:g}", d0, d1, g))
rows = []
for name, d0, d1, g in bands:
    m = mask(g); L = land[m].sum()
    for p in U: rows.append(dict(product=p, band=name, d0=d0, d1=d1, d_mid=None if d0 is None else (d0 + d1) / 2, land_km2=L, urban_pct=100 * (U[p][m] * land[m]).sum() / L))
out = pd.DataFrame(rows).sort_values(["product"], kind="stable"); os.makedirs(TAB_MS, exist_ok=True)
order = {"GLC30": 0, "GHSL100": 1, "ESA300": 2}; out["o"] = out["product"].map(order); out = out.sort_values("o", kind="stable").drop(columns="o")
out.to_csv(os.path.join(TAB_MS, "pa_distance_profile_2020.csv"), index=False); print(out.groupby("product").urban_pct.describe().round(2))
