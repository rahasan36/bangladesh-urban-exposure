# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

import os, sys, time, numpy as np, pandas as pd, geopandas as gpd, rasterio
from rasterio.features import rasterize
from common import *
SP = os.path.join(OUT, "Species"); t0 = time.time()
K = NCLASS + 1
with rasterio.open(os.path.join(HARM, "mask_esa300.tif")) as r:
    esa_mask = r.read(1).astype(bool); et = r.transform; EH, EW = r.height, r.width
    EPROF = dict(driver="GTiff", crs=r.crs, transform=et, height=EH, width=EW, compress="LZW", tiled=True, blockxsize=256, blockysize=256, count=1, dtype="float32", nodata=-9999)
with rasterio.open(os.path.join(HARM, "mask_glc30.tif")) as r:
    gt = r.transform; GH, GW = r.height, r.width
lat_c = gt.f + (np.arange(GH) + 0.5) * gt.e; lon_c = gt.c + (np.arange(GW) + 0.5) * gt.a
rowmap = np.floor((et.f - lat_c) / -et.e).astype(np.int64); colmap = np.floor((lon_c - et.c) / et.a).astype(np.int64)
valid_rows = (rowmap >= 0) & (rowmap < EH); rowmap = np.where(valid_rows, rowmap, 0)
gra = row_areas_km2(gt, GH, GW)
era = np.repeat(row_areas_km2(et, EH, EW), EW).reshape(EH, EW)

def wopen(path):
    return rasterio.open(path, "r+") if os.path.exists(path) else rasterio.open(path, "w", **EPROF)

land_glc = os.path.join(HARM, "LandArea300_GLC30.tif")
for yr in GLC_YEARS:
    fo = os.path.join(HARM, f"UrbanFrac300_GLC30_{yr}.tif")
    if os.path.exists(fo) and os.path.exists(land_glc): continue
    if time.time() - t0 > float(os.environ.get("BGD_BUDGET", 1e12)): print("CONTINUE"); sys.exit(0)
    urb = np.zeros(EH * EW + 1); land = np.zeros(EH * EW + 1)
    with rasterio.open(os.path.join(HARM, f"GLC30_{yr}.tif")) as r:
        for row0 in range(0, GH, 2048):
            n = min(2048, GH - row0)
            a = r.read(1, window=rasterio.windows.Window(0, row0, GW, n))
            key = (rowmap[row0:row0 + n][:, None] * EW + colmap[None, :]); key[~valid_rows[row0:row0 + n], :] = EH * EW
            key = key.ravel(); w = np.repeat(gra[row0:row0 + n], GW); af = a.ravel()
            land += np.bincount(key[af > 0], weights=w[af > 0], minlength=EH * EW + 1)
            urb += np.bincount(key[af == 1], weights=w[af == 1], minlength=EH * EW + 1)
    land = land[:-1].reshape(EH, EW); urb = urb[:-1].reshape(EH, EW)
    frac = np.where(land > 0, urb / np.maximum(land, 1e-12), 0).astype(np.float32)
    with wopen(fo) as d: d.write(frac, 1)
    if not os.path.exists(land_glc):
        with wopen(land_glc) as d: d.write(land.astype(np.float32), 1)
    log(f"UrbanFrac300 GLC30 {yr}: urban {urb.sum():.0f} km2 on {land.sum():.0f} km2 land")
for yr in ESA_YEARS:
    fo = os.path.join(HARM, f"UrbanFrac300_ESA300_{yr}.tif")
    if os.path.exists(fo): continue
    with rasterio.open(os.path.join(HARM, f"ESA300_{yr}.tif")) as r: a = r.read(1)
    with wopen(fo) as d: d.write((a == 1).astype(np.float32), 1)
land_esa = os.path.join(HARM, "LandArea300_ESA300.tif")
if not os.path.exists(land_esa):
    with wopen(land_esa) as d: d.write((era * esa_mask).astype(np.float32), 1)

sp = gpd.read_file(os.path.join(SP, "species_ranges_bgd.geojson"))
stack = {}
for ds, yrs in [("GLC30", GLC_YEARS), ("ESA300", ESA_YEARS)]:
    with rasterio.open(os.path.join(HARM, f"LandArea300_{ds}.tif")) as r: land = r.read(1)
    for yr in yrs:
        with rasterio.open(os.path.join(HARM, f"UrbanFrac300_{ds}_{yr}.tif")) as r: stack[(ds, yr)] = r.read(1) * land
    stack[(ds, "land")] = land
recs = []
for row in sp.itertuples():
    m = rasterize([(row.geometry, 1)], out_shape=(EH, EW), transform=et, fill=0, dtype="uint8", all_touched=False).astype(bool)
    touched = False
    if m.sum() < 10:
        m = rasterize([(row.geometry, 1)], out_shape=(EH, EW), transform=et, fill=0, dtype="uint8", all_touched=True).astype(bool); touched = True
    for ds, yrs in [("GLC30", GLC_YEARS), ("ESA300", ESA_YEARS)]:
        land = float(stack[(ds, "land")][m].sum())
        for yr in yrs:
            u = float(stack[(ds, yr)][m].sum())
            recs.append(dict(id_no=row.id_no, dataset=ds, year=yr, n_cells=int(m.sum()), all_touched=touched, range_land_km2=land, urban_km2=u, urban_pct=100 * u / land if land > 0 else np.nan))
df = pd.DataFrame(recs).merge(sp.drop(columns="geometry"), on="id_no")
df.to_csv(os.path.join(SP, "species_urban_exposure.csv"), index=False)
log(f"species exposure: {sp.shape[0]} species, {len(df)} rows ({time.time()-t0:.0f}s)")
s = df[(df.year.isin([2000, 2020]))].groupby(["dataset", "year", "taxon", "status_group"]).urban_pct.agg(["median", "mean", "count"]).round(2)
print(s.to_string()); print("DONE")
