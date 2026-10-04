# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

import os, sys, numpy as np, pandas as pd, rasterio
from rasterio.warp import reproject, Resampling
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, REPO)
from config import HARM, MR_HARM, SDM_RASTERS, SDM_WORK

ref = rasterio.open(os.path.join(SDM_RASTERS, "irreplacecable_areas_up2.tif")); T, C, SH = ref.transform, ref.crs, ref.shape
irr = ref.read(1); pri = rasterio.open(os.path.join(SDM_RASTERS, "spatial_prioritization_up.tif")).read(1)

def to_grid(arr, src):
    o = np.zeros(SH, np.float64)
    reproject(arr.astype(np.float64), o, src_transform=src.transform, src_crs=src.crs, dst_transform=T, dst_crs=C, resampling=Resampling.sum)
    return o

rl = rasterio.open(os.path.join(HARM, "LandArea300_GLC30.tif")); land = np.where(rl.read(1) < 0, 0, rl.read(1))
out = {"L": to_grid(land, rl)}
R = 6371.0072; lat = T.f + T.e * np.arange(SH[0] + 1); dlon = np.radians(T.a)
band = R ** 2 * dlon * np.abs(np.sin(np.radians(lat[:-1])) - np.sin(np.radians(lat[1:])))
out["CA"] = np.repeat(band[:, None], SH[1], axis=1)
for p, d in [("GLC30", HARM), ("GHSL100", MR_HARM), ("ESA300", HARM)]:
    for y in (2000, 2020):
        with rasterio.open(os.path.join(d, f"UrbanFrac300_{p}_{y}.tif")) as r:
            out[f"{p}_{y}"] = to_grid(np.clip(r.read(1), 0, None) * land, r)
np.savez_compressed(os.path.join(SDM_WORK, "urban_sdm_grid.npz"), **out)

meta = pd.read_csv(os.path.join(SDM_RASTERS, "targetRepresentation_taxonomicDetails_up3.csv"))
nval = np.zeros(SH, np.int16); rich_t = np.zeros(SH, np.int16); rich_a = np.zeros(SH, np.int16)
for sp, ts in zip(meta.species, meta.threatStatus):
    with rasterio.open(os.path.join(SDM_RASTERS, f"{sp}.tif")) as r:
        a = r.read(1); ok = a != r.nodata; pres = ok & (a == 1)
    nval += ok; rich_a += pres
    if ts == "Threatened": rich_t += pres
for k, v in dict(nval=nval, rich_t=rich_t, rich_a=rich_a, irr=irr, pri=pri).items():
    np.save(os.path.join(SDM_WORK, f"{k}.npy"), v)
A = (nval == len(meta)) & (irr > -1e30) & (out["L"] / out["CA"] > 0.5)
print(f"models {len(meta)}, domain {out['L'][A].sum():,.0f} km2 of {out['L'].sum():,.0f} km2 land, GLC30 2020 urban {out['GLC30_2020'].sum():,.0f} km2")

np.savez_compressed(os.path.join(SDM_WORK, "urban_on_sdm_grid.npz"), irr=irr, pri=pri, cellkm=out["CA"], GLC30_2020_l=out["L"],
                    **{f"{p}_{y}_u": out[f"{p}_{y}"] for p in ("GLC30", "GHSL100", "ESA300") for y in (2000, 2020)})
np.savez_compressed(os.path.join(SDM_WORK, "urban_sdm_layers.npz"), **{k: v for k, v in out.items() if k not in ("L", "CA")})
np.save(os.path.join(SDM_WORK, "valid.npy"), nval == len(meta))
np.save(os.path.join(SDM_WORK, "rich_threatened.npy"), rich_t); np.save(os.path.join(SDM_WORK, "rich_all.npy"), rich_a)
