# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

import os, numpy as np, pandas as pd, geopandas as gpd, rasterio
from rasterio.features import rasterize
from scipy.ndimage import uniform_filter
from common import *
SP = os.path.join(OUT, "Species")
with rasterio.open(os.path.join(HARM, "mask_esa300.tif")) as r:
    et = r.transform; EH, EW = r.height, r.width; crs = r.crs
F = 3
H, W = EH // F, EW // F
kt = rasterio.Affine(et.a * F, 0, et.c, 0, et.e * F, et.f)
def agg(a): return a[:H * F, :W * F].reshape(H, F, W, F).sum(axis=(1, 3))
def gistar(x, mask):
    x = np.where(mask, x, 0.0); n = mask.sum(); xbar = x[mask].mean(); s = x[mask].std()
    wsum = uniform_filter(x, 3, mode="constant") * 9
    wn = uniform_filter(mask.astype(float), 3, mode="constant") * 9
    denom = s * np.sqrt((n * wn - wn ** 2) / (n - 1))
    z = np.where(denom > 0, (wsum - xbar * wn) / np.where(denom > 0, denom, 1), 0)
    return np.where(mask, z, np.nan)
sp = gpd.read_file(os.path.join(SP, "species_ranges_bgd.geojson"))
rich = {}
for grp_name, sel in [("threatened", sp.status_group == "Threatened"), ("all", sp.status_group != "Other")]:
    acc = np.zeros((EH, EW), dtype=np.int16)
    for g in sp[sel].geometry:
        acc += rasterize([(g, 1)], out_shape=(EH, EW), transform=et, fill=0, dtype="int16")
    rich[grp_name] = agg(acc.astype(float)) / (F * F)
    log(f"richness {grp_name}: {int(sel.sum())} species, max per cell {rich[grp_name].max():.1f}")
adm = gpd.read_file(ADM1_SHP).to_crs(4326)
adm_r = rasterize([(g, i + 1) for i, g in enumerate(adm.geometry)], out_shape=(H, W), transform=kt, fill=0, dtype="uint8")
rows = []
for ds in ["GLC30", "ESA300"]:
    with rasterio.open(os.path.join(HARM, f"LandArea300_{ds}.tif")) as r: land = agg(r.read(1).astype(float))
    u = {}
    for yr in [2000, 2020]:
        with rasterio.open(os.path.join(HARM, f"UrbanFrac300_{ds}_{yr}.tif")) as r, rasterio.open(os.path.join(HARM, f"LandArea300_{ds}.tif")) as l:
            u[yr] = agg(r.read(1).astype(float) * l.read(1).astype(float))
    mask = land > 0.5 * np.nanmax(land) * 0.2
    mask = land > 0.1
    gain = u[2020] - u[2000]
    z_gain = gistar(gain, mask)
    expo = gain * rich["threatened"]
    z_expo = gistar(expo, mask)
    hs_urban = (z_gain > 2.58) & mask; hs_expo = (z_expo > 2.58) & mask
    cls = np.zeros((H, W), dtype=np.uint8); cls[hs_urban] = 1; cls[hs_expo] = 2; cls[hs_urban & hs_expo] = 3
    prof = dict(driver="GTiff", crs=crs, transform=kt, height=H, width=W, compress="LZW")
    fo = os.path.join(HARM, f"Hotspots_1km_{ds}.tif")
    with (rasterio.open(fo, "r+") if os.path.exists(fo) else rasterio.open(fo, "w", count=6, dtype="float32", nodata=-9999, **prof)) as d:
        out = np.stack([gain, z_gain, rich["threatened"], expo, z_expo, cls.astype(float)]).astype(np.float32)
        out[:, ~mask] = -9999; d.write(out)
        for i, n in enumerate(["urban_gain_km2_2000_2020", "gistar_z_urban_gain", "threatened_richness", "exposure_gain_x_richness", "gistar_z_exposure", "hotspot_class_1urban_2exposure_3both"]): d.set_band_description(i + 1, n)
    for i, name in enumerate(list(adm.shapeName) + ["Bangladesh"]):
        sel = (adm_r == i + 1) if i < len(adm) else mask
        sel = sel & mask
        rows.append(dict(dataset=ds, division=name, land_km2=land[sel].sum(), urban_gain_km2=gain[sel].sum(),
                         urban_hotspot_km2=land[sel & hs_urban].sum(), exposure_hotspot_km2=land[sel & hs_expo].sum(), coincident_km2=land[sel & hs_urban & hs_expo].sum(),
                         urban_gain_in_hotspots_km2=gain[sel & hs_urban].sum(), mean_threatened_richness=rich["threatened"][sel].mean()))
    log(f"{ds}: urban hotspot cells {hs_urban.sum()}, exposure hotspot cells {hs_expo.sum()}, coincident {(hs_urban & hs_expo).sum()} of {mask.sum()} land cells")
df = pd.DataFrame(rows); df.to_csv(os.path.join(TAB, "hotspots_by_division.csv"), index=False)
print(df.round(1).to_string(index=False)); print("DONE")
