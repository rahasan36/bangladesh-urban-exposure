# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

import os, sys, json, numpy as np, pandas as pd, rasterio, geopandas as gpd
from rasterio.features import rasterize
from rasterio.transform import Affine
from shapely.geometry import box
_R = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, os.path.join(_R, "pipeline")); sys.path.insert(0, _R)
from common import row_areas_km2

from config import HARM, OUT as P26, MR as OUT, GHSL_DIR
FACTORS = {1: "30 m", 3: "90 m", 10: "300 m", 17: "510 m", 33: "990 m"}
YEARS_GLC = [2000, 2005, 2010, 2015, 2020]
YEARS_GHSL = [2000, 2005, 2010, 2015, 2020]
UTM = 32646
BUFFERS = [("core", 0, 0), ("r0_5", 0, 5000), ("r5_10", 5000, 10000), ("r10_25", 10000, 25000)]

def blocksum(a, f):
    H, W = a.shape; Hp, Wp = -(-H // f) * f, -(-W // f) * f
    if (Hp, Wp) != (H, W):
        b = np.zeros((Hp, Wp), dtype=np.uint32); b[:H, :W] = a; a = b
    return a.reshape(Hp // f, f, Wp // f, f).sum(axis=(1, 3), dtype=np.uint32)

def block_transform(t, f): return Affine(t.a * f, t.b, t.c, t.d, t.e * f, t.f)

def aggregate_glc(arr, f):
    land = blocksum((arr > 0).astype(np.uint8), f)
    urb = blocksum((arr == 1).astype(np.uint8), f)
    best = urb.copy(); arg = np.where(urb > 0, 1, 0).astype(np.uint8)
    for c in range(2, 11):
        cnt = blocksum((arr == c).astype(np.uint8), f)
        sel = cnt > best; best[sel] = cnt[sel]; arg[sel] = c
    arg[land == 0] = 0
    return arg, urb, land

def zones_gdf():
    pa = gpd.read_file(f"{P26}/BGD_PA_Jan2026.geojson").to_crs(UTM)
    kb = gpd.read_file(f"{P26}/BGD_KBA_2021.geojson").to_crs(UTM)
    out = []
    for t, g in [("PA", pa.geometry.union_all()), ("KBA", kb.geometry.union_all())]:
        for n, d0, d1 in BUFFERS:
            z = g if n == "core" else g.buffer(d1).difference(g.buffer(d0) if d0 > 0 else g)
            out.append(dict(feature_type=t, zone=n, geometry=z))
    return gpd.GeoDataFrame(out, crs=UTM)

def zonal_share(mask_geom_gdf, transform, shape, crs, urban_bool, land_bool, cell_area):
    g = mask_geom_gdf.to_crs(crs)
    rows = []
    if cell_area.ndim == 1: cell_area = np.repeat(cell_area[:, None], shape[1], axis=1)
    for i, r in g.iterrows():
        m = rasterize([(r.geometry, 1)], out_shape=shape, transform=transform, fill=0, dtype="uint8").astype(bool)
        L = (cell_area * (m & land_bool)).sum(); U = (cell_area * (m & land_bool & urban_bool)).sum()
        rows.append(dict(**{k: v for k, v in r.items() if k != "geometry"}, land_km2=L, urban_km2=U, urban_pct=100 * U / L if L > 0 else np.nan))
    return pd.DataFrame(rows)

def species_share(sp, transform, shape, crs, urban_val, land_bool, cell_area, frac=None):
    g = sp.to_crs(crs); H, W = shape
    if cell_area.ndim == 1: cell_area = np.repeat(cell_area[:, None], W, axis=1)
    rows = []
    for _, r in g.iterrows():
        minx, miny, maxx, maxy = r.geometry.bounds
        c0, r1 = ~transform * (minx, miny); c1, r0 = ~transform * (maxx, maxy)
        r0, r1, c0, c1 = max(int(np.floor(r0)), 0), min(int(np.ceil(r1)) + 1, H), max(int(np.floor(c0)), 0), min(int(np.ceil(c1)) + 1, W)
        if r1 <= r0 or c1 <= c0: rows.append(dict(id_no=r.id_no, land_km2=0, urban_maj_km2=0, urban_frac_km2=np.nan)); continue
        wt = transform * Affine.translation(c0, r0)
        m = rasterize([(r.geometry, 1)], out_shape=(r1 - r0, c1 - c0), transform=wt, fill=0, dtype="uint8", all_touched=False).astype(bool)
        if not m.any():
            m = rasterize([(r.geometry, 1)], out_shape=(r1 - r0, c1 - c0), transform=wt, fill=0, dtype="uint8", all_touched=True).astype(bool)
        A = cell_area[r0:r1, c0:c1]; Lm = m & land_bool[r0:r1, c0:c1]
        L = (A * Lm).sum(); U = (A * (Lm & urban_val[r0:r1, c0:c1])).sum()
        Uf = (A * Lm * frac[r0:r1, c0:c1]).sum() if frac is not None else np.nan
        rows.append(dict(id_no=r.id_no, land_km2=L, urban_maj_km2=U, urban_frac_km2=Uf))
    d = pd.DataFrame(rows)
    d["urban_maj_pct"] = 100 * d.urban_maj_km2 / d.land_km2; d["urban_frac_pct"] = 100 * d.urban_frac_km2 / d.land_km2
    return d
