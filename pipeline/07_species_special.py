# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

import os, sys, time, json, numpy as np, pandas as pd, geopandas as gpd, rasterio
from rasterio.features import rasterize, shapes as rshapes
from rasterio.windows import from_bounds, Window
from shapely.geometry import shape, box
from common import *
OUTD = os.path.join(OUT, "SpeciesSpecial"); os.makedirs(OUTD, exist_ok=True); t0 = time.time()
K = NCLASS + 1; UTM = 32646
def save_vec(gdf, path):
    if not os.path.exists(path): gdf.to_file(path, driver="GeoJSON")
sp = gpd.read_file(os.path.join(OUT, "Species", "species_ranges_bgd.geojson"))
pa = gpd.read_file(PA_SHP).to_crs(4326); kba = gpd.read_file(os.path.join(OUT, "BGD_KBA_2021.geojson"))
log(f"loaded {len(sp)} species ({time.time()-t0:.0f}s)")
sund = pa[pa.NAME_ENG.str.contains("Sundarban", case=False)].geometry.union_all()
hill_names = ["Sangu Matamuhari", "Pablakhali", "Kaptai", "Rampahar", "Hazarikhil", "Teknaf", "Himchari", "Inani", "Chunati", "Fasiakhali", "Dudpukuria", "Baraiyadhala", "Sangu", "Lawachara", "Rema-Kalenga", "Rajkandi", "Satchari", "Khadimnagar", "Medhakachhapia", "Hajarikhil"]
hill = pd.concat([pa[pa.NAME_ENG.str.contains("|".join(hill_names), case=False)].geometry, kba[kba.IntName.str.contains("|".join(hill_names), case=False)].geometry]).union_all()
haor = kba[kba.IntName.str.contains("Haor|Beel", case=False)].geometry.union_all()
river = kba[kba.IntName.str.contains("Jamuna|Ganges", case=False)].geometry.union_all()
regions = gpd.GeoDataFrame(dict(region=["Sundarbans", "Hill forest", "Haor wetlands", "Riverine"]), geometry=[sund, hill, haor, river], crs=4326)
save_vec(regions, os.path.join(OUTD, "regions.geojson"))
ru = regions.to_crs(UTM); env = {r.region: r.geometry.buffer(0).buffer(10000) for r in ru.itertuples()}
spu = sp.to_crs(UTM); spu["geometry"] = spu.geometry.make_valid()
rows = []
for r in spu.itertuples():
    a = r.geometry.area
    rows.append(dict(id_no=r.id_no, **{f"share_{k.split()[0].lower()}": (r.geometry.intersection(v).area / a if a > 0 else np.nan) for k, v in env.items()}))
pd.DataFrame(rows).to_csv(os.path.join(OUTD, "species_region_share.csv"), index=False)
log(f"region shares done ({time.time()-t0:.0f}s)")
sund_u = ru[ru.region == "Sundarbans"].geometry.iloc[0].buffer(0); allcores = pa.to_crs(UTM).geometry.make_valid().union_all()
zones = [("core", sund_u)] + [(n, sund_u.buffer(d1).difference(sund_u.buffer(d0) if d0 > 0 else sund_u).difference(allcores)) for n, d0, d1 in [("r0_5", 0, 5000), ("r5_10", 5000, 10000), ("r10_25", 10000, 25000)]]
zdf = gpd.GeoDataFrame(dict(zone=[z[0] for z in zones]), geometry=[z[1] for z in zones], crs=UTM).to_crs(4326)
save_vec(zdf, os.path.join(OUTD, "sundarbans_zones.geojson"))
recs = []
with rasterio.open(os.path.join(HARM, "GLC30_2020.tif")) as r:
    b = zdf.total_bounds; win = from_bounds(*b, transform=r.transform).round_offsets().round_lengths(); wt = r.window_transform(win)
    zid = rasterize([(g, i + 1) for i, g in enumerate(zdf.geometry)], out_shape=(int(win.height), int(win.width)), transform=wt, fill=0, dtype="uint8")
    ra = row_areas_km2(wt, int(win.height), int(win.width)); Wg = np.repeat(ra, int(win.width)).reshape(int(win.height), int(win.width))
for yr in GLC_YEARS:
    with rasterio.open(os.path.join(HARM, f"GLC30_{yr}.tif")) as r: a = r.read(1, window=win).astype(np.int64)
    acc = np.bincount((zid.astype(np.int64) * K + a).ravel(), weights=Wg.ravel(), minlength=5 * K)
    for zi, zn in enumerate(zdf.zone):
        row = acc[(zi + 1) * K:(zi + 2) * K]; recs.append(dict(zone=zn, year=yr, land_km2=row[1:].sum(), urban_km2=row[1], urban_pct=row[1] / max(row[1:].sum(), 1e-9) * 100))
pd.DataFrame(recs).to_csv(os.path.join(OUTD, "sundarbans_buffers_GLC30.csv"), index=False)
def urban_window(bb, name):
    with rasterio.open(os.path.join(HARM, "GLC30_2000.tif")) as r0, rasterio.open(os.path.join(HARM, "GLC30_2020.tif")) as r1:
        win = from_bounds(*bb, transform=r0.transform).round_offsets().round_lengths(); wt = r0.window_transform(win)
        a0 = r0.read(1, window=win); a1 = r1.read(1, window=win)
        out = np.zeros(a0.shape, dtype=np.uint8); out[a0 > 0] = 1; out[a1 == 1] = 3; out[(a0 == 1)] = 2
        out[(a0 == 1) & (a1 == 1)] = 2; out[(a0 != 1) & (a1 == 1)] = 3
        prof = dict(driver="GTiff", dtype="uint8", count=1, nodata=0, compress="LZW", crs=r0.crs, transform=wt, height=int(win.height), width=int(win.width))
        p_ = os.path.join(OUTD, f"{name}.tif")
        with (rasterio.open(p_, "r+") if os.path.exists(p_) else rasterio.open(p_, "w", **prof)) as d: d.write(out, 1)
urban_window((88.9, 21.5, 90.2, 23.1), "urban30_sundarbans_window")
log(f"sundarbans done ({time.time()-t0:.0f}s)")
narrow = sp[sp.range_bgd_km2 < 200].copy()
meta = []
for r in narrow.itertuples():
    b = r.geometry.bounds; pad = max((b[2] - b[0]), (b[3] - b[1])) * 0.25 + 0.01
    bb = (b[0] - pad, b[1] - pad, b[2] + pad, b[3] + pad)
    if not os.path.exists(os.path.join(OUTD, f"narrow_{r.id_no}.tif")): urban_window(bb, f"narrow_{r.id_no}")
    if time.time() - t0 > float(os.environ.get("BGD_BUDGET", 1e12)): print("CONTINUE"); sys.exit(0)
    with rasterio.open(os.path.join(OUTD, f"narrow_{r.id_no}.tif")) as d:
        a = d.read(1); m = rasterize([(r.geometry, 1)], out_shape=a.shape, transform=d.transform, fill=0, dtype="uint8", all_touched=True).astype(bool)
        ra = row_areas_km2(d.transform, d.height, d.width)
        land = (m & (a > 0)).sum(1) @ ra; u00 = (m & (a == 2)).sum(1) @ ra; u20 = (m & ((a == 2) | (a == 3))).sum(1) @ ra
    meta.append(dict(id_no=r.id_no, sci_name=r.sci_name, taxon=r.taxon, category=r.category, range_bgd_km2=r.range_bgd_km2, range_bgd_share=r.range_bgd_share, land_km2=land, urban2000_pct=u00 / max(land, 1e-9) * 100, urban2020_pct=u20 / max(land, 1e-9) * 100))
save_vec(narrow, os.path.join(OUTD, "narrow_ranges.geojson"))
pd.DataFrame(meta).to_csv(os.path.join(OUTD, "narrow_species_30m.csv"), index=False)
log(f"narrow-range: {len(meta)} species ({time.time()-t0:.0f}s)")
with rasterio.open(os.path.join(HARM, "Hotspots_1km_GLC30.tif")) as r:
    cls = r.read(6); land = r.read(1) > -9999; tr = r.transform; H, Wd = cls.shape
    ra = row_areas_km2(tr, H, Wd); Wg = np.repeat(ra, Wd).reshape(H, Wd)
rows = []
if os.path.exists(os.path.join(OUTD, "species_hotspot_overlap.csv")): print("DONE"); sys.exit(0)
for r in sp.itertuples():
    m = rasterize([(r.geometry, 1)], out_shape=(H, Wd), transform=tr, fill=0, dtype="uint8", all_touched=(r.range_bgd_km2 < 50)).astype(bool) & land
    tot = Wg[m].sum(); rows.append(dict(id_no=r.id_no, range_cells_km2=tot, coincident_hotspot_km2=Wg[m & (cls == 3)].sum(), urban_hotspot_km2=Wg[m & ((cls == 1) | (cls == 3))].sum()))
pd.DataFrame(rows).to_csv(os.path.join(OUTD, "species_hotspot_overlap.csv"), index=False)
log(f"hotspot overlap done ({time.time()-t0:.0f}s)"); print("DONE")
