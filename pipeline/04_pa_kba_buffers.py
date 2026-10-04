# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

import os, sys, argparse, json
import numpy as np, rasterio, pandas as pd, geopandas as gpd
from rasterio.features import rasterize
from rasterio.windows import from_bounds, Window
from shapely.geometry import box
from common import *

ap = argparse.ArgumentParser(); ap.add_argument("--budget", type=float, default=float(os.environ.get("BGD_BUDGET", 1e12))); args = ap.parse_args()
B = Budget(args.budget)
K = NCLASS + 1
UTM = 32646
BUFFERS = [("core", 0, 0), ("r0_5", 0, 5000), ("r5_10", 5000, 10000), ("r10_25", 10000, 25000)]
YEARS = {"GLC30": GLC_YEARS, "ESA300": ESA_YEARS}

def load_features():
    pa = gpd.read_file(PA_SHP).to_crs(UTM)
    pa = pa.rename(columns={"SITE_ID": "feature_id", "NAME_ENG": "name", "IUCN_CAT": "iucn_cat", "DESIG_ENG": "designation", "REALM": "realm", "STATUS_YR": "status_yr", "GIS_AREA": "gis_area_km2"})
    pa["feature_type"] = "PA"
    kb = gpd.read_file(KBA_SHP, bbox=(88, 20.5, 92.8, 26.7))
    kb = kb[kb.ISO3 == "BGD"].to_crs(UTM)
    kb = kb.rename(columns={"SitRecID": "feature_id", "IntName": "name", "SitArea": "gis_area_ha"})
    kb["iucn_cat"] = "KBA"; kb["designation"] = "KBA"; kb["realm"] = ""; kb["status_yr"] = pd.to_datetime(kb.AddedDate).dt.year
    kb["gis_area_km2"] = kb.geometry.area / 1e6; kb["feature_type"] = "KBA"
    cols = ["feature_type", "feature_id", "name", "iucn_cat", "designation", "realm", "status_yr", "gis_area_km2", "geometry"]
    feats = pd.concat([pa[cols], kb[cols]], ignore_index=True)
    feats["geometry"] = feats.geometry.make_valid()
    feats["feature_id"] = feats.feature_id.astype(str)
    return gpd.GeoDataFrame(feats, crs=UTM)

feats = load_features()
meta_path = os.path.join(TAB, "pa_kba_features.csv")
if not os.path.exists(meta_path):
    feats.drop(columns="geometry").to_csv(meta_path, index=False)
    log(f"features: {int((feats.feature_type=='PA').sum())} PAs, {int((feats.feature_type=='KBA').sum())} KBAs")
cores = {t: feats[feats.feature_type == t].geometry.union_all() for t in ["PA", "KBA"]}

def zones_for(row):
    g = row.geometry; allcores = cores[row.feature_type]
    out = [("core", g)]
    for name, d0, d1 in BUFFERS[1:]:
        buf = g.buffer(d1).difference(g.buffer(d0) if d0 > 0 else g).difference(allcores)
        out.append((name, buf))
    return out

def estate_rows():
    rows = []
    for t in ["PA", "KBA"]:
        g = cores[t]
        rows.append(dict(feature_type=t, feature_id="ALL", name=f"All {t}s (dissolved)", zones=[("core", g)] +
                         [(n, g.buffer(d1).difference(g.buffer(d0) if d0 > 0 else g)) for n, d0, d1 in BUFFERS[1:]]))
    return rows

def zonal(ds, zones_utm):
    res = []
    zdf = gpd.GeoDataFrame(dict(zone=[z[0] for z in zones_utm]), geometry=[z[1] for z in zones_utm], crs=UTM).to_crs(4326)
    zdf = zdf[~zdf.geometry.is_empty]
    if len(zdf) == 0: return res
    with rasterio.open(os.path.join(HARM, f"{ds}_{YEARS[ds][0]}.tif")) as r:
        b = zdf.total_bounds
        b = [max(b[0], r.bounds.left), max(b[1], r.bounds.bottom), min(b[2], r.bounds.right), min(b[3], r.bounds.top)]
        if b[0] >= b[2] or b[1] >= b[3]: return res
        win = from_bounds(*b, transform=r.transform).round_offsets().round_lengths()
        win = Window(max(0, win.col_off), max(0, win.row_off), min(win.width, r.width - max(0, win.col_off)), min(win.height, r.height - max(0, win.row_off)))
        wt = r.window_transform(win)
        ra = row_areas_km2(wt, int(win.height), int(win.width))
        shapes = [(geom, i + 1) for i, geom in enumerate(zdf.geometry)]
        zid = rasterize(shapes, out_shape=(int(win.height), int(win.width)), transform=wt, fill=0, dtype="uint8", all_touched=False)
    if zid.max() == 0: return res
    W = np.repeat(ra, int(win.width)).reshape(int(win.height), int(win.width))
    names = list(zdf.zone)
    for yr in YEARS[ds]:
        with rasterio.open(os.path.join(HARM, f"{ds}_{yr}.tif")) as r:
            a = r.read(1, window=win).astype(np.int64)
        acc = np.bincount((zid.astype(np.int64) * K + a).ravel(), weights=W.ravel(), minlength=(len(names) + 1) * K)
        for zi, zn in enumerate(names):
            row = acc[(zi + 1) * K:(zi + 2) * K]
            for c in range(1, K):
                res.append(dict(zone=zn, year=yr, class_id=c, class_name=CLASS_NAMES[c], area_km2=row[c]))
    return res

for ds in ["GLC30", "ESA300"]:
    out = os.path.join(TAB, f"pa_kba_zonal_{ds}.csv")
    done = set(pd.read_csv(out).feature_key.unique()) if os.path.exists(out) else set()
    units = [dict(feature_type=r.feature_type, feature_id=r.feature_id, name=r.name, zones=zones_for(r)) for r in feats.itertuples()]
    for u in units:
        key = f"{u['feature_type']}:{u['feature_id']}"
        if key in done: continue
        if not B.ok(30): print("CONTINUE"); sys.exit(0)
        recs = zonal(ds, u["zones"])
        df = pd.DataFrame(recs)
        if len(df) == 0:
            df = pd.DataFrame([dict(zone="core", year=YEARS[ds][0], class_id=0, class_name="none", area_km2=0.0)])
        df.insert(0, "feature_key", key); df.insert(1, "feature_type", u["feature_type"]); df.insert(2, "feature_id", u["feature_id"]); df.insert(3, "name", u["name"]); df.insert(4, "dataset", ds)
        df.to_csv(out, mode="a", header=not os.path.exists(out), index=False)
        done.add(key)
        urb = df[(df.zone == "core") & (df.class_id == 1) & (df.year == YEARS[ds][-1])].area_km2.sum()
        tot = df[(df.zone == "core") & (df.year == YEARS[ds][-1])].area_km2.sum()
        log(f"{ds} {key} {u['name'][:30]}: core land {tot:.1f} km2, urban {urb:.2f} km2 in {YEARS[ds][-1]}")

BLOCK = 2048
for ds in ["GLC30", "ESA300"]:
    out = os.path.join(TAB, f"pa_kba_zonal_{ds}.csv")
    existing = pd.read_csv(out) if os.path.exists(out) else pd.DataFrame(columns=["feature_key", "year"])
    for u in estate_rows():
        key = f"{u['feature_type']}:ALL"
        zpath = os.path.join(HARM, f"zones_{u['feature_type']}_{ds}.tif")
        if not os.path.exists(zpath):
            if not B.ok(60): print("CONTINUE"); sys.exit(0)
            zdf = gpd.GeoDataFrame(dict(zone=[z[0] for z in u["zones"]]), geometry=[z[1] for z in u["zones"]], crs=UTM).to_crs(4326)
            with rasterio.open(os.path.join(HARM, f"{ds}_{YEARS[ds][0]}.tif")) as r:
                zid = rasterize([(g, i + 1) for i, g in enumerate(zdf.geometry)], out_shape=(r.height, r.width), transform=r.transform, fill=0, dtype="uint8")
                prof = r.profile; prof.update(dtype="uint8", nodata=0, compress="LZW")
            with rasterio.open(zpath, "w", **prof) as d: d.write(zid, 1)
            log(f"estate zone raster written {os.path.basename(zpath)}")
        done_years = set(existing[existing.feature_key == key].year.unique())
        names = [z[0] for z in u["zones"]]
        for yr in YEARS[ds]:
            if yr in done_years: continue
            if not B.ok(45): print("CONTINUE"); sys.exit(0)
            acc = np.zeros((len(names) + 1) * K)
            with rasterio.open(os.path.join(HARM, f"{ds}_{yr}.tif")) as r, rasterio.open(zpath) as z:
                ra = row_areas_km2(r.transform, r.height, r.width)
                for row0 in range(0, r.height, BLOCK):
                    n = min(BLOCK, r.height - row0); win = Window(0, row0, r.width, n)
                    a = r.read(1, window=win).astype(np.int64); zz = z.read(1, window=win).astype(np.int64)
                    acc += np.bincount((zz * K + a).ravel(), weights=np.repeat(ra[row0:row0 + n], r.width), minlength=(len(names) + 1) * K)[:(len(names) + 1) * K]
            recs = []
            for zi, zn in enumerate(names):
                row = acc[(zi + 1) * K:(zi + 2) * K]
                for c in range(1, K):
                    recs.append(dict(feature_key=key, feature_type=u["feature_type"], feature_id="ALL", name=u["name"], dataset=ds, zone=zn, year=yr, class_id=c, class_name=CLASS_NAMES[c], area_km2=row[c]))
            pd.DataFrame(recs).to_csv(out, mode="a", header=not os.path.exists(out), index=False)
            log(f"{ds} {key} {yr}: core urban {acc[K+1]:.1f} km2 of {acc[K+1:2*K].sum():.0f} km2")
print("DONE")
