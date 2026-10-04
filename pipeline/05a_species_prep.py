# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

import os, sys, time, numpy as np, pandas as pd, geopandas as gpd
from common import *
SP = os.path.join(OUT, "Species")
t0 = time.time()
bnd = gpd.read_file(BOUNDARY_SHP).to_crs(4326).geometry.union_all()
files = IUCN_FILES
BB = (88, 20.5, 92.8, 26.7)
parts = {}
for taxon, fl in files.items():
    for f in fl:
        cache = os.path.join(SP, "part_" + os.path.basename(f).split(".")[0] + ".geojson")
        if os.path.exists(cache):
            parts.setdefault(taxon, []).append(gpd.read_file(cache)); continue
        if time.time() - t0 > float(os.environ.get("BGD_BUDGET", 1e12)): print("CONTINUE"); sys.exit(0)
        src = f if os.path.isabs(f) else os.path.join(ROOT, f)
        import struct
        with open(src, "rb") as fh: declared = struct.unpack(">i", fh.read(28)[24:28])[0] * 2
        assert os.path.getsize(src) == declared, f"{src}: declared {declared} bytes, on disk {os.path.getsize(src)} (truncated download?)"
        g = gpd.read_file(src, bbox=BB)
        g["n_polys_bbox"] = 1
        g = g[g.presence.isin([1, 2]) & g.origin.isin([1, 2]) & g.seasonal.isin([1, 2, 3])].copy()
        g["poly_global_km2"] = g.geometry.to_crs(6933).area / 1e6
        g["geometry"] = g.geometry.clip_by_rect(*BB).make_valid()
        g = g[~g.geometry.is_empty]
        g = g[g.intersects(bnd)]
        g["geometry"] = g.geometry.intersection(bnd)
        g = g[~g.geometry.is_empty]
        g["taxon"] = taxon
        g = g[["id_no", "sci_name", "category", "order_", "family", "marine", "freshwater", "taxon", "poly_global_km2", "geometry"]]
        g.to_file(cache, driver="GeoJSON")
        log(f"{os.path.basename(f)}: {len(g)} polygons, {g.id_no.nunique()} species after filters, clipped to Bangladesh ({time.time()-t0:.0f}s)")
        parts.setdefault(taxon, []).append(g)
allsp = []
for taxon in files:
    g = pd.concat(parts[taxon], ignore_index=True)
    attrs = g.groupby("id_no").agg(sci_name=("sci_name", "first"), category=("category", "first"), order_=("order_", "first"), family=("family", "first"),
                                   marine=("marine", "first"), freshwater=("freshwater", "first"), range_global_km2=("poly_global_km2", "sum")).reset_index()
    dis = g[["id_no", "geometry"]].dissolve(by="id_no").reset_index().merge(attrs, on="id_no")
    dis["taxon"] = taxon
    dis["range_bgd_km2"] = dis.geometry.to_crs(6933).area / 1e6
    allsp.append(dis)
sp = gpd.GeoDataFrame(pd.concat(allsp, ignore_index=True), crs=4326)
grp = {"CR": "Threatened", "EN": "Threatened", "VU": "Threatened", "NT": "Near Threatened", "LC": "Least Concern", "DD": "Data Deficient"}
sp["status_group"] = sp.category.map(grp).fillna("Other")
sp["range_bgd_share"] = sp.range_bgd_km2 / sp.range_global_km2
sp.drop(columns="geometry").to_csv(os.path.join(SP, "species_table.csv"), index=False)
sp.to_file(os.path.join(SP, "species_ranges_bgd.geojson"), driver="GeoJSON")
print(sp.groupby(["taxon", "status_group"]).size().unstack(fill_value=0))
print(sp.groupby(["taxon", "category"]).size().unstack(fill_value=0))
print("DONE")
