# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

import os, sys, numpy as np, netCDF4
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from config import ESA_DIR as OUT
GLOBAL = os.environ.get("BGD_ESA_GLOBAL", "").split(os.pathsep)
def find(stem):
    for d in GLOBAL:
        if d and os.path.exists(os.path.join(d, stem + ".nc")): return d
    raise SystemExit(f"{stem}.nc not found in BGD_ESA_GLOBAL")
SRC = {1992: (None, "ESACCI-LC-L4-LCCS-Map-300m-P1Y-1992-v2.0.7cds"), 1996: (None, "ESACCI-LC-L4-LCCS-Map-300m-P1Y-1996-v2.0.7cds"),
       2000: (None, "ESACCI-LC-L4-LCCS-Map-300m-P1Y-2000-v2.0.7cds"), 2004: (None, "ESACCI-LC-L4-LCCS-Map-300m-P1Y-2004-v2.0.7cds"),
       2008: (None, "ESACCI-LC-L4-LCCS-Map-300m-P1Y-2008-v2.0.7cds"), 2012: (None, "ESACCI-LC-L4-LCCS-Map-300m-P1Y-2012-v2.0.7cds"),
       2016: (None, "C3S-LC-L4-LCCS-Map-300m-P1Y-2016-v2.1.1"), 2020: (None, "C3S-LC-L4-LCCS-Map-300m-P1Y-2020-v2.1.1")}
N, S, W, E = 27.0, 20.5, 88.0, 93.0
for yr in map(int, sys.argv[1:]):
    stem = SRC[yr][1]; d = find(stem); os.makedirs(OUT, exist_ok=True); out = os.path.join(OUT, f"{stem}.area-subset.27.93.20.5.88.nc")
    if os.path.exists(out): print(yr, "exists"); continue
    g = netCDF4.Dataset(os.path.join(d, stem + ".nc"))
    r0, r1 = int(round((90 - N) * 360)), int(round((90 - S) * 360)); c0, c1 = int(round((W + 180) * 360)), int(round((E + 180) * 360))
    lat = g["lat"][r0:r1]; lon = g["lon"][c0:c1]; a = g["lccs_class"][0, r0:r1, c0:c1]
    assert abs(lat[0] - (N - 1 / 720)) < 1e-6 and abs(lon[0] - (W + 1 / 720)) < 1e-6, (lat[0], lon[0])
    o = netCDF4.Dataset(out + ".tmp", "w")
    o.createDimension("time", 1); o.createDimension("lat", len(lat)); o.createDimension("lon", len(lon))
    v = o.createVariable("lat", "f8", ("lat",)); v[:] = lat; v.units = "degrees_north"
    v = o.createVariable("lon", "f8", ("lon",)); v[:] = lon; v.units = "degrees_east"
    v = o.createVariable("lccs_class", "u1", ("time", "lat", "lon"), zlib=True, complevel=4); v[0] = np.asarray(a)
    o.source = f"subset of {stem}.nc, rows {r0}:{r1}, cols {c0}:{c1}"; o.close(); g.close()
    os.replace(out + ".tmp", out); print(yr, "written", a.shape)
