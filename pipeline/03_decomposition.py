# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

import os, sys, argparse
import numpy as np, rasterio, pandas as pd
from common import *

ap = argparse.ArgumentParser(); ap.add_argument("--budget", type=float, default=float(os.environ.get("BGD_BUDGET", 1e12))); args = ap.parse_args()
B = Budget(args.budget)
K = NCLASS + 1
YEARS = COMMON_YEARS
BLOCK = 2048

with rasterio.open(os.path.join(HARM, "mask_esa300.tif")) as r:
    esa_mask = r.read(1).astype(bool); et = r.transform; EH, EW = r.height, r.width
    EPROF = dict(driver="GTiff", crs=r.crs, transform=et, height=EH, width=EW, compress="LZW", tiled=True, blockxsize=256, blockysize=256)
with rasterio.open(os.path.join(HARM, "mask_glc30.tif")) as r:
    gt = r.transform; GH, GW = r.height, r.width
lat_c = gt.f + (np.arange(GH) + 0.5) * gt.e
lon_c = gt.c + (np.arange(GW) + 0.5) * gt.a
rowmap = np.floor((et.f - lat_c) / -et.e).astype(np.int64)
colmap = np.floor((lon_c - et.c) / et.a).astype(np.int64)
valid_rows = (rowmap >= 0) & (rowmap < EH); valid_cols = (colmap >= 0) & (colmap < EW)
DUMMY = EH * EW
rowmap = np.where(valid_rows, rowmap, 0); colmap = np.where(valid_cols, colmap, 0)
log(f"GLC rows outside ESA subset: {(~valid_rows).sum()} of {GH}; cols: {(~valid_cols).sum()} of {GW}")

gra = row_areas_km2(gt, GH, GW)
emask_flat = np.append(esa_mask.ravel(), False)
def aggregate(yr):
    counts = np.zeros((K, EH * EW), dtype=np.int32)
    a30c = np.zeros(K)
    with rasterio.open(os.path.join(HARM, f"GLC30_{yr}.tif")) as r:
        for row0 in range(0, GH, BLOCK):
            n = min(BLOCK, GH - row0)
            a = r.read(1, window=rasterio.windows.Window(0, row0, GW, n))
            key = (rowmap[row0:row0 + n][:, None] * EW + colmap[None, :])
            key[~valid_rows[row0:row0 + n], :] = DUMMY; key[:, ~valid_cols] = DUMMY
            key = key.ravel(); af = a.ravel()
            inm = emask_flat[key]
            a30c += np.bincount(af[inm], weights=np.repeat(gra[row0:row0 + n], GW)[inm], minlength=K)[:K]
            for c in range(1, K):
                sel = af == c
                if sel.any(): counts[c] += np.bincount(key[sel], minlength=EH * EW + 1)[:EH * EW].astype(np.int32)
    tot = counts[1:].sum(0)
    frac = np.zeros((K - 1, EH * EW), dtype=np.float32)
    nz = tot > 0
    frac[:, nz] = counts[1:, nz] / tot[nz]
    maj = (np.argmax(counts[1:], axis=0) + 1).astype(np.uint8); maj[~nz] = 0
    maj = maj.reshape(EH, EW) * esa_mask
    frac = frac.reshape(K - 1, EH, EW) * esa_mask
    pd.DataFrame(dict(year=yr, class_id=range(1, K), glc30_common_km2=a30c[1:])).to_csv(os.path.join(TAB, f"glc30_commonmask_{yr}.csv"), index=False)
    return maj, frac, tot.reshape(EH, EW)

for yr in YEARS:
    fm = os.path.join(HARM, f"GLC300maj_{yr}.tif"); ff = os.path.join(HARM, f"GLC300frac_{yr}.tif")
    if os.path.exists(os.path.join(TAB, f"glc30_commonmask_{yr}.csv")): continue
    if not B.ok(70): print("CONTINUE"); sys.exit(0)
    maj, frac, tot = aggregate(yr)
    with (rasterio.open(fm, "r+") if os.path.exists(fm) else rasterio.open(fm, "w", dtype="uint8", count=1, nodata=0, **EPROF)) as d: d.write(maj, 1)
    with (rasterio.open(ff, "r+") if os.path.exists(ff) else rasterio.open(ff, "w", dtype="uint8", count=K - 1, nodata=255, **EPROF)) as d:
        pct = np.round(frac * 100).astype(np.uint8); pct[:, ~esa_mask] = 255
        d.write(pct)
        for c in range(1, K): d.set_band_description(c, CLASS_NAMES[c])
    log(f"GLC30 -> 300 m aggregation {yr}: GLC pixels per cell median {np.median(tot[esa_mask]):.0f}")

ra = row_areas_km2(et, EH, EW)
W = np.repeat(ra, EW).reshape(EH, EW)
g30 = pd.read_csv(os.path.join(TAB, "class_areas_GLC30.csv"))
e300 = pd.read_csv(os.path.join(TAB, "class_areas_ESA300.csv"))
recs = []; agree = []
for yr in YEARS:
    with rasterio.open(os.path.join(HARM, f"GLC300maj_{yr}.tif")) as r: maj = r.read(1)
    with rasterio.open(os.path.join(HARM, f"GLC300frac_{yr}.tif")) as r: frac = r.read().astype(np.float32)
    with rasterio.open(os.path.join(HARM, f"ESA300_{yr}.tif")) as r: esa = r.read(1)
    frac[frac == 255] = 0; frac /= 100.0
    a30c = pd.read_csv(os.path.join(TAB, f"glc30_commonmask_{yr}.csv")).set_index("class_id").glc30_common_km2
    for c in range(1, K):
        a30 = float(g30[(g30.year == yr) & (g30.class_id == c)].area_km2.iloc[0])
        a300f = float((frac[c - 1] * W).sum())
        a300m = float(((maj == c) * W).sum())
        ae = float(e300[(e300.year == yr) & (e300.class_id == c)].area_km2.iloc[0])
        recs.append(dict(year=yr, class_id=c, class_name=CLASS_NAMES[c], glc30_km2=float(a30c[c]), glc30_fullcountry_km2=a30, glc300_fraction_km2=a300f,
                         glc300_majority_km2=a300m, esa300_km2=ae,
                         resolution_effect_km2=a300m - float(a30c[c]), product_effect_km2=ae - a300m, total_difference_km2=ae - float(a30c[c])))
    m = np.bincount((maj.astype(np.int64) * K + esa.astype(np.int64)).ravel(), weights=W.ravel(), minlength=K * K).reshape(K, K)[1:, 1:]
    oa = np.trace(m) / m.sum()
    for c in range(1, K):
        i = c - 1
        agree.append(dict(year=yr, class_id=c, class_name=CLASS_NAMES[c], overall_agreement=oa,
                          glc300_area=m[i].sum(), esa300_area=m[:, i].sum(), agreeing_area=m[i, i],
                          producers_acc_glc_as_ref=m[i, i] / max(m[i].sum(), 1e-9), users_acc_glc_as_ref=m[i, i] / max(m[:, i].sum(), 1e-9)))
    fu = frac[0]; mu = (maj == 1).astype(np.float32); eu = (esa == 1).astype(np.float32)
    fd = os.path.join(HARM, f"UrbanDiff300_{yr}.tif")
    with (rasterio.open(fd, "r+") if os.path.exists(fd) else rasterio.open(fd, "w", dtype="float32", count=4, nodata=-9999, **EPROF)) as d:
        out = np.stack([fu, mu - fu, eu - mu, eu - fu]); out[:, ~esa_mask] = -9999
        d.write(out)
        for i, n in enumerate(["glc30_urban_fraction", "resolution_effect", "product_effect", "total_difference"]): d.set_band_description(i + 1, n)
pd.DataFrame(recs).to_csv(os.path.join(TAB, "decomposition_by_class.csv"), index=False)
pd.DataFrame(agree).to_csv(os.path.join(TAB, "agreement_GLC300maj_vs_ESA300.csv"), index=False)
d = pd.DataFrame(recs); u = d[d.class_id == 1]
print(u[["year", "glc30_km2", "glc300_majority_km2", "esa300_km2", "resolution_effect_km2", "product_effect_km2", "total_difference_km2"]].round(1).to_string(index=False))
print("DONE")
