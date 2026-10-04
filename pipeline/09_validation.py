# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

import sys, os, csv, numpy as np, rasterio
from rasterio.warp import reproject, Resampling, transform_bounds
from rasterio.windows import from_bounds, Window
from common import HARM, TAB, Budget, log, GHSL_DIR, GAUD_DIR
GHSL = GHSL_DIR; GAUD = GAUD_DIR
b = Budget(float(sys.argv[1]) if len(sys.argv) > 1 else float(os.environ.get("BGD_BUDGET", 1e12)))
with rasterio.open(os.path.join(HARM, "LandArea300_GLC30.tif")) as r:
    DST_T, DST_CRS, H, W = r.transform, r.crs, r.height, r.width
    LAND = r.read(1); LAND = np.where(LAND < 0, 0, LAND)
BNDS = (88.0, 21.0, 93.0, 27.0)

def write(path, arr):
    with rasterio.open(path, "r+" if os.path.exists(path) else "w", driver="GTiff", height=H, width=W, count=1, dtype="float32",
                       crs=DST_CRS, transform=DST_T, nodata=-9999, compress="deflate", tiled=True) as d:
        d.write(arr.astype("float32"), 1)

for yr in (2000, 2015, 2020):
    out = os.path.join(HARM, f"GHSL_frac300_{yr}.tif")
    if os.path.exists(out) or not b.ok(60): continue
    src_path = os.path.join(GHSL, f"GHS_BUILT_S_E{yr}_GLOBE_R2023A_54009_100_V1_0.tif")
    with rasterio.open(src_path) as s:
        sb = transform_bounds(DST_CRS, s.crs, *BNDS, densify_pts=21)
        win = from_bounds(*sb, transform=s.transform).round_offsets().round_lengths()
        a = s.read(1, window=win).astype("float32"); a[a == s.nodata] = 0.0
        st = s.window_transform(win)
        dst = np.zeros((H, W), dtype="float32")
        reproject(a, dst, src_transform=st, src_crs=s.crs, dst_transform=DST_T, dst_crs=DST_CRS, resampling=Resampling.average, src_nodata=None, dst_nodata=-9999)
        frac = np.clip(dst / 10000.0, 0, 1)
        frac[dst < 0] = 0.0
        write(out, frac); log(f"GHSL {yr}: built-up {(frac * LAND).sum():.0f} km2 on common grid")

GAUD_YEARS = (2000, 2015)
for Y in GAUD_YEARS:
    out = os.path.join(HARM, f"GAUD_frac300_{Y}.tif"); prog = out + ".rows"
    if os.path.exists(prog) and open(prog).read().strip() == "DONE": continue
    dst = rasterio.open(out).read(1) if os.path.exists(out) else np.full((H, W), -9999, dtype="float32")
    done = int(open(prog).read()) if os.path.exists(prog) else 0
    tiles = [rasterio.open(os.path.join(GAUD, f"urban_grid_{i}.tif")) for i in (116, 117)]
    STRIP = 90
    while done < H and b.ok(25):
        r0, r1 = done, min(done + STRIP, H)
        top = DST_T.f + r0 * DST_T.e; bottom = DST_T.f + r1 * DST_T.e
        strip = np.zeros((r1 - r0, W), dtype="float32"); wsum = np.zeros_like(strip)
        for t in tiles:
            tb = t.bounds
            x0, x1 = max(BNDS[0], tb.left), min(BNDS[2], tb.right)
            if x1 <= x0: continue
            win = from_bounds(x0, bottom, x1, top, transform=t.transform).round_offsets().round_lengths()
            if win.width <= 0 or win.height <= 0: continue
            a = t.read(1, window=win)
            binar = ((a > 0) & (a <= Y)).astype("float32")
            part = np.full(strip.shape, -9999, dtype="float32")
            reproject(binar, part, src_transform=t.window_transform(win), src_crs=t.crs,
                      dst_transform=rasterio.Affine(DST_T.a, 0, DST_T.c, 0, DST_T.e, top), dst_crs=DST_CRS,
                      resampling=Resampling.average, dst_nodata=-9999)
            ok = part >= 0
            strip[ok] += part[ok]; wsum[ok] += 1
        strip = np.where(wsum > 0, strip / np.maximum(wsum, 1), 0.0)
        dst[r0:r1] = strip; done = r1
        write(out, dst); open(prog, "w").write(str(done))
    if done >= H:
        open(prog, "w").write("DONE"); log(f"GAUD {Y}: urban {(np.clip(dst, 0, 1) * LAND).sum():.0f} km2 on common grid")
    for t in tiles: t.close()

outcsv = os.path.join(TAB, "validation_metrics.csv")
need = all(os.path.exists(os.path.join(HARM, f)) for f in [f"GHSL_frac300_{y}.tif" for y in (2000, 2015, 2020)]) and all(os.path.exists(os.path.join(HARM, f"GAUD_frac300_{y}.tif.rows")) and open(os.path.join(HARM, f"GAUD_frac300_{y}.tif.rows")).read().strip() == "DONE" for y in GAUD_YEARS)
if need and not os.path.exists(outcsv) and b.ok(60):
    def rd(f): return rasterio.open(os.path.join(HARM, f)).read(1)
    with rasterio.open(os.path.join(HARM, "Hotspots_1km_GLC30.tif")) as hs: hot9 = hs.read(6)
    hot = np.kron(hot9, np.ones((3, 3)))[:H, :W] == 3
    mask = LAND > 0
    def product(name, yr):
        if name == "GLC30_fraction":
            a = rd(f"UrbanFrac300_GLC30_{yr}.tif"); a[a < 0] = 0; return a.astype("float32")
        if name == "GLC30_majority": return (rasterio.open(os.path.join(HARM, f"GLC300maj_{yr}.tif")).read(1) == 1).astype("float32")
        if name == "ESA300": return (rd(f"ESA300_{yr}.tif") == 1).astype("float32")
    rows = []
    pairs = [("GHSL", 2000, 2000), ("GHSL", 2020, 2020), ("GHSL", 2015, 2016), ("GAUD", 2000, 2000), ("GAUD", 2015, 2016)]
    for ref, ry, py in pairs:
        R = np.clip(rd(f"{ref}_frac300_{ry}.tif"), 0, 1)
        for prod in ("GLC30_fraction", "GLC30_majority", "ESA300"):
            P = product(prod, py)
            for stratum, m in (("national", mask), ("hotspot", mask & hot), ("outside_hotspot", mask & ~hot)):
                p, r, l = P[m], R[m], LAND[m]
                ref_km2 = float((r * l).sum()); prod_km2 = float((p * l).sum())
                D = float((np.abs(p - r) * l).sum()); Q = abs(prod_km2 - ref_km2); A = D - Q
                pr = float(np.corrcoef(p, r)[0, 1])
                row = dict(reference=ref, ref_year=ry, product=prod, prod_year=py, stratum=stratum, land_km2=float(l.sum()),
                           ref_km2=ref_km2, prod_km2=prod_km2, ratio_prod_ref=prod_km2 / ref_km2 if ref_km2 else np.nan,
                           mae_frac=float(np.abs(p - r).mean()), bias_frac=float((p - r).mean()), pearson_r=pr,
                           total_disagreement_km2=D, quantity_km2=Q, allocation_km2=A, agreement_km2=float((np.minimum(p, r) * l).sum()))
                row["capture_any"] = float((r * l)[p > 0].sum()) / ref_km2 if ref_km2 else np.nan
                row["capture_maj"] = float((r * l)[p >= 0.5].sum()) / ref_km2 if ref_km2 else np.nan
                row["cells_any"] = int((p > 0).sum()); row["cells_maj"] = int((p >= 0.5).sum())
                for t in (0.05, 0.25, 0.5):
                    rb = r >= t
                    pb = (p >= t) if prod == "GLC30_fraction" else (p > 0.5)
                    tp = float((rb & pb).sum()); fp = float((~rb & pb).sum()); fn = float((rb & ~pb).sum())
                    row[f"ref_cells_ge{t}"] = int(rb.sum()); row[f"prod_cells_t{t}"] = int(pb.sum())
                    row[f"F1_t{t}"] = 2 * tp / (2 * tp + fp + fn) if tp + fp + fn else np.nan
                    row[f"producers_t{t}"] = tp / (tp + fn) if tp + fn else np.nan
                    row[f"users_t{t}"] = tp / (tp + fp) if tp + fp else np.nan
                rows.append(row)
        log(f"validation {ref} {ry} done")
    with open(outcsv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print("metrics written")
print("stage end; left %.0fs" % b.left())
