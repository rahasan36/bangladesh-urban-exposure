# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

import sys, os, csv, numpy as np, rasterio
from common import HARM, TAB, row_areas_km2, Budget, log
OUT = os.path.join(TAB, "tile_test_90E.csv")
YEARS = [1985, 1990, 1995, 2000, 2004, 2005, 2008, 2010, 2012, 2015, 2016, 2020]
done = {}
if os.path.exists(OUT):
    for r in csv.DictReader(open(OUT)): done[int(r["year"])] = r
b = Budget(float(sys.argv[1]) if len(sys.argv) > 1 else float(os.environ.get("BGD_BUDGET", 1e12)))
for yr in YEARS:
    if yr in done: continue
    if not b.ok(40): break
    with rasterio.open(os.path.join(HARM, f"GLC30_{yr}.tif")) as src:
        t = src.transform; H, W = src.height, src.width
        areas = row_areas_km2(t, H, W)
        col90 = int((90.0 - t.c) / t.a)
        west = east = 0.0
        for r0 in range(0, H, 2048):
            win = rasterio.windows.Window(0, r0, W, min(2048, H - r0))
            a = src.read(1, window=win) == 1
            ar = areas[r0:r0 + a.shape[0]][:, None]
            west += float((a[:, :col90] * ar).sum()); east += float((a[:, col90:] * ar).sum())
    done[yr] = dict(year=yr, west_km2=round(west, 1), east_km2=round(east, 1))
    log(f"tile test {yr}: west {west:.1f} east {east:.1f}")
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["year", "west_km2", "east_km2"]); w.writeheader()
        for y in sorted(done): w.writerow(done[y])
print("done", sorted(done))
