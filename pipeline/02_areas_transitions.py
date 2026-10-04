# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

import os, sys, json, argparse
import numpy as np, rasterio, pandas as pd
from common import *

ap = argparse.ArgumentParser(); ap.add_argument("--budget", type=float, default=float(os.environ.get("BGD_BUDGET", 1e12))); args = ap.parse_args()
B = Budget(args.budget)
BLOCK = 2048
K = NCLASS + 1

def path(ds, yr): return os.path.join(HARM, f"{ds}_{yr}.tif")

def class_areas(ds, yr):
    with rasterio.open(path(ds, yr)) as r:
        ra = row_areas_km2(r.transform, r.height, r.width)
        acc = np.zeros(K)
        for row0 in range(0, r.height, BLOCK):
            n = min(BLOCK, r.height - row0)
            a = r.read(1, window=rasterio.windows.Window(0, row0, r.width, n))
            w = np.repeat(ra[row0:row0 + n], r.width)
            acc += np.bincount(a.ravel(), weights=w, minlength=K)[:K]
    return acc

def transition(ds, y1, y2):
    with rasterio.open(path(ds, y1)) as r1, rasterio.open(path(ds, y2)) as r2:
        ra = row_areas_km2(r1.transform, r1.height, r1.width)
        acc = np.zeros(K * K)
        for row0 in range(0, r1.height, BLOCK):
            n = min(BLOCK, r1.height - row0)
            win = rasterio.windows.Window(0, row0, r1.width, n)
            a = r1.read(1, window=win).astype(np.int32); b = r2.read(1, window=win).astype(np.int32)
            w = np.repeat(ra[row0:row0 + n], r1.width)
            acc += np.bincount((a * K + b).ravel(), weights=w, minlength=K * K)[:K * K]
    return acc.reshape(K, K)

PAIRS = {
    "GLC30": [(1985, 1990), (1990, 1995), (1995, 2000), (2000, 2005), (2005, 2010), (2010, 2015), (2015, 2020),
              (1985, 2020), (2000, 2020), (2000, 2004), (2004, 2008), (2008, 2012), (2012, 2016), (2016, 2020)],
    "ESA300": [(1992, 1996), (1996, 2000), (2000, 2004), (2004, 2008), (2008, 2012), (2012, 2016), (2016, 2020),
               (1992, 2020), (2000, 2020)],
}
YEARS = {"GLC30": GLC_YEARS, "ESA300": ESA_YEARS}

for ds in ["ESA300", "GLC30"]:
    f = os.path.join(TAB, f"class_areas_{ds}.csv")
    rows = pd.read_csv(f).to_dict("records") if os.path.exists(f) else []
    done_years = {r["year"] for r in rows}
    for yr in YEARS[ds]:
        if yr in done_years: continue
        if not B.ok(40): print("CONTINUE"); sys.exit(0)
        acc = class_areas(ds, yr)
        for c in range(1, K):
            rows.append(dict(dataset=ds, year=yr, class_id=c, class_name=CLASS_NAMES[c], area_km2=acc[c]))
        pd.DataFrame(rows).to_csv(f, index=False)
        log(f"{ds} {yr} areas: total {acc[1:].sum():.0f} km2, urban {acc[1]:.1f} km2")
    for (y1, y2) in PAIRS[ds]:
        f = os.path.join(TAB, f"transition_{ds}_{y1}_{y2}.csv")
        if os.path.exists(f): continue
        if not B.ok(60): print("CONTINUE"); sys.exit(0)
        m = transition(ds, y1, y2)
        recs = [dict(dataset=ds, year_from=y1, year_to=y2, from_id=i, from_name=CLASS_NAMES[i], to_id=j, to_name=CLASS_NAMES[j], area_km2=m[i, j])
                for i in range(1, K) for j in range(1, K)]
        pd.DataFrame(recs).to_csv(f, index=False)
        newurb = m[2:, 1].sum(); log(f"{ds} {y1}->{y2}: new urban {newurb:.1f} km2, from cropland {m[2,1]/max(newurb,1e-9)*100:.1f}%")
print("DONE")
