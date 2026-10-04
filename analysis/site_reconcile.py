# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

import os, sys
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path[:0] = [REPO, os.path.join(REPO, "figures"), os.path.join(REPO, "sdm")]
from config import HARM, MR_TAB, TAB, OUT_FIG, SDM_RASTERS, SDM_WORK, SDM_DATA, REV, SUPP
import pandas as pd, numpy as np
from scipy.stats import spearmanr
def load(f):
    z = pd.read_csv(f); z = z[z.feature_id.astype(str) != "ALL"]
    t = z.groupby(["feature_type", "name", "zone", "year"]).area_km2.sum()
    u = z[z.class_id == 1].groupby(["feature_type", "name", "zone", "year"]).area_km2.sum().reindex(t.index).fillna(0)
    return pd.DataFrame({"land": t, "urban": u}).reset_index()
src = {"GLC30": os.path.join(TAB, "pa_kba_zonal_GLC30.csv"), "ESA300": os.path.join(TAB, "pa_kba_zonal_ESA300.csv"), "GHSL100": os.path.join(MR_TAB, "pa_kba_zonal_GHSL100.csv")}
D = {k: load(v) for k, v in src.items()}
rows = None
for k, d in D.items():
    c = d[(d.zone == "core") & (d.year == 2020)].set_index(["feature_type", "name"])
    x = pd.DataFrame({f"core_land_{k}": c.land, f"core_urban_{k}": c.urban, f"core_pct_{k}": 100 * c.urban / c.land.where(c.land > 0)})
    rows = x if rows is None else rows.join(x, how="outer")
rows = rows.reset_index()
rows["incl_core_land_30m_ge_half_km2"] = rows.core_land_GLC30 >= 0.5
rows["incl_all_three_grids"] = rows["incl_core_land_30m_ge_half_km2"] & rows.core_pct_GHSL100.notna() & rows.core_pct_ESA300.notna()
print("sites:", len(rows), "| core land >=0.5 km2 at 30 m:", rows["incl_core_land_30m_ge_half_km2"].sum(), "| and land on all three grids:", rows.incl_all_three_grids.sum())
a = rows[rows.incl_all_three_grids]
print("rho 30-100 %.2f, 30-300 %.2f, 100-300 %.2f" % (spearmanr(a.core_pct_GLC30, a.core_pct_GHSL100)[0], spearmanr(a.core_pct_GLC30, a.core_pct_ESA300)[0], spearmanr(a.core_pct_GHSL100, a.core_pct_ESA300)[0]))
print("zero cores among these: 30m %d, 100m %d, 300m %d" % ((a.core_pct_GLC30 == 0).sum(), (a.core_pct_GHSL100 == 0).sum(), (a.core_pct_ESA300 == 0).sum()))
b = rows[rows["incl_core_land_30m_ge_half_km2"]]
print("zero cores among 30 m >=0.5 km2 set: 30m", (b.core_pct_GLC30 == 0).sum())
ex = rows[~rows.incl_all_three_grids]
print("excluded:\n", ex[["feature_type", "name", "core_land_GLC30", "core_land_GHSL100", "core_land_ESA300"]].round(3).to_string(index=False))
rows.to_csv(os.path.join(REV, "sites.csv"), index=False)
