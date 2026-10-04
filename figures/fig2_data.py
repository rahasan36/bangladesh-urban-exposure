# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

import os, sys, pandas as pd
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path[:0] = [REPO]
from config import TAB, MR_TAB, TAB_MS
rows = []
for p, f in [("GLC30", os.path.join(TAB, "pa_kba_zonal_GLC30.csv")), ("GHSL100", os.path.join(MR_TAB, "pa_kba_zonal_GHSL100.csv")), ("ESA300", os.path.join(TAB, "pa_kba_zonal_ESA300.csv"))]:
    z = pd.read_csv(f); z = z[z.feature_id.astype(str) == "ALL"]
    t = z.groupby(["feature_type", "zone", "year"]).area_km2.sum(); u = z[z.class_id == 1].groupby(["feature_type", "zone", "year"]).area_km2.sum().reindex(t.index).fillna(0)
    d = pd.DataFrame({"land": t, "urban": u}).reset_index()
    r10 = d[d.zone.isin(["r0_5", "r5_10"])].groupby(["feature_type", "year"])[["land", "urban"]].sum().reset_index(); r10["zone"] = "r0_10"
    d = pd.concat([d, r10]); d["urban_pct"] = 100 * d.urban / d.land; d["product"] = p
    rows.append(d[["product", "feature_type", "zone", "year", "urban_pct"]])
out = pd.concat(rows).sort_values(["product", "feature_type", "zone", "year"]); os.makedirs(TAB_MS, exist_ok=True)
out.to_csv(os.path.join(TAB_MS, "pa_kba_core_buffer.csv"), index=False); print(out[out.year == 2020].round(3).to_string(index=False))
