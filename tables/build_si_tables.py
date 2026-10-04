# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

import os, sys
import numpy as np, pandas as pd
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path[:0] = [REPO, os.path.join(REPO, "pipeline")]
from config import TAB, MR, MR_TAB, MR_REVIEW, TAB_MS, SDM_DATA, REV, OUT_TABLES

PN = {"GLC30": "GLC_FCS30D 30 m", "GHSL100": "GHSL 100 m", "ESA300": "ESA CCI 300 m"}
P3 = ["GLC30", "GHSL100", "ESA300"]
esa_nat = pd.read_csv(f"{TAB}/class_areas_ESA300.csv"); esa_nat = esa_nat[esa_nat.class_id == 1].set_index("year").area_km2

def decomposition():
    d = pd.read_csv(f"{TAB}/decomposition_by_class.csv"); d = d[d.class_id == 1].drop_duplicates("year").sort_values("year")
    rows = [["GLC_FCS30D 30 m", int(r.year), r.glc30_km2, r.glc300_majority_km2, r.esa300_km2, r.resolution_effect_km2, r.product_effect_km2,
             r.total_difference_km2] for _, r in d.iterrows()]
    for _, r in pd.read_csv(f"{MR_TAB}/decomposition_three_products.csv").iterrows():
        res = r.GHSL100_maj300_km2 - r.GHSL100_frac_km2; tot = r.ESA300_km2 - r.GHSL100_frac_km2
        rows.append(["GHSL 100 m", int(r.year), r.GHSL100_frac_km2, r.GHSL100_maj300_km2, r.ESA300_km2, res, tot - res, tot])
    out = pd.DataFrame(rows, columns=["product", "year", "fine resolution km²", "aggregated to 300 m km²", "ESA CCI km²", "resolution component km²",
                                      "product component km²", "total difference km²"])
    out["resolution share (%)"] = 100 * out["resolution component km²"] / out["total difference km²"]
    return out.round(1)

def resolution_ladder():
    g = pd.read_csv(f"{MR}/T_ladder_glc_national.csv"); rows = []
    for yr, grp in g.groupby("year"):
        nat = grp[grp.factor == 1].urban_maj_km2.iloc[0]
        for _, r in grp.sort_values("factor").iterrows():
            rows.append(["GLC_FCS30D, plurality rule", int(yr), r.cell, r.urban_maj_km2, 100 * r.urban_maj_km2 / nat])
    t = pd.read_csv(f"{MR}/T_ladder_ghsl_thresholds.csv")
    for thr in [0.05, 0.10, 0.20, 0.30, 0.50]:
        for yr, grp in t[np.isclose(t.threshold, thr)].groupby("year"):
            nat = grp[grp.factor == 1].urban_km2.iloc[0]
            for _, r in grp.sort_values("factor").iterrows():
                rows.append([f"GHSL-derived, built ≥ {int(round(thr * 100))}%", int(yr), r.cell, r.urban_km2, 100 * r.urban_km2 / nat if nat > 0 else np.nan])
    return pd.DataFrame(rows, columns=["representation", "year", "cell", "urban_km2", "share_of_native_pct"]).round(1)

def aggregation_rules():
    r = pd.read_csv(f"{MR_REVIEW}/threshold_sensitivity_300m.csv")
    return r.rename(columns={"fraction_km2": "fraction km²", "plurality_km2": "plurality km²", "esa_km2": "ESA CCI km²", "thr10_km2": "≥10% km²",
                             "thr10_resshare_pct": "≥10% resolution share (%)", "thr25_km2": "≥25% km²", "thr25_resshare_pct": "≥25% resolution share (%)",
                             "thr50_km2": "≥50% km²", "thr50_resshare_pct": "≥50% resolution share (%)",
                             "plurality_resshare_pct": "plurality resolution share (%)"}).round(1)

def table_s1(t1, s3):
    d3 = pd.read_csv(f"{MR_TAB}/decomposition_three_products.csv").set_index("year"); s3 = s3.set_index("year"); rows = []
    for _, r in t1.iterrows():
        y = int(r.year); nat = r["fine resolution km²"]; agg = r["aggregated to 300 m km²"]; res = r["resolution component km²"]
        if r["product"].startswith("GLC"):
            frac, esa, prod, tot = s3.loc[y, "fraction km²"], r["ESA CCI km²"], r["product component km²"], r["total difference km²"]
        else:
            frac, esa = d3.loc[y, "GHSL100_frac_km2"], esa_nat[y]; tot = esa - nat; prod = esa - agg
        rows.append([r["product"], y, nat, frac, agg, esa, res, prod, tot, 100 * res / r["total difference km²"]])
    return pd.DataFrame(rows, columns=["Map", "Year", "Native", "Fraction", "Aggregated", "ESA CCI", "Resolution", "Product", "Total",
                                       "Resolution share (%)"])

def table_s2():
    g = pd.read_csv(f"{MR}/T_ladder_glc_national.csv"); h = pd.read_csv(f"{MR}/T_ladder_ghsl_thresholds.csv"); h = h[np.isclose(h.threshold, 0.2)]
    z = pd.read_csv(f"{MR}/T_zones_ladder.csv"); z = z[(z.feature_type == "PA") & (z.year == 2020)].set_index(["product", "zone"]).urban_pct
    e = pd.read_csv(f"{TAB_MS}/pa_kba_core_buffer.csv"); e = e[(e["product"] == "ESA300") & (e.feature_type == "PA") & (e.year == 2020)].set_index("zone").urban_pct
    base = g[(g.factor == 1) & (g.year == 2020)].urban_maj_km2.iloc[0]; rows = []
    for cell in ["30 m", "90 m", "300 m", "510 m", "990 m"]:
        a = g[g.cell == cell].set_index("year").urban_maj_km2; k = "GLC 30 m" if cell == "30 m" else f"GLC {cell} majority"
        rows.append(["GLC_FCS30D", cell, a[2000], a[2020], 100 * a[2020] / base, z[(k, "core")], z[(k, "r5_10")]])
    hb = h[(h.factor == 1) & (h.year == 2020)].urban_km2.iloc[0]
    for cell, k in [("100 m", "GHSL 100 m built>=0.2"), ("300 m", "GHSL 300 m majority"), ("500 m", None), ("1000 m", None)]:
        a = h[h.cell == cell].set_index("year").urban_km2
        rows.append(["GHSL-derived", cell.replace("1000 m", "1,000 m"), a[2000], a[2020], 100 * a[2020] / hb,
                     z[(k, "core")] if k else np.nan, z[(k, "r5_10")] if k else np.nan])
    rows.append(["ESA CCI", "300 m", esa_nat[2000], esa_nat[2020], 100 * esa_nat[2020] / base, e["core"], e["r5_10"]])
    out = pd.DataFrame(rows, columns=["Product", "Cell", "Urban 2000 (km²)", "Urban 2020 (km²)", "Retained 2020 (%)", "Core (%)", "Buffer 5 to 10 km (%)"])
    out["Buffer to core ratio"] = out["Buffer 5 to 10 km (%)"] / out["Core (%)"]
    return out

def table_s3(s3):
    s = s3.set_index(s3.year.astype(int)); rows = []
    for y, r in s.iterrows():
        rows.append([int(y), r["fraction km²"], r["plurality km²"], esa_nat[y], r["≥10% km²"], r["≥25% km²"], r["≥50% km²"],
                     r["plurality resolution share (%)"], r["≥25% resolution share (%)"], r["≥50% resolution share (%)"]])
    return pd.DataFrame(rows, columns=["Year", "Fraction", "Plurality", "ESA CCI", "≥10%", "≥25%", "≥50%", "Share, plurality", "Share, ≥25%", "Share, ≥50%"])

def table_s4():
    ts = pd.read_csv(f"{SDM_DATA}/sdm_tests.csv"); gg = pd.read_csv(f"{REV}/group_gain_tests.csv"); sp = pd.read_csv(f"{SDM_DATA}/sdm_species.csv")
    cols = ["Comparison", "n threatened", "n non-threatened"] + [f"{lab} {k}" for lab in ["30 m", "100 m", "300 m"] for k in ["thr.", "non", "P or q"]]
    rows = []
    for lab, meas in [("All species, share 2020", "2020"), ("All species, gain 2000 to 2020", "gain")]:
        row = [lab, int((sp.thr == 1).sum()), int((sp.threatStatus == "Non-threatened").sum())]
        for p in P3:
            x = ts[(ts.rep == p) & (ts.measure == meas)].iloc[0]; c = f"{p}_{meas}"
            row += [sp[sp.thr == 1][c].median(), sp[sp.threatStatus == "Non-threatened"][c].median(), x.MWU_P]
        rows.append(row)
    for grp in ["Birds", "Butterflies", "Freshwater fishes", "Mammals", "Reptiles", "Crustaceans"]:
        a = gg[gg.group == grp]; r30 = a[a["product"] == "30 m"].iloc[0]
        row = [f"{grp}, gain (r = {r30.rank_biserial:.2f})", int(r30.n_thr), int(r30.n_non)]
        for pr in ["30 m", "100 m", "300 m"]:
            x = a[a["product"] == pr].iloc[0]; row += [x.median_thr, x.median_non, x.q_BH]
        rows.append(row)
    return pd.DataFrame(rows, columns=cols)

def table_s5():
    sp = pd.read_csv(f"{SDM_DATA}/sdm_species.csv"); rows = []
    for c in ["CR", "EN", "VU", "NT", "LC", "DD"]:
        g = sp[sp.nationalStatus == c]; row = [c, len(g)]
        for p in P3: row += [g[f"{p}_2020"].median(), g[f"{p}_gain"].median(), 100 * (g[f"{p}_2020"] > 1).mean()]
        rows.append(row)
    a = pd.DataFrame(rows, columns=["Category", "n"] + [f"{lab} {k}" for lab in ["30 m", "100 m", "300 m"] for k in ["share", "gain", ">1%"]])
    fm = pd.read_csv(f"{REV}/family_mixed_model.csv").set_index("product"); sh = pd.read_csv(f"{REV}/share2020_models.csv").set_index("product")
    mo = pd.read_csv(f"{SDM_DATA}/sdm_models.csv"); rows = []
    for lab, src, est, lo, hi in [("Gain, least squares (HC3): threatened", fm, "ols_thr", "ols_lo", "ols_hi"),
                                  ("Gain, family random intercept: threatened", fm, "mixed_thr", "mixed_lo", "mixed_hi"),
                                  ("Share 2020, least squares (HC3): threatened", sh, "ols", "lo", "hi"),
                                  ("Share 2020, family random intercept: threatened", sh, "mixed", "mlo", "mhi")]:
        for pr, p in zip(["30 m", "100 m", "300 m"], P3): rows.append([lab, PN[p], src.loc[pr, est], src.loc[pr, lo], src.loc[pr, hi]])
    for p in P3:
        x = mo[(mo.rep == p) & (mo.model == "OLS gain, HC3") & (mo.term == "log_area")].iloc[0]
        rows.append(["Gain, least squares (HC3): tenfold habitat area", PN[p], x.est, x.lo, x.hi])
    b = pd.DataFrame(rows, columns=["Estimate (percentage points)", "Product", "Estimate", "Lower 95%", "Upper 95%"])
    return a, b

def table_s6():
    SN = {"threatened richness": "threatened richness", "irreplaceability>=1": "irreplaceability ≥ 1"}
    hs = pd.read_csv(f"{SDM_DATA}/hotspot_sensitivity_sdm.csv")
    al = pd.read_csv(f"{REV}/hotspot_allland_sensitivity.csv"); al = al[al.domain == "all land"]
    rows = [[r.neighborhood.replace("x", " by "), r["rule"], PN[r["product"]], SN[r.surface], r.urban_hotspot_km2, r.obs_pct, r.exp_pct, r.ratio] for _, r in hs.iterrows()]
    rows += [[r.neighborhood, "z>2.58, all land", PN[r["product"]], SN[r.surface], r.hotspot_km2, r.observed_pct, r.expected_pct, r.ratio] for _, r in al.iterrows()]
    a = pd.DataFrame(rows, columns=["Neighborhood", "Rule", "Product", "Surface", "Urban hotspot (km²)", "Observed (%)", "Expected (%)", "Ratio"])
    ic = pd.read_csv(f"{SDM_DATA}/irreplaceability_class_gain.csv")
    b = pd.DataFrame({"Irreplaceability class": ic.irreplaceability_class, "Share of domain land (%)": ic.land_share_pct,
                      **{f"{PN[p]} ratio": ic[f"{p}_ratio"] for p in P3}})
    return a, b

def table_s7():
    from common import ESA_XWALK, GLC_XWALK, CLASS_NAMES
    LABEL = {3: "Mosaic cropland and natural vegetation", 9: "Bare or sparse"}; rows = []
    for prod, xw in [("GLC_FCS30D", GLC_XWALK), ("ESA CCI", ESA_XWALK)]:
        inv = {}
        for code, cls in xw.items(): inv.setdefault(cls, []).append(code)
        for cls, codes in sorted(inv.items(), key=lambda kv: min(kv[1])):
            rows.append([prod, ", ".join(str(c) for c in sorted(codes)), LABEL.get(cls, CLASS_NAMES[cls])])
    return pd.DataFrame(rows, columns=["Product", "Source codes", "Common class"])

if __name__ == "__main__":
    t1, lad, s3 = decomposition(), resolution_ladder(), aggregation_rules()
    t1.to_csv(f"{TAB_MS}/decomposition_300m.csv", index=False); lad.to_csv(f"{TAB_MS}/resolution_ladder.csv", index=False)
    s3.to_csv(f"{TAB_MS}/aggregation_rules_300m.csv", index=False)
    s5a, s5b = table_s5(); s6a, s6b = table_s6()
    out = {"TableS1_decomposition": table_s1(t1, s3).round(2), "TableS2_resolution_ladder": table_s2().round(2), "TableS3_aggregation_rules": table_s3(s3).round(2),
           "TableS4_threatened_species": table_s4(), "TableS5_red_list_categories": s5a, "TableS5_model_estimates": s5b,
           "TableS6_hotspot_overlap": s6a, "TableS6_irreplaceability_classes": s6b, "TableS7_class_crosswalk": table_s7()}
    for k, v in out.items():
        v.to_csv(f"{OUT_TABLES}/{k}.csv", index=False, encoding="utf-8-sig"); print(k, len(v))
