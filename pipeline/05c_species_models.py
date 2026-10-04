# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

import os, numpy as np, pandas as pd, statsmodels.formula.api as smf
from scipy import stats
from common import *
SP = os.path.join(OUT, "Species")
df = pd.read_csv(os.path.join(SP, "species_urban_exposure.csv"))
df = df[df.range_land_km2 > 0].copy()
df["status_group"] = pd.Categorical(df.status_group, ["Least Concern", "Near Threatened", "Threatened", "Data Deficient"])
df["taxon"] = pd.Categorical(df.taxon, ["Mammal", "Reptile", "Amphibian"])
df["year_c"] = df.year - 2000
df["log_urban"] = np.log(df.urban_pct + 0.01)
summ = df.groupby(["dataset", "year", "taxon", "status_group"], observed=True).urban_pct.agg(n="count", mean="mean", median="median", q25=lambda x: x.quantile(.25), q75=lambda x: x.quantile(.75)).reset_index()
summ.to_csv(os.path.join(SP, "summary_by_status_taxon_year.csv"), index=False)
thr = []
for (ds, yr, tx, sg), g in df.groupby(["dataset", "year", "taxon", "status_group"], observed=True):
    thr.append(dict(dataset=ds, year=yr, taxon=tx, status_group=sg, n=len(g), **{f"n_over_{t}pct": int((g.urban_pct > t).sum()) for t in [1, 5, 10, 20]}))
pd.DataFrame(thr).to_csv(os.path.join(SP, "threshold_crossings.csv"), index=False)
slopes = []
for (ds, sid), g in df[df.year >= 2000].groupby(["dataset", "id_no"]):
    if len(g) < 3: continue
    res = stats.linregress(g.year, g.urban_pct)
    r = g.iloc[0]
    slopes.append(dict(dataset=ds, id_no=sid, sci_name=r.sci_name, taxon=r.taxon, category=r.category, status_group=r.status_group, range_bgd_km2=r.range_bgd_km2, range_bgd_share=r.range_bgd_share,
                       urban_pct_2000=float(g[g.year == 2000].urban_pct.iloc[0]), urban_pct_2020=float(g[g.year == 2020].urban_pct.iloc[0]), slope_pct_per_yr=res.slope,
                       change_pp=float(g[g.year == 2020].urban_pct.iloc[0] - g[g.year == 2000].urban_pct.iloc[0])))
sl = pd.DataFrame(slopes); sl.to_csv(os.path.join(SP, "species_slopes_ranking.csv"), index=False)
tests = []
for ds, g in sl[sl.status_group != "Data Deficient"].groupby("dataset"):
    grp = [x.slope_pct_per_yr.values for _, x in g.groupby("status_group") if len(x) > 0]
    grp2 = [x.urban_pct_2020.values for _, x in g.groupby("status_group") if len(x) > 0]
    tests.append(dict(dataset=ds, test="Kruskal-Wallis slope by status", H=stats.kruskal(*grp).statistic, p=stats.kruskal(*grp).pvalue))
    tests.append(dict(dataset=ds, test="Kruskal-Wallis 2020 exposure by status", H=stats.kruskal(*grp2).statistic, p=stats.kruskal(*grp2).pvalue))
    a = g[g.status_group == "Threatened"]; b = g[g.status_group == "Least Concern"]
    mw = stats.mannwhitneyu(a.urban_pct_2020, b.urban_pct_2020, alternative="two-sided")
    tests.append(dict(dataset=ds, test="Mann-Whitney 2020 exposure Threatened vs LC", H=mw.statistic, p=mw.pvalue, median_threatened=a.urban_pct_2020.median(), median_LC=b.urban_pct_2020.median()))
    mw = stats.mannwhitneyu(a.slope_pct_per_yr, b.slope_pct_per_yr, alternative="two-sided")
    tests.append(dict(dataset=ds, test="Mann-Whitney slope Threatened vs LC", H=mw.statistic, p=mw.pvalue, median_threatened=a.slope_pct_per_yr.median(), median_LC=b.slope_pct_per_yr.median()))
pd.DataFrame(tests).to_csv(os.path.join(SP, "status_group_tests.csv"), index=False)
print(pd.DataFrame(tests).round(4).to_string())
out = []
for ds in ["GLC30", "ESA300"]:
    d = df[(df.dataset == ds) & (df.year >= 2000) & (df.status_group != "Data Deficient")].copy()
    d["status_group"] = d.status_group.cat.remove_unused_categories()
    for resp in ["urban_pct", "log_urban"]:
        m = smf.mixedlm(f"{resp} ~ status_group * year_c + taxon * year_c", d, groups=d["id_no"]).fit(reml=True)
        for k, v in m.params.items():
            out.append(dict(dataset=ds, response=resp, term=k, estimate=v, se=m.bse[k], p=m.pvalues[k]))
        log(f"MixedLM {ds} {resp}: n_obs {m.nobs}, n_species {d.id_no.nunique()}, converged {m.converged}")
pd.DataFrame(out).to_csv(os.path.join(SP, "mixed_effects_results.csv"), index=False)
o = pd.DataFrame(out); print(o[o.term.str.contains("year_c")].round(4).to_string())
print("DONE")
