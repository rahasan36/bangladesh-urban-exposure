# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

import os, sys, numpy as np, pandas as pd, warnings; warnings.filterwarnings("ignore")
from scipy import stats
import statsmodels.formula.api as smf
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, REPO)
from config import SDM_DATA
src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(SDM_DATA, "sdm_species.csv"); sfx = sys.argv[2] if len(sys.argv) > 2 else ""
d = pd.read_csv(src); x = d[d.threatStatus != "Data Deficient"]
t, n = x[x.thr == 1], x[x.thr == 0]
tests, models, counts = [], [], []
for p in ["GLC30", "GHSL100", "ESA300"]:
    for meas, col in [("2020", f"{p}_2020"), ("gain", f"{p}_gain")]:
        tests.append(dict(rep=p, measure=meas, med_thr=round(t[col].median(), 3), med_non=round(n[col].median(), 3), MWU_P=stats.mannwhitneyu(t[col], n[col]).pvalue, n_thr=float(len(t)), n_non=float(len(n))))
    cats = ["CR", "EN", "VU", "NT", "LC", "DD"]
    tests.append(dict(rep=p, measure="KW_category_2020", MWU_P=stats.kruskal(*[d[d.nationalStatus == c][f"{p}_2020"] for c in cats]).pvalue))
    for g, gg in x.groupby("group"):
        a, b = gg[gg.thr == 1][f"{p}_2020"], gg[gg.thr == 0][f"{p}_2020"]
        if len(a) < 5: continue
        tests.append(dict(rep=p, measure=f"group:{g}", med_thr=round(a.median(), 3), med_non=round(b.median(), 3), MWU_P=stats.mannwhitneyu(a, b).pvalue, n_thr=float(len(a)), n_non=float(len(b))))
    ys = sorted(int(c[-4:]) for c in d.columns if c.startswith(f"traj_{p}_"))
    lg = x.melt(id_vars=["species", "thr", "group"], value_vars=[f"traj_{p}_{y}" for y in ys], var_name="v", value_name="share")
    lg["yc"] = lg.v.str[-4:].astype(int) - 2000
    for nm, re in [("mixed intercept", None), ("mixed intercept+slope", "~yc")]:
        m = smf.mixedlm("share ~ thr*yc + C(group)", lg, groups=lg.species, re_formula=re).fit(reml=True); ci = m.conf_int()
        for term in ["thr:yc", "thr"]: models.append(dict(rep=p, model=nm, term=term, est=m.params[term], lo=ci.loc[term, 0], hi=ci.loc[term, 1], P=m.pvalues[term]))
    m = smf.ols(f"{p}_gain ~ thr + log_area + C(group)", x).fit(cov_type="HC3"); ci = m.conf_int()
    for term in ["thr", "log_area"]: models.append(dict(rep=p, model="OLS gain, HC3", term=term, est=m.params[term], lo=ci.loc[term, 0], hi=ci.loc[term, 1], P=m.pvalues[term]))
    r = stats.spearmanr(x[f"{p}_gain"], x.log_area); models.append(dict(rep=p, model="Spearman gain vs log area", term="rho", est=r.statistic, P=r.pvalue))
CELL = {"GLC30": 30, "GLCp90": 90, "GLCp300": 300, "GLCp510": 510, "GLCp990": 990, "GHSL100": 100, "GHSLm300": 300, "GHSLm500": 500, "GHSLm1000": 1000, "ESA300": 300}
for p, c in CELL.items():
    s = d[f"{p}_2020"]; st = d[d.thr == 1][f"{p}_2020"]
    counts.append(dict(rep=p, cell=c, n_all_gt1=int((s > 1).sum()), n_thr_gt1=int((st > 1).sum()), n_all_gt2=int((s > 2).sum()), n_thr_gt2=int((st > 2).sum()), median_all=round(s.median(), 3), median_thr=round(st.median(), 3)))
out = os.path.dirname(src)
pd.DataFrame(tests).to_csv(os.path.join(out, f"sdm_tests{'_' + sfx if sfx else ''}.csv"), index=False)
pd.DataFrame(models).to_csv(os.path.join(out, f"sdm_models{'_' + sfx if sfx else ''}.csv"), index=False)
pd.DataFrame(counts).to_csv(os.path.join(out, f"sdm_counts_ladder{'_' + sfx if sfx else ''}.csv"), index=False)
print(pd.DataFrame(counts).to_string(index=False))
