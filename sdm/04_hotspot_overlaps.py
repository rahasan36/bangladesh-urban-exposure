# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

import os, sys, numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path[:0] = [REPO, os.path.join(REPO, "sdm")]
from config import SDM_WORK, SDM_DATA
from gistar import gistar
z = np.load(os.path.join(SDM_WORK, "urban_sdm_grid.npz")); L, CA = z["L"], z["CA"]
irr = np.load(os.path.join(SDM_WORK, "irr.npy")); nval = np.load(os.path.join(SDM_WORK, "nval.npy")); rt = np.load(os.path.join(SDM_WORK, "rich_t.npy")).astype(float)
A = (nval == 1098) & (irr > -1e30) & (L / CA > 0.5)
P = ["GLC30", "GHSL100", "ESA300"]
def hot(zz, rule):
    if rule == "z>2.58": return (zz > 2.58) & A
    p = 1 - stats.norm.cdf(np.nan_to_num(zz[A], nan=-10)); out = np.zeros_like(A); out[A] = multipletests(p, alpha=0.01, method="fdr_bh")[0]; return out
rows = []
for k in (3, 5):
    for rule in ("z>2.58", "FDR q=0.01"):
        RH = hot(gistar(np.where(A, rt, 0.0), A, k), rule); IH = A & (irr >= 1)
        for p in P:
            UH = hot(gistar(z[f"{p}_2020"] - z[f"{p}_2000"], A, k), rule)
            for surf, S in [("threatened richness", RH), ("irreplaceability>=1", IH)]:
                obs, exp = 100 * L[UH & S].sum() / L[UH].sum(), 100 * L[S].sum() / L[A].sum()
                rows.append(dict(neighborhood=f"{k}x{k}", rule=rule, product=p, surface=surf, urban_hotspot_km2=round(L[UH].sum()), obs_pct=round(obs, 1), exp_pct=round(exp, 1), ratio=round(obs / exp, 2)))
pd.DataFrame(rows).to_csv(os.path.join(SDM_DATA, "hotspot_sensitivity_sdm.csv"), index=False)
cls = [("0", (irr <= 0)), ("0 to 1", (irr > 0) & (irr < 1)), ("1 to 2", (irr >= 1) & (irr < 2)), ("2 to 5", (irr >= 2) & (irr < 5)), ("5 to 15", (irr >= 5))]
out = []
for nm, m in cls:
    m = m & A; r = dict(irreplaceability_class=nm, land_share_pct=100 * L[m].sum() / L[A].sum())
    for p in P:
        g = (z[f"{p}_2020"] - z[f"{p}_2000"]); r[f"{p}_gain_share_pct"] = 100 * g[m].sum() / g[A].sum(); r[f"{p}_ratio"] = r[f"{p}_gain_share_pct"] / r["land_share_pct"]
    out.append(r)
pd.DataFrame(out).to_csv(os.path.join(SDM_DATA, "irreplaceability_class_gain.csv"), index=False)
print(pd.DataFrame(rows).to_string(index=False)); print(pd.DataFrame(out).round(2).to_string(index=False))
