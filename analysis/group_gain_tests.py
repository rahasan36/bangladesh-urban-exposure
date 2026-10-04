# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

import os, sys
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path[:0] = [REPO, os.path.join(REPO, "figures"), os.path.join(REPO, "sdm")]
from config import HARM, MR_TAB, TAB, OUT_FIG, SDM_RASTERS, SDM_WORK, SDM_DATA, REV, SUPP
import pandas as pd, numpy as np
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
raw = pd.read_excel(os.path.join(SUPP, "Data_S2_species_habitat_models.xlsx"), sheet_name="Data S2", header=None)
g = raw.iloc[0].ffill(); d = raw.iloc[2:].reset_index(drop=True)
d.columns = [f"{a}|{str(b).replace('.0','')}" if j >= 8 else str(raw.iloc[1, j]) for j, (a, b) in enumerate(zip(g, raw.iloc[1]))]
P = {"GLC_FCS30D 30 m": "30 m", "GHSL-derived 100 m": "100 m", "ESA CCI 300 m": "300 m", "GLC_FCS30D aggregated to 300 m": "30 m aggregated to 300 m"}
rows = []
for col, lab in P.items():
    gain = pd.to_numeric(d[f"{col}|Change 2000 to 2020"])
    for grp in ["All species"] + sorted(d["Taxonomic group"].unique()):
        m = np.ones(len(d), bool) if grp == "All species" else (d["Taxonomic group"] == grp).values
        a = gain[m & (d["National threat status"] == "Threatened").values]; b = gain[m & (d["National threat status"] == "Non-threatened").values]
        if len(a) < 3 or len(b) < 3:
            rows.append(dict(product=lab, group=grp, n_thr=len(a), n_non=len(b))); continue
        u, p = mannwhitneyu(a, b, alternative="two-sided")
        rows.append(dict(product=lab, group=grp, n_thr=len(a), n_non=len(b), median_thr=a.median(), median_non=b.median(),
                         diff=a.median() - b.median(), rank_biserial=1 - 2 * u / (len(a) * len(b)), p=p))
t = pd.DataFrame(rows)
for lab in t["product"].unique():
    m = (t["product"] == lab) & (t.group != "All species") & t.p.notna()
    t.loc[m, "q_BH"] = multipletests(t.loc[m, "p"], method="fdr_bh")[1]
t.to_csv(os.path.join(REV, "group_gain_tests.csv"), index=False)
pd.set_option("display.width", 200)
print(t.round(4).to_string(index=False))
