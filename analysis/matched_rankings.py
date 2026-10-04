# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

import os, sys
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path[:0] = [REPO, os.path.join(REPO, "figures"), os.path.join(REPO, "sdm")]
from config import HARM, MR_TAB, TAB, OUT_FIG, SDM_RASTERS, SDM_WORK, SDM_DATA, REV, SUPP
import pandas as pd, numpy as np
from scipy.stats import spearmanr, kendalltau
def load(f, sh, nfix):
    raw = pd.read_excel(f, sheet_name=sh, header=None); g = raw.iloc[0].ffill(); d = raw.iloc[2:].reset_index(drop=True)
    d.columns = [f"{a}|{str(b).replace('.0','')}" if j >= nfix else str(raw.iloc[1, j]) for j, (a, b) in enumerate(zip(g, raw.iloc[1]))]
    return d
s2 = load(os.path.join(SUPP, "Data_S2_species_habitat_models.xlsx"), "Data S2", 8)
s3 = load(os.path.join(SUPP, "Data_S3_species_IUCN_ranges.xlsx"), "Data S3", 6)
key = lambda s: " ".join(str(s).strip().split()[:2]).lower()
s2["k"] = s2.Species.map(key); s3["k"] = s3.Species.map(key)
c2 = "GLC_FCS30D 30 m|2020"; c3 = [c for c in s3.columns if c.endswith("GLC_FCS30D 30 m|2020")][0]
m = s2[["k", "Species", "Taxonomic group", "National Red List category", "Modeled habitat (km²)", c2]].merge(
    s3[["k", "Class", "IUCN Red List category", "Range in Bangladesh (km²)", c3]], on="k")
m = m.rename(columns={c2: "share_habitat_2020", c3: "share_range_2020"}).drop_duplicates("k")
for c in ["share_habitat_2020", "share_range_2020", "Modeled habitat (km²)", "Range in Bangladesh (km²)"]: m[c] = pd.to_numeric(m[c])
n = len(m); print("matched species:", n, m["Taxonomic group"].value_counts().to_dict())
rho = spearmanr(m.share_habitat_2020, m.share_range_2020); tau = kendalltau(m.share_habitat_2020, m.share_range_2020)
print(f"Spearman rho {rho[0]:.2f} (p={rho[1]:.1e}), Kendall tau {tau[0]:.2f}")
print("median share habitat %.2f, range %.2f; habitat/range area median %.2f" % (m.share_habitat_2020.median(), m.share_range_2020.median(), (m['Modeled habitat (km²)'] / m['Range in Bangladesh (km²)']).median()))
m["rank_habitat"] = m.share_habitat_2020.rank(ascending=False, method="average"); m["rank_range"] = m.share_range_2020.rank(ascending=False, method="average")
for k in (12, 20, int(round(0.1 * n)), int(round(0.25 * n))):
    a = set(m.nsmallest(k, "rank_habitat").k); b = set(m.nsmallest(k, "rank_range").k)
    print(f"top {k:3d}: shared {len(a & b):3d} ({100*len(a & b)/k:.0f}%), expected under independence {k*k/n:.1f}")
m["abs_rank_shift"] = (m.rank_habitat - m.rank_range).abs()
print("median absolute rank shift %.1f of %d" % (m.abs_rank_shift.median(), n))
m.sort_values("rank_habitat").to_csv(os.path.join(REV, "matched_species.csv"), index=False)
print(m.sort_values("rank_habitat").head(12)[["Species", "Taxonomic group", "share_habitat_2020", "share_range_2020", "rank_range", "Modeled habitat (km²)", "Range in Bangladesh (km²)"]].round(2).to_string(index=False))
print(m.sort_values("rank_range").head(12)[["Species", "share_range_2020", "share_habitat_2020", "rank_habitat", "Range in Bangladesh (km²)"]].round(2).to_string(index=False))
