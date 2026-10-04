# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

import os, sys
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path[:0] = [REPO, os.path.join(REPO, "figures"), os.path.join(REPO, "sdm")]
from config import HARM, MR_TAB, TAB, OUT_FIG, SDM_RASTERS, SDM_WORK, SDM_DATA, REV, SUPP
import pandas as pd, numpy as np, statsmodels.formula.api as smf, warnings
raw = pd.read_excel(os.path.join(SUPP, "Data_S2_species_habitat_models.xlsx"), sheet_name="Data S2", header=None)
g = raw.iloc[0].ffill(); d = raw.iloc[2:].reset_index(drop=True)
d.columns = [f"{a}|{str(b).replace('.0','')}" if j >= 8 else str(raw.iloc[1, j]) for j, (a, b) in enumerate(zip(g, raw.iloc[1]))]
d = d[d["National threat status"] != "Data Deficient"].copy()
d["thr"] = (d["National threat status"] == "Threatened").astype(int); d["la"] = np.log10(pd.to_numeric(d["Modeled habitat (km²)"]))
d["grp"] = d["Taxonomic group"]; d["fam"] = d["Family"]
rows = []
for col, lab in [("GLC_FCS30D 30 m", "30 m"), ("GHSL-derived 100 m", "100 m"), ("ESA CCI 300 m", "300 m")]:
    d["y"] = pd.to_numeric(d[f"{col}|Change 2000 to 2020"])
    ols = smf.ols("y ~ thr + la + C(grp)", d).fit(cov_type="HC3")
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        mm = smf.mixedlm("y ~ thr + la + C(grp)", d, groups=d["fam"]).fit(reml=True)
    ci = mm.conf_int().loc["thr"]; oci = ols.conf_int().loc["thr"]
    rows.append(dict(product=lab, n=len(d), families=d.fam.nunique(), ols_thr=ols.params["thr"], ols_lo=oci[0], ols_hi=oci[1],
                     mixed_thr=mm.params["thr"], mixed_lo=ci[0], mixed_hi=ci[1], family_var=float(mm.cov_re.iloc[0, 0]), converged=mm.converged, warnings=len(w)))
t = pd.DataFrame(rows); t.to_csv(os.path.join(REV, "family_mixed_model.csv"), index=False); print(t.round(3).to_string(index=False))
