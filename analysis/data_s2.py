# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

import sys, pandas as pd, openpyxl
src, tpl, out = sys.argv[1:4]
d = pd.read_csv(src); d["key"] = d.species.str.replace("_", " ")
d = d.set_index("key")
REPS = ["GLC30", "GLCp90", "GLCp300", "GLCp510", "GLCp990", "GHSL100", "GHSLm300", "GHSLm500", "GHSLm1000", "ESA300"]
wb = openpyxl.load_workbook(tpl); ws = wb["Data S2"]; n = 0; miss = []
for r in range(3, ws.max_row + 1):
    k = ws.cell(r, 1).value
    if k is None: continue
    if k not in d.index: miss.append(k); continue
    row = d.loc[k]; ws.cell(r, 8).value = round(float(row.habitat_km2), 1); c = 9
    for p in REPS:
        for s in ("2000", "2020", "gain"): ws.cell(r, c).value = round(float(row[f"{p}_{s}"]), 4); c += 1
    n += 1
wb.save(out); print(f"{n} rows written, {len(miss)} not found {miss[:5]}")
