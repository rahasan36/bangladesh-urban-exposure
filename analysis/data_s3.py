# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

import sys, pandas as pd, openpyxl
src, tpl, out = sys.argv[1:4]
s = pd.read_csv(src); s = s[s.year.isin([2000, 2020])].pivot_table(index="id_no", columns=["dataset", "year"], values="urban_pct")
wb = openpyxl.load_workbook(tpl); ws = wb["Data S3"]; n = 0; miss = []
for r in range(3, ws.max_row + 1):
    k = ws.cell(r, 1).value
    if k is None: continue
    if k not in s.index: miss.append(k); continue
    for ds, c0 in [("GLC30", 7), ("ESA300", 13)]:
        a, b = float(s.loc[k, (ds, 2000)]), float(s.loc[k, (ds, 2020)])
        ws.cell(r, c0).value = round(a, 4); ws.cell(r, c0 + 1).value = round(b, 4); ws.cell(r, c0 + 2).value = round(round(b, 4) - round(a, 4), 4)
    n += 1
wb.save(out); print(f"{n} rows written, {len(miss)} not found")
