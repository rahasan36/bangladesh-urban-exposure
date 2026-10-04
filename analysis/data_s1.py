# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

import os, sys
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path[:0] = [REPO, os.path.join(REPO, "figures"), os.path.join(REPO, "sdm")]
from config import HARM, MR_TAB, TAB, OUT_FIG, SDM_RASTERS, SDM_WORK, SDM_DATA, REV, SUPP
import pandas as pd, numpy as np, os
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment, PatternFill, Border, Side
from openpyxl.utils import get_column_letter
src = {"GLC30": os.path.join(TAB, "pa_kba_zonal_GLC30.csv"), "GHSL100": os.path.join(MR_TAB, "pa_kba_zonal_GHSL100.csv"), "ESA300": os.path.join(TAB, "pa_kba_zonal_ESA300.csv")}
lab = {"GLC30": "GLC_FCS30D 30 m", "GHSL100": "GHSL-derived 100 m", "ESA300": "ESA CCI 300 m"}
ids = None; tabs = {}
for k, f in src.items():
    z = pd.read_csv(f); z = z[z.feature_id.astype(str) != "ALL"]; z = z[z.year.isin([2000, 2020])]
    if ids is None: ids = z.groupby(["feature_type", "name"]).feature_id.apply(lambda s: ", ".join(sorted(set(map(str, s)))))
    t = z.groupby(["feature_type", "name", "zone", "year"]).area_km2.sum(); u = z[z.class_id == 1].groupby(["feature_type", "name", "zone", "year"]).area_km2.sum().reindex(t.index).fillna(0)
    tabs[k] = pd.DataFrame({"land": t, "urban": u}).reset_index()
keys = pd.concat([t[["feature_type", "name", "zone"]] for t in tabs.values()]).drop_duplicates()
out = keys.copy()
for k, t in tabs.items():
    for y in (2000, 2020):
        s = t[t.year == y].set_index(["feature_type", "name", "zone"])
        out = out.merge((100 * s.urban / s.land.where(s.land > 0)).rename(f"{k}_{y}").reset_index(), how="left", on=["feature_type", "name", "zone"])
    s = t[(t.year == 2020)].set_index(["feature_type", "name", "zone"]).land.rename(f"{k}_land").reset_index()
    out = out.merge(s, how="left", on=["feature_type", "name", "zone"])
core = out[out.zone == "core"].set_index(["feature_type", "name"])
reason = {}
for (ft, nm), r in core.iterrows():
    if not (r.GLC30_land >= 0.5): reason[(ft, nm)] = "core land below 0.5 km2 at 30 m"
    elif pd.isna(r.ESA300_2020) or not (r.ESA300_land > 0): reason[(ft, nm)] = "no core land on the ESA CCI grid" + (" (marine site)" if "Marine" in nm else "")
    elif pd.isna(r.GHSL100_2020): reason[(ft, nm)] = "no core land on the GHSL grid"
for ft, nm in out[["feature_type", "name"]].drop_duplicates().itertuples(index=False):
    if (ft, nm) not in core.index: reason[(ft, nm)] = "no core land on any grid"
out["in_comparison"] = [((ft, nm) not in reason) for ft, nm in zip(out.feature_type, out.name)]
out["exclusion"] = [reason.get((ft, nm), "") for ft, nm in zip(out.feature_type, out.name)]
out["ids"] = [ids.get((ft, nm), "") for ft, nm in zip(out.feature_type, out.name)]
zo = {"core": 0, "r0_5": 1, "r5_10": 2, "r10_25": 3}; zl = {"core": "Core", "r0_5": "Buffer 0 to 5 km", "r5_10": "Buffer 5 to 10 km", "r10_25": "Buffer 10 to 25 km"}
out = out.sort_values(["feature_type", "name", "zone"], key=lambda c: c.map(zo) if c.name == "zone" else (c.map({"PA": 0, "KBA": 1}) if c.name == "feature_type" else c)).reset_index(drop=True)
n_sites = out[["feature_type", "name"]].drop_duplicates().shape[0]; n_in = out[out.in_comparison][["feature_type", "name"]].drop_duplicates().shape[0]
print("rows", len(out), "sites", n_sites, "in comparison", n_in); print(out[~out.in_comparison][["name", "exclusion"]].drop_duplicates().to_string(index=False))
FN = "Arial"; HF = PatternFill("solid", fgColor="D9D9D9"); GF = PatternFill("solid", fgColor="EFEFEF"); BOT = Border(bottom=Side(style="thin", color="808080"))
wb = Workbook(); ws = wb.active; ws.title = "README"
rd = [("Supplemental Data S1", None), ("Mapping hidden urban exposure around Bangladesh's protected areas and species habitats", None), ("Rakibul Ahasan, Palash Chandra Das, Shawan Chowdhury", None),
      ("One Earth (submitted). Corresponding author: Rakibul Ahasan (rahasan@esri.com)", None), (None, None),
      ("Description", f"Urban share of land (percent) in the core and in the 0 to 5, 5 to 10, and 10 to 25 km buffers of the {n_sites} protected areas and Key Biodiversity Areas of Bangladesh in 2000 and 2020, and the 2000 to 2020 change, under three urban maps. Sites are grouped by type and name, the unit used in the manuscript, so a site made of several polygons appears once."),
      ("Urban maps", "GLC_FCS30D 30 m land cover (Zhang et al. 2024), GHS-BUILT-S R2023A 100 m built-up surface labeled urban where at least 20 percent of a cell is built (Pesaresi and Politis 2023), and ESA CCI and C3S 300 m land cover (ESA 2017)."),
      ("Sites", "Protected areas from the World Database on Protected Areas (January 2026) and Key Biodiversity Areas from the World Database of Key Biodiversity Areas (March 2021). Site IDs list the WDPA or KBA identifiers of the polygons that make up the site."),
      ("Zones", "Core is the site polygon. Buffers are distance bands from the site boundary built in UTM zone 46N and exclude the cores of all sites of the same type."),
      ("Calculation", "Urban share is urban area divided by land area in the zone, with cell areas computed on each product's native grid. Change is the 2020 share minus the 2000 share, in percentage points. Core land is the 2020 land area of the core on each product's grid."),
      ("Site comparison", f"In comparison marks the {n_in} sites used for the rank correlations and zero-core counts in the main text: core land of at least 0.5 km2 at 30 m and core land on all three grids. Exclusion gives the reason for the other sites."),
      ("Missing values", "Blank cells mark zones with no land on that product's grid."),
      (None, None), ("Column dictionary", None)]
r = 1
for k, v in rd:
    if k: ws.cell(r, 1, k)
    if v: ws.cell(r, 2, v)
    r += 1
cols = [("Site type", "Protected area or Key Biodiversity Area", "Text"), ("Site name", "Name in the source database", "Text"), ("Site IDs", "Identifiers of the polygons in the site", "Text"),
        ("Zone", "Core or buffer", "Text"), ("2000, 2020", "Urban share of land in the zone under the product named in the group header", "Percent of land"), ("Change", "2020 share minus 2000 share", "Percentage points"),
        ("Core land", "2020 land area of the zone on the product's grid", "km2"), ("In comparison", "Site used in the rank correlations and zero-core counts (n = %d)" % n_in, "TRUE or FALSE"), ("Exclusion", "Reason a site is outside the comparison", "Text")]
for j, h in enumerate(["Sheet column", "Description", "Units"], 1):
    c = ws.cell(r, j, h); c.font = Font(name=FN, bold=True, size=10); c.fill = HF
for i, row in enumerate(cols, 1):
    for j, v in enumerate(row, 1): ws.cell(r + i, j, v)
for row in ws.iter_rows():
    for c in row:
        c.alignment = Alignment(wrap_text=True, vertical="top")
        if c.font is None or not c.font.bold: c.font = Font(name=FN, size=10, bold=(c.column == 1 and c.row > 5))
ws["A1"].font = Font(name=FN, size=12, bold=True); ws["A2"].font = Font(name=FN, size=11, bold=True); ws.column_dimensions["A"].width = 22; ws.column_dimensions["B"].width = 100; ws.column_dimensions["C"].width = 18
ds = wb.create_sheet("Data S1"); hdr1 = ["", "", "", ""]; hdr2 = ["Site type", "Site name", "Site IDs", "Zone"]
for k in ["GLC30", "GHSL100", "ESA300"]: hdr1 += [lab[k]] * 4; hdr2 += ["2000", "2020", "Change", "Core land (km2)"]
hdr1 += ["", ""]; hdr2 += ["In comparison", "Exclusion"]
for j, (a, b) in enumerate(zip(hdr1, hdr2), 1):
    for rr, v in ((1, a), (2, b)):
        c = ds.cell(rr, j, v if v else None); c.font = Font(name=FN, bold=True, size=10); c.fill = GF if rr == 1 else HF; c.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
        if rr == 2: c.border = BOT
for g in range(3): ds.merge_cells(start_row=1, start_column=5 + 4 * g, end_row=1, end_column=8 + 4 * g)
for i, x in enumerate(out.itertuples(index=False), 3):
    vals = [{"PA": "Protected area", "KBA": "Key Biodiversity Area"}[x.feature_type], x.name, x.ids, zl[x.zone]]
    for k in ["GLC30", "GHSL100", "ESA300"]:
        a, b = getattr(x, f"{k}_2000"), getattr(x, f"{k}_2020"); land = getattr(x, f"{k}_land")
        vals += [None if pd.isna(a) else round(a, 3), None if pd.isna(b) else round(b, 3), None if (pd.isna(a) or pd.isna(b)) else round(b - a, 3), None if (pd.isna(land) or x.zone != "core") else round(land, 3)]
    vals += [bool(x.in_comparison), x.exclusion or None]
    for j, v in enumerate(vals, 1):
        c = ds.cell(i, j, v); c.font = Font(name=FN, size=10)
        if isinstance(v, float): c.number_format = "0.000"
ds.freeze_panes = "C3"; ds.auto_filter.ref = f"A2:{get_column_letter(len(hdr2))}{len(out) + 2}"; ds.row_dimensions[2].height = 30
for j, w in enumerate([20, 34, 22, 16] + [10] * 12 + [12, 40], 1): ds.column_dimensions[get_column_letter(j)].width = w
wb.save(os.path.join(os.environ.get("BGD_DATAS1_OUT", SUPP), "Data_S1_protected_areas_KBAs.xlsx"))
from scipy.stats import spearmanr
c = out[(out.zone == "core") & out.in_comparison]
print("from workbook: n", len(c), "rho 30-100 %.2f, 30-300 %.2f" % (spearmanr(c.GLC30_2020, c.GHSL100_2020)[0], spearmanr(c.GLC30_2020, c.ESA300_2020)[0]), "zeros", (c.GLC30_2020 == 0).sum(), (c.GHSL100_2020 == 0).sum(), (c.ESA300_2020 == 0).sum())
