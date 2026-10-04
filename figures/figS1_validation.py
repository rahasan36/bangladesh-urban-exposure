# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

import os, sys
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path[:0] = [REPO, os.path.join(REPO, "figures"), os.path.join(REPO, "sdm")]
from config import HARM, MR_TAB, TAB, TAB_MS, OUT_FIG, SDM_RASTERS, SDM_WORK, SDM_DATA, REV, SUPP
from figstyle import *

import pandas as pd
from matplotlib.colors import LogNorm
d = pd.read_csv(f"{P}/Tables/validation_metrics.csv")
nat = d[d.stratum == "national"].set_index(["reference", "ref_year", "product"])

fig, axs = plt.subplots(2, 2, figsize=(W, 6.6), gridspec_kw=dict(wspace=0.42, hspace=0.42, left=0.08, right=0.98, top=0.95, bottom=0.07))
(a, b), (c, e) = axs
GRAYD = "#4d4d4d"; LB = "#9dbbe8"
title(a, "a", "Area on the common 300 m grid")
labels = ["GHSL built-up surface", "GHSL-derived 100 m map (≥20% built-up)", "GAUD urban extent (2000 and 2015)", "GLC_FCS30D 30 m", "GLC_FCS30D 300 m plurality", "ESA CCI 300 m"]
cols = [GRAYD, "#a6a6a6", GREEN, BLUE, LB, ORANGE]
d3 = pd.read_csv(os.path.join(MR_TAB, "decomposition_three_products.csv")).set_index("year")
t1 = pd.read_csv(os.path.join(TAB_MS, "decomposition_300m.csv")); t1 = t1[t1["product"].str.startswith("GLC")].set_index("year")
esa = pd.read_csv(f"{P}/Tables/class_areas_ESA300.csv"); esa = esa[esa.class_id == 1].set_index("year").area_km2
for k, yr in enumerate((2000, 2020)):
    gy = nat.loc[("GAUD", 2000 if yr == 2000 else 2015, "ESA300"), "ref_km2"]
    vals = [nat.loc[("GHSL", yr, "GLC30_fraction"), "ref_km2"], d3.loc[yr, "GHSL100_frac_km2"], gy, nat.loc[("GHSL", yr, "GLC30_fraction"), "prod_km2"],
            t1.loc[yr, "aggregated to 300 m km²"], esa[yr]]
    x = np.arange(6) + k * 7
    a.bar(x, vals, color=cols, width=0.8, zorder=3)
    for xi, v in zip(x, vals): a.text(xi, v + 60, f"{v:,.0f}", ha="center", fontsize=6, color=INK)
    a.text(x.mean(), -350, str(yr), ha="center", fontsize=8, color=INK)
a.set_xticks([]); a.legend(handles=[Patch(color=cc, label=l) for cc, l in zip(cols, labels)], fontsize=6, frameon=False, loc="upper left")
a.set_ylabel("Area (km²)"); a.set_ylim(0, 5600); a.set_xlim(-0.7, 13.2); a.grid(False)
title(b, "b", "GHSL built-up surface captured, 2020")
strata = [("national", "All land"), ("hotspot", "Inside hotspots"), ("outside_hotspot", "Outside hotspots")]
vh = pd.read_csv(os.path.join(SDM_DATA, "validation_sdm_hotspots_2020.csv"))
prods = [("GLC30_fraction", "GLC_FCS30D 30 m, any urban", BLUE), ("GLC30_plurality", "GLC_FCS30D 300 m plurality", LB), ("ESA300", "ESA CCI 300 m", ORANGE)]
for j, (pk, lab, cc) in enumerate(prods):
    vals = [100 * vh[(vh["product"] == pk) & (vh.stratum == s_)].capture_any.iloc[0] for s_, _ in strata]
    x = np.arange(3) * 4 + j
    b.bar(x, vals, color=cc, width=0.85, label=lab, zorder=3)
    for xi, v in zip(x, vals): b.text(xi, v + 1.5, f"{v:.0f}%", ha="center", fontsize=6, color=INK)
b.set_xticks(np.arange(3) * 4 + 1); b.set_xticklabels(["All land", "Inside hotspots", "Outside hotspots"], fontsize=7); b.grid(False)
b.set_ylabel("Built-up surface captured (%)"); b.set_ylim(0, 125); b.legend(fontsize=6, frameon=False, loc="upper right", ncol=1)
title(c, "c", "GLC_FCS30D against GHSL, 2020")
g, _ = read300(f"{H}/GHSL_frac300_2020.tif"); u, _ = read300(f"{H}/UrbanFrac300_GLC30_2020.tif"); land, _ = read300(f"{H}/LandArea300_GLC30.tif")
m = (land > 0) & ~np.isnan(g) & ~np.isnan(u) & ((g > 0) | (u > 0))
hb = c.hexbin(g[m] * 100, u[m] * 100, gridsize=45, norm=LogNorm(), cmap="Blues", mincnt=1, linewidths=0.1, extent=(0, 100, 0, 100))
c.plot([0, 100], [0, 100], color=INK, lw=0.6, ls="--")
r = nat.loc[("GHSL", 2020, "GLC30_fraction"), "pearson_r"]
c.text(0.04, 0.95, f"Pearson r = {r:.2f}\n{m.sum():,} cells with built-up or urban land", transform=c.transAxes, fontsize=6.5, va="top", color=INK)
c.set_xlabel("GHSL built-up surface (% of cell)"); c.set_ylabel("GLC_FCS30D urban land (% of cell)")
c.grid(False); cax = c.inset_axes([0.64, 0.16, 0.30, 0.035]); cb = fig.colorbar(hb, cax=cax, orientation="horizontal"); cb.set_label("Cells per hexagon", fontsize=6); cb.ax.tick_params(labelsize=6, length=2); cb.outline.set_visible(False)
title(e, "d", "Detection of built-up cells, 2020")
rows = [("GLC30_fraction", "GLC 30 m", BLUE), ("GLC30_majority", "GLC 300 m plurality", LB), ("ESA300", "ESA CCI 300 m", ORANGE)]
y = 0
yt, yl = [], []
for pk, lab, cc in rows:
    for t, mk in ((0.05, "o"), (0.25, "s")):
        pr = nat.loc[("GHSL", 2020, pk), f"producers_t{t}"] * 100; us = nat.loc[("GHSL", 2020, pk), f"users_t{t}"] * 100
        e.plot([pr], [y], marker=mk, color=cc, ms=6, mec=INK, mew=0.4, ls="none")
        e.plot([us], [y], marker=mk, color="white", ms=6, mec=cc, mew=1.2, ls="none")
        e.hlines(y, min(pr, us), max(pr, us), color=cc, lw=1)
        yt.append(y); yl.append(f"{lab}\nBuilt-up ≥ {int(t*100)}%"); y -= 1
    y -= 0.5
e.set_yticks(yt); e.set_yticklabels(yl, fontsize=6); e.set_xlim(0, 105); e.set_ylim(min(yt) - 1.6, 0.6); e.set_xlabel("Percent of cells"); e.grid(False)
e.legend(handles=[Line2D([], [], marker="o", color=GRAY, ls="none", label="Producer's: built-up cells mapped urban"),
                  Line2D([], [], marker="o", color="white", mec=GRAY, ls="none", label="User's: urban cells that are built-up")],
         fontsize=6, frameon=False, loc="lower left", bbox_to_anchor=(0.0, -0.02))
for ax in (a, b, e):
    ax.spines[["top", "right"]].set_visible(False)
save(fig, "FigureS1_validation", supp=True)
