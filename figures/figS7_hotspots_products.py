# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

import os, sys
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path[:0] = [REPO, os.path.join(REPO, "figures"), os.path.join(REPO, "sdm")]
from config import HARM, MR_TAB, TAB, OUT_FIG, SDM_RASTERS, SDM_WORK, SDM_DATA, REV, SUPP
os.chdir(SDM_WORK)
import os, sys, numpy as np, pandas as pd, geopandas as gpd, rasterio
import figstyle as fs; from figstyle import *
from matplotlib.colors import ListedColormap; from matplotlib.patches import Patch
z = np.load("urban_on_sdm_grid.npz"); valid = np.load("valid.npy"); U = np.load("urban_sdm_layers.npz"); irr = z["irr"]; CA = z["cellkm"]; L = z["GLC30_2020_l"]
A = valid & (irr > -1e30) & (L / CA > 0.5)
def _nat(f, ds):
    c = pd.read_csv(f); c = c[(c.class_id == 1) & (c.dataset == ds)].set_index("year").area_km2; return c[2020] - c[2000]
NAT = {"GLC30": _nat(os.path.join(TAB, "class_areas_GLC30.csv"), "GLC30"), "GHSL100": _nat(os.path.join(MR_TAB, "class_areas_GHSL100.csv"), "GHSL100"),
       "ESA300": _nat(os.path.join(TAB, "class_areas_ESA300.csv"), "ESA300")}
P = ["GLC30", "GHSL100", "ESA300"]; C3 = {"GLC30": BLUE, "GHSL100": GREEN, "ESA300": ORANGE}; N3 = {"GLC30": "GLC_FCS30D 30 m", "GHSL100": "GHSL-derived 100 m", "ESA300": "ESA CCI 300 m"}
bnd = gpd.read_file(f"{REPO}/data/pipeline/BGD_Boundary.geojson"); from rasterio.features import rasterize
with rasterio.open("hotspots_sdm_GLC30.tif") as r0: land = rasterize(bnd.to_crs(r0.crs).geometry, out_shape=r0.shape, transform=r0.transform) == 1
adm = gpd.read_file(f"{REPO}/data/pipeline/ADM1/geoBoundaries-BGD-ADM1.geojson")
fig = plt.figure(figsize=(W, 3.6)); gs = fig.add_gridspec(1, 3, wspace=0.04, left=0.0, right=0.72, top=0.95, bottom=0.2); gd = fig.add_gridspec(1, 1, left=0.81, right=0.995, top=0.95, bottom=0.2)
rows = []
for i, p in enumerate(P):
    with rasterio.open(f"hotspots_sdm_{p}.tif") as r: h = r.read(1) == 1; ext = [r.bounds.left, r.bounds.right, r.bounds.bottom, r.bounds.top]
    g = U[f"{p}_2020"] - U[f"{p}_2000"]
    tot = NAT[p]; inh = g[h].sum(); dom = g[A].sum()
    rows.append(dict(product=p, hotspot_km2=L[h].sum(), gain_total_km2=tot, gain_hotspots_km2=inh, gain_domain_other_km2=dom - inh, gain_outside_domain_km2=tot - dom))
    ax = fig.add_subplot(gs[i]); mapax(ax); ax.set_xlim(88.0, 92.69); ax.set_ylim(20.72, 26.64); title(ax, "abc"[i], "")
    k = np.full(h.shape, np.nan); k[land & ~A] = 1; k[A] = 0; k[h] = 2
    ax.imshow(np.ma.masked_invalid(k), extent=ext, cmap=ListedColormap(["#ececea", "#bdbdb8", C3[p]]), vmin=-0.5, vmax=2.5, interpolation="nearest")
    adm.boundary.plot(ax=ax, color="#7f7f7f", lw=0.3); bnd.boundary.plot(ax=ax, color=INK, lw=0.5)
    ax.legend(handles=[Patch(color=C3[p], label=f"Hotspot, {L[h].sum():,.0f} km²"), Patch(color="#ececea", label="Analysis domain"), Patch(color="#bdbdb8", label="Outside domain")],
              fontsize=6, loc="upper left", bbox_to_anchor=(0.03, 0.0), handlelength=1.0, labelspacing=0.2, borderaxespad=0.1, frameon=False, title=N3[p], title_fontsize=6.5, alignment="left")
    ax.set_xlim(88.0, 92.69); ax.set_ylim(20.72, 26.64); ax.set_aspect(1 / np.cos(np.deg2rad(23.7)), adjustable="box"); ax.set_anchor("N")
t = pd.DataFrame(rows); t.to_csv(os.path.join(SDM_DATA, "hotspots_products_sdm.csv"), index=False); print(t.round(1).to_string(index=False))
ax = fig.add_subplot(gd[0]); title(ax, "d", "")
for i, r in t.iterrows():
    c = C3[r["product"]]; parts = [(r.gain_hotspots_km2, 1.0, "white"), (r.gain_domain_other_km2, 0.35, INK), (r.gain_outside_domain_km2, 0.12, INK)]; b0 = 0
    for v, al, tc in parts:
        ax.bar(i, v, bottom=b0, color=c, alpha=al, width=0.7, lw=0)
        if v / r.gain_total_km2 > 0.06: ax.text(i, b0 + v / 2, f"{100 * v / r.gain_total_km2:.0f}%", ha="center", va="center", fontsize=6, color=tc)
        b0 += v
    ax.text(i, b0 + 15, f"{r.gain_total_km2:,.0f}", ha="center", va="bottom", fontsize=6, color=INK)
ax.set_xticks([0, 1, 2]); ax.set_xticklabels(["30 m", "100 m", "300 m"]); ax.set_xlim(-0.6, 2.6); ax.set_ylabel("Urban gain, 2000 to 2020 (km²)"); ax.grid(False)
ax.set_ylim(0, 1750)
ax.legend(handles=[Patch(color=GRAY, alpha=1.0, label="Inside hotspots"), Patch(color=GRAY, alpha=0.35, label="Elsewhere in domain"), Patch(color=GRAY, alpha=0.12, label="Outside domain")],
          fontsize=6, loc="upper right", frameon=False, handlelength=1.0)
fs.save(fig, "FigureS7_hotspots_products", supp=True)
