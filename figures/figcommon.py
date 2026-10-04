# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import *
from config import SDM_DATA, SDM_WORK, SDM_RASTERS, REV, TAB_MS, MR_HARM, MR_TAB, SUPP
from matplotlib.patches import Patch
plt.rcParams.update({"lines.linewidth": 1.2, "lines.markersize": 4, "axes.grid": False})
P3 = ["GLC30", "GHSL100", "ESA300"]
PC = {"GLC30": BLUE, "GHSL100": GREEN, "ESA300": ORANGE}
PN = {"GLC30": "GLC_FCS30D 30 m", "GHSL100": "GHSL 100 m", "ESA300": "ESA CCI 300 m"}
PM = {"GLC30": "o", "GHSL100": "s", "ESA300": "^"}
BINS = [-1e-9, 1e-9, 1, 2.5, 5, 10, 25, 1e9]
BLAB = ["0", "0 to 1", "1 to 2.5", "2.5 to 5", "5 to 10", "10 to 25", "> 25"]
def class_cmap(name):
    cols = ["#f0f0f0"] + [plt.get_cmap(name)(v) for v in np.linspace(0.18, 1.0, 6)]
    return ListedColormap(cols), BoundaryNorm(BINS, 7), cols
def class_legend(ax, cols, title_):
    h = [Patch(facecolor=c, edgecolor="#bbbbbb", linewidth=0.3) for c in cols]
    leg = ax.legend(h, BLAB, title=title_, ncol=4, loc="upper center", bbox_to_anchor=(0.5, 0.0), fontsize=7, title_fontsize=7.5,
                    handlelength=1.4, handleheight=0.9, columnspacing=0.8, labelspacing=0.25, borderaxespad=0.2)
    return leg
def bd_layers():
    return bnd, adm
def draw_bounds(ax, color="#9a9a9a", lw=0.45):
    adm.boundary.plot(ax=ax, color=color, linewidth=lw); bnd.boundary.plot(ax=ax, color=color, linewidth=lw)
def sp3():
    return pd.read_csv(os.path.join(SDM_DATA, "sdm_species.csv"))
def lab(v, nd=0): return f"{v:,.{nd}f}"
def agg_km(a, land, f=3):
    h, w = a.shape[0] // f, a.shape[1] // f
    L = land[:h * f, :w * f].reshape(h, f, w, f).sum(axis=(1, 3)); U = (a * land)[:h * f, :w * f].reshape(h, f, w, f).sum(axis=(1, 3))
    return np.where(L > 0, U / np.where(L > 0, L, 1), np.nan)
