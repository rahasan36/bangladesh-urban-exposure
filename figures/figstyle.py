# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

import os, sys, numpy as np, pandas as pd, geopandas as gpd, rasterio
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from config import OUT as P, TAB as T, SP as S, HARM as H, SPS as SS, OUT_FIG, ADM1_OUT, MR, MR_HARM, MR_TAB, MR_SP, MR_REVIEW, TAB_MS
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm, LinearSegmentedColormap, TwoSlopeNorm, LogNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter
from matplotlib.gridspec import GridSpec
import matplotlib.patheffects as PEF
from matplotlib import font_manager as _fm
if not any("Arial" in f.name for f in _fm.fontManager.ttflist):
    print("NOTE: Arial not found; falling back to Liberation Sans / DejaVu Sans. Install Arial (or ttf-mscorefonts) for the exact house style.")
    plt.rcParams["font.family"] = "sans-serif"; plt.rcParams["font.sans-serif"] = ["Liberation Sans", "DejaVu Sans"]
    ARIAL = plt.rcParams["font.sans-serif"][0]
else:
    ARIAL = "Arial"
FIG_MAIN = os.path.join(OUT_FIG, "main"); FIG_SUPP = os.path.join(OUT_FIG, "supplementary")
os.makedirs(FIG_MAIN, exist_ok=True); os.makedirs(FIG_SUPP, exist_ok=True)
BG = "#ffffff"; INK = "#000000"; INK2 = "#000000"; GRID = "#e6e6e6"
BLUE, ORANGE, PURPLE, YELLOW, GRAY, NAVY, PINK, PEACH, GREEN, TEAL = "#3b7dd8", "#e8632b", "#7b2c8c", "#e6a100", "#8c8c8c", "#1f1f5e", "#d94a6a", "#f4a26a", "#2a9d5c", "#1a9e9a"
C = {"GLC30": BLUE, "ESA300": ORANGE}; DSN = {"GLC30": "GLC_FCS30D, 30 m", "ESA300": "ESA CCI, 300 m"}
SCOL = {"Least Concern": PEACH, "Near Threatened": PINK, "Threatened": PURPLE, "Data Deficient": GRAY}
W = 183 / 25.4
plt.rcParams.update({"font.family": ARIAL, "font.size": 8.5, "axes.titlesize": 10, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": INK, "axes.linewidth": 0.8, "axes.labelcolor": INK,
    "xtick.color": INK, "ytick.color": INK, "xtick.labelsize": 8, "ytick.labelsize": 8, "axes.labelsize": 8.5,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "legend.frameon": False, "legend.fontsize": 8,
    "pdf.fonttype": 42, "ps.fonttype": 42, "figure.facecolor": BG, "axes.facecolor": BG, "savefig.facecolor": BG, "savefig.dpi": 300, "lines.linewidth": 2.6, "lines.markersize": 7})
fmt = FuncFormatter(lambda v, p: f"{v:g}")
def title(ax, letter, text): ax.set_title(f"{letter}.  {text}", pad=8, color=INK)
def fit_titles(fig, min_size=8.0):
    fig.canvas.draw(); r = fig.canvas.get_renderer()
    for ax in fig.axes:
        for t in (ax.title, getattr(ax, "_left_title", None), getattr(ax, "_right_title", None)):
            if t is None or not t.get_text(): continue
            axw = ax.get_window_extent(r).width; size = t.get_fontsize()
            while t.get_window_extent(r).width > axw and size > min_size:
                size -= 0.5; t.set_fontsize(size); fig.canvas.draw()
            if t.get_window_extent(r).width > axw or size < plt.rcParams["axes.titlesize"]:
                print(f"  title fitted: '{t.get_text()}' -> {size} pt" + ("" if t.get_window_extent(r).width <= axw else " (STILL TOO WIDE)"))

def save(fig, name, supp=False):
    d = FIG_SUPP if supp else FIG_MAIN
    fit_titles(fig)
    fig.savefig(f"{d}/{name}.png", bbox_inches="tight", pad_inches=0.04); plt.close(fig); print("saved", name)
def mapax(ax):
    ax.grid(False); ax.set_axis_off(); ax.set_xlim(88.0, 92.7); ax.set_ylim(20.6, 26.65); ax.set_aspect("equal"); ax.set_anchor("N")
def read300(path, band=1):
    with rasterio.open(path) as r:
        a = r.read(band).astype(float); a[a == r.nodata] = np.nan; ext = [r.bounds.left, r.bounds.right, r.bounds.bottom, r.bounds.top]
    return a, ext
def block_mean(a, f):
    h, w = a.shape[0] // f, a.shape[1] // f
    return np.nanmean(a[:h * f, :w * f].reshape(h, f, w, f), axis=(1, 3))
bnd = gpd.read_file(f"{P}/BGD_Boundary.geojson"); adm = gpd.read_file(f"{ADM1_OUT}/geoBoundaries-BGD-ADM1.geojson")
adm["shapeName"] = adm.shapeName.replace({"Rajshani": "Rajshahi"})
pa = gpd.read_file(f"{P}/BGD_PA_Jan2026.geojson"); kba = gpd.read_file(f"{P}/BGD_KBA_2021.geojson")
g = pd.read_csv(f"{T}/class_areas_GLC30.csv"); e = pd.read_csv(f"{T}/class_areas_ESA300.csv")
gu = g[g.class_id == 1].sort_values("year").set_index("year").area_km2; eu = e[e.class_id == 1].sort_values("year").set_index("year").area_km2
sp = pd.read_csv(f"{S}/species_table.csv"); sx = pd.read_csv(f"{S}/species_urban_exposure.csv"); sl = pd.read_csv(f"{S}/species_slopes_ranking.csv")
reg = pd.read_csv(f"{SS}/species_region_share.csv"); hs = pd.read_csv(f"{SS}/species_hotspot_overlap.csv")
g30 = sl[sl.dataset == "GLC30"].merge(reg, on="id_no").merge(hs, on="id_no")
SUND = ["Panthera tigris", "Batagur baska", "Chelonia mydas", "Caretta caretta", "Dermochelys coriacea", "Eretmochelys imbricata", "Lepidochelys olivacea"]
RIVER = ["Batagur dhongoka", "Batagur kachuga", "Nilssonia gangetica", "Nilssonia hurum", "Nilssonia nigricans", "Pangshura tecta", "Pangshura sylhetensis", "Hardella thurjii",
         "Geoclemys hamiltonii", "Morenia petersi", "Gavialis gangeticus", "Pelochelys cantorii", "Lissemys punctata", "Crocodylus palustris", "Varanus flavescens", "Melanochelys tricarinata", "Cuora amboinensis", "Caprolagus hispidus", "Xenochrophis cerasogaster"]
def assemblage(n): return "Sundarbans and coast" if n in SUND else ("Rivers and floodplain" if n in RIVER else "Hill forest")
ACOL = {"Sundarbans and coast": TEAL, "Rivers and floodplain": BLUE, "Hill forest": GREEN}
