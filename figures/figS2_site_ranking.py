# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

from figstyle import *

from matplotlib.colors import LogNorm
z = pd.read_csv(f"{T}/pa_kba_zonal_GLC30.csv"); z = z[z.feature_id != "ALL"]
tot = z.groupby(["feature_type", "name", "zone", "year"]).area_km2.sum(); urb = z[z.class_id == 1].groupby(["feature_type", "name", "zone", "year"]).area_km2.sum().reindex(tot.index).fillna(0)
pct = (urb / tot * 100).rename("pct").reset_index(); landc = tot.rename("land").reset_index()
keep = landc[(landc.zone == "core") & (landc.year == 2020) & (landc.land >= 0.5)][["feature_type", "name"]]
pct = pct.merge(keep, on=["feature_type", "name"])
core20 = pct[(pct.zone == "core") & (pct.year == 2020)].set_index(["feature_type", "name"]).pct.sort_values(ascending=False)
buf = pct[pct.zone == "r0_5"].pivot_table(index=["feature_type", "name"], columns="year", values="pct").reindex(core20.index)
from rasterio.features import shapes as rshapes
from shapely.geometry import shape as shp
with rasterio.open(f"{H}/Hotspots_1km_GLC30.tif") as r: cls = r.read(6); tr_ = r.transform
m = (cls == 3).astype(np.uint8); hs = gpd.GeoDataFrame(geometry=[shp(g) for g, v in rshapes(m, mask=m == 1, transform=tr_)], crs=4326).to_crs(32646).union_all()
near = set()
for ft, df_, col in [("PA", pa.to_crs(32646), "NAME_ENG"), ("KBA", kba.to_crs(32646), "IntName")]:
    for _, r_ in df_.iterrows():
        if r_.geometry.distance(hs) <= 5000: near.add((ft, r_[col]))
n = len(core20); fig = plt.figure(figsize=(W, 8.2)); gs = fig.add_gridspec(1, 2, width_ratios=[1.25, 1], wspace=0.06)
ax = fig.add_subplot(gs[0]); title(ax, "a", "Site exposure ranking, 2020")
yy = np.arange(n)
def dots(a, items, ms=6.5, buf=12):
    xs = []; ys = []
    for i, ((ft, nm), v) in enumerate(items):
        col = BLUE if ft == "PA" else "#e0452b"
        a.plot(v, i, "o", color=col, ms=ms, zorder=3, mew=0); xs.append(v); ys.append(i)
        if (ft, nm) in near: a.plot(v, i, "o", mfc="none", mec=INK, mew=1.3, ms=buf, zorder=4)
    a.plot(xs, ys, color="#bdbdbd", lw=0.9, ls=(0, (1, 2)), zorder=2)
dots(ax, list(core20.items()))
ax.set_xlim(-0.3, core20.max() * 1.25); ax.set_yticks(yy); ax.set_yticklabels([nm[:26] for _, nm in core20.index], fontsize=6.4); ax.invert_yaxis(); ax.set_ylim(n - 0.5, -0.5)
ax.set_xlabel("Urban share of core land (%)"); ax.grid(axis="y", color="#e9e9e5", lw=0.6); ax.grid(axis="x", visible=False); ax.tick_params(axis="y", length=0)
for k in [("PA", "Bhawal"), ("PA", "National Botanical Garden"), ("PA", "Purbachal"), ("KBA", "Hail Haor")]:
    if k in core20.index:
        left = k[1] == "Bhawal"
        ax.annotate(k[1], (core20[k], list(core20.index).index(k)), xytext=(-10 if left else 9, 0), textcoords="offset points", ha="right" if left else "left", va="center", fontsize=7, color=INK, fontweight="bold")
NLOW = 24; low = core20.iloc[-NLOW:]
ins = ax.inset_axes([0.47, 0.09, 0.51, 0.40]); ins.set_facecolor("#ffffff")
dots(ins, list(low.items()), ms=4.5, buf=8.5)
ins.set_yticks(range(NLOW)); ins.set_yticklabels([nm[:24] for _, nm in low.index], fontsize=6, color=INK2); ins.invert_yaxis(); ins.set_ylim(NLOW - 0.5, -0.5); ins.tick_params(axis="y", length=0)
ins.set_xlim(-0.001, max(low.max() * 1.15, 0.02)); ins.tick_params(axis="x", labelsize=6); ins.set_xlabel("Urban share (%)", fontsize=6.5, labelpad=2); ins.grid(axis="y", color="#e2e2de", lw=0.5); ins.grid(axis="x", visible=False)
ins.set_title(f"Low-exposure sites  (n = {NLOW})", fontsize=6.8, fontweight="bold", color=INK2, loc="left", pad=3)
for sp_ in ins.spines.values(): sp_.set_color("#cfcfca")
ax.legend(handles=[Line2D([], [], color=BLUE, marker="o", lw=0, ms=6, label="Protected area"), Line2D([], [], color="#e0452b", marker="o", lw=0, ms=6, label="Key Biodiversity Area"),
                   Line2D([], [], mfc="none", mec=INK, mew=1.3, marker="o", lw=0, ms=10, label="Urban expansion hotspot within 5 km")], loc="center right", bbox_to_anchor=(1.0, 0.62), fontsize=7)
ax2 = fig.add_subplot(gs[1], sharey=ax); title(ax2, "b", "0 to 5 km buffer by year")
yrs = list(buf.columns); vals = np.clip(buf.values, 0.05, None)
edges = np.array([yrs[0] - 2.5] + [(a + b) / 2 for a, b in zip(yrs[:-1], yrs[1:])] + [yrs[-1] + 2.5])
im = ax2.pcolormesh(edges, np.arange(n + 1) - 0.5, vals, cmap=LinearSegmentedColormap.from_list("h", ["#ffffff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"]), norm=LogNorm(0.1, 50), shading="flat")
ax2.set_xticks([1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020]); ax2.tick_params(axis="x", labelsize=7); ax2.grid(False); plt.setp(ax2.get_yticklabels(), visible=False); ax2.tick_params(axis="y", length=0)
ax2.axvline(1997.5, color=INK, ls="--", lw=1); ax2.text(1997.5, n - 0.9, "five-year | annual maps", ha="center", va="center", fontsize=6.5, color=INK, path_effects=[PEF.withStroke(linewidth=2.5, foreground=BG)])
cb = fig.colorbar(im, ax=ax2, fraction=0.03, pad=0.02); cb.set_label("Urban share of 0 to 5 km buffer (%)", fontsize=7.5); cb.ax.tick_params(labelsize=7); cb.outline.set_visible(False)
cb.ax.yaxis.set_major_formatter(FuncFormatter(lambda v, p: f"{v:g}"))
ax2.set_xlabel("Year")
save(fig, "FigureS2_site_ranking", supp=True)
