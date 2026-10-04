# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

from figcommon import *
t = pd.read_csv(os.path.join(TAB_MS, "decomposition_300m.csv"))
g20 = t[(t["product"] == "GLC_FCS30D 30 m") & (t.year == 2020)].iloc[0]; h20 = t[(t["product"] == "GHSL 100 m") & (t.year == 2020)].iloc[0]
fig = plt.figure(figsize=(W, 8.7))
gs = GridSpec(3, 3, figure=fig, height_ratios=[1, 1.45, 1.45], hspace=0.38, wspace=0.42, left=0.07, right=0.985, top=0.975, bottom=0.07)
top = gs[0, :].subgridspec(1, 3, width_ratios=[1.4, 1, 1], wspace=0.62)
ax = fig.add_subplot(top[0, 0])
vals = [g20["fine resolution km²"], g20["aggregated to 300 m km²"], h20["fine resolution km²"], h20["aggregated to 300 m km²"], g20["ESA CCI km²"]]
cols = [BLUE, "#9dbbe8", GREEN, "#82c99f", ORANGE]; x = [0, 1, 2, 3, 4]
ax.bar(x, vals, color=cols, width=0.62)
for xi, v in zip(x, vals): ax.text(xi, v + 50, f"{v:,.0f}", ha="center", va="bottom", fontsize=7)
for (x0, x1, v0, v1, c, rc) in [(0, 1, vals[0], vals[1], BLUE, g20["resolution component km²"]), (2, 3, vals[2], vals[3], GREEN, h20["resolution component km²"])]:
    ax.annotate("", xy=(x1 - 0.12, v1 + 1150), xytext=(x0 + 0.33, v0 - 40), arrowprops=dict(arrowstyle="->", color=c, lw=0.9))
    ax.text(x1 + 0.08, v1 + 330, "Resolution\n\u2212" + f"{abs(rc):,.0f} km\u00b2", ha="center", va="bottom", fontsize=7, color=c, fontweight="bold", linespacing=1.0)
ax.text(4.35, 4350, f"Product against ESA CCI:\nGLC_FCS30D −{abs(g20['product component km²']):,.0f} km²\nGHSL −{abs(h20['product component km²']):,.0f} km²",
        ha="right", va="top", fontsize=7, linespacing=1.15)
ax.set_xticks(x); ax.set_xticklabels(["30 m", "300 m", "100 m", "300 m", "300 m"], fontsize=7)
for xc, nm in [(0.5, "GLC_FCS30D"), (2.5, "GHSL"), (4, "ESA CCI")]: ax.text(xc, -0.12, nm, transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=7)
ax.set_ylim(0, 4400); ax.set_ylabel("Urban land (km²)"); ax.tick_params(axis="x", length=0)
ax = fig.add_subplot(top[0, 1])
ca = {"GLC30": pd.read_csv(os.path.join(T, "class_areas_GLC30.csv")), "ESA300": pd.read_csv(os.path.join(T, "class_areas_ESA300.csv")),
      "GHSL100": pd.read_csv(os.path.join(MR_TAB, "class_areas_GHSL100.csv"))}
for p in P3:
    s = ca[p]; s = s[(s.class_id == 1) & (s.year >= 2000)].sort_values("year")
    ax.plot(s.year, s.area_km2, color=PC[p], marker=PM[p], markersize=4)
    nm = {"GLC30": "GLC_FCS30D\n30 m", "GHSL100": "GHSL\n100 m", "ESA300": "ESA CCI\n300 m"}[p]
    ax.text(2021.2, s.area_km2.iloc[-1], nm, color=PC[p], fontsize=7.5, va="center", ha="left")
ax.set_xlim(1998, 2029); ax.set_ylim(0, 4400); ax.set_xticks([2000, 2010, 2020]); ax.set_xlabel("Year"); ax.set_ylabel("Urban land (km²)")
ax = fig.add_subplot(top[0, 2])
for prod, p, ls, mk in [("GLC_FCS30D 30 m", "GLC30", "-", "o"), ("GHSL 100 m", "GHSL100", "--", "s")]:
    s = t[t["product"] == prod].sort_values("year")
    ax.plot(s.year, s["resolution share (%)"], ls, color=PC[p], marker=mk, markersize=4, label={"GLC30": "30 m vs 300 m", "GHSL100": "100 m vs 300 m"}[p])
ax.set_ylim(0, 100); ax.set_xticks([2000, 2010, 2020]); ax.set_xlabel("Year"); ax.set_ylabel("Resolution share of the gap (%)")
ax.legend(loc="lower left", fontsize=7.5)
for a_, l in zip(fig.axes, "abc"): a_.text(-0.02, 1.04, l, transform=a_.transAxes, fontsize=10, fontweight="bold", ha="right", va="bottom")
with rasterio.open(os.path.join(H, "LandArea300_GLC30.tif")) as r: land = r.read(1); ext = [r.bounds.left, r.bounds.right, r.bounds.bottom, r.bounds.top]
land = np.where(land > 0, land, 0.0)
def rd(p):
    with rasterio.open(p) as r: return np.clip(r.read(1).astype(float), 0, 1)
fg = rd(os.path.join(H, "UrbanFrac300_GLC30_2020.tif")) * 100
with rasterio.open(os.path.join(H, "GLC300maj_2020.tif")) as r: cg = (r.read(1) == 1) * 100.0
fh = rd(os.path.join(MR_HARM, "UrbanFrac300_GHSL100_2020.tif")) * 100; ch = (fh > 50) * 100.0
rows = [("GLC_FCS30D, 30 m", fg, "Blues", "Urban land (% of land)"), ("GLC_FCS30D, aggregated\nto 300 m", cg, "Blues", "Urban land (% of land)"),
        ("Removed by\naggregation", np.clip(fg - cg, 0, None), "Reds", "Urban land removed (% of land)"),
        ("GHSL, 100 m", fh, "Greens", "Urban land (% of land)"), ("GHSL, aggregated\nto 300 m", ch, "Greens", "Urban land (% of land)"),
        ("Removed by\naggregation", np.clip(fh - ch, 0, None), "Reds", "Urban land removed (% of land)")]
for i, (lbl, a, cm, lt) in enumerate(rows):
    ax = fig.add_subplot(gs[1 + i // 3, i % 3]); cmap, norm, cols = class_cmap(cm)
    ax.imshow(agg_km(a, land), extent=ext, cmap=cmap, norm=norm, interpolation="nearest")
    draw_bounds(ax); mapax(ax); ax.set_ylim(20.55, 26.7)
    ax.text(0.98, 0.98, lbl, transform=ax.transAxes, ha="right", va="top", fontsize=7.5)
    class_legend(ax, cols, lt)
    ax.text(0.02, 1.0, "defghi"[i], transform=ax.transAxes, fontsize=10, fontweight="bold", ha="left", va="bottom")
save(fig, "Figure1_resolution")
