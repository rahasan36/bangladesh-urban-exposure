# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

from figcommon import *
d = pd.read_csv(os.path.join(TAB_MS, "pa_kba_core_buffer.csv")); prof = pd.read_csv(os.path.join(TAB_MS, "pa_distance_profile_2020.csv"))
YRS = {"GLC30": [2000, 2004, 2005, 2008, 2010, 2012, 2015, 2016, 2020], "GHSL100": [2000, 2005, 2010, 2015, 2020], "ESA300": [2000, 2004, 2008, 2012, 2016, 2020]}
def fan(ax, x, items, side, ymin, ymax, gap):
    items = sorted(items, key=lambda t: t[0]); pos = []
    for y, txt, c in items:
        p = max(y, pos[-1] + gap) if pos else max(y, gap); pos.append(p)
    over = pos[-1] - ymax if pos and pos[-1] > ymax else 0
    pos = [p - over for p in pos]
    for (y, txt, c), p in zip(items, pos):
        dx = -1.0 if side == "left" else 1.0
        ax.plot([x, x + dx * 0.9, x + dx * 1.5], [y, p, p], color=c, lw=0.6); ax.text(x + dx * 1.65, p, txt, color=c, fontsize=7, va="center", ha="right" if side == "left" else "left")
fig = plt.figure(figsize=(W, 5.6))
gs = GridSpec(2, 2, figure=fig, hspace=0.38, wspace=0.30, left=0.09, right=0.98, top=0.97, bottom=0.08)
for k, (ft, ttl, ymax, gap) in enumerate([("PA", "Protected areas", 5, 0.21), ("KBA", "Key Biodiversity Areas", 3, 0.125)]):
    ax = fig.add_subplot(gs[0, k]); L, R = [], []
    for p in P3:
        for zone, ls, mfc in [("r0_10", "-", PC[p]), ("core", "--", "white")]:
            s = d[(d["product"] == p) & (d.feature_type == ft) & (d.zone == zone) & (d.year.isin(YRS[p]))].sort_values("year")
            ax.plot(s.year, s.urban_pct, ls, color=PC[p], marker="o" if zone == "r0_10" else "s", markersize=3.5 if zone == "r0_10" else 3.2, mfc=mfc, mew=1.0)
            L.append((s.urban_pct.iloc[0], f"{s.urban_pct.iloc[0]:.2f}", PC[p])); R.append((s.urban_pct.iloc[-1], f"{s.urban_pct.iloc[-1]:.2f}", PC[p]))
    fan(ax, 2000, L, "left", 0, ymax, gap); fan(ax, 2020, R, "right", 0, ymax, gap)
    ax.set_xlim(1996, 2024); ax.set_ylim(0, ymax); ax.set_xticks([2000, 2005, 2010, 2015, 2020]); ax.set_xlabel("Year"); ax.set_ylabel("Urban share of land (%)")
    ax.grid(axis="y", color=GRID, lw=0.6); ax.text(0.02, 0.98, ttl, transform=ax.transAxes, va="top", fontsize=8)
    if ft == "KBA":
        h1 = [Line2D([], [], color=PC[p], lw=1.4) for p in P3]; h2 = [Line2D([], [], color="k", marker="o", ms=3.5, lw=1.2), Line2D([], [], color="k", ls="--", marker="s", mfc="white", ms=3.2, lw=1.2)]
        lg = ax.legend(h1, [PN[p] for p in P3], loc="upper left", bbox_to_anchor=(0.0, 0.93), fontsize=7.5, handlelength=1.8); ax.add_artist(lg)
        ax.legend(h2, ["0 to 10 km buffer", "Core"], loc="upper left", bbox_to_anchor=(0.5, 0.93), fontsize=7.5, handlelength=1.8)
    ax.text(-0.02, 1.02, "ab"[k], transform=ax.transAxes, fontsize=10, fontweight="bold", ha="right", va="bottom")
ax = fig.add_subplot(gs[1, 0]); rows = []
ZN = {"core": "core", "r0_5": "0 to 5 km", "r5_10": "5 to 10 km", "r10_25": "10 to 25 km"}; yl, yt = [], []; y = 0
for ft in ["PA", "KBA"]:
    for z in ["core", "r0_5", "r5_10", "r10_25"]:
        for j, p in enumerate(P3):
            s = d[(d["product"] == p) & (d.feature_type == ft) & (d.zone == z)].set_index("year").urban_pct
            yy = y + (j - 1) * 0.25; ax.plot([s[2000], s[2020]], [yy, yy], color=PC[p], lw=1.4)
            ax.plot(s[2000], yy, "o", mfc="white", mec=PC[p], ms=4, mew=1.1); ax.plot(s[2020], yy, "o", color=PC[p], ms=4.5)
        yt.append(y); yl.append(f"{ft} {ZN[z]}"); y += 1
    y += 0.5
ax.set_yticks(yt); ax.set_yticklabels(yl); ax.invert_yaxis(); ax.set_xlim(0, 5); ax.set_xlabel("Urban share of land (%)"); ax.grid(axis="x", color=GRID, lw=0.6)
h = [Line2D([], [], color=PC[p], marker="o", ms=4.5, lw=1.4) for p in P3] + [Line2D([], [], ls="", marker="o", mfc="white", mec="k", ms=4), Line2D([], [], ls="", marker="o", color="k", ms=4.5)]
ax.legend(h, [PN[p] for p in P3] + ["2000", "2020"], loc="lower right", fontsize=7.5)
ax.text(-0.02, 1.02, "c", transform=ax.transAxes, fontsize=10, fontweight="bold", ha="right", va="bottom")
ax = fig.add_subplot(gs[1, 1]); XC = 0.12
ax.axvspan(17.5, 22.5, color="#efefef", lw=0); ax.axvline(5, color="#999999", lw=0.6, ls=":"); ax.axvline(10, color="#999999", lw=0.6, ls=":")
for k, p in enumerate(P3):
    g = prof[prof["product"] == p]; r = g[g.band != "core"].sort_values("d_mid"); core = g[g.band == "core"].urban_pct.iloc[0]
    ax.plot(r.d_mid, r.urban_pct, color=PC[p], marker="o", ms=1.8, lw=0.9, label=PN[p])
    ax.plot(XC, 0.33 / (1.6 ** k), "D", color=PC[p], ms=4); ax.text(XC * 1.25, 0.33 / (1.6 ** k), f"{core:.2f}%", color=PC[p], fontsize=7, va="center")
    pk = r.loc[r.urban_pct.idxmax()]; ax.plot(pk.d_mid, pk.urban_pct, "o", mfc="none", mec=PC[p], ms=7, mew=1.0)
    plat = r[(r.d0 >= 17.5) & (r.d1 <= 22.5)].urban_pct.mean()
    ax.text(21.5, plat * {0: 1.22, 1: 0.74, 2: 0.6}[k], f"{plat:.1f}%", color=PC[p], fontsize=7, ha="center")
    yt_ = {0: 6.7, 1: 4.7, 2: 0.6}[k]; xt = {0: 21.0, 1: 1.05, 2: 0.9}[k]
    ax.annotate(f"Peak, {pk.d_mid:g} km: {pk.urban_pct:.1f}%", xy=(pk.d_mid, pk.urban_pct), xytext=(xt, yt_), fontsize=7, color=PC[p], ha="right" if k == 0 else "left",
                arrowprops=dict(arrowstyle="-", color=PC[p], lw=0.6))
ax.text(0.99, 0.99, "Lowest plateau, 17.5 to 22.5 km", transform=ax.transAxes, ha="right", va="top", fontsize=7)
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(0.09, 27); ax.set_ylim(0.1, 10)
ax.set_xticks([XC, 0.25, 0.5, 1, 2, 5, 10, 25]); ax.set_xticklabels(["Core", "0.25", "0.5", "1", "2", "5", "10", "25"]); ax.minorticks_off()
ax.set_yticks([0.1, 0.2, 0.5, 1, 2, 5, 10]); ax.set_yticklabels(["0.1", "0.2", "0.5", "1", "2", "5", "10"])
ax.set_xlabel("Distance from protected area boundary (km)"); ax.set_ylabel("Urban share of land, 2020 (%)")
h = [Line2D([], [], color=PC[p], marker="o", ms=2, lw=0.9) for p in P3] + [Line2D([], [], ls="", marker="D", color="k", ms=4)]
ax.legend(h, [PN[p] for p in P3] + ["Core"], loc="lower right", fontsize=7.5)
ax.text(-0.02, 1.02, "d", transform=ax.transAxes, fontsize=10, fontweight="bold", ha="right", va="bottom")
save(fig, "Figure2_pa_kba")
