# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

from figcommon import *
from scipy.stats import gaussian_kde
sp = sp3(); FLOOR = 0.003
SC = {"Threatened": PURPLE, "Non-threatened": PEACH, "Data Deficient": GRAY}
fig = plt.figure(figsize=(W, 5.7))
gs = GridSpec(2, 3, figure=fig, hspace=0.42, wspace=0.42, left=0.08, right=0.985, top=0.97, bottom=0.1)
rng = np.random.default_rng(7)
for k, p in enumerate(P3):
    ax = fig.add_subplot(gs[0, k])
    for i, st in enumerate(["Threatened", "Non-threatened", "Data Deficient"]):
        v = np.maximum(sp[sp.threatStatus == st][f"{p}_2020"].values, FLOOR); lv = np.log10(v)
        ax.scatter(i - 0.08 + rng.uniform(-0.17, 0.17, len(v)), v, s=4, color=SC[st], alpha=0.45, lw=0, zorder=2)
        q1, md, q3 = np.percentile(lv, [25, 50, 75]); iqr = q3 - q1
        lo = lv[lv >= q1 - 1.5 * iqr].min(); hi = lv[lv <= q3 + 1.5 * iqr].max()
        ax.add_patch(plt.Rectangle((i - 0.2, 10 ** q1), 0.24, 10 ** q3 - 10 ** q1, facecolor=(1, 1, 1, 0.55), edgecolor="k", lw=0.9, zorder=3))
        ax.plot([i - 0.2, i + 0.04], [10 ** md] * 2, color="k", lw=1.1, zorder=4)
        for a_, b_ in [(lo, q1), (q3, hi)]: ax.plot([i - 0.08] * 2, [10 ** a_, 10 ** b_], color="k", lw=0.9, zorder=3)
        for e_ in (lo, hi): ax.plot([i - 0.13, i - 0.03], [10 ** e_] * 2, color="k", lw=0.9, zorder=3)
        grid = np.linspace(lv.min(), lv.max(), 200); dens = gaussian_kde(lv)(grid); dens = dens / dens.max() * 0.3
        ax.fill_betweenx(10 ** grid, i + 0.1, i + 0.1 + dens, color=SC[st], alpha=0.6, lw=0); ax.plot(i + 0.1 + dens, 10 ** grid, color=SC[st], lw=0.8)
    ax.set_yscale("log"); ax.set_ylim(0.002, 40); ax.set_xlim(-0.5, 2.55)
    ax.set_yticks([0.01, 0.1, 1, 10]); ax.set_yticklabels(["0.01", "0.1", "1", "10"]); ax.grid(axis="y", color=GRID, lw=0.6)
    ax.set_xticks([0, 1, 2]); ax.set_xticklabels(["Threatened", "Non-\nthreatened", "Data\nDeficient"], fontsize=7.5)
    ax.text(0.03, 0.97, PN[p], color=PC[p], transform=ax.transAxes, va="top", fontsize=7.5)
    if k == 0: ax.set_ylabel("Urban share of modeled habitat, 2020 (%)")
    ax.text(-0.02, 1.02, "abc"[k], transform=ax.transAxes, fontsize=10, fontweight="bold", ha="right", va="bottom")
ax = fig.add_subplot(gs[1, 0])
s = pd.read_csv(os.path.join(REV, "share2020_models.csv")).set_index("product"); g = pd.read_csv(os.path.join(REV, "family_mixed_model.csv")).set_index("product")
rows = [("Share 2020,\nleast squares", s, "ols", "lo", "hi"), ("Share 2020,\nfamily intercept", s, "mixed", "mlo", "mhi"),
        ("Gain,\nleast squares", g, "ols_thr", "ols_lo", "ols_hi"), ("Gain,\nfamily intercept", g, "mixed_thr", "mixed_lo", "mixed_hi")]
LAB = {"30 m": "GLC30", "100 m": "GHSL100", "300 m": "ESA300"}
for i, (lab_, t_, e, lo, hi) in enumerate(rows):
    for j, pr in enumerate(["30 m", "100 m", "300 m"]):
        y = i + (j - 1) * 0.22; p = LAB[pr]
        ax.plot([t_.loc[pr, lo], t_.loc[pr, hi]], [y, y], color=PC[p], lw=1.2); ax.plot(t_.loc[pr, e], y, "o", color=PC[p], ms=4)
ax.axvline(0, color="k", ls=":", lw=0.9); ax.set_yticks(range(4)); ax.set_yticklabels([r[0] for r in rows], fontsize=7.5); ax.invert_yaxis()
ax.set_ylim(5.3, -0.6); ax.set_xlabel("Threatened minus\nnon-threatened (pp)")
ax.legend([Line2D([], [], color=PC[p], marker="o", ms=4, lw=1.2) for p in P3], [PN[p] for p in P3], loc="lower left", fontsize=7)
ax.text(-0.02, 1.02, "d", transform=ax.transAxes, fontsize=10, fontweight="bold", ha="right", va="bottom")
ax = fig.add_subplot(gs[1, 1]); c = pd.read_csv(os.path.join(SDM_DATA, "sdm_counts_ladder.csv")).set_index("rep"); z = lambda v: np.maximum(v, 0.5)
glc = ["GLC30", "GLCp90", "GLCp300", "GLCp510", "GLCp990"]; gh = ["GHSL100", "GHSLm300", "GHSLm500", "GHSLm1000"]
for col, colr in [("n_all_gt1", "k"), ("n_thr_gt1", PURPLE)]:
    ax.plot(c.loc[glc, "cell"], z(c.loc[glc, col]), "-o", color=colr, ms=4, lw=1.2)
    ax.plot(c.loc[gh, "cell"], z(c.loc[gh, col]), "--s", color=colr, ms=4, lw=1.0, mfc="white")
    ax.plot(300, z(c.loc["ESA300", col]), "^", color=ORANGE, mec=colr, mew=1.0, ms=6, zorder=5)
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(25, 1150); ax.set_ylim(0.4, 2500); ax.minorticks_off()
ax.set_xticks([30, 100, 300, 1000]); ax.set_xticklabels(["30", "100", "300", "1000"]); ax.set_yticks([1, 10, 100, 1000]); ax.set_yticklabels(["1", "10", "100", "1000"])
ax.grid(color=GRID, lw=0.6); ax.set_xlabel("Cell size (m)"); ax.set_ylabel("Species above 1% urban")
h = [Line2D([], [], color="k", lw=1.2), Line2D([], [], color=PURPLE, lw=1.2), Line2D([], [], color="k", marker="o", ms=4, lw=1.2),
     Line2D([], [], color="k", ls="--", marker="s", mfc="white", ms=4, lw=1.0), Line2D([], [], ls="", marker="^", color=ORANGE, mec="k", ms=6)]
ax.legend(h, ["All species", "Threatened", "GLC ladder", "GHSL ladder", "ESA CCI"], loc="lower left", fontsize=7, labelspacing=0.25)
ax.text(-0.02, 1.02, "e", transform=ax.transAxes, fontsize=10, fontweight="bold", ha="right", va="bottom")
ax = fig.add_subplot(gs[1, 2]); acc = pd.read_csv(os.path.join(REV, "domain_accounting.csv")).set_index("zone")
segs = [("analysis domain", "#4d4d4d", "In modeled habitat"), ("removed as built-up by the models", "#a6a6a6", "Removed as built-up"), ("outside model extent", "#dcdcdc", "Outside model extent")]
for i, (p, nm) in enumerate([("GLC30", "30 m"), ("GHSL100", "100 m"), ("ESA300", "300 m")]):
    left = 0
    for z_, colr, _ in segs:
        v = acc.loc[z_, f"{p}_urban2020_km2"]; ax.barh(i, v, left=left, color=colr, height=0.6, edgecolor="white", lw=0.5); left += v
    ax.text(left + 50, i, f"{100 * acc.loc['analysis domain', f'{p}_urban2020_km2'] / acc.loc['national land', f'{p}_urban2020_km2']:.0f}%", va="center", fontsize=7.5)
ax.set_yticks([0, 1, 2]); ax.set_yticklabels(["30 m", "100 m", "300 m"]); ax.invert_yaxis(); ax.set_xlim(0, 4800); ax.set_ylim(3.6, -0.6)
ax.set_xticks([0, 1000, 2000, 3000, 4000]); ax.set_xlabel("Urban land in 2020 (km²)")
ax.legend([Patch(color=c_) for _, c_, _ in segs], [l_ for _, _, l_ in segs], loc="lower right", fontsize=7, labelspacing=0.25)
ax.text(-0.02, 1.02, "f", transform=ax.transAxes, fontsize=10, fontweight="bold", ha="right", va="bottom")
save(fig, "Figure3_species_status")
