# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

from figcommon import *
sp = sp3(); mo = pd.read_csv(os.path.join(SDM_DATA, "sdm_models.csv"))
GO = ["Birds", "Butterflies", "Freshwater fishes", "Reptiles", "Mammals", "Crustaceans", "Amphibians"]
GN = {"Freshwater fishes": "Fishes"}
GC = {"Birds": BLUE, "Butterflies": YELLOW, "Freshwater fishes": TEAL, "Reptiles": GREEN, "Mammals": "#b8213a", "Crustaceans": PURPLE, "Amphibians": ORANGE}
CATS = ["CR", "EN", "VU", "NT", "LC", "DD"]
fig = plt.figure(figsize=(W, 6.3))
gs = GridSpec(2, 3, figure=fig, hspace=0.45, wspace=0.42, left=0.08, right=0.985, top=0.97, bottom=0.08)
ax = fig.add_subplot(gs[0, 0]); fl = lambda v: np.maximum(v, 0.001)
nt = sp[sp.thr == 0]; ax.scatter(nt.habitat_km2, fl(nt.GLC30_gain), s=4, color="#cfcfcf", lw=0, label="Non-threatened")
for g_ in GO:
    t_ = sp[(sp.thr == 1) & (sp.group == g_)]; ax.scatter(t_.habitat_km2, fl(t_.GLC30_gain), s=6, color=GC[g_], lw=0, label=GN.get(g_, g_))
rho = mo[(mo.rep == "GLC30") & (mo.model == "Spearman gain vs log area")].est.iloc[0]
ax.text(0.96, 0.04, "\u03c1 = " + f"{rho:.2f}", transform=ax.transAxes, va="bottom", ha="right", fontsize=7.5)
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(40, 150000); ax.set_ylim(0.0007, 60)
ax.set_xticks([100, 1000, 10000, 100000]); ax.set_xticklabels(["100", "1000", "10000", "100000"]); ax.set_yticks([0.001, 0.01, 0.1, 1, 10]); ax.set_yticklabels(["0.001", "0.01", "0.1", "1", "10"])
ax.set_xlabel("Modeled habitat (km²)"); ax.set_ylabel("Gain in urban share, 2000 to 2020 (pp)")
h, l = ax.get_legend_handles_labels(); order = list(range(1, 8)) + [0]
ax.legend([h[i] for i in order], [l[i] for i in order], ncol=2, loc="upper right", fontsize=6.5, columnspacing=0.6, handletextpad=0.1, labelspacing=0.15, markerscale=1.6)
ax = fig.add_subplot(gs[0, 1])
for i, g_ in enumerate(GO):
    for st, dx, c in [(1, -0.15, PURPLE), (0, 0.15, GRAY)]:
        v = sp[(sp.group == g_) & (sp.thr == st) & (sp.threatStatus != "Data Deficient")].GLC30_gain
        if st == 1 and len(v) < 5: continue
        q1, md, q3 = v.quantile([0.25, 0.5, 0.75]); ax.plot([i + dx] * 2, [q1, q3], color=c, lw=1.4); ax.plot(i + dx, md, "o", color=c, ms=5)
ax.set_xticks(range(7)); ax.set_xticklabels([GN.get(g_, g_) for g_ in GO], rotation=35, ha="right", fontsize=7.5); ax.set_ylim(0, 1.6)
ax.set_ylabel("Gain at 30 m (pp)"); ax.legend([Line2D([], [], color=PURPLE, marker="o", ms=5, lw=1.4), Line2D([], [], color=GRAY, marker="o", ms=5, lw=1.4)], ["Threatened", "Non-threatened"], loc="upper right", fontsize=7.5)
ax = fig.add_subplot(gs[0, 2])
for i, term in enumerate(["thr", "log_area"]):
    for j, p in enumerate(P3):
        r = mo[(mo.rep == p) & (mo.model == "OLS gain, HC3") & (mo.term == term)].iloc[0]; y = i + (j - 1) * 0.22
        ax.errorbar(r.est, y, xerr=[[r.est - r.lo], [r.hi - r.est]], fmt="o", color=PC[p], ms=5, lw=1.2, capsize=2.5)
ax.axvline(0, color="k", ls=":", lw=1.0); ax.set_yticks([0, 1]); ax.set_yticklabels(["Threatened", "Habitat\narea, ×10"], fontsize=7.5); ax.set_ylim(2.6, -0.6)
ax.grid(axis="x", color=GRID, lw=0.6); ax.set_xlabel("Effect on gain (pp, 95% CI)")
ax.legend([Line2D([], [], color=PC[p], marker="o", ms=5, lw=1.2) for p in P3], [PN[p] for p in P3], loc="lower left", fontsize=7)
for k, (col, ylab, logy) in enumerate([("2020", "Urban share, 2020 (%)", True), ("gain", "Gain, 2000 to 2020 (pp)", False)]):
    ax = fig.add_subplot(gs[1, k])
    for i, cat in enumerate(CATS):
        for j, p in enumerate(P3):
            v = sp[sp.nationalStatus == cat][f"{p}_{col}"]; q1, md, q3 = v.quantile([0.25, 0.5, 0.75]); x = i + (j - 1) * 0.22
            ax.plot([x, x], [q1, q3], color=PC[p], lw=1.4); ax.plot(x, md, "o", color=PC[p], ms=5)
        if k == 0: ax.text(i, 0.023, f"{(sp.nationalStatus == cat).sum()}", ha="center", fontsize=7)
    ax.set_xticks(range(6)); ax.set_xticklabels(CATS); ax.set_xlabel("National Red List category"); ax.set_ylabel(ylab)
    if logy: ax.set_yscale("log"); ax.set_ylim(0.015, 10); ax.set_yticks([0.1, 1, 10]); ax.set_yticklabels(["0.1", "1", "10"])
    else:
        ax.set_ylim(0, 1.3); ax.legend([Line2D([], [], color=PC[p], marker="o", ms=5, lw=1.4) for p in P3], [PN[p] for p in P3], loc="upper left", fontsize=7.5)
ax = fig.add_subplot(gs[1, 2])
for i, cat in enumerate(CATS):
    for j, p in enumerate(P3):
        v = 100 * (sp[sp.nationalStatus == cat][f"{p}_2020"] > 1).mean(); x = i + (j - 1) * 0.27
        ax.bar(x, v, width=0.27, color=PC[p]); ax.text(x, v + 1.5, f"{v:.0f}", ha="center", fontsize=6)
ax.set_xticks(range(6)); ax.set_xticklabels(CATS); ax.set_ylim(0, 105); ax.set_xlabel("National Red List category"); ax.set_ylabel("Above 1% urban (% of category)")
for a_, l_ in zip(fig.axes, "abcdef"): a_.text(-0.02, 1.02, l_, transform=a_.transAxes, fontsize=10, fontweight="bold", ha="right", va="bottom")
save(fig, "Figure4_species_groups")
