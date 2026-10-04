# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

from figcommon import *
sp = sp3(); th = sp[sp.thr == 1]
fig, axs = plt.subplots(1, 3, figsize=(W, 3.1), gridspec_kw=dict(wspace=0.38, left=0.08, right=0.95, top=0.92, bottom=0.17))
CM = {"GLC30": "Blues", "GHSL100": "Greens", "ESA300": "Oranges"}
for k, p in enumerate(P3):
    ax = axs[k]; ys = sorted(int(c[-4:]) for c in sp.columns if c.startswith(f"traj_{p}_"))
    ax.axhline(288, color="k", ls=":", lw=1.0); ax.text(2000, 293, "All 288 threatened species", fontsize=7, va="bottom")
    ends = []
    for j, t in enumerate([0.5, 1, 2, 3]):
        n = [int((th[f"traj_{p}_{y}"] > t).sum()) for y in ys]; c = plt.get_cmap(CM[p])(0.3 + 0.23 * j)
        ax.plot(ys, n, marker=PM[p], color=c, ms=3.5, lw=1.2); ends.append((n[-1], f"> {t:g}%: {n[-1]}", c, t))
    ends.sort(key=lambda e: (e[0], -e[3])); pos = []
    for v, _, _, _ in ends: pos.append(max(v, pos[-1] + 13) if pos else v)
    for (v, txt, c, _), py in zip(ends, pos): ax.text(2021, py, txt, color=c, fontsize=7, va="center")
    ax.set_ylim(-40, 320); ax.set_xlim(1998, 2027); ax.set_xticks([2000, 2005, 2010, 2015, 2020]); ax.tick_params(axis="x", labelsize=7)
    ax.set_xlabel("Year"); ax.set_ylabel("Number of threatened species"); ax.text(0.03, 0.03, {"GLC30": "GLC_FCS30D, 30 m", "GHSL100": "GHSL, 100 m", "ESA300": "ESA CCI, 300 m"}[p], transform=ax.transAxes, fontsize=7.5)
    ax.text(-0.02, 1.02, "abc"[k], transform=ax.transAxes, fontsize=10, fontweight="bold", ha="right", va="bottom")
save(fig, "FigureS4_threshold_counts", supp=True)
t5 = th.sort_values("GLC30_2020", ascending=False).reset_index(drop=True); FL = 0.01; per = int(np.ceil(len(t5) / 3))
fig, axs = plt.subplots(1, 3, figsize=(W, 10.2), gridspec_kw=dict(wspace=1.25, left=0.2, right=0.985, top=0.995, bottom=0.07))
for k in range(3):
    ax = axs[k]; s = t5.iloc[k * per:(k + 1) * per].reset_index(drop=True)
    for i, r in s.iterrows():
        v = [max(r[f"{p}_2020"], FL) for p in P3]; ax.plot([min(v), max(v)], [i, i], color="#d9d9d9", lw=0.6, zorder=1)
        for p, x in zip(P3, v): ax.plot(x, i, PM[p], color=PC[p], ms=3.2 if p != "ESA300" else 3.0, zorder=2, lw=0)
    ax.axvline(1, color="k", ls="--", lw=0.9); ax.set_xscale("log"); ax.set_xlim(0.007, 15); ax.minorticks_off()
    ax.set_xticks([FL, 0.1, 1, 10]); ax.set_xticklabels(["None", "0.1", "1", "10"], fontsize=7)
    ax.set_yticks(range(len(s))); ax.set_yticklabels([f"{' '.join(n.split('_')[:2])} ({c})" for n, c in zip(s.species, s.nationalStatus)], fontsize=6.3, fontstyle="italic")
    ax.set_ylim(len(s) - 0.4, -0.6); ax.tick_params(axis="y", length=0, pad=1)
    if k == 1: ax.set_xlabel("Urban share of modeled habitat, 2020 (%)")
h = [Line2D([], [], ls="", marker=PM[p], color=PC[p], ms=5) for p in P3] + [Line2D([], [], color="k", ls="--", lw=0.9)]
fig.legend(h, ["GLC_FCS30D, 30 m", "GHSL, 100 m", "ESA CCI, 300 m", "1% threshold"], loc="lower center", ncol=4, fontsize=7.5, bbox_to_anchor=(0.55, 0.0))
save(fig, "FigureS5_threatened_species", supp=True)
raw = pd.read_excel(os.path.join(SUPP, "Data_S3_species_IUCN_ranges.xlsx"), sheet_name="Data S3", header=None).iloc[2:]
iu = pd.DataFrame({"name": raw[1], "cat": raw[3], "range": raw[4].astype(float), "a": raw[6].astype(float), "b": raw[7].astype(float)})
GRP = lambda c: {"CR": "Threatened", "EN": "Threatened", "VU": "Threatened", "NT": "Near Threatened", "DD": "Data Deficient"}.get(c, "Least Concern")
fig, axs = plt.subplots(1, 2, figsize=(W, 3.6), gridspec_kw=dict(wspace=1.05, left=0.17, right=0.97, top=0.95, bottom=0.27))
rng_ = lambda v: f"{v:,.1f}" if v < 1 else f"{v:,.0f}"
for k, (d, nmc, ac, bc, catc, area, xl) in enumerate([(iu, "name", "a", "b", "cat", "range", "Urban share of range, 2020 (%)\nIUCN range polygons, 420 vertebrates"),
                                                    (sp.assign(name=sp.species.map(lambda s: " ".join(s.split("_")[:2]))), "name", "GLC30_2000", "GLC30_2020", "nationalStatus", "habitat_km2", "Urban share of habitat, 2020 (%)\nHabitat models, 1,097 species")]):
    ax = axs[k]; top = d.sort_values(bc, ascending=False).head(12).reset_index(drop=True)
    for i, r in top.iterrows():
        ax.barh(i, r[bc], color=SCOL[GRP(r[catc])], height=0.62); ax.plot([r[ac]] * 2, [i - 0.31, i + 0.31], color="k", lw=1.6)
        ax.text(r[bc] + 1.5, i, f"{r[bc]:.1f}%  ·  {rng_(r[area])} km²  ·  {r[catc]}", va="center", fontsize=7.5)
    ax.set_yticks(range(12)); ax.set_yticklabels(top[nmc], fontstyle="italic", fontsize=7.5); ax.set_ylim(11.6, -0.6); ax.set_xlim(0, 120); ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel(xl); ax.text(-0.02, 1.02, "ab"[k], transform=ax.transAxes, fontsize=10, fontweight="bold", ha="right", va="bottom")
h = [Patch(color=SCOL[g]) for g in ["Threatened", "Near Threatened", "Least Concern", "Data Deficient"]] + [Line2D([], [], color="k", lw=1.6)]
fig.legend(h, ["Threatened", "Near Threatened", "Least Concern", "Data Deficient", "Value in 2000"], loc="lower center", ncol=5, fontsize=7.5, bbox_to_anchor=(0.55, 0.0))
save(fig, "FigureS6_most_exposed", supp=True)
