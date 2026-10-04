# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

from figcommon import *
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "sdm")); from gistar import gistar
z = np.load(os.path.join(SDM_WORK, "urban_sdm_grid.npz")); L, CA = z["L"], z["CA"]
irr = np.load(os.path.join(SDM_WORK, "irr.npy")); nval = np.load(os.path.join(SDM_WORK, "nval.npy")); rt = np.load(os.path.join(SDM_WORK, "rich_t.npy")).astype(float)
A = (nval == 1098) & (irr > -1e30) & (L / CA > 0.5)
RH = (gistar(np.where(A, rt, 0.0), A) > 2.58) & A; IH = A & (irr >= 1)
with rasterio.open(os.path.join(SDM_WORK, "hotspots_sdm_GLC30.tif")) as r: UH = (r.read(1) == 1) & A; ext = [r.bounds.left, r.bounds.right, r.bounds.bottom, r.bounds.top]
fig = plt.figure(figsize=(W, 6.9))
gs = GridSpec(2, 2, figure=fig, height_ratios=[1.55, 1], hspace=0.12, wspace=0.18, left=0.07, right=0.985, top=0.975, bottom=0.08)
for k, (S, nm, col) in enumerate([(RH, "Threatened richness hotspot", "#a98fd0"), (IH, "Irreplaceability ≥ 1", "#5ab4c8")]):
    ax = fig.add_subplot(gs[0, k]); cls = np.where(A, 1, 0); cls[S] = 2; cls[UH] = 3; cls[S & UH] = 4
    ax.imshow(np.ma.masked_equal(cls, 0), extent=ext, cmap=ListedColormap(["#e6e6e6", col, "#c0392b", "#ffe000"]), vmin=0.5, vmax=4.5, interpolation="nearest")
    adm.boundary.plot(ax=ax, color="k", linewidth=0.6); bnd.boundary.plot(ax=ax, color="k", linewidth=0.6); mapax(ax); ax.set_ylim(20.55, 26.7); ax.set_anchor("C")
    ax.legend([Patch(color=col), Patch(color="#c0392b"), Patch(color="#ffe000")], [nm, "Urban expansion hotspot, 30 m", "Both"], loc="lower left", fontsize=7.5, handlelength=1.4, labelspacing=0.3)
    ax.text(0.0, 1.0, "ab"[k], transform=ax.transAxes, fontsize=10, fontweight="bold", ha="left", va="top")
hs = pd.read_csv(os.path.join(SDM_DATA, "hotspot_sensitivity_sdm.csv")); hs = hs[(hs.neighborhood == "3x3") & (hs["rule"] == "z>2.58")]
ax = fig.add_subplot(gs[1, 0])
for i, surf in enumerate(["threatened richness", "irreplaceability>=1"]):
    for j, p in enumerate(P3):
        v = hs[(hs.surface == surf) & (hs["product"] == p)].ratio.iloc[0]; x = i + (j - 1) * 0.26
        ax.bar(x, v, width=0.25, color=PC[p], label=PN[p] if i == 0 else None); ax.text(x, v + 0.04, f"{v:.2f}", ha="center", fontsize=6.5)
ax.axhline(1, color="k", ls="--", lw=0.9); ax.set_xticks([0, 1]); ax.set_xticklabels(["Threatened richness", "Irreplaceability ≥ 1"]); ax.set_ylim(0, 3.0)
ax.set_ylabel("Observed ÷ expected overlap"); ax.legend(loc="upper left", fontsize=7.5)
ax.text(-0.02, 1.02, "c", transform=ax.transAxes, fontsize=10, fontweight="bold", ha="right", va="bottom")
ic = pd.read_csv(os.path.join(SDM_DATA, "irreplaceability_class_gain.csv"))
ax = fig.add_subplot(gs[1, 1])
for i, r in ic.iterrows():
    for j, p in enumerate(P3):
        v = r[f"{p}_ratio"]; x = i + (j - 1) * 0.26; ax.bar(x, v, width=0.25, color=PC[p]); ax.text(x, v + 0.06, f"{v:.1f}", ha="center", fontsize=6.5)
ax.axhline(1, color="k", ls="--", lw=0.9); ax.set_xticks(range(len(ic))); ax.set_xticklabels(ic.irreplaceability_class); ax.set_xlabel("Irreplaceability score")
ax.set_ylabel("Share of urban gain ÷ share of land"); ax.set_ylim(0, max(4.6, ic[[f"{p}_ratio" for p in P3]].values.max() * 1.1))
ax.text(-0.02, 1.02, "d", transform=ax.transAxes, fontsize=10, fontweight="bold", ha="right", va="bottom")
save(fig, "Figure5_hotspots")
