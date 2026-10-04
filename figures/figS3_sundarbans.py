# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

from figcommon import *
from PIL import Image
from rasterio.features import rasterize
base = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "figure_inputs", "FigureS3_panels_ab.png"); b = sp3()
z = np.load(os.path.join(SDM_WORK, "urban_sdm_grid.npz")); Lg, CA = z["L"], z["CA"]
A = (np.load(os.path.join(SDM_WORK, "nval.npy")) == 1098) & (np.load(os.path.join(SDM_WORK, "irr.npy")) > -1e30) & (Lg / CA > 0.5)
ref = rasterio.open(os.path.join(SDM_RASTERS, "irreplacecable_areas_up2.tif"))
sund = rasterize([(pa[pa.NAME_ENG.str.contains("Sundarban", case=False)].geometry.union_all(), 1)], out_shape=ref.shape, transform=ref.transform, fill=0, dtype="uint8").astype(bool)
share = []
for s in b.species:
    with rasterio.open(os.path.join(SDM_RASTERS, f"{s}.tif")) as r: h = (r.read(1) == 1) & A
    share.append(Lg[h & sund].sum() / Lg[h].sum())
b["share_sundarbans"] = share; o = b[b.share_sundarbans > 0.5].copy()
o = o.sort_values("GLC30_2020").reset_index(drop=True)
STC = {"CR": PURPLE, "EN": PURPLE, "VU": PURPLE, "NT": PINK, "LC": PEACH, "DD": GRAY}
im = Image.open(base); Wpx, Hpx = im.size; CUT = int(round(Hpx * 0.535))
fig = plt.figure(figsize=(Wpx / 300, (Hpx - CUT) / 300), dpi=300)
ax = fig.add_axes([0.235, 0.115, 0.74, 0.84])
for i, r in o.iterrows():
    c = STC[r.nationalStatus]
    ax.plot([r.GLC30_2000, r.GLC30_2020], [i, i], color=c, lw=1.4, zorder=2)
    ax.plot(r.GHSL100_2020, i, "s", color=c, ms=4.5, zorder=3); ax.plot(r.ESA300_2020, i, "^", color=c, ms=5, zorder=3)
    ax.plot(r.GLC30_2000, i, "o", mfc="white", mec=c, mew=1.3, ms=5, zorder=4); ax.plot(r.GLC30_2020, i, "o", color=c, ms=5.5, zorder=4)
    xm = max(r.GLC30_2020, r.GLC30_2000, r.GHSL100_2020, r.ESA300_2020)
    ax.text(xm + 0.03, i, f"{r.nationalStatus}, {100 * r.share_sundarbans:.0f}% in Sundarbans", va="center", fontsize=7)
ax.set_yticks(range(len(o))); ax.set_yticklabels([" ".join(s.split("_")[:2]) for s in o.species], fontstyle="italic", fontsize=7.5)
ax.set_xlim(0, 1.3); ax.set_ylim(-0.6, len(o) - 0.4); ax.set_xlabel("Urban share of modeled habitat (%)")
h = [Line2D([], [], ls="", marker="o", mfc="white", mec="k", mew=1.3, ms=5), Line2D([], [], ls="", marker="o", color="k", ms=5.5),
     Line2D([], [], ls="", marker="s", color="k", ms=4.5), Line2D([], [], ls="", marker="^", color="k", ms=5)] + [Patch(color=c) for c in [PURPLE, PINK, PEACH, GRAY]]
ax.legend(h, ["GLC_FCS30D 30 m, 2000", "GLC_FCS30D 30 m, 2020", "GHSL 100 m, 2020", "ESA CCI 300 m, 2020", "Threatened", "Near Threatened", "Least Concern", "Data Deficient"],
          loc="lower right", fontsize=7, labelspacing=0.25, handlelength=1.4)
ax.text(-0.02, 1.0, "c", transform=ax.transAxes, fontsize=10, fontweight="bold", ha="right", va="bottom")
tmp = os.path.join(os.path.expanduser("~"), "_s3c_tmp.png"); fig.savefig(tmp, dpi=300, facecolor="white"); plt.close(fig)
c_ = Image.open(tmp).convert("RGB").resize((Wpx, Hpx - CUT)); out = im.convert("RGB").copy(); out.paste(c_, (0, CUT))
out.save(os.path.join(FIG_SUPP, "FigureS3_sundarbans.png")); os.remove(tmp)
o.to_csv(os.path.join(SDM_DATA, "figS3_sundarbans_species.csv"), index=False); print(o[["species", "share_sundarbans", "GLC30_2020"]].round(3).to_string())
