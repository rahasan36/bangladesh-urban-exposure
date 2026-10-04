# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

from mr_common import *
from scipy.ndimage import uniform_filter
from scipy import stats
from statsmodels.stats.multitest import multipletests
OUTD = f"{OUT}/Review"; os.makedirs(OUTD, exist_ok=True)

if os.path.exists(f"{OUTD}/pa_designation_year.csv"): print("designation-year table present, skipped", flush=True)
else:
    pa = gpd.read_file(f"{P26}/BGD_PA_Jan2026.geojson").to_crs(UTM); pa["geometry"] = pa.geometry.make_valid()
    groups = {"designated by 2000 (n=%d)" % (pa.STATUS_YR <= 2000).sum(): pa[pa.STATUS_YR <= 2000], "designated 2001-2010 (n=%d)" % ((pa.STATUS_YR > 2000) & (pa.STATUS_YR <= 2010)).sum(): pa[(pa.STATUS_YR > 2000) & (pa.STATUS_YR <= 2010)],
              "designated after 2010 (n=%d)" % (pa.STATUS_YR > 2010).sum(): pa[pa.STATUS_YR > 2010], "all PAs (n=%d)" % len(pa): pa}
    allcores = pa.geometry.union_all()
    zl = []
    for gname, g in groups.items():
        u = g.geometry.union_all()
        zl.append(dict(group=gname, zone="core", geometry=u))
        for n, d0, d1 in BUFFERS[1:]: zl.append(dict(group=gname, zone=n, geometry=u.buffer(d1).difference(u.buffer(d0) if d0 > 0 else u).difference(allcores)))
    ZG = gpd.GeoDataFrame(zl, crs=UTM)
    rows = []
    for yr in [2000, 2020]:
        with rasterio.open(f"{HARM}/GLC30_{yr}.tif") as r: arr = r.read(1); T = r.transform; crs = r.crs
        with rasterio.open(f"{HARM}/mask_glc30.tif") as r: msk = r.read(1) > 0
        arr[~msk] = 0; px = row_areas_km2(T, *arr.shape).astype(np.float32)
        z = zonal_share(ZG, T, arr.shape, crs, arr == 1, arr > 0, px); z["dataset"] = "GLC30"; z["year"] = yr; rows.append(z); del arr, msk
        with rasterio.open(f"{OUT}/GHSL_f1_{yr}_builtfrac.tif") as r: fr = r.read(1); Tg = r.transform; crsg = r.crs
        with rasterio.open(f"{OUT}/GHSL_land100.tif") as r: land = r.read(1).astype(bool)
        z = zonal_share(ZG, Tg, fr.shape, crsg, fr >= 0.2, land, np.full(fr.shape, 0.01, np.float32)); z["dataset"] = "GHSL100"; z["year"] = yr; rows.append(z)
    DES = pd.concat(rows); DES.to_csv(f"{OUTD}/pa_designation_year.csv", index=False)
    print(DES.pivot_table(index=["group", "zone"], columns=["dataset", "year"], values="urban_pct").round(3).to_string(), flush=True)

with rasterio.open(f"{HARM}/LandArea300_GLC30.tif") as r: land300 = r.read(1).astype(float)
dec = pd.read_csv(f"{P26}/Tables/decomposition_by_class.csv"); du = dec[dec.class_id == 1].set_index("year")
rows = []
for yr in [2000, 2004, 2008, 2012, 2016, 2020]:
    with rasterio.open(f"{HARM}/UrbanFrac300_GLC30_{yr}.tif") as r: fl = r.read(1).astype(float)
    with rasterio.open(f"{HARM}/UrbanFrac300_ESA300_{yr}.tif") as r: fe = r.read(1).astype(float)
    A_frac = (fl * land300).sum(); A_esa = (fe * land300).sum()
    d = dict(year=yr, fraction_km2=A_frac, plurality_km2=du.loc[yr, "glc300_majority_km2"], esa_km2=A_esa)
    for t in [0.10, 0.25, 0.50]:
        A = ((fl >= t) * land300).sum(); d[f"thr{int(t*100)}_km2"] = A; d[f"thr{int(t*100)}_resshare_pct"] = 100 * (A_frac - A) / (A_frac - A_esa)
    d["plurality_resshare_pct"] = 100 * (A_frac - d["plurality_km2"]) / (A_frac - A_esa)
    rows.append(d)
THR = pd.DataFrame(rows); THR.to_csv(f"{OUTD}/threshold_sensitivity_300m.csv", index=False); print(THR.round(1).to_string(), flush=True)

with rasterio.open(f"{HARM}/mask_esa300.tif") as r: et = r.transform; EH, EW = r.height, r.width; ecrs = r.crs
F = 3; H, W_ = EH // F, EW // F; kt = Affine(et.a * F, 0, et.c, 0, et.e * F, et.f)
def agg(a): return a[:H * F, :W_ * F].reshape(H, F, W_, F).sum(axis=(1, 3))
def gistar(x, mask, k=3):
    x = np.where(mask, x, 0.0); n = mask.sum(); xbar = x[mask].mean(); s = x[mask].std()
    wsum = uniform_filter(x, k, mode="constant") * k * k; wn = uniform_filter(mask.astype(float), k, mode="constant") * k * k
    denom = s * np.sqrt((n * wn - wn ** 2) / (n - 1)); z = np.where(denom > 0, (wsum - xbar * wn) / np.where(denom > 0, denom, 1), 0)
    return np.where(mask, z, np.nan)
def fdr_hot(z, mask, q=0.01):
    p = 1 - stats.norm.cdf(np.nan_to_num(z[mask], nan=-10)); rej = multipletests(p, alpha=q, method="fdr_bh")[0]
    out = np.zeros_like(mask); out[mask] = rej; return out
sp = gpd.read_file(f"{P26}/Species/species_ranges_bgd.geojson"); thr = sp[sp.status_group == "Threatened"]
rich = np.zeros((H, W_)); wrich = np.zeros((H, W_))
for _, r in thr.iterrows():
    m = rasterize([(r.geometry, 1)], out_shape=(H, W_), transform=kt, fill=0, dtype="uint8").astype(float); rich += m; wrich += m * r.range_bgd_share
landk = agg(land300); mask = landk > 0.1
gains = {}
for ds, hp in [("GLC30", f"{HARM}/UrbanFrac300_GLC30_{{yr}}.tif"), ("GHSL100", f"{OUT}/Harmonized/UrbanFrac300_GHSL100_{{yr}}.tif")]:
    u = {}
    for yr in [2000, 2020]:
        with rasterio.open(hp.format(yr=yr)) as r: u[yr] = agg(r.read(1).astype(float) * land300)
    gains[ds] = u[2020] - u[2000]
rows = []
for k in [3, 5]:
    z_rich = gistar(rich, mask, k); z_wr = gistar(wrich, mask, k)
    for ds, gain in gains.items():
        z_g = gistar(gain, mask, k); expo = gain * rich; z_e = gistar(expo, mask, k)
        for rule in ["z>2.58", "FDR q=0.01"]:
            if rule == "z>2.58": hg, hr, hw, he = (z_g > 2.58) & mask, (z_rich > 2.58) & mask, (z_wr > 2.58) & mask, (z_e > 2.58) & mask
            else: hg, hr, hw, he = fdr_hot(z_g, mask), fdr_hot(z_rich, mask), fdr_hot(z_wr, mask), fdr_hot(z_e, mask)
            A = lambda b: landk[b].sum()
            rows.append(dict(dataset=ds, neighbourhood=f"{k}x{k}", rule=rule, urban_hotspot_km2=A(hg), richness_hotspot_km2=A(hr), coincident_urban_x_richness_km2=A(hg & hr),
                             pct_of_urban_hotspot_in_richness_hotspot=100 * A(hg & hr) / A(hg), expected_pct_if_independent=100 * A(hr) / A(mask),
                             weighted_richness_hotspot_km2=A(hw), coincident_urban_x_wrichness_km2=A(hg & hw), pct_urban_in_wrichness=100 * A(hg & hw) / A(hg), expected_pct_wrichness=100 * A(hw) / A(mask),
                             exposure_product_hotspot_km2=A(he), coincident_product_construction_km2=A(hg & he), pct_product_construction=100 * A(hg & he) / A(he)))
        if k == 3 and ds == "GLC30":
            prof = dict(driver="GTiff", crs=ecrs, transform=kt, height=H, width=W_, compress="LZW")
            with rasterio.open(f"{OUTD}/Hotspots_independent_1km.tif", "w", count=5, dtype="float32", nodata=-9999, **prof) as dd:
                o = np.stack([rich, wrich, z_rich, z_wr, (((z_g > 2.58) & mask).astype(float) + 2 * ((z_rich > 2.58) & mask).astype(float))]).astype(np.float32); o[:, ~mask] = -9999; dd.write(o)
                for i, n in enumerate(["threatened_richness", "rangeshare_weighted_richness", "gistar_z_richness", "gistar_z_wrichness", "class_1urbanGLC_2richness_3both"]): dd.set_band_description(i + 1, n)
IND = pd.DataFrame(rows); IND.to_csv(f"{OUTD}/hotspots_independent_surface.csv", index=False); print(IND.round(1).to_string(), flush=True)

sx = pd.read_csv(f"{P26}/Species/species_urban_exposure.csv"); sx = sx[(sx.dataset == "GLC30") & (sx.year.isin([2000, 2020]))]
small = sp[(sp.range_bgd_km2 < 200) & (sp.range_bgd_km2 > 0)]
rows = []
for yr in [2000, 2020]:
    with rasterio.open(f"{HARM}/GLC30_{yr}.tif") as r: arr = r.read(1); T = r.transform; crs = r.crs
    with rasterio.open(f"{HARM}/mask_glc30.tif") as r: msk = r.read(1) > 0
    arr[~msk] = 0; px = row_areas_km2(T, *arr.shape).astype(np.float32)
    s = species_share(small, T, arr.shape, crs, arr == 1, arr > 0, px); s["year"] = yr; rows.append(s); del arr, msk
d30 = pd.concat(rows).rename(columns={"urban_maj_pct": "urban_pct_30m_direct"})[["id_no", "year", "land_km2", "urban_pct_30m_direct"]]
cmp_ = d30.merge(sx[["id_no", "sci_name", "status_group", "range_bgd_km2", "year", "urban_pct"]].rename(columns={"urban_pct": "urban_pct_300m_fraction"}), on=["id_no", "year"])
cmp_["diff_pp"] = cmp_.urban_pct_30m_direct - cmp_.urban_pct_300m_fraction; cmp_.to_csv(f"{OUTD}/narrow_range_concordance.csv", index=False)
for yr in [2000, 2020]:
    c = cmp_[cmp_.year == yr].dropna()
    print(f"narrow ranges (<200 km2) {yr}: n={len(c)}  Spearman rho={stats.spearmanr(c.urban_pct_30m_direct, c.urban_pct_300m_fraction).correlation:.3f}  median |diff|={c.diff_pp.abs().median():.2f} pp  "
          f"90th pct |diff|={c.diff_pp.abs().quantile(0.9):.2f} pp  agree on >5% threshold: {((c.urban_pct_30m_direct>5)==(c.urban_pct_300m_fraction>5)).mean()*100:.0f}%  agree on >1%: {((c.urban_pct_30m_direct>1)==(c.urban_pct_300m_fraction>1)).mean()*100:.0f}%")
print("DONE")
