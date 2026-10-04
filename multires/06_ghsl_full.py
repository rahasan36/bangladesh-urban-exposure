# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

from mr_common import *
from rasterio.warp import reproject, Resampling
from scipy import stats
from scipy.ndimage import uniform_filter
import statsmodels.formula.api as smf, time, warnings; warnings.filterwarnings("ignore")
T0 = 0.2; DS = "GHSL100"; YEARS = [2000, 2005, 2010, 2015, 2020]
TABD = f"{OUT}/Tables"; SPD = f"{OUT}/Species"; HD = f"{OUT}/Harmonized"
for d in (TABD, SPD, HD): os.makedirs(d, exist_ok=True)
with rasterio.open(f"{OUT}/GHSL_land100.tif") as r: land = r.read(1).astype(bool); T = r.transform; crs = r.crs
cat = {}
for yr in YEARS:
    with rasterio.open(f"{OUT}/GHSL_f1_{yr}_builtfrac.tif") as r: fr = r.read(1)
    cat[yr] = (fr >= T0) & land
    prof = dict(driver="GTiff", height=land.shape[0], width=land.shape[1], count=1, crs=crs, transform=T, compress="deflate")
    with rasterio.open(f"{HD}/{DS}_{yr}.tif", "w", dtype="uint8", nodata=0, **prof) as d: d.write(np.where(land, np.where(cat[yr], 1, 2), 0).astype(np.uint8), 1)
CELL = 0.01
if not os.path.exists(f"{TABD}/class_areas_{DS}.csv"): pd.DataFrame([dict(dataset=DS, year=yr, class_id=1, class_name="Urban", area_km2=cat[yr].sum() * CELL) for yr in YEARS] +
             [dict(dataset=DS, year=yr, class_id=2, class_name="Other land", area_km2=(land & ~cat[yr]).sum() * CELL) for yr in YEARS]).to_csv(f"{TABD}/class_areas_{DS}.csv", index=False)
print("national urban km2:", {yr: round(cat[yr].sum() * CELL) for yr in YEARS}, flush=True)

SKIP_ZONES = os.path.exists(f"{TABD}/pa_kba_zonal_{DS}.csv")
pa = gpd.read_file(f"{P26}/BGD_PA_Jan2026.geojson").to_crs(UTM).rename(columns={"SITE_ID": "feature_id", "NAME_ENG": "name"}); pa["feature_type"] = "PA"
kb = gpd.read_file(f"{P26}/BGD_KBA_2021.geojson").to_crs(UTM).rename(columns={"SitRecID": "feature_id", "IntName": "name"}); kb["feature_type"] = "KBA"
feats = gpd.GeoDataFrame(pd.concat([pa[["feature_type", "feature_id", "name", "geometry"]], kb[["feature_type", "feature_id", "name", "geometry"]]], ignore_index=True), crs=UTM)
feats["geometry"] = feats.geometry.make_valid(); feats["feature_id"] = feats.feature_id.astype(str)
cores = {t: feats[feats.feature_type == t].geometry.union_all() for t in ["PA", "KBA"]}
zone_list = []
for _, row in feats.iterrows():
    g = row.geometry; allc = cores[row.feature_type]
    zone_list.append((f"{row.feature_type}:{row.feature_id}", row.feature_type, row.feature_id, row["name"], "core", g))
    for n, d0, d1 in BUFFERS[1:]:
        zone_list.append((f"{row.feature_type}:{row.feature_id}", row.feature_type, row.feature_id, row["name"], n, g.buffer(d1).difference(g.buffer(d0) if d0 > 0 else g).difference(allc)))
for t in ["PA", "KBA"]:
    g = cores[t]; zone_list.append((f"{t}:ALL", t, "ALL", f"All {t}s (dissolved)", "core", g))
    for n, d0, d1 in BUFFERS[1:]: zone_list.append((f"{t}:ALL", t, "ALL", f"All {t}s (dissolved)", n, g.buffer(d1).difference(g.buffer(d0) if d0 > 0 else g)))
zg = None if SKIP_ZONES else gpd.GeoDataFrame([dict(feature_key=a, feature_type=b, feature_id=c, name=d, zone=e, geometry=f) for a, b, c, d, e, f in zone_list], crs=UTM).to_crs(crs)
t0 = time.time(); rows = []
for _, z in (zg.iterrows() if not SKIP_ZONES else []):
    if z.geometry.is_empty: continue
    m = rasterize([(z.geometry, 1)], out_shape=land.shape, transform=T, fill=0, dtype="uint8").astype(bool) & land
    if not m.any():
        m = rasterize([(z.geometry, 1)], out_shape=land.shape, transform=T, fill=0, dtype="uint8", all_touched=True).astype(bool) & land
    n_land = m.sum()
    for yr in YEARS:
        nu = (m & cat[yr]).sum()
        rows.append(dict(feature_key=z.feature_key, feature_type=z.feature_type, feature_id=z.feature_id, name=z["name"], dataset=DS, zone=z.zone, year=yr, class_id=1, class_name="Urban", area_km2=nu * CELL))
        rows.append(dict(feature_key=z.feature_key, feature_type=z.feature_type, feature_id=z.feature_id, name=z["name"], dataset=DS, zone=z.zone, year=yr, class_id=2, class_name="Other land", area_km2=(n_land - nu) * CELL))
zon = pd.read_csv(f"{TABD}/pa_kba_zonal_{DS}.csv") if SKIP_ZONES else pd.DataFrame(rows)
if not SKIP_ZONES: zon.to_csv(f"{TABD}/pa_kba_zonal_{DS}.csv", index=False)
zon["feature_id"] = zon.feature_id.astype(str)
e = zon[(zon.feature_id == "ALL")].pivot_table(index=["feature_type", "zone", "year"], columns="class_id", values="area_km2").reset_index()
e["urban_pct"] = 100 * e[1] / (e[1] + e[2]); print(e.pivot_table(index=["feature_type", "zone"], columns="year", values="urban_pct").round(3).to_string(), f"\n  zones {time.time()-t0:.0f}s", flush=True)

if os.path.exists(f"{SPD}/mixed_effects_results_{DS}.csv"): print("species outputs present, skipped", flush=True)
else:
    sp = gpd.read_file(f"{P26}/Species/species_ranges_bgd.geojson"); tab = pd.read_csv(f"{P26}/Species/species_table.csv")
    ca = np.full(land.shape, CELL, dtype=np.float32); srows = []
    for yr in YEARS:
        t0 = time.time()
        s = species_share(sp, T, land.shape, crs, cat[yr], land, ca); s["dataset"] = DS; s["year"] = yr; srows.append(s)
        print(f"  species {yr}: median {s.urban_maj_pct.median():.2f}  {time.time()-t0:.0f}s", flush=True)
    sx = pd.concat(srows).rename(columns={"land_km2": "range_land_km2", "urban_maj_km2": "urban_km2", "urban_maj_pct": "urban_pct"}).drop(columns=["urban_frac_km2", "urban_frac_pct"])
    sx = sx.merge(tab, on="id_no"); sx.to_csv(f"{SPD}/species_urban_exposure_{DS}.csv", index=False)
    sl = []
    for sid, g in sx.groupby("id_no"):
        g = g.sort_values("year"); r = tab[tab.id_no == sid].iloc[0]
        if g.range_land_km2.iloc[0] <= 0: continue
        res = stats.linregress(g.year, g.urban_pct)
        sl.append(dict(dataset=DS, id_no=sid, sci_name=r.sci_name, taxon=r.taxon, category=r.category, status_group=r.status_group, range_bgd_km2=r.range_bgd_km2, range_bgd_share=r.range_bgd_share,
                       urban_pct_2000=float(g[g.year == 2000].urban_pct.iloc[0]), urban_pct_2020=float(g[g.year == 2020].urban_pct.iloc[0]), slope_pct_per_yr=res.slope, change_pp=float(g[g.year == 2020].urban_pct.iloc[0] - g[g.year == 2000].urban_pct.iloc[0])))
    sl = pd.DataFrame(sl); sl.to_csv(f"{SPD}/species_slopes_ranking_{DS}.csv", index=False)
    thr = [dict(dataset=DS, year=yr, threshold=t, n_threatened=int(((sx.year == yr) & (sx.status_group == "Threatened") & (sx.urban_pct > t)).sum()), n_all=int(((sx.year == yr) & (sx.urban_pct > t)).sum())) for yr in YEARS for t in [0.5, 1, 2, 3, 5]]
    pd.DataFrame(thr).to_csv(f"{SPD}/threshold_crossings_{DS}.csv", index=False)
    th, lc = sl[sl.status_group == "Threatened"], sl[sl.status_group == "Least Concern"]
    tests = []
    for what, a, b in [("2020", th.urban_pct_2020, lc.urban_pct_2020), ("change", th.change_pp, lc.change_pp), ("slope", th.slope_pct_per_yr, lc.slope_pct_per_yr)]:
        mw = stats.mannwhitneyu(a, b); tests.append(dict(dataset=DS, test=f"Mann-Whitney {what} Threatened vs LC", U=mw.statistic, p=mw.pvalue, median_threatened=a.median(), median_LC=b.median()))
    g3 = sl[sl.status_group.isin(["Threatened", "Near Threatened", "Least Concern"])]
    for what, col in [("2020", "urban_pct_2020"), ("slope", "slope_pct_per_yr")]:
        k = stats.kruskal(*[x[col].values for _, x in g3.groupby("status_group")]); tests.append(dict(dataset=DS, test=f"Kruskal-Wallis {what} by status", U=k.statistic, p=k.pvalue))
    pd.DataFrame(tests).to_csv(f"{SPD}/status_group_tests_{DS}.csv", index=False); print(pd.DataFrame(tests).round(4).to_string(), flush=True)
    d = sx[(sx.status_group != "Data Deficient") & (sx.range_land_km2 > 0)].copy(); d["year_c"] = d.year - 2000
    d["status_group"] = pd.Categorical(d.status_group, ["Least Concern", "Near Threatened", "Threatened"]); d["taxon"] = pd.Categorical(d.taxon, ["Mammal", "Reptile", "Amphibian"])
    out = []
    for resp, col in [("urban_pct", "urban_pct"), ("log_urban_pct", None)]:
        if col is None: d["log_urban_pct"] = np.log(d.urban_pct + 0.01)
        m = smf.mixedlm(f"{resp} ~ status_group * year_c + taxon * year_c", d, groups=d["id_no"]).fit(reml=True)
        for term in m.params.index: out.append(dict(dataset=DS, response=resp, term=term, estimate=m.params[term], se=m.bse[term], p=m.pvalues[term], n_obs=m.nobs, n_species=d.id_no.nunique()))
    mx = pd.DataFrame(out); mx.to_csv(f"{SPD}/mixed_effects_results_{DS}.csv", index=False); print(mx[mx.term.str.contains("year_c")].round(4).to_string(), flush=True)

with rasterio.open(f"{HARM}/LandArea300_GLC30.tif") as r: land300 = r.read(1).astype(float); et = r.transform; EH, EW = r.height, r.width; ecrs = r.crs
for yr in YEARS:
    dst = np.zeros((EH, EW), dtype=np.float32)
    reproject(cat[yr].astype(np.float32), dst, src_transform=T, src_crs=crs, dst_transform=et, dst_crs=ecrs, resampling=Resampling.average)
    prof = dict(driver="GTiff", height=EH, width=EW, count=1, crs=ecrs, transform=et, compress="deflate")
    with rasterio.open(f"{HD}/UrbanFrac300_{DS}_{yr}.tif", "w", dtype="float32", **prof) as dd: dd.write(dst, 1)
F = 3; H, W = EH // F, EW // F; kt = Affine(et.a * F, 0, et.c, 0, et.e * F, et.f)
def agg(a): return a[:H * F, :W * F].reshape(H, F, W, F).sum(axis=(1, 3))
def gistar(x, mask):
    x = np.where(mask, x, 0.0); n = mask.sum(); xbar = x[mask].mean(); s = x[mask].std()
    wsum = uniform_filter(x, 3, mode="constant") * 9; wn = uniform_filter(mask.astype(float), 3, mode="constant") * 9
    denom = s * np.sqrt((n * wn - wn ** 2) / (n - 1)); z = np.where(denom > 0, (wsum - xbar * wn) / np.where(denom > 0, denom, 1), 0)
    return np.where(mask, z, np.nan)
with rasterio.open(f"{HARM}/Hotspots_1km_GLC30.tif") as r: rich = r.read(3); glc_cls = r.read(6); glc_mask = glc_cls != -9999
rich = np.where(glc_mask, rich, 0.0)
landk = agg(land300); u = {}
for yr in [2000, 2020]:
    with rasterio.open(f"{HD}/UrbanFrac300_{DS}_{yr}.tif") as r: u[yr] = agg(r.read(1).astype(float) * land300)
mask = landk > 0.1; gain = u[2020] - u[2000]; z_gain = gistar(gain, mask); expo = gain * rich; z_expo = gistar(expo, mask)
hs_u = (z_gain > 2.58) & mask; hs_e = (z_expo > 2.58) & mask
cls = np.zeros((H, W), dtype=np.uint8); cls[hs_u] = 1; cls[hs_e] = 2; cls[hs_u & hs_e] = 3
prof = dict(driver="GTiff", crs=ecrs, transform=kt, height=H, width=W, compress="LZW")
with rasterio.open(f"{HD}/Hotspots_1km_{DS}.tif", "w", count=6, dtype="float32", nodata=-9999, **prof) as dd:
    o = np.stack([gain, z_gain, rich, expo, z_expo, cls.astype(float)]).astype(np.float32); o[:, ~mask] = -9999; dd.write(o)
adm = gpd.read_file(f"{P26}/ADM1/geoBoundaries-BGD-ADM1.geojson").to_crs(4326)
adm_r = rasterize([(g, i + 1) for i, g in enumerate(adm.geometry)], out_shape=(H, W), transform=kt, fill=0, dtype="uint8")
rows = []
for i, name in enumerate(list(adm.shapeName) + ["Bangladesh"]):
    sel = ((adm_r == i + 1) if i < len(adm) else mask) & mask
    rows.append(dict(dataset=DS, division=name, land_km2=landk[sel].sum(), urban_gain_km2=gain[sel].sum(), urban_hotspot_km2=landk[sel & hs_u].sum(), exposure_hotspot_km2=landk[sel & hs_e].sum(),
                     coincident_km2=landk[sel & hs_u & hs_e].sum(), urban_gain_in_hotspots_km2=gain[sel & hs_u].sum(), mean_threatened_richness=rich[sel].mean()))
hd = pd.DataFrame(rows); hd.to_csv(f"{TABD}/hotspots_by_division_{DS}.csv", index=False); print(hd.round(1).to_string(index=False), flush=True)
with rasterio.open(f"{HARM}/Hotspots_1km_ESA300.tif") as r: esa_cls = r.read(6)
co = []
for nm, other in [("GLC30", glc_cls), ("ESA300", esa_cls)]:
    a = (other == 1) | (other == 3); b = hs_u
    co.append(dict(pair=f"{DS} vs {nm}", ghsl_hotspot_km2=landk[b].sum(), other_hotspot_km2=landk[a & mask].sum(), overlap_km2=landk[a & b].sum(), jaccard=landk[a & b].sum() / landk[(a | b) & mask].sum()))
pd.DataFrame(co).to_csv(f"{TABD}/hotspot_overlap_{DS}.csv", index=False); print(pd.DataFrame(co).round(3).to_string(index=False))
dec = []
for yr in [2000, 2020]:
    with rasterio.open(f"{HD}/UrbanFrac300_{DS}_{yr}.tif") as r: fg = r.read(1).astype(float)
    with rasterio.open(f"{HARM}/UrbanFrac300_GLC30_{yr}.tif") as r: fl = r.read(1).astype(float)
    with rasterio.open(f"{HARM}/UrbanFrac300_ESA300_{yr}.tif") as r: fe = r.read(1).astype(float)
    dec.append(dict(year=yr, GLC30_frac_km2=(fl * land300).sum(), GLC30_maj300_km2=((fl > 0.5) * land300).sum(), GHSL100_frac_km2=(fg * land300).sum(), GHSL100_maj300_km2=((fg > 0.5) * land300).sum(), ESA300_km2=(fe * land300).sum()))
pd.DataFrame(dec).to_csv(f"{TABD}/decomposition_three_products.csv", index=False); print(pd.DataFrame(dec).round(0).to_string(index=False)); print("DONE")
