# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

from mr_common import *
from scipy.stats import mannwhitneyu, kruskal
import time
Z = zones_gdf()
sp = gpd.read_file(f"{P26}/Species/species_ranges_bgd.geojson")
tab = pd.read_csv(f"{P26}/Species/species_table.csv")[["id_no", "sci_name", "status_group", "taxon", "range_bgd_share"]]
zrows, srows = [], []

def run(label, yr, urban_bool, land_bool, T, shape, crs, cell_area, frac=None, do_species=True):
    t0 = time.time()
    z = zonal_share(Z, T, shape, crs, urban_bool, land_bool, cell_area); z["product"] = label; z["year"] = yr; zrows.append(z)
    msg = f"  {label} {yr}: PA core {z[(z.feature_type=='PA')&(z.zone=='core')].urban_pct.iloc[0]:.3f}  PA 5-10 {z[(z.feature_type=='PA')&(z.zone=='r5_10')].urban_pct.iloc[0]:.2f}"
    if do_species:
        s = species_share(sp, T, shape, crs, urban_bool, land_bool, cell_area, frac); s["product"] = label; s["year"] = yr; srows.append(s)
        msg += f"  species median {s.urban_maj_pct.median():.2f}"
        pd.concat(srows).to_csv(f"{OUT}/T_species_ladder_partial.csv", index=False)
    pd.concat(zrows).to_csv(f"{OUT}/T_zones_ladder_partial.csv", index=False)
    print(msg + f"  {time.time()-t0:.0f}s", flush=True)

for yr in [2000, 2020]:
    with rasterio.open(f"{HARM}/GLC30_{yr}.tif") as r: arr = r.read(1); T = r.transform; crs = r.crs
    with rasterio.open(f"{HARM}/mask_glc30.tif") as r: msk = r.read(1) > 0
    arr[~msk] = 0; px = row_areas_km2(T, *arr.shape).astype(np.float32)
    run("GLC 30 m", yr, arr == 1, arr > 0, T, arr.shape, crs, px, do_species=False)
    del arr, msk
    for f in [3, 10, 17, 33]:
        with rasterio.open(f"{OUT}/GLC_f{f}_{yr}_class.tif") as r: cls = r.read(1); T = r.transform; crs = r.crs
        with rasterio.open(f"{OUT}/GLC_f{f}_{yr}_urbfrac.tif") as r: fr = r.read(1)
        with rasterio.open(f"{OUT}/GLC_f{f}_landfrac.tif") as r: lf = r.read(1)
        px = row_areas_km2(T, *cls.shape).astype(np.float32)
        cell_area = (px[:, None] * lf).astype(np.float32)
        run(f"GLC {FACTORS[f]} majority", yr, cls == 1, lf > 0, T, cls.shape, crs, cell_area, frac=fr)
    with rasterio.open(f"{OUT}/GHSL_f1_{yr}_builtfrac.tif") as r: fr = r.read(1); T = r.transform; crs = r.crs
    with rasterio.open(f"{OUT}/GHSL_land100.tif") as r: land = r.read(1).astype(bool)
    ca = np.full(fr.shape, 0.01, dtype=np.float32)
    run("GHSL 100 m built>=0.2", yr, fr >= 0.2, land, T, fr.shape, crs, ca, frac=fr)
    with rasterio.open(f"{OUT}/GHSL_f3_{yr}_majcat20.tif") as r: mj = r.read(1).astype(bool); T3 = r.transform
    with rasterio.open(f"{OUT}/GHSL_f3_landfrac.tif") as r: lf3 = r.read(1)
    with rasterio.open(f"{OUT}/GHSL_f3_{yr}_meanfrac.tif") as r: mf3 = r.read(1)
    run("GHSL 300 m majority", yr, mj, lf3 > 0, T3, mj.shape, crs, (0.09 * lf3).astype(np.float32), frac=mf3)

sx = pd.read_csv(f"{P26}/Species/species_urban_exposure.csv")
sx = sx[(sx.dataset == "GLC30") & (sx.year.isin([2000, 2020]))]
srows.append(pd.DataFrame(dict(id_no=sx.id_no, land_km2=sx.range_land_km2, urban_maj_km2=np.nan, urban_frac_km2=np.nan, urban_maj_pct=sx.urban_pct, urban_frac_pct=sx.urban_pct, product="GLC 30 m", year=sx.year)))
zz = pd.concat(zrows); zz.to_csv(f"{OUT}/T_zones_ladder.csv", index=False)
ss = pd.concat(srows).merge(tab, on="id_no"); ss.to_csv(f"{OUT}/T_species_ladder.csv", index=False)

res = []
for prod, g in ss.groupby("product"):
    w = g.pivot_table(index=["id_no", "status_group"], columns="year", values="urban_maj_pct").reset_index()
    w = w[w[2020].notna()]; w["gain"] = w[2020] - w[2000]
    th, lc = w[w.status_group == "Threatened"], w[w.status_group == "Least Concern"]
    u1 = mannwhitneyu(th[2020], lc[2020]); u2 = mannwhitneyu(th.gain, lc.gain)
    k = kruskal(*[x[2020].values for _, x in w[w.status_group.isin(["Threatened", "Near Threatened", "Least Concern"])].groupby("status_group")])
    res.append(dict(product=prod, n_species=len(w), median_all_2000=w[2000].median(), median_all_2020=w[2020].median(),
                    thr_median_2020=th[2020].median(), lc_median_2020=lc[2020].median(), MWU_2020_p=u1.pvalue,
                    thr_median_gain=th.gain.median(), lc_median_gain=lc.gain.median(), MWU_gain_p=u2.pvalue, KW_2020_p=k.pvalue,
                    n_thr_above1pct_2020=(th[2020] > 1).sum(), n_all_above1pct_2020=(w[2020] > 1).sum(),
                    n_thr_zero_2020=(th[2020] == 0).sum()))
R = pd.DataFrame(res); R.to_csv(f"{OUT}/T_species_tests_ladder.csv", index=False)
print(R.round(3).to_string())
print(zz.pivot_table(index=["feature_type", "zone"], columns=["product", "year"], values="urban_pct").round(2).to_string())
