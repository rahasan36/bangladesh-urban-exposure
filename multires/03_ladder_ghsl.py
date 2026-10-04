# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

from mr_common import *
from rasterio.merge import merge
import time

bnd = gpd.read_file(f"{P26}/BGD_Boundary.geojson")
PART = f"{OUT}/T_ladder_ghsl_national_partial.csv"
rows = pd.read_csv(PART).to_dict("records") if os.path.exists(PART) else []
done = {r["year"] for r in rows}; first = not os.path.exists(f"{OUT}/GHSL_land100.tif")
if not first:
    with rasterio.open(f"{OUT}/GHSL_land100.tif") as r: land = r.read(1).astype(bool); prof = dict(r.profile); prof.pop("nodata", None); prof.pop("dtype", None)
for yr in YEARS_GHSL:
    if yr in done: continue
    t0 = time.time()
    tiles = [f"{GHSL_DIR}/GHS_BUILT_S_E{yr}_GLOBE_R2023A_54009_100_V1_0_{t}.tif" for t in ["R6_C27", "R7_C27"]]
    if not all(os.path.exists(p) for p in tiles): tiles = [f"{GHSL_DIR}/GHS_BUILT_S_E{yr}_GLOBE_R2023A_54009_100_V1_0.tif"]
    srcs = [rasterio.open(p) for p in tiles]
    crs = srcs[0].crs; b = bnd.to_crs(crs); minx, miny, maxx, maxy = b.total_bounds
    mos, T = merge(srcs, bounds=(minx - 200, miny - 200, maxx + 200, maxy + 200))
    nod = srcs[0].nodata; [s.close() for s in srcs]
    a = mos[0].astype(np.float32); a[a == nod] = 0
    if first:
        land = rasterize([(g, 1) for g in b.geometry], out_shape=a.shape, transform=T, fill=0, dtype="uint8").astype(bool)
        prof = dict(driver="GTiff", height=a.shape[0], width=a.shape[1], count=1, crs=crs, transform=T, compress="deflate")
        with rasterio.open(f"{OUT}/GHSL_land100.tif", "w", dtype="uint8", **prof) as d: d.write(land.astype(np.uint8), 1)
        first = False
    frac = np.clip(a / 10000.0, 0, 1); frac[~land] = 0
    with rasterio.open(f"{OUT}/GHSL_f1_{yr}_builtfrac.tif", "w", dtype="float32", **prof) as d: d.write(frac, 1)
    cell_km2 = 0.01
    A_land = land.sum() * cell_km2; A_built = frac.sum() * cell_km2
    u100 = frac >= 0.5
    rows.append(dict(year=yr, factor=1, cell="100 m", land_km2=A_land, built_km2=A_built, urban_cat_maj_km2=u100.sum() * cell_km2,
                     urban_meanfrac_km2=u100.sum() * cell_km2, urban_cat20_km2=(frac >= 0.2).sum() * cell_km2))
    print(yr, "100 m: land", round(A_land), "built", round(A_built), "urban(>=0.5)", round(u100.sum() * cell_km2), "urban(>=0.2)", round((frac >= 0.2).sum() * cell_km2), f"{time.time()-t0:.0f}s", flush=True)
    for f, nm in [(3, "300 m"), (5, "500 m"), (10, "1000 m")]:
        L = blocksum(land.astype(np.uint8), f); S = blocksum((frac * 10000).astype(np.uint32), f) / 10000.0
        U = blocksum(u100.astype(np.uint8), f)
        landfrac = L / float(f * f); meanfrac = np.where(L > 0, S / np.maximum(L, 1), 0)
        Tc = block_transform(T, f); ca = cell_km2 * f * f
        A_maj_cat = ((U > L / 2.0) & (L > 0)).astype(float) * landfrac * ca
        A_maj_mean = ((meanfrac >= 0.5) & (L > 0)).astype(float) * landfrac * ca
        rows.append(dict(year=yr, factor=f, cell=nm, land_km2=(landfrac * ca).sum(), built_km2=S.sum() * cell_km2,
                         urban_cat_maj_km2=A_maj_cat.sum(), urban_meanfrac_km2=A_maj_mean.sum(), urban_cat20_km2=np.nan))
        profc = dict(prof, height=L.shape[0], width=L.shape[1], transform=Tc)
        with rasterio.open(f"{OUT}/GHSL_f{f}_{yr}_meanfrac.tif", "w", dtype="float32", **profc) as d: d.write(meanfrac.astype(np.float32), 1)
        with rasterio.open(f"{OUT}/GHSL_f{f}_{yr}_majcat.tif", "w", dtype="uint8", **profc) as d: d.write(((U > L / 2.0) & (L > 0)).astype(np.uint8), 1)
        if yr == YEARS_GHSL[0]:
            with rasterio.open(f"{OUT}/GHSL_f{f}_landfrac.tif", "w", dtype="float32", **profc) as d: d.write(landfrac.astype(np.float32), 1)
        print(f"  {nm}: built {S.sum()*cell_km2:.0f}  maj(cat) {A_maj_cat.sum():.0f}  maj(meanfrac) {A_maj_mean.sum():.0f}", flush=True)
    pd.DataFrame(rows).to_csv(PART, index=False)
pd.DataFrame(rows).to_csv(f"{OUT}/T_ladder_ghsl_national.csv", index=False)
print(pd.DataFrame(rows).round(1).to_string())
