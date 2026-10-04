# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

from mr_common import *
with rasterio.open(f"{OUT}/GHSL_land100.tif") as r: land = r.read(1).astype(bool); T = r.transform; crs = r.crs
rows = []
for yr in [2000, 2010, 2020]:
    with rasterio.open(f"{OUT}/GHSL_f1_{yr}_builtfrac.tif") as r: frac = r.read(1)
    for t in [0.05, 0.1, 0.2, 0.3, 0.5]:
        u = (frac >= t) & land
        rows.append(dict(year=yr, threshold=t, factor=1, cell="100 m", urban_km2=u.sum() * 0.01))
        for f, nm in [(3, "300 m"), (5, "500 m"), (10, "1000 m")]:
            L = blocksum(land.astype(np.uint8), f); U = blocksum(u.astype(np.uint8), f)
            A = (((U > L / 2.0) & (L > 0)) * (L / float(f * f)) * 0.01 * f * f).sum()
            rows.append(dict(year=yr, threshold=t, factor=f, cell=nm, urban_km2=A))
            if t == 0.2:
                profc = dict(driver="GTiff", height=L.shape[0], width=L.shape[1], count=1, crs=crs, transform=block_transform(T, f), compress="deflate")
                with rasterio.open(f"{OUT}/GHSL_f{f}_{yr}_majcat20.tif", "w", dtype="uint8", **profc) as d: d.write(((U > L / 2.0) & (L > 0)).astype(np.uint8), 1)
        if t == 0.2:
            prof = dict(driver="GTiff", height=land.shape[0], width=land.shape[1], count=1, crs=crs, transform=T, compress="deflate")
            with rasterio.open(f"{OUT}/GHSL_f1_{yr}_cat20.tif", "w", dtype="uint8", **prof) as d: d.write(u.astype(np.uint8), 1)
d = pd.DataFrame(rows); d.to_csv(f"{OUT}/T_ladder_ghsl_thresholds.csv", index=False)
print(d.pivot_table(index=["year", "threshold"], columns="cell", values="urban_km2").round(0)[["100 m", "300 m", "500 m", "1000 m"]].to_string())
