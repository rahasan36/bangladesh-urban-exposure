# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

from mr_common import *
import time

rows = []
for yr in YEARS_GLC:
    t0 = time.time()
    with rasterio.open(f"{HARM}/GLC30_{yr}.tif") as r:
        arr = r.read(1); T = r.transform; crs = r.crs
    with rasterio.open(f"{HARM}/mask_glc30.tif") as r: msk = r.read(1) > 0
    arr[~msk] = 0
    px = row_areas_km2(T, *arr.shape)
    land30 = ((arr > 0) * px[:, None]).sum(); urb30 = ((arr == 1) * px[:, None]).sum()
    rows.append(dict(year=yr, factor=1, cell=FACTORS[1], land_km2=land30, urban_frac_km2=urb30, urban_maj_km2=urb30))
    print(yr, "30 m: land", round(land30), "urban", round(urb30), f"{time.time()-t0:.0f}s", flush=True)
    for f in [3, 10, 17, 33]:
        arg, urb, land = aggregate_glc(arr, f)
        Tc = block_transform(T, f); pxc = row_areas_km2(Tc, *arg.shape)
        frac = np.where(land > 0, urb / np.maximum(land, 1), 0).astype(np.float32)
        landfrac = (land / float(f * f)).astype(np.float32)
        A_frac = (frac * landfrac * pxc[:, None]).sum()
        A_maj = ((arg == 1) * landfrac * pxc[:, None]).sum()
        A_land = (landfrac * pxc[:, None]).sum()
        rows.append(dict(year=yr, factor=f, cell=FACTORS[f], land_km2=A_land, urban_frac_km2=A_frac, urban_maj_km2=A_maj))
        prof = dict(driver="GTiff", height=arg.shape[0], width=arg.shape[1], count=1, crs=crs, transform=Tc, compress="deflate")
        with rasterio.open(f"{OUT}/GLC_f{f}_{yr}_class.tif", "w", dtype="uint8", nodata=0, **prof) as d: d.write(arg, 1)
        with rasterio.open(f"{OUT}/GLC_f{f}_{yr}_urbfrac.tif", "w", dtype="float32", **prof) as d: d.write(frac, 1)
        if yr == YEARS_GLC[0]:
            with rasterio.open(f"{OUT}/GLC_f{f}_landfrac.tif", "w", dtype="float32", **prof) as d: d.write(landfrac, 1)
        print(f"  f={f} ({FACTORS[f]}): frac {A_frac:.0f}  maj {A_maj:.0f}  ratio {A_maj/A_frac:.3f}", flush=True)
    del arr, msk
pd.DataFrame(rows).to_csv(f"{OUT}/T_ladder_glc_national.csv", index=False)
print(pd.DataFrame(rows).round(1).to_string())
