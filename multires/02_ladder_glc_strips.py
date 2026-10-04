# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

from mr_common import *
import time
BUDGET = float(os.environ.get("BGD_BUDGET", 1e12)); t_start = time.time()
YRS = [2000, 2020]; rows = []
PART = f"{OUT}/T_ladder_glc_national_strips.csv"
done = pd.read_csv(PART).to_dict("records") if os.path.exists(PART) else []
have = {(r["year"], r["factor"]) for r in done}
for yr in YRS:
    for f in [3, 10, 17, 33]:
        if (yr, f) in have: continue
        if time.time() - t_start > BUDGET - 60: print("CONTINUE"); sys.exit(0)
        t0 = time.time()
        with rasterio.open(f"{HARM}/GLC30_{yr}.tif") as r, rasterio.open(f"{HARM}/mask_glc30.tif") as m:
            T = r.transform; crs = r.crs; H, W = r.height, r.width
            Hc, Wc = -(-H // f), -(-W // f); Tc = block_transform(T, f); pxc = row_areas_km2(Tc, Hc, Wc)
            arg = np.zeros((Hc, Wc), np.uint8); frac = np.zeros((Hc, Wc), np.float32); landfrac = np.zeros((Hc, Wc), np.float32)
            step = f * max(1, 2048 // f)
            for r0 in range(0, H, step):
                n = min(step, H - r0); win = rasterio.windows.Window(0, r0, W, n)
                a = r.read(1, window=win); a[m.read(1, window=win) == 0] = 0
                g, u, l = aggregate_glc(a, f); c0 = r0 // f
                arg[c0:c0 + g.shape[0]] = g; frac[c0:c0 + g.shape[0]] = np.where(l > 0, u / np.maximum(l, 1), 0); landfrac[c0:c0 + g.shape[0]] = l / float(f * f)
        A_frac = (frac * landfrac * pxc[:, None]).sum(); A_maj = ((arg == 1) * landfrac * pxc[:, None]).sum(); A_land = (landfrac * pxc[:, None]).sum()
        prof = dict(driver="GTiff", height=Hc, width=Wc, count=1, crs=crs, transform=Tc, compress="deflate")
        with rasterio.open(f"{OUT}/GLC_f{f}_{yr}_class.tif", "w", dtype="uint8", nodata=0, **prof) as d: d.write(arg, 1)
        with rasterio.open(f"{OUT}/GLC_f{f}_{yr}_urbfrac.tif", "w", dtype="float32", **prof) as d: d.write(frac, 1)
        if not os.path.exists(f"{OUT}/GLC_f{f}_landfrac.tif"):
            with rasterio.open(f"{OUT}/GLC_f{f}_landfrac.tif", "w", dtype="float32", **prof) as d: d.write(landfrac, 1)
        done.append(dict(year=yr, factor=f, cell=FACTORS[f], land_km2=A_land, urban_frac_km2=A_frac, urban_maj_km2=A_maj)); pd.DataFrame(done).to_csv(PART, index=False)
        print(yr, FACTORS[f], f"frac {A_frac:.1f} maj {A_maj:.1f} land {A_land:.1f}  {time.time()-t0:.0f}s", flush=True)
print("DONE")
