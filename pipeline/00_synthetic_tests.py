# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

import numpy as np, rasterio, geopandas as gpd
from rasterio.transform import from_origin
from rasterio.features import rasterize
from shapely.geometry import box
from common import *
ok = True
t = from_origin(88, 27, 1/360, 1/360); ra = row_areas_km2(t, 2160, 1800)
tot = ra.sum() * 1800
R = R_AUTHALIC / 1e3; exact = R**2 * np.deg2rad(5) * (np.sin(np.deg2rad(27)) - np.sin(np.deg2rad(21)))
print(f"area check: grid {tot:.1f} km2 vs analytic {exact:.1f} km2 -> {'OK' if abs(tot-exact)/exact < 1e-9 else 'FAIL'}"); ok &= abs(tot-exact)/exact < 1e-9
lut = make_lut(GLC_XWALK); codes = np.array([190, 10, 185, 182, 210, 51, 120, 130, 200, 0, 250], dtype=np.uint8)
exp = np.array([1, 2, 5, 8, 10, 4, 6, 7, 9, 0, 0]); print("GLC LUT:", "OK" if (lut[codes] == exp).all() else "FAIL"); ok &= (lut[codes] == exp).all()
lut = make_lut(ESA_XWALK); codes = np.array([190, 30, 40, 170, 180, 210, 110, 100, 0], dtype=np.uint8)
exp = np.array([1, 3, 3, 5, 8, 10, 7, 4, 0]); print("ESA LUT:", "OK" if (lut[codes] == exp).all() else "FAIL"); ok &= (lut[codes] == exp).all()
K = 11; a = np.array([[1, 2, 2], [4, 4, 10], [2, 2, 1]]); b = np.array([[1, 1, 2], [4, 1, 10], [2, 1, 1]])
w = np.ones(9); m = np.bincount((a * K + b).ravel(), weights=w, minlength=K * K).reshape(K, K)
print("transition:", "OK" if m[2, 1] == 2 and m[2, 2] == 2 and m[4, 1] == 1 and m[1, 1] == 2 and m[10, 10] == 1 else "FAIL"); ok &= m[2, 1] == 2 and m[4, 1] == 1
fine = np.zeros((30, 30), dtype=np.uint8) + 2; fine[:10, :10][:6, :] = 1
fine[10:20, 10:20][:3, :] = 1
EH = EW = 3; rowmap = np.repeat(np.arange(3), 10); colmap = np.repeat(np.arange(3), 10)
counts = np.zeros((K, 9), dtype=np.int64); key = (rowmap[:, None] * EW + colmap[None, :]).ravel(); af = fine.ravel()
for c in range(1, K):
    sel = af == c
    if sel.any(): counts[c] += np.bincount(key[sel], minlength=9)
frac = counts[1] / counts[1:].sum(0); maj = np.argmax(counts[1:], axis=0) + 1
print("aggregation:", "OK" if abs(frac[0] - .6) < 1e-9 and abs(frac[4] - .3) < 1e-9 and maj[0] == 1 and maj[4] == 2 and abs(frac.sum() - .9) < 1e-9 else "FAIL")
ok &= abs(frac[0] - .6) < 1e-9 and maj[4] == 2
t = from_origin(0, 10, 1, 1); poly = box(2, 3, 6, 7); z = rasterize([(poly, 1)], out_shape=(10, 10), transform=t, dtype="uint8")
img = np.zeros((10, 10), dtype=np.uint8) + 2; img[3:7, 2:4] = 1
acc = np.bincount((z.astype(int) * K + img).ravel(), minlength=2 * K)[K:]
print("zonal:", "OK" if acc[1] == 8 and acc[2] == 8 and z.sum() == 16 else "FAIL"); ok &= acc[1] == 8 and z.sum() == 16
print("ALL OK" if ok else "SOME FAILED")
