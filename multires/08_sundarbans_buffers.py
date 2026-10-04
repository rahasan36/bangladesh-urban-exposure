# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

import rasterio, numpy as np, pandas as pd, geopandas as gpd
from rasterio.mask import mask
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); from config import OUT as UP, MR_HARM as MR, MR_TAB
zones = gpd.read_file(f"{UP}/SpeciesSpecial/sundarbans_zones.geojson")
rows = []
with rasterio.open(f"{UP}/Harmonized/LandArea300_ESA300.tif") as r: land_prof = r
for yr in [2000, 2004, 2008, 2012, 2016, 2020]:
    for _, z in zones.iterrows():
        with rasterio.open(f"{UP}/Harmonized/LandArea300_ESA300.tif") as r: la, _ = mask(r, [z.geometry], crop=True, nodata=0, all_touched=False)
        with rasterio.open(f"{UP}/Harmonized/UrbanFrac300_ESA300_{yr}.tif") as r: uf, _ = mask(r, [z.geometry], crop=True, nodata=0, all_touched=False)
        la = np.where(la[0] > 0, la[0], 0); uf = np.clip(uf[0], 0, 1)
        L = la.sum(); U = (la * uf).sum(); rows.append(dict(dataset="ESA300", zone=z.zone, year=yr, land_km2=L, urban_km2=U, urban_pct=100 * U / L))
with rasterio.open(f"{MR}/GHSL100_2020.tif") as r: crs = r.crs
zm = zones.to_crs(crs)
for yr in [2000, 2005, 2010, 2015, 2020]:
    for _, z in zm.iterrows():
        with rasterio.open(f"{MR}/GHSL100_{yr}.tif") as r: a, _ = mask(r, [z.geometry], crop=True, nodata=0)
        a = a[0]; L = (a > 0).sum() * 0.01; U = (a == 1).sum() * 0.01
        rows.append(dict(dataset="GHSL100", zone=z.zone, year=yr, land_km2=L, urban_km2=U, urban_pct=100 * U / L))
df = pd.DataFrame(rows); df.to_csv(f"{MR_TAB}/sundarbans_buffers_ESA300_GHSL100.csv", index=False)
print(df.pivot_table(index=["dataset", "zone"], columns="year", values="urban_pct").round(3))
