# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

import os, geopandas as gpd
from common import *

def write(gdf, path):
    if os.path.exists(path): print("exists", os.path.basename(path)); return
    gdf.to_file(path, driver="GeoJSON"); print("wrote", os.path.basename(path), len(gdf), "features")

adm = gpd.read_file(ADM1_SHP).to_crs(4326)
write(adm, os.path.join(ADM1_OUT, "geoBoundaries-BGD-ADM1.geojson"))
write(adm, os.path.join(OUT, "BGD_Boundary.geojson"))
pa = gpd.read_file(PA_SHP).to_crs(4326)
write(pa, os.path.join(OUT, "BGD_PA_Jan2026.geojson"))
kb = gpd.read_file(KBA_SHP, bbox=(88, 20.5, 92.8, 26.7)); kb = kb[kb.ISO3 == "BGD"].to_crs(4326)
write(kb, os.path.join(OUT, "BGD_KBA_2021.geojson"))
if not os.path.exists(BOUNDARY_SHP):
    adm.dissolve().to_file(BOUNDARY_SHP); print("wrote", BOUNDARY_SHP)
print("DONE")
