# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

import os, time, json
import numpy as np

import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from config import (ROOT, OUT, HARM, TAB, LOG, SP, SPS, ESA_DIR, GLC_EARLY, GLC_LATE, BOUNDARY_SHP, PA_SHP, KBA_SHP, ADM1_SHP,
                    IUCN_FILES, GHSL_DIR, GAUD_DIR, ADM1_OUT)

CLASS_NAMES = {1: "Urban", 2: "Cropland", 3: "Mosaic cropland-natural", 4: "Forest",
               5: "Flooded forest", 6: "Shrubland", 7: "Grassland", 8: "Wetland",
               9: "Bare/sparse", 10: "Water"}
NCLASS = 10

ESA_XWALK = {
    10: 2, 11: 2, 12: 2, 20: 2,
    30: 3, 40: 3,
    50: 4, 60: 4, 61: 4, 62: 4, 70: 4, 71: 4, 72: 4, 80: 4, 81: 4, 82: 4, 90: 4, 100: 4,
    110: 7, 130: 7,
    120: 6, 121: 6, 122: 6,
    140: 9, 150: 9, 151: 9, 152: 9, 153: 9, 200: 9, 201: 9, 202: 9, 220: 9,
    160: 5, 170: 5,
    180: 8,
    190: 1,
    210: 10,
}
GLC_XWALK = {
    10: 2, 11: 2, 12: 2, 20: 2,
    51: 4, 52: 4, 61: 4, 62: 4, 71: 4, 72: 4, 81: 4, 82: 4, 91: 4, 92: 4,
    120: 6, 121: 6, 122: 6,
    130: 7,
    140: 9, 150: 9, 152: 9, 153: 9, 200: 9, 201: 9, 202: 9, 220: 9,
    181: 5, 185: 5,
    180: 8, 182: 8, 183: 8, 184: 8, 186: 8, 187: 8,
    190: 1,
    210: 10,
}

def make_lut(xwalk, size=256):
    lut = np.zeros(size, dtype=np.uint8)
    for k, v in xwalk.items():
        lut[k] = v
    return lut

ESA_YEARS = [1992, 1996, 2000, 2004, 2008, 2012, 2016, 2020]
GLC_YEARS = [1985, 1990, 1995, 2000, 2004, 2005, 2008, 2010, 2012, 2015, 2016, 2020]
COMMON_YEARS = [2000, 2004, 2008, 2012, 2016, 2020]

ESA_FILES = {
    1992: "ESACCI-LC-L4-LCCS-Map-300m-P1Y-1992-v2.0.7cds.area-subset.27.93.20.5.88.nc",
    1996: "ESACCI-LC-L4-LCCS-Map-300m-P1Y-1996-v2.0.7cds.area-subset.27.93.20.5.88.nc",
    2000: "ESACCI-LC-L4-LCCS-Map-300m-P1Y-2000-v2.0.7cds.area-subset.27.93.20.5.88.nc",
    2004: "ESACCI-LC-L4-LCCS-Map-300m-P1Y-2004-v2.0.7cds.area-subset.27.93.20.5.88.nc",
    2008: "ESACCI-LC-L4-LCCS-Map-300m-P1Y-2008-v2.0.7cds.area-subset.27.93.20.5.88.nc",
    2012: "ESACCI-LC-L4-LCCS-Map-300m-P1Y-2012-v2.0.7cds.area-subset.27.93.20.5.88.nc",
    2016: "C3S-LC-L4-LCCS-Map-300m-P1Y-2016-v2.1.1.area-subset.27.93.20.5.88.nc",
    2020: "C3S-LC-L4-LCCS-Map-300m-P1Y-2020-v2.1.1.area-subset.27.93.20.5.88.nc",
}

R_AUTHALIC = 6371007.181

def row_areas_km2(transform, height, width):
    dlon = abs(transform.a)
    dlat = abs(transform.e)
    top = transform.f
    lat_top = top - np.arange(height) * dlat
    lat_bot = lat_top - dlat
    a = (R_AUTHALIC ** 2) * np.deg2rad(dlon) * (np.sin(np.deg2rad(lat_top)) - np.sin(np.deg2rad(lat_bot)))
    return a / 1e6

class Budget:
    def __init__(self, seconds):
        self.t0 = time.time(); self.s = seconds
    def left(self):
        return self.s - (time.time() - self.t0)
    def ok(self, need=0):
        return self.left() > need

def log(msg, name="pipeline.log"):
    line = time.strftime("%H:%M:%S ") + msg
    print(line, flush=True)
    with open(os.path.join(LOG, name), "a") as f:
        f.write(line + "\n")
