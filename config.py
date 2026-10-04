# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

import os

PKG = os.path.dirname(os.path.abspath(__file__))

ROOT = os.environ.get("BGD_DATA_ROOT", os.path.join(PKG, "data", "raw"))
ESA_DIR = os.environ.get("BGD_ESA_DIR", os.path.join(ROOT, "ESACCI"))
GLC_TILE_DIR = os.path.join(ROOT, "GLC_FCS30D_tiles")
GLC_EARLY = os.path.join(ROOT, "Data", "Processed", "GLC_Early_BGD_30m.tif")
GLC_LATE = os.path.join(ROOT, "Data", "Processed", "GLC_Late_BGD_30m.tif")
BOUNDARY_SHP = os.path.join(ROOT, "BGD_Boundary.shp")
PA_SHP = os.path.join(ROOT, "BGD_PA.shp")
KBA_SHP = os.path.join(ROOT, "KBAsGlobal_2021_March_01-20210718T193202Z-001", "KBAsGlobal_2021_March_01", "KBAsGlobal_2021_March_01_POL.shp")
ADM1_SHP = os.path.join(ROOT, "geoBoundaries-BGD-ADM1-all", "geoBoundaries-BGD-ADM1.shp")
IUCN_FILES = {"Mammal": ["MAMMALS_TERRESTRIAL_ONLY.shp"], "Reptile": ["REPTILES_PART1.shp", "REPTILES_PART2.shp"],
              "Amphibian": ["AMPHIBIANS_PART1.shp", "AMPHIBIANS_PART2.shp"]}
GHSL_DIR = os.environ.get("BGD_GHSL_DIR", os.path.join(ROOT, "GHSL"))
GAUD_DIR = os.environ.get("BGD_GAUD_DIR", os.path.join(ROOT, "GAUD"))
SDM_RASTERS = os.environ.get("BGD_SDM_RASTERS", os.path.join(ROOT, "BinaryRasters_Built-upAreasRemoved"))

OUT = os.environ.get("BGD_OUT", os.path.join(PKG, "data", "pipeline"))
HARM = os.path.join(OUT, "Harmonized")
TAB = os.path.join(OUT, "Tables")
SP = os.path.join(OUT, "Species")
SPS = os.path.join(OUT, "SpeciesSpecial")
LOG = os.path.join(OUT, "logs")
ADM1_OUT = os.path.join(OUT, "ADM1")

MR = os.environ.get("BGD_MR", os.path.join(PKG, "data", "multires"))
MR_HARM = os.path.join(MR, "Harmonized")
MR_TAB = os.path.join(MR, "Tables")
MR_SP = os.path.join(MR, "Species")
MR_REVIEW = os.path.join(MR, "Sensitivity")

TAB_MS = os.environ.get("BGD_TAB_MS", os.path.join(PKG, "data", "tables"))
SDM_WORK = os.environ.get("BGD_SDM_WORK", os.path.join(PKG, "data", "sdm_work"))
SDM_DATA = os.environ.get("BGD_SDM_DATA", os.path.join(PKG, "data", "sdm"))
REV = os.environ.get("BGD_REV", os.path.join(PKG, "data", "analysis"))
SUPP = os.environ.get("BGD_SUPP", os.path.join(PKG, "data", "supplemental"))

OUT_FIG = os.environ.get("BGD_OUT_FIG", os.path.join(PKG, "outputs", "figures"))
OUT_TABLES = os.environ.get("BGD_OUT_TABLES", os.path.join(PKG, "outputs", "tables"))

for d in (OUT, HARM, TAB, SP, SPS, LOG, ADM1_OUT, MR, MR_HARM, MR_TAB, MR_SP, MR_REVIEW, TAB_MS, SDM_WORK, SDM_DATA, REV, SUPP,
          OUT_FIG, os.path.join(OUT_FIG, "main"), os.path.join(OUT_FIG, "supplementary"), OUT_TABLES):
    os.makedirs(d, exist_ok=True)
