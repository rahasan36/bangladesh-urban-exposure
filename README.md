# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats

Rakibul Ahasan, Palash Chandra Das, Shawan Chowdhury

This repository holds the code, derived data, figures, and tables for a study that maps the urban land of Bangladesh at 30, 100, and 300 m (GLC_FCS30D, a GHSL-derived map, and ESA CCI) and measures how map resolution changes the urban exposure of protected areas, Key Biodiversity Areas, and the modeled habitat of 1,097 species from the national conservation prioritization. Aggregating the finer maps to one label per 300 m cell reproduces most of the difference in urban area among the three products, and that loss carries into protected area exposure, threshold counts of exposed species, and the overlap of urban expansion hotspots with irreplaceable land.

## Folder structure

```
config.py            paths and environment variables
run_all.py           python run_all.py <stage>: pipeline, multires, sdm, analysis, tables, figures
requirements.txt     Python packages (environment.yml for conda)
pipeline/            download and processing of the land cover, protected area, KBA, and IUCN inputs
multires/            resolution ladders, GHSL-derived map, and sensitivity analyses
sdm/                 urban share of modeled habitat, species tests and models, hotspot overlaps
analysis/            group tests, family mixed models, hotspot domain accounting, Data S1 to S3
tables/              Supplementary Tables S1 to S7
figures/             one script per figure
data/                derived data; tables and figures rebuild from it without the raw inputs
  pipeline/          harmonized 300 m rasters, zonal and species tables
  multires/          GHSL-derived rasters, ladder and sensitivity tables
  sdm/               habitat-model species, test, model, and hotspot tables
  analysis/          outputs of analysis/
  tables/            tables behind Figures 1 and 2 and Tables S1 to S3
  supplemental/      Data S1 to S3
  figure_inputs/     panels a and b of Figure S3
outputs/
  figures/           Figures 1 to 5 and S1 to S7
  tables/            Tables S1 to S7
```

## License

MIT License. See `LICENSE`.
