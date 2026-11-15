# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats

This repository accompanies the paper above, which maps the urban land of Bangladesh at 30, 100 and 300 m with GLC_FCS30D, a GHSL-derived map and ESA CCI, and measures how map resolution changes the urban exposure of 55 protected areas, 20 Key Biodiversity Areas and the modeled habitat of 1,097 species from the national conservation prioritization. It includes the Python code from the source products to the tables of the paper, the scripts that draw the figures of the paper and its Supplemental Information with the figures they draw, the derived data the table and figure scripts read, and Data S1 to S3.

The source products are not redistributed here and are available from their providers (see Sources).

**Keywords:** urban expansion, spatial resolution, land cover products, protected areas, Key Biodiversity Areas, species habitat, Bangladesh

## Analyses

| Folder | Contents |
|---|---|
| pipeline | From the source products to the 300 m common grid and the protected area, KBA and IUCN range measurements. 00a_download_inputs.py downloads geoBoundaries and GHS-BUILT-S and lists the other inputs with their sources (--check reports which are present). 00b_prepare_vectors.py, 00c_mosaic_glc.py (GLC_FCS30D mosaic) and 00d_esa_subsets.py (ESA CCI and C3S subsets, 20.5 to 27°N and 88 to 93°E, cut from the global files) prepare the inputs. 00_synthetic_tests.py checks the cell areas, the aggregation and the zonal extraction on synthetic grids. 01_harmonize.py (ten common classes and the 300 m urban fractions), 02_areas_transitions.py (class areas and transitions), 03_decomposition.py (resolution and product components), 04_pa_kba_buffers.py (cores and 0 to 5, 5 to 10 and 10 to 25 km buffers), 05a_species_prep.py, 05b_species_zonal.py and 05c_species_models.py (IUCN ranges), 06_hotspots.py (Getis-Ord Gi*), 07_species_special.py (Sundarbans), 08_tile_test.py (check across the 90°E boundary of the GLC_FCS30D tiles) and 09_validation.py (comparison with GHS-BUILT-S and GAUD). |
| multires | Resolution ladders and the GHSL-derived map. 01_ladder_glc.py and 02_ladder_glc_strips.py (GLC_FCS30D at 90 to 990 m), 03_ladder_ghsl.py and 04_ghsl_thresholds.py (GHSL-derived map at 100 to 1,000 m and its built-up thresholds), 05_zones_species.py (protected area zones and species along the ladders), 06_ghsl_full.py (all measurements under the GHSL-derived map), 07_sensitivities.py (designation year, labeling rules at 300 m, an independent hotspot surface and narrow-range species), 08_sundarbans_buffers.py. |
| sdm | Urban land in the habitat models of the national prioritization. 00_build_sdm_grid.py (1 km grid, urban area of each product, richness, irreplaceability and priority areas), 01_species_urban_share.py (urban share of each species' modeled habitat along the aggregation ladders), 02_species_tests_models.py (threat-status tests, models and threshold counts), 03_hotspot_masks.py and 04_hotspot_overlaps.py (urban expansion hotspots and their overlap with threatened richness and irreplaceable land), 05_validation_hotspots.py (validation inside and outside the hotspots). |
| analysis | group_gain_tests.py (tests within taxonomic groups), family_mixed_model.py and share2020_models.py (least squares and family random-intercept models), hotspot_domain.py (analysis-domain accounting and the all-land hotspots), matched_rankings.py (habitat models and IUCN ranges for the species in both), site_reconcile.py (sites in the rank comparison), data_s1.py, data_s2.py and data_s3.py (Data S1 to S3). |
| tables | build_si_tables.py writes Tables S1 to S7 to outputs/tables and the tables behind Figures 1 and 2 to data/tables. |
| figures | The script behind each figure (see Figures). figstyle.py and figcommon.py hold the shared style and helpers. |

run_all.py runs the stages in order. From the repository folder,

    python run_all.py pipeline
    python run_all.py multires
    python run_all.py sdm
    python run_all.py analysis
    python run_all.py tables
    python run_all.py figures

or `python run_all.py all` for every stage, and `python run_all.py --check` lists the source products that are present. One script of a stage runs with `python run_all.py <stage> <script>`, for example `python run_all.py figures fig3`. The pipeline and multires stages read the source products and resume from their last completed step when they are started again. The tables and figures stages run from the derived data in data and the protected area and KBA polygons (see Data). The paths are set in config.py, and the environment variables below point the scripts to other copies.

| Variable | Default | Contents |
|---|---|---|
| BGD_DATA_ROOT | data/raw | Source products (see Sources) |
| BGD_ESA_DIR | BGD_DATA_ROOT/ESACCI | ESA CCI and C3S subsets |
| BGD_ESA_GLOBAL | none | Folders with the global ESA CCI and C3S files, for 00d_esa_subsets.py |
| BGD_GHSL_DIR, BGD_GAUD_DIR | BGD_DATA_ROOT/GHSL, BGD_DATA_ROOT/GAUD | GHS-BUILT-S and GAUD |
| BGD_SDM_RASTERS | BGD_DATA_ROOT/BinaryRasters_Built-upAreasRemoved | Habitat models, irreplaceability surface and priority areas of the national prioritization |
| BGD_SDM_WORK | data/sdm_work | 1 km grids written by the sdm stage |
| BGD_OUT_FIG, BGD_OUT_TABLES | outputs/figures, outputs/tables | Figures and tables |

## Figures

The folder figures holds the script behind each figure of the paper and its Supplemental Information. Each script writes its figure as PNG to outputs/figures/main or outputs/figures/supplementary. Run a script from the repository folder, for example

    python figures/fig3_species_status.py

Figure 5 and Supplementary Figs 3 and 7 also read the 1 km grids in data/sdm_work, which the sdm stage writes from the habitat models.

| Figure | Content | Script (in figures/) | Inputs (in data unless noted) |
|---|---|---|---|
| Figure 1 | Urban land under three products before and after aggregation | fig1_resolution.py | tables/decomposition_300m.csv, class areas in pipeline/Tables and multires/Tables, 300 m urban fractions in pipeline/Harmonized and multires/Harmonized |
| Figure 2 | Urban land within and around protected areas and Key Biodiversity Areas | fig2_data.py, fig2_distance_profile.py, fig2_pa_kba.py | Zonal tables in pipeline/Tables and multires/Tables, tables/pa_kba_core_buffer.csv, tables/pa_distance_profile_2020.csv |
| Figure 3 | Urban land in the modeled habitat of 1,097 species by national threat status | fig3_species_status.py | sdm/sdm_species.csv, sdm/sdm_counts_ladder.csv, analysis |
| Figure 4 | Urban land in modeled habitat by taxonomic group, habitat area and national Red List category | fig4_species_groups.py | sdm/sdm_species.csv, sdm/sdm_models.csv |
| Figure 5 | Urban expansion hotspots against threatened species richness and irreplaceability | fig5_hotspots.py | sdm/hotspot_sensitivity_sdm.csv, sdm/irreplaceability_class_gain.csv, sdm_work |
| Supplementary Fig. 1 | Urban land of three products against independent built-up references | figS1_validation.py | pipeline/Tables/validation_metrics.csv, tables/decomposition_300m.csv, sdm/validation_sdm_hotspots_2020.csv, pipeline/Harmonized |
| Supplementary Fig. 2 | Urban land in the cores and 0 to 5 km buffers of protected areas and Key Biodiversity Areas | figS2_site_ranking.py | pipeline/Tables/pa_kba_zonal_GLC30.csv, pipeline/Harmonized/Hotspots_1km_GLC30.tif |
| Supplementary Fig. 3 | Urban land around the Sundarbans and in the habitat of species centered on it | figS3_sundarbans.py | figure_inputs/FigureS3_panels_ab.png, sdm/sdm_species.csv, sdm_work, habitat models |
| Supplementary Figs 4 to 6 | Threatened species above urban share thresholds, urban share of each threatened species, and the most exposed species | figS4_S5_S6.py | sdm/sdm_species.csv, supplemental/Data_S3_species_IUCN_ranges.xlsx |
| Supplementary Fig. 7 | Urban expansion hotspots under the three products | figS7_hotspots_products.py | sdm_work, class areas in pipeline/Tables and multires/Tables |

## Data

| Folder in data | Contents |
|---|---|
| pipeline | Boundaries, the 300 m rasters of the common grid (urban fractions, land area, masks, hotspots, GHSL and GAUD fractions), class areas, transitions, decomposition, zonal and validation tables, and the IUCN range and Sundarbans tables, written by the pipeline stage. |
| multires | GHSL-derived rasters at 100 and 300 m, the ladder tables, the zonal, species and decomposition tables under the GHSL-derived map, and the sensitivity tables, written by the multires stage. |
| sdm | Urban share of each species' modeled habitat, the threat-status tests, models and threshold counts, and the hotspot overlap and validation tables, written by the sdm stage. |
| analysis | Tables written by the analysis scripts. |
| tables | Tables behind Figures 1 and 2 and Tables S1 to S3. |
| supplemental | Data S1 to S3 as Excel workbooks. |
| figure_inputs | Panels a and b of Supplementary Fig. 3. |

The protected area, KBA and IUCN range polygons and the layers built from them are not redistributed. After downloading the World Database on Protected Areas and the World Database of Key Biodiversity Areas (see Sources), `python run_all.py pipeline 00b` writes data/pipeline/BGD_PA_Jan2026.geojson and data/pipeline/BGD_KBA_2021.geojson, which the figure scripts read. With the IUCN ranges added, pipeline/07_species_special.py writes regions.geojson, sundarbans_zones.geojson and narrow_ranges.geojson to data/pipeline/SpeciesSpecial.

## Environment

The scripts were run with Python 3.10, numpy 2.2, pandas 2.3, scipy 1.15, statsmodels 0.15, rasterio 1.4 (GDAL 3.10), geopandas 1.1, shapely 2.1, pyproj 3.7, pyogrio 0.13, matplotlib 3.10, openpyxl 3.1, netCDF4 1.7 and pillow 12.3. environment.yml creates a conda environment (`conda env create -f environment.yml`), and requirements.txt lists the packages for pip.

## Sources

Land cover from GLC_FCS30D (Zhang et al. 2024, Earth System Science Data 16, 1353-1381, https://doi.org/10.5281/zenodo.8239305) and from the Copernicus Climate Change Service, land cover classification gridded maps from 1992 to present derived from satellite observations, versions 2.0.7 and 2.1.1 (https://doi.org/10.24381/cds.006f2c9a). Built-up surface from GHS-BUILT-S R2023A, European Commission, Joint Research Centre (Pesaresi and Politis 2023, https://doi.org/10.2905/9F06F36F-4B11-47EC-ABB0-4F8B7B1D72EA). Urban extent from GAUD (Liu et al. 2020, Nature Sustainability 3, 564-570, https://doi.org/10.6084/m9.figshare.11513178). Protected areas from the World Database on Protected Areas, January 2026 (UNEP-WCMC and IUCN, https://www.protectedplanet.net). Key Biodiversity Areas from the World Database of Key Biodiversity Areas, March 2021 (KBA Partnership, https://www.keybiodiversityareas.org). Species ranges from the IUCN Red List of Threatened Species, spatial data version 2025-1 (https://www.iucnredlist.org). Habitat models, irreplaceability surface and priority areas from the national conservation prioritization of Bangladesh (Chowdhury et al. 2023, One Earth 6, 1315-1325), provided by its authors, with the code at https://github.com/ShawanChowdhury/bd_cons_plan. Administrative boundaries from geoBoundaries (https://www.geoboundaries.org). The data in data keep the licenses of their sources.

## License

CC BY-NC-ND 4.0. See LICENSE.
