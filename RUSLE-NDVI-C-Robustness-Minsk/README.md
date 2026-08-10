# RUSLE NDVI-C parameterisation robustness in Minsk Region, Belarus

This repository contains the Google Earth Engine (GEE) analysis code and the principal derived tabular outputs supporting the manuscript **“Scale-dependent robustness of RUSLE erosion-hazard prioritisation to NDVI-C parameterisation in Minsk Region, Belarus.”**

## Scientific design

The study is a controlled parameterisation-sensitivity experiment. The rainfall erosivity (R), soil erodibility (K), topographic (LS), and support-practice (P) factors are held identical between scenarios; only the empirical transformation from annual Sentinel-2 NDVI to the RUSLE cover-management factor (C) changes.

Two literature-supported formulations are compared:

- **VDK:** `C = exp[-2 * NDVI / (1 - NDVI)]`
- **Linear:** `C = (1 - NDVI) / 2`

The principal question is how this C-factor choice propagates from continuous erosion magnitude to fixed erosion classes, pixel-level hotspot geography, and district-level prioritisation.

## Repository structure

```text
scripts/
  01_build_C_and_RUSLE_scenarios.js
  02_fixed_threshold_agreement.js
  03_matched_percentile_hotspots.js
  04_district_prioritisation.js
  05_domain_coverage_audit.js
  06_figure3_exports.js
  99_full_pipeline_paper_version.js

docs/
  input_data.md
  workflow.md
  output_dictionary.md
  reproducibility_notes.md

data/
  README.md

results/
  principal CSV outputs used in the manuscript

figures/
  figure3_symbology.md
```

## Important reproducibility note

The repository provides the analytical workflow used in the study, but it does **not** redistribute restricted or institution-specific Belarusian geospatial source data. Exact reproduction of the published Minsk Region outputs requires the prepared regional district, soil-derived K, topographic L/S or LS, rainfall erosivity, and land-cover inputs described in the manuscript.

The hydrological conditioning of the DEM and construction of L and S were performed outside GEE and are treated here as prepared inputs.

## Running the scripts

1. Create or select an Earth Engine project.
2. Upload/provide the required regional inputs described in `docs/input_data.md`.
3. Replace all `projects/YOUR_EE_PROJECT/assets/...` placeholders.
4. Run the scripts in numerical order.
5. Use `99_full_pipeline_paper_version.js` as the archival, paper-version workflow corresponding to the final analysis.

## Public data used directly in GEE

- `COPERNICUS/S2_SR_HARMONIZED`
- `USGS/SRTMGL1_003` (used in the modelling-domain slope audit)

## What is not included

- restricted/institutional Belarusian soil-map source data;
- large GeoTIFF rasters;
- private Earth Engine asset identifiers;
- exploratory experiments not used in the manuscript.

## License

Code is released under the MIT License. This license does not grant redistribution rights for third-party or restricted datasets.

## Citation

Please cite the associated article and this repository when reusing the workflow. A `CITATION.cff` file is included for GitHub citation support.
