# Minsk Region RUSLE analysis, 2024

This software accompanies the manuscript *Scale-dependent propagation of NDVI–C parameterisation uncertainty in regional RUSLE soil-erosion modelling: evidence from Minsk Region, Belarus*. It compares the exponential Van der Knijff (VDK) and Linear NDVI-to-C scenarios on their common 30 m analysis grid. The revised common modelling domain covers 14,753.21 km², or 36.89% of Minsk Region.

## Analysis scripts

1. [Regional comparisons and map exports](scripts/01_main_rusle_si2024_pipeline.js) report continuous distributions, five-class agreement, fixed-threshold and matched-percentile hotspot overlap, and district rankings.
2. [Land-type and environmental audits](scripts/02_exp2_landcover_and_ndvi_slope_audit_si2024.js) compare valid and excluded land, summarise disagreement by land type, and calculate NDVI and slope differences.
3. [Supplementary district rankings](scripts/03_exp3_extra_district_rankings_si2024.js) compare district medians and area fractions above the 5 and 10 thresholds.

The three files in `scripts/` are byte-for-byte copies of the JavaScript files supplied by the authors for this release. Read the [run notes](docs/RUN_NOTES.md) before running the scripts. They identify the six matched-percentile cutoffs in the revised output and the settings present in the archived scripts.

## Required inputs

Run the scripts in the [Google Earth Engine Code Editor](https://code.earthengine.google.com/) using an Earth Engine project with access to the required assets. Set the district boundary, VDK and Linear scenario images, and official land-type asset identifiers to the corresponding assets in that project. The official land-type field is `LANDTYPE`. The image names used for the revised analysis are `RUSLE_VDK_SI_2024` and `RUSLE_LINEAR_SI_2024`.

The analysis starts from two prepared erosion rasters. Producing those rasters from the national soil data is a separate upstream step. The scripts access Sentinel-2 and SRTM through the public Earth Engine catalog. The Belarusian national digital soil map and official land-type layer are subject to institutional access conditions. The manuscript's supplementary material contains the principal derived tables.

## Data, citation and licence

The [availability statement](docs/DATA_AVAILABILITY_STATEMENT.md) gives public dataset identifiers and describes access to the regional inputs. Software attribution is supplied in [CITATION.cff](CITATION.cff). The code is distributed under the [MIT License](LICENSE).

The supplied software package acknowledges support from the Fundamental Research Funds for the Central Universities (Project No. `2024CDJYDYL006`).
