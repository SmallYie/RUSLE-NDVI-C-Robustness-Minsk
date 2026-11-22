# Minsk Region RUSLE analysis code

## Revised SI 2024 manuscript workflow

The [fixed manuscript version](https://github.com/SmallYie/RUSLE-NDVI-C-Robustness-Minsk/tree/2b811c2d6023f0440e1f044cae98d340d646ffd7) contains the three author-supplied Earth Engine analysis scripts and the settings for interpreting their outputs. The [downloadable ZIP](https://codeload.github.com/SmallYie/RUSLE-NDVI-C-Robustness-Minsk/zip/2b811c2d6023f0440e1f044cae98d340d646ffd7) contains this version alone. The [run notes](https://github.com/SmallYie/RUSLE-NDVI-C-Robustness-Minsk/tree/2b811c2d6023f0440e1f044cae98d340d646ffd7/docs) cover input assets and matched-percentile thresholds.

## Earlier analysis

The following documentation and the original files on this branch describe the earlier analysis.

This repository contains the Google Earth Engine code and derived tabular
outputs supporting the manuscript:

"Scale-dependent robustness of RUSLE erosion-hazard prioritisation to
NDVI–C parameterisation in Minsk Region, Belarus"

## Study design

The study compares two literature-supported NDVI–C parameterisations:

C_VDK = exp[-2 NDVI / (1 - NDVI)]

C_Linear = (1 - NDVI) / 2

R, K, LS and P are held identical between scenarios.

## Repository structure

scripts/
docs/
data/
results/

## Main analyses

1. RUSLE scenario construction
2. Fixed-threshold class agreement
3. Matched-percentile hotspot robustness
4. District prioritisation robustness
5. Modelling-domain coverage audit

## Required data

Sentinel-2 and SRTM are publicly available.
The Belarusian soil-map and official land-type datasets are not redistributed.

## Reproducibility

Users should replace placeholder Earth Engine asset IDs with their own
equivalent datasets.

## Citation

Please cite the associated article when using this code.
