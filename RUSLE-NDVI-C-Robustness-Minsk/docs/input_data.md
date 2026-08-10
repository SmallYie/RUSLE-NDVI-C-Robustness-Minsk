# Input data

| Input | GEE type | Required content | Redistribution |
|---|---|---|---|
| District boundaries | FeatureCollection | 22 district polygons and district-name field | user/institution supplied |
| R factor | Image | annual rainfall erosivity | source-derived input |
| K factor | Image | regional soil-erodibility parameterisation | not redistributed |
| L factor | Image | slope-length component | derived input |
| S factor | Image | slope-steepness component | derived input |
| Sentinel-2 | ImageCollection | B4, B8, SCL | public in GEE |
| Land-cover / land-type code | Image | official regional class code | not redistributed |

## Spatial standard

- CRS: EPSG:32635
- nominal analysis scale: 30 m
- year of the principal RUSLE comparison: 2024

## Sentinel-2 preprocessing

The paper uses `COPERNICUS/S2_SR_HARMONIZED`, the full 2024 calendar year, scene-level `CLOUDY_PIXEL_PERCENTAGE < 40`, and masks SCL classes 3, 8, 9, 10 and 11. NDVI is computed from B8 and B4 after the 0.0001 reflectance scaling and summarised as a pixel-wise annual median.

## Restricted inputs

The exact Belarusian soil-map and official land-type inputs are not redistributed in this repository. Users can apply the workflow to equivalent local datasets after matching the required raster structure and projection.
