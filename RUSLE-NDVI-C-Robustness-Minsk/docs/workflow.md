# Workflow

1. Prepare R, K, L/S (or LS), district boundaries and land-cover inputs.
2. Build the 2024 annual-median Sentinel-2 NDVI composite.
3. Construct the VDK and Linear C-factor surfaces.
4. Generate `A_VDK` and `A_LINEAR` with identical R, K, LS and P=1.
5. Restrict comparisons to the common valid non-negative modelling domain.
6. Apply the shared five-class scheme: `<5`, `5-<10`, `10-<20`, `20-<50`, `>=50` t ha-1 yr-1.
7. Quantify fixed-threshold agreement (transition matrix, kappa, weighted kappa, Jaccard and Dice).
8. Estimate robust log-space percentiles and compare matched Top 20%, Top 10% and Top 5% hotspot footprints.
9. Rank districts by mean A and by the percentages at or above 20 and 50 t ha-1 yr-1; compare ranks and top-six membership.
10. Audit modelling-domain selectivity using NDVI, slope, land-cover composition, district coverage, and coverage-vs-rank-instability diagnostics.
11. Export the three classified rasters used to prepare Figure 3.
