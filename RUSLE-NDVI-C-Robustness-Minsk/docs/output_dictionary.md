# Output dictionary

## Core regional outputs
- mean / SD / min / max of `A_VDK` and `A_LINEAR`
- robust P50/P80/P90/P95/P99 thresholds
- valid modelling area and coverage

## Five-class agreement
- 5 x 5 VDK-to-Linear transition matrix
- exact agreement
- Cohen's kappa
- quadratic weighted kappa
- mean absolute class difference

## Fixed thresholds
For 5, 10, 20 and 50 t ha-1 yr-1:
- VDK area
- Linear area
- intersection
- union
- VDK-only
- Linear-only
- Jaccard
- Dice

## Matched-percentile hotspots
For nominal Top 20%, Top 10% and Top 5% sets:
- scenario-specific thresholds
- actual selected area percentages
- intersection / union / scenario-only areas
- Jaccard / Dice
- shared fraction of each scenario's selected set

## District prioritisation
For mean A, area >=20 and area >=50:
- VDK and Linear values
- VDK and Linear ranks
- signed and absolute rank change
- stable top-six membership
- scenario-sensitive top-six membership
- Spearman rank correlation

## Domain audit
- exact and raster QA areas
- valid vs excluded NDVI statistics
- valid vs excluded slope statistics
- standardised mean differences
- land-cover composition
- district coverage associations
- coverage vs absolute rank-change associations
