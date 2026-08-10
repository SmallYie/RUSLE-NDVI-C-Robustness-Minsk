# Reproducibility notes

## Scope
The public repository is intended to reproduce the **analysis logic** and derived statistics. It is not a complete raw-data redistribution package.

## Prepared RUSLE inputs
The final paper-version analysis script starts from prepared `A_VDK` and `A_LINEAR` rasters because the full regional R, K and LS preparation chain includes external and restricted inputs. Script 01 documents the GEE Sentinel-2/NDVI/C/RUSLE multiplication step when prepared factor rasters are available.

## Percentile calculation
RUSLE distributions are strongly right-skewed. Percentiles are therefore estimated after the strictly monotonic transformation `logA = ln(1 + A)` and back-transformed as `A = exp(logA) - 1`. This preserves pixel rank.

## Priority group
The manuscript reports an explicit **top-six** district priority group. Public scripts therefore use a cutoff rank of 6 so the public scripts and copied result-column names use explicit `top6` terminology.

## Terminology
Public scripts use `erosion class` / `erosion-hazard class` rather than the older internal `risk class` wording, because exposure and vulnerability are not modelled.
