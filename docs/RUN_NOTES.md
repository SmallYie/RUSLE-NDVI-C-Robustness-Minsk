# Settings for the revised 2024 results

The source archive designated by the authors is `minsk-rusle-ndvi-c-robustness-github.zip` (SHA-256 `e9601be55291be6d8bf570377df62eca5affdaeaf48d682d72e9d50b819dfdd2`). The JavaScript, licence, software citation and ignore file retain their supplied contents. The repository overview and availability text were updated for this public deposit, and these notes were added. The supplied citation metadata identifies software contributors.

## Input paths

Replace the full `projects/YOUR_EE_PROJECT/assets/...` paths with authorised asset identifiers. The two default scenario names in the main script repeat the suffix (`_SI_SI_2024`). The revised input names end in `_SI_2024`. The two modular scripts already use that spelling. Use the configured output destination in your own project for any optional exports.

## Matched-percentile hotspots

Section 19A of the main script computes the scenario-specific percentiles from the input images. Section 19.10.1 separately assigns six numeric constants to the hotspot calculation. The static values in the supplied file and the thresholds recorded in the authors' revised EXP1B export differ.

| Priority set | VDK default | Linear default | VDK in revised EXP1B | Linear in revised EXP1B |
| --- | ---: | ---: | ---: | ---: |
| Top 20% | 13.39958024 | 33.27425248 | 1.8066835252070508 | 4.221570814139686 |
| Top 10% | 30.13698445 | 58.05030540 | 3.9611039542280633 | 7.357909009689854 |
| Top 5% | 54.92359721 | 91.12832498 | 7.146722059512486 | 11.588918459094655 |

To reproduce the revised table with its prepared images, link the variables in Section 19.10.1 to the Section 19A audit results before running the EXP1B export:

```javascript
var EXP1B_VDK_P80 = auditVDKP80;
var EXP1B_VDK_P90 = auditVDKP90;
var EXP1B_VDK_P95 = auditVDKP95;
var EXP1B_LINEAR_P80 = auditLinearP80;
var EXP1B_LINEAR_P90 = auditLinearP90;
var EXP1B_LINEAR_P95 = auditLinearP95;
```

The revised EXP1B table reports Jaccard values of 0.477, 0.439 and 0.412 for these three priority sets. The original case-specific flag `QA_LINEAR_P50_gt_10` in the main script evaluates to false for the revised regional median of 1.438. This flag is a reported diagnostic.

## Environmental and district summaries

The environmental audit in the main script uses the Sentinel-2 scene classification layer (SCL classes 3, 8, 9, 10 and 11 masked) and scene cloudiness below 40%. The stand-alone land-type and environmental audit uses QA60 cloud bits 10 and 11 and scene cloudiness below 20%. Specify the script and masking rule when reporting NDVI summaries.

The main regional percentile calculation uses log(1 + A) followed by back-transformation. The supplied stand-alone district-ranking script applies the Earth Engine median reducer directly to the district A values. Regional and district summaries should be described with their respective spatial units and reducers.

The land-type audit in script 02 reports area-weighted mean A, signed mean differences and the proportions assigned to the same or different classes. Exact reproduction of these outputs requires the prepared scenario rasters and access to the supporting regional layers. The source archive documents the saved settings and does not contain a history of individual Code Editor sessions.
