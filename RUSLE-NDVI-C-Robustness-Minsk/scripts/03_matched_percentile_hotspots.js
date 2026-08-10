// ============================================================================
// PUBLIC REPRODUCIBILITY VERSION
// Manuscript: Scale-dependent robustness of RUSLE erosion-hazard prioritisation
//             to NDVI-C parameterisation in Minsk Region, Belarus
//
// IMPORTANT:
// - Replace placeholder Earth Engine asset IDs before running.
// - Restricted/institutional Belarusian source datasets are not redistributed.
// - Exact published outputs require the prepared R, K, L/S or LS, district and
//   land-cover inputs described in the manuscript and repository documentation.
// ============================================================================

// ============================================================================
// COMMON INPUTS FOR THE PAPER ANALYSES
// ============================================================================
var DISTRICTS_ASSET = 'projects/YOUR_EE_PROJECT/assets/MINSK_DISTRICTS';
var VDK_RUSLE_ASSET = 'projects/YOUR_EE_PROJECT/assets/RUSLE_VDK_2024';
var LINEAR_RUSLE_ASSET = 'projects/YOUR_EE_PROJECT/assets/RUSLE_LINEAR_2024';
var DISTRICT_NAME_FIELD = 'F09';
var ANALYSIS_CRS = 'EPSG:32635';
var ANALYSIS_SCALE = 30;

var districts = ee.FeatureCollection(DISTRICTS_ASSET);
var region = districts.geometry();
var rusleVDK0 = ee.Image(VDK_RUSLE_ASSET).select(0).toFloat().rename('A_VDK');
var rusleLinear0 = ee.Image(LINEAR_RUSLE_ASSET).select(0).toFloat().rename('A_LINEAR');

function numberOrZero(value) {
  return ee.Number(ee.Algorithms.If(value, value, 0));
}
function safeDivide(numerator, denominator) {
  numerator = ee.Number(numerator);
  denominator = ee.Number(denominator);
  return ee.Number(ee.Algorithms.If(denominator.neq(0), numerator.divide(denominator), 0));
}
function safePercent(numerator, denominator) {
  return safeDivide(numerator, denominator).multiply(100);
}

var commonMask = rusleVDK0.mask()
  .and(rusleLinear0.mask())
  .and(rusleVDK0.gte(0))
  .and(rusleLinear0.gte(0));
var rusleVDK = rusleVDK0.updateMask(commonMask).clip(region);
var rusleLinear = rusleLinear0.updateMask(commonMask).clip(region);
var pixelAreaHa = ee.Image.pixelArea().divide(10000).rename('area_ha');
var regionValidHa = numberOrZero(pixelAreaHa.updateMask(commonMask).reduceRegion({
  reducer: ee.Reducer.sum(), geometry: region, scale: ANALYSIS_SCALE,
  crs: ANALYSIS_CRS, maxPixels: 1e13, tileScale: 16
}).get('area_ha'));

// 19. EXPERIMENT 1A - FINAL ROBUST VERSION
// LOG-SPACE PERCENTILE THRESHOLD EXTRACTION
//
// Problem addressed:
// Large, strongly right-skewed RUSLE distributions can lead to coarse
// histogram-based percentile estimates in Earth Engine.
//
// Solution:
// 1. Transform A -> log(1 + A)
// 2. Calculate P50, P80, P90, P95 and P99 in log space
// 3. Back-transform threshold:
//        A = exp(logA) - 1
//
// The transformation is strictly monotonic and therefore preserves
// pixel ranking.
// ============================================================================


// ============================================================================
// 19.1 SETTINGS
// ============================================================================

// Histogram settings.
// 4096 buckets provide high resolution in compressed log space
// without creating an excessively heavy reducer.

var EXP1_MAX_BUCKETS = 4096;

var EXP1_MIN_BUCKET_WIDTH = 0.0001;

var EXP1_MAX_RAW = 10000;


// ============================================================================
// 19.2 PERCENTILE REDUCER
// ============================================================================

var exp1LogPercentileReducer =
  ee.Reducer.percentile(

    [
      50,
      80,
      90,
      95,
      99
    ],

    [
      'P50',
      'P80',
      'P90',
      'P95',
      'P99'
    ],

    EXP1_MAX_BUCKETS,

    EXP1_MIN_BUCKET_WIDTH,

    EXP1_MAX_RAW
  );


// ============================================================================
// 19.3 FUNCTION:
// CALCULATE LOG-SPACE PERCENTILES FOR ONE RUSLE SCENARIO
// ============================================================================

function exp1CalculateLogPercentiles(
  rusleImage,
  scenarioName
) {

  // ------------------------------------------------------------
  // Input
  // ------------------------------------------------------------

  var A =
    ee.Image(rusleImage)

      .updateMask(
        commonMask
      )

      .max(0)

      .toFloat();


  // ------------------------------------------------------------
  // Monotonic transformation:
  //
  // logA = ln(1 + A)
  // ------------------------------------------------------------

  var logA =
    A.add(1)

      .log()

      .rename(
        'LOG_A'
      );


  // ------------------------------------------------------------
  // Calculate percentiles in log space.
  // ------------------------------------------------------------

  var stats =
    logA.reduceRegion({

      reducer:
        exp1LogPercentileReducer,

      geometry:
        region,

      scale:
        ANALYSIS_SCALE,

      crs:
        ANALYSIS_CRS,

      maxPixels:
        1e13,

      tileScale:
        16
    });


  stats =
    ee.Dictionary(
      stats
    );


  // ------------------------------------------------------------
  // Explicitly retrieve log-space percentile values.
  //
  // Expected keys:
  //
  // LOG_A_P50
  // LOG_A_P80
  // LOG_A_P90
  // LOG_A_P95
  // LOG_A_P99
  // ------------------------------------------------------------

  var logP50 =
    ee.Number(
      stats.get(
        'LOG_A_P50'
      )
    );


  var logP80 =
    ee.Number(
      stats.get(
        'LOG_A_P80'
      )
    );


  var logP90 =
    ee.Number(
      stats.get(
        'LOG_A_P90'
      )
    );


  var logP95 =
    ee.Number(
      stats.get(
        'LOG_A_P95'
      )
    );


  var logP99 =
    ee.Number(
      stats.get(
        'LOG_A_P99'
      )
    );


  // ------------------------------------------------------------
  // Back-transform:
  //
  // A = exp(logA) - 1
  // ------------------------------------------------------------

  var p50 =
    logP50
      .exp()
      .subtract(1);


  var p80 =
    logP80
      .exp()
      .subtract(1);


  var p90 =
    logP90
      .exp()
      .subtract(1);


  var p95 =
    logP95
      .exp()
      .subtract(1);


  var p99 =
    logP99
      .exp()
      .subtract(1);


  // ------------------------------------------------------------
  // QA tests
  // ------------------------------------------------------------

  var qaP50P80 =
    p50.lt(
      p80
    );


  var qaP80P90 =
    p80.lt(
      p90
    );


  var qaP90P95 =
    p90.lt(
      p95
    );


  var qaP95P99 =
    p95.lt(
      p99
    );


  // ------------------------------------------------------------
  // Output Feature
  // ------------------------------------------------------------

  return ee.Feature(null, {

    Experiment:
      'Equal-area hotspot percentile extraction',

    Scenario:
      scenarioName,

    Spatial_unit:
      'Minsk Region',

    Analysis_year:
      2024,

    Analysis_scale_m:
      ANALYSIS_SCALE,

    Analysis_CRS:
      ANALYSIS_CRS,

    Domain:
      'Common VDK-Linear valid modelling domain',


    // ----------------------------------------------------------
    // ORIGINAL-SCALE RUSLE THRESHOLDS
    // ----------------------------------------------------------

    A_P50:
      p50,

    A_P80:
      p80,

    A_P90:
      p90,

    A_P95:
      p95,

    A_P99:
      p99,


    // ----------------------------------------------------------
    // LOG-SPACE VALUES
    // ----------------------------------------------------------

    LOG_P50:
      logP50,

    LOG_P80:
      logP80,

    LOG_P90:
      logP90,

    LOG_P95:
      logP95,

    LOG_P99:
      logP99,


    // ----------------------------------------------------------
    // QA
    // ----------------------------------------------------------

    QA_P50_lt_P80:
      qaP50P80,

    QA_P80_lt_P90:
      qaP80P90,

    QA_P90_lt_P95:
      qaP90P95,

    QA_P95_lt_P99:
      qaP95P99,


    Histogram_maxBuckets:
      EXP1_MAX_BUCKETS,

    Histogram_minBucketWidth:
      EXP1_MIN_BUCKET_WIDTH,

    Histogram_maxRaw:
      EXP1_MAX_RAW
  });
}


// ============================================================================
// 19.4 VDK
// ============================================================================

var exp1VDKThresholdFeature =
  exp1CalculateLogPercentiles(

    rusleVDK,

    'VDK'
  );


// ============================================================================
// 19.5 LINEAR
// ============================================================================

var exp1LinearThresholdFeature =
  exp1CalculateLogPercentiles(

    rusleLinear,

    'Linear'
  );


// ============================================================================
// 19.6 EXPORT VDK SEPARATELY
//
// Run as one independent Batch Task.
// ============================================================================

Export.table.toDrive({

  collection:
    ee.FeatureCollection([
      exp1VDKThresholdFeature
    ]),

  description:
    'EXP1A_VDK_LogSpace_Percentiles_2024',

  fileFormat:
    'CSV'
});


// ============================================================================
// 19.7 EXPORT LINEAR SEPARATELY
//
// Run as another independent Batch Task.
// ============================================================================

Export.table.toDrive({

  collection:
    ee.FeatureCollection([
      exp1LinearThresholdFeature
    ]),

  description:
    'EXP1A_LINEAR_LogSpace_Percentiles_2024',

  fileFormat:
    'CSV'
});


// ============================================================================
// 19.8 NO HEAVY CONSOLE EVALUATION
//
// Do NOT print the full Features here.
// Do NOT create hotspot rasters yet.
// Do NOT add Map layers yet.
//
// Only lightweight confirmation.
// ============================================================================

print(
  'EXP1A prepared:',
  'Run VDK and Linear percentile exports separately.'
);


// ============================================================================
// 19.9 REGIONAL RUSLE STATISTICS AUDIT
//
// Purpose:
// 1. Verify regional mean / SD on original A scale.
// 2. Use robust log-space P50/P80/P90/P95/P99.
// 3. Cross-check percentiles against the official five-class area statistics.
// ============================================================================


// ---------------------------------------------------------------------------
// Original-scale mean / SD / min / max
// ---------------------------------------------------------------------------

var auditVDKRaw =
  rusleVDK
    .reduceRegion({
      reducer:
        ee.Reducer.mean()
          .combine({
            reducer2: ee.Reducer.stdDev(),
            sharedInputs: true
          })
          .combine({
            reducer2: ee.Reducer.minMax(),
            sharedInputs: true
          }),

      geometry: region,
      scale: ANALYSIS_SCALE,
      crs: ANALYSIS_CRS,
      maxPixels: 1e13,
      tileScale: 16
    });


var auditLinearRaw =
  rusleLinear
    .reduceRegion({
      reducer:
        ee.Reducer.mean()
          .combine({
            reducer2: ee.Reducer.stdDev(),
            sharedInputs: true
          })
          .combine({
            reducer2: ee.Reducer.minMax(),
            sharedInputs: true
          }),

      geometry: region,
      scale: ANALYSIS_SCALE,
      crs: ANALYSIS_CRS,
      maxPixels: 1e13,
      tileScale: 16
    });


// ---------------------------------------------------------------------------
// Robust percentiles already defined in Section 19A
// ---------------------------------------------------------------------------

var auditVDKP50 =
  ee.Number(exp1VDKThresholdFeature.get('A_P50'));

var auditVDKP80 =
  ee.Number(exp1VDKThresholdFeature.get('A_P80'));

var auditVDKP90 =
  ee.Number(exp1VDKThresholdFeature.get('A_P90'));

var auditVDKP95 =
  ee.Number(exp1VDKThresholdFeature.get('A_P95'));

var auditVDKP99 =
  ee.Number(exp1VDKThresholdFeature.get('A_P99'));


var auditLinearP50 =
  ee.Number(exp1LinearThresholdFeature.get('A_P50'));

var auditLinearP80 =
  ee.Number(exp1LinearThresholdFeature.get('A_P80'));

var auditLinearP90 =
  ee.Number(exp1LinearThresholdFeature.get('A_P90'));

var auditLinearP95 =
  ee.Number(exp1LinearThresholdFeature.get('A_P95'));

var auditLinearP99 =
  ee.Number(exp1LinearThresholdFeature.get('A_P99'));


// ---------------------------------------------------------------------------
// Cumulative official-threshold areas
// ---------------------------------------------------------------------------

var auditVDKGE5 =
  ee.Number(regionSummary.get('T5_VDK_pct_valid'));

var auditVDKGE10 =
  ee.Number(regionSummary.get('T10_VDK_pct_valid'));

var auditVDKGE20 =
  ee.Number(regionSummary.get('T20_VDK_pct_valid'));

var auditVDKGE50 =
  ee.Number(regionSummary.get('T50_VDK_pct_valid'));


var auditLinearGE5 =
  ee.Number(regionSummary.get('T5_LINEAR_pct_valid'));

var auditLinearGE10 =
  ee.Number(regionSummary.get('T10_LINEAR_pct_valid'));

var auditLinearGE20 =
  ee.Number(regionSummary.get('T20_LINEAR_pct_valid'));

var auditLinearGE50 =
  ee.Number(regionSummary.get('T50_LINEAR_pct_valid'));


// ---------------------------------------------------------------------------
// Final audit feature
// ---------------------------------------------------------------------------

var regionalStatisticsAudit =
  ee.Feature(null, {

    Spatial_unit:
      'Minsk Region',

    Valid_area_ha:
      regionValidHa,

    Valid_pct_region:
      safePercent(
        regionValidHa,
        totalRegionAreaHa
      ),


    // ============================================================
    // VDK
    // ============================================================

    VDK_mean:
      auditVDKRaw.get('A_VDK_mean'),

    VDK_stdDev:
      auditVDKRaw.get('A_VDK_stdDev'),

    VDK_min:
      auditVDKRaw.get('A_VDK_min'),

    VDK_max:
      auditVDKRaw.get('A_VDK_max'),

    VDK_P50:
      auditVDKP50,

    VDK_P80:
      auditVDKP80,

    VDK_P90:
      auditVDKP90,

    VDK_P95:
      auditVDKP95,

    VDK_P99:
      auditVDKP99,

    VDK_pct_below_5:
      ee.Number(100).subtract(auditVDKGE5),

    VDK_pct_below_10:
      ee.Number(100).subtract(auditVDKGE10),

    VDK_pct_below_20:
      ee.Number(100).subtract(auditVDKGE20),

    VDK_pct_below_50:
      ee.Number(100).subtract(auditVDKGE50),

    VDK_pct_GE5:
      auditVDKGE5,

    VDK_pct_GE10:
      auditVDKGE10,

    VDK_pct_GE20:
      auditVDKGE20,

    VDK_pct_GE50:
      auditVDKGE50,


    // ============================================================
    // LINEAR
    // ============================================================

    LINEAR_mean:
      auditLinearRaw.get('A_LINEAR_mean'),

    LINEAR_stdDev:
      auditLinearRaw.get('A_LINEAR_stdDev'),

    LINEAR_min:
      auditLinearRaw.get('A_LINEAR_min'),

    LINEAR_max:
      auditLinearRaw.get('A_LINEAR_max'),

    LINEAR_P50:
      auditLinearP50,

    LINEAR_P80:
      auditLinearP80,

    LINEAR_P90:
      auditLinearP90,

    LINEAR_P95:
      auditLinearP95,

    LINEAR_P99:
      auditLinearP99,

    LINEAR_pct_below_5:
      ee.Number(100).subtract(auditLinearGE5),

    LINEAR_pct_below_10:
      ee.Number(100).subtract(auditLinearGE10),

    LINEAR_pct_below_20:
      ee.Number(100).subtract(auditLinearGE20),

    LINEAR_pct_below_50:
      ee.Number(100).subtract(auditLinearGE50),

    LINEAR_pct_GE5:
      auditLinearGE5,

    LINEAR_pct_GE10:
      auditLinearGE10,

    LINEAR_pct_GE20:
      auditLinearGE20,

    LINEAR_pct_GE50:
      auditLinearGE50,


    // ============================================================
    // INTERNAL QA
    // ============================================================

    QA_VDK_P50_lt_5:
      auditVDKP50.lt(5),

    QA_LINEAR_P50_gt_10:
      auditLinearP50.gt(10),

    QA_LINEAR_P50_lt_20:
      auditLinearP50.lt(20),

    QA_VDK_P80_lt_P90:
      auditVDKP80.lt(auditVDKP90),

    QA_VDK_P90_lt_P95:
      auditVDKP90.lt(auditVDKP95),

    QA_LINEAR_P80_lt_P90:
      auditLinearP80.lt(auditLinearP90),

    QA_LINEAR_P90_lt_P95:
      auditLinearP90.lt(auditLinearP95)
  });


Export.table.toDrive({

  collection:
    ee.FeatureCollection([
      regionalStatisticsAudit
    ]),

  description:
    'RUSLE_Regional_Statistics_Audit_2024',

  fileFormat:
    'CSV'
});

// ============================================================================
// 19.10 EXPERIMENT 1B
// MATCHED-PERCENTILE / NOMINAL EQUAL-AREA HOTSPOT ROBUSTNESS
//
// Percentile thresholds were obtained in Section 19A.
// No percentile reducer is called in this section.
// ============================================================================


// ============================================================================
// 19.10.1 FIXED THRESHOLDS
// ============================================================================

var EXP1B_VDK_P80 = 13.39958024;
var EXP1B_VDK_P90 = 30.13698445;
var EXP1B_VDK_P95 = 54.92359721;

var EXP1B_LINEAR_P80 = 33.27425248;
var EXP1B_LINEAR_P90 = 58.05030540;
var EXP1B_LINEAR_P95 = 91.12832498;


// ============================================================================
// 19.10.2 BUILD AREA BANDS
// ============================================================================

function exp1BHotspotBands(
  tag,
  vdkThreshold,
  linearThreshold
) {

  var vdkHot =
    rusleVDK
      .gte(vdkThreshold)
      .updateMask(commonMask);

  var linearHot =
    rusleLinear
      .gte(linearThreshold)
      .updateMask(commonMask);

  var intersection =
    vdkHot.and(linearHot);

  var union =
    vdkHot.or(linearHot);

  var vdkOnly =
    vdkHot.and(linearHot.not());

  var linearOnly =
    linearHot.and(vdkHot.not());


  return pixelAreaHa

    .updateMask(vdkHot)
    .rename(tag + '_VDK_ha')

    .addBands(
      pixelAreaHa
        .updateMask(linearHot)
        .rename(tag + '_LINEAR_ha')
    )

    .addBands(
      pixelAreaHa
        .updateMask(intersection)
        .rename(tag + '_Intersection_ha')
    )

    .addBands(
      pixelAreaHa
        .updateMask(union)
        .rename(tag + '_Union_ha')
    )

    .addBands(
      pixelAreaHa
        .updateMask(vdkOnly)
        .rename(tag + '_VDKOnly_ha')
    )

    .addBands(
      pixelAreaHa
        .updateMask(linearOnly)
        .rename(tag + '_LinearOnly_ha')
    );
}


// Build all three percentile levels.

var exp1BAreaStack =
  exp1BHotspotBands(
    'Top20',
    EXP1B_VDK_P80,
    EXP1B_LINEAR_P80
  )

  .addBands(
    exp1BHotspotBands(
      'Top10',
      EXP1B_VDK_P90,
      EXP1B_LINEAR_P90
    )
  )

  .addBands(
    exp1BHotspotBands(
      'Top5',
      EXP1B_VDK_P95,
      EXP1B_LINEAR_P95
    )
  );


// ONE area reduction only.

var exp1BAreaDictionary =
  exp1BAreaStack
    .reduceRegion({

      reducer:
        ee.Reducer.sum(),

      geometry:
        region,

      scale:
        ANALYSIS_SCALE,

      crs:
        ANALYSIS_CRS,

      maxPixels:
        1e13,

      tileScale:
        16
    });


// ============================================================================
// 19.10.3 RESULT FEATURE
// ============================================================================

function exp1BResultFeature(
  tag,
  percentile,
  targetTopPct,
  vdkThreshold,
  linearThreshold
) {

  var d =
    ee.Dictionary(
      exp1BAreaDictionary
    );

  var vdkArea =
    numberOrZero(
      d.get(tag + '_VDK_ha')
    );

  var linearArea =
    numberOrZero(
      d.get(tag + '_LINEAR_ha')
    );

  var intersectionArea =
    numberOrZero(
      d.get(tag + '_Intersection_ha')
    );

  var unionArea =
    numberOrZero(
      d.get(tag + '_Union_ha')
    );

  var vdkOnlyArea =
    numberOrZero(
      d.get(tag + '_VDKOnly_ha')
    );

  var linearOnlyArea =
    numberOrZero(
      d.get(tag + '_LinearOnly_ha')
    );

  var vdkPct =
    safePercent(
      vdkArea,
      regionValidHa
    );

  var linearPct =
    safePercent(
      linearArea,
      regionValidHa
    );


  return ee.Feature(null, {

    Priority_set:
      tag,

    Percentile:
      percentile,

    Target_top_pct:
      targetTopPct,

    VDK_threshold:
      vdkThreshold,

    LINEAR_threshold:
      linearThreshold,

    Valid_area_ha:
      regionValidHa,

    VDK_area_ha:
      vdkArea,

    LINEAR_area_ha:
      linearArea,

    VDK_actual_pct_valid:
      vdkPct,

    LINEAR_actual_pct_valid:
      linearPct,

    VDK_deviation_from_target_pp:
      vdkPct.subtract(targetTopPct),

    LINEAR_deviation_from_target_pp:
      linearPct.subtract(targetTopPct),

    Area_difference_ha:
      linearArea.subtract(vdkArea),

    Area_difference_pct_points:
      linearPct.subtract(vdkPct),

    Intersection_ha:
      intersectionArea,

    Union_ha:
      unionArea,

    VDK_only_ha:
      vdkOnlyArea,

    LINEAR_only_ha:
      linearOnlyArea,

    Jaccard:
      safeDivide(
        intersectionArea,
        unionArea
      ),

    Dice:
      safeDivide(
        intersectionArea.multiply(2),
        vdkArea.add(linearArea)
      ),

    Intersection_fraction_of_VDK:
      safeDivide(
        intersectionArea,
        vdkArea
      ),

    Intersection_fraction_of_LINEAR:
      safeDivide(
        intersectionArea,
        linearArea
      )
  });
}


// ============================================================================
// 19.10.4 FINAL EXP1B TABLE
// ============================================================================

var exp1BHotspotTable =
  ee.FeatureCollection([

    exp1BResultFeature(
      'Top20',
      80,
      20,
      EXP1B_VDK_P80,
      EXP1B_LINEAR_P80
    ),

    exp1BResultFeature(
      'Top10',
      90,
      10,
      EXP1B_VDK_P90,
      EXP1B_LINEAR_P90
    ),

    exp1BResultFeature(
      'Top5',
      95,
      5,
      EXP1B_VDK_P95,
      EXP1B_LINEAR_P95
    )
  ]);


Export.table.toDrive({

  collection:
    exp1BHotspotTable,

  description:
    'EXP1B_RUSLE_MatchedPercentile_Hotspot_Robustness_2024',

  fileFormat:
    'CSV'
});

// ============================================================================