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
// RUSLE VDK–Linear scenario comparison
// Five-class erosion scheme used in the manuscript
// Minsk Region, Belarus, 2024
//
// Official class boundaries:
// Class 1: A < 5
// Class 2: 5 <= A < 10
// Class 3: 10 <= A < 20
// Class 4: 20 <= A < 50
// Class 5: A >= 50
//
// Main outputs:
// 1. Regional continuous statistics
// 2. District-level continuous and agreement statistics
// 3. Five-class transition matrix
// 4. Class-area table
// 5. Threshold-overlap metrics at 5, 10, 20 and 50
// 6. Cohen's kappa and quadratic weighted kappa
// 7. District-ranking robustness
// 8. Comparison raster
// ============================================================================


// ============================================================================
// 0. USER SETTINGS
// ============================================================================

// Administrative districts.
var districts = ee.FeatureCollection(
  'projects/YOUR_EE_PROJECT/assets/MINSK_DISTRICTS'
);

// District-name field.
// Change only if the final district layer uses another field.
var DISTRICT_NAME_FIELD = 'F09';

// Final RUSLE products.
var rusleVDK0 = ee.Image(
  'projects/YOUR_EE_PROJECT/assets/RUSLE_VDK_2024'
)
  .select(0)
  .toFloat()
  .rename('A_VDK');

var rusleLinear0 = ee.Image(
  'projects/YOUR_EE_PROJECT/assets/RUSLE_LINEAR_2024'
)
  .select(0)
  .toFloat()
  .rename('A_LINEAR');

// Analysis projection and scale.
var ANALYSIS_CRS = 'EPSG:32635';
var ANALYSIS_SCALE = 30;

// Ratio is unstable when the denominator is close to zero.
// This value can be changed, but it must be reported in the manuscript.
var MIN_VDK_FOR_RATIO = 0.1;

// Output asset.
var OUTPUT_COMPARISON_ASSET =
  'projects/YOUR_EE_PROJECT/assets/' +
  'RUSLE_VDK_LINEAR_Comparison_2024_Official5Class';

// Region geometry.
var region = districts.geometry();

Map.centerObject(region, 7);


// ============================================================================
// 1. HELPER FUNCTIONS
// ============================================================================

function numberOrZero(value) {
  return ee.Number(
    ee.Algorithms.If(value, value, 0)
  );
}

function safeDivide(numerator, denominator) {
  numerator = ee.Number(numerator);
  denominator = ee.Number(denominator);

  return ee.Number(
    ee.Algorithms.If(
      denominator.neq(0),
      numerator.divide(denominator),
      0
    )
  );
}

function safePercent(numerator, denominator) {
  return safeDivide(numerator, denominator)
    .multiply(100);
}


// ============================================================================
// 2. COMMON COMPARISON MASK
// ============================================================================

// Compare scenarios only where both RUSLE products are valid.
// Zero values are retained.
var commonMask = rusleVDK0.mask()
  .and(rusleLinear0.mask())
  .and(rusleVDK0.gte(0))
  .and(rusleLinear0.gte(0));

var rusleVDK = rusleVDK0
  .updateMask(commonMask)
  .clip(region);

var rusleLinear = rusleLinear0
  .updateMask(commonMask)
  .clip(region);


// ============================================================================
// 3. CONTINUOUS SCENARIO-DIFFERENCE PRODUCTS
// ============================================================================

// Positive delta means Linear > VDK.
var delta = rusleLinear
  .subtract(rusleVDK)
  .rename('Delta_LINEAR_minus_VDK');

// Absolute magnitude of disagreement.
var absDelta = delta
  .abs()
  .rename('Absolute_Delta');

// Ratio mask.
// Linear must also be positive to avoid log(0).
var ratioMask = commonMask
  .and(rusleVDK.gt(MIN_VDK_FOR_RATIO))
  .and(rusleLinear.gt(0));

var ratio = rusleLinear
  .divide(rusleVDK)
  .updateMask(ratioMask)
  .rename('Ratio_LINEAR_to_VDK');

var logRatio = ratio
  .log()
  .updateMask(ratio.gt(0))
  .rename('Log_Ratio_LINEAR_to_VDK');

// Continuous comparison stack.
var comparisonContinuous = rusleVDK
  .addBands(rusleLinear)
  .addBands(delta)
  .addBands(absDelta)
  .addBands(ratio)
  .addBands(logRatio);


// ============================================================================
// 4. OFFICIAL FIVE-CLASS EROSION SCHEME
// ============================================================================

function classifyErosionClass5(img) {
  return ee.Image.constant(1)

    // Class 2: 5 <= A < 10
    .where(
      img.gte(5).and(img.lt(10)),
      2
    )

    // Class 3: 10 <= A < 20
    .where(
      img.gte(10).and(img.lt(20)),
      3
    )

    // Class 4: 20 <= A < 50
    .where(
      img.gte(20).and(img.lt(50)),
      4
    )

    // Class 5: A >= 50
    .where(
      img.gte(50),
      5
    )

    .updateMask(commonMask)
    .clip(region)
    .toInt8();
}

var classVDK = classifyErosionClass5(rusleVDK)
  .rename('ErosionClass_VDK');

var classLinear = classifyErosionClass5(rusleLinear)
  .rename('ErosionClass_LINEAR');

// Positive values mean Linear is assigned to a higher class.
var classDifference = classLinear
  .subtract(classVDK)
  .rename('ErosionClass_Difference');

// Binary agreement and direction layers.
var sameClass = classVDK.eq(classLinear);

var linearHigher = classLinear.gt(classVDK);
var vdkHigher = classVDK.gt(classLinear);

var linearUp1 = classDifference.eq(1);
var linearUp2Plus = classDifference.gte(2);

var vdkUp1 = classDifference.eq(-1);
var vdkUp2Plus = classDifference.lte(-2);


// ============================================================================
// 5. AREA IMAGES
// ============================================================================

var pixelAreaHa = ee.Image.pixelArea()
  .divide(10000);

var pixelAreaKm2 = ee.Image.pixelArea()
  .divide(1e6)
  .rename('area_km2');


// ----------------------------------------------------------------------------
// 5.1 Threshold-specific overlap bands
// ----------------------------------------------------------------------------

function thresholdComparisonBands(threshold) {

  var tag = 'T' + threshold;

  var vdkHigh = rusleVDK.gte(threshold);
  var linearHigh = rusleLinear.gte(threshold);

  var intersection = vdkHigh.and(linearHigh);
  var union = vdkHigh.or(linearHigh);

  var linearOnly = linearHigh.and(vdkHigh.not());
  var vdkOnly = vdkHigh.and(linearHigh.not());

  return pixelAreaHa
    .updateMask(vdkHigh)
    .rename(tag + '_VDK_ha')

    .addBands(
      pixelAreaHa
        .updateMask(linearHigh)
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
        .updateMask(linearOnly)
        .rename(tag + '_LinearOnly_ha')
    )

    .addBands(
      pixelAreaHa
        .updateMask(vdkOnly)
        .rename(tag + '_VDKOnly_ha')
    );
}

var thresholdBands = thresholdComparisonBands(5)
  .addBands(thresholdComparisonBands(10))
  .addBands(thresholdComparisonBands(20))
  .addBands(thresholdComparisonBands(50));


// ----------------------------------------------------------------------------
// 5.2 General agreement-area bands
// ----------------------------------------------------------------------------

// Area-weighted absolute class difference.
// Dividing its sum by valid area gives the mean absolute class difference.
var weightedAbsoluteClassDifference = pixelAreaHa
  .multiply(classDifference.abs())
  .updateMask(commonMask)
  .rename('Absolute_Class_Difference_x_Area_ha');

var areaBands = pixelAreaHa
  .updateMask(commonMask)
  .rename('Valid_RUSLE_Area_ha')

  .addBands(
    pixelAreaHa
      .updateMask(ratioMask)
      .rename('Ratio_Valid_Area_ha')
  )

  .addBands(
    pixelAreaHa
      .updateMask(sameClass)
      .rename('Same_Class_Area_ha')
  )

  .addBands(
    pixelAreaHa
      .updateMask(linearHigher)
      .rename('Linear_Higher_Class_Area_ha')
  )

  .addBands(
    pixelAreaHa
      .updateMask(vdkHigher)
      .rename('VDK_Higher_Class_Area_ha')
  )

  .addBands(
    pixelAreaHa
      .updateMask(linearUp1)
      .rename('Linear_Up1_Area_ha')
  )

  .addBands(
    pixelAreaHa
      .updateMask(linearUp2Plus)
      .rename('Linear_Up2Plus_Area_ha')
  )

  .addBands(
    pixelAreaHa
      .updateMask(vdkUp1)
      .rename('VDK_Up1_Area_ha')
  )

  .addBands(
    pixelAreaHa
      .updateMask(vdkUp2Plus)
      .rename('VDK_Up2Plus_Area_ha')
  )

  .addBands(weightedAbsoluteClassDifference)
  .addBands(thresholdBands);

// ============================================================================
// 6. STATISTICAL REDUCERS
// ============================================================================
//
// IMPORTANT:
// Mean / SD / min / max are calculated on the original RUSLE scale.
//
// Median and upper percentiles are NOT calculated here.
// Robust percentile statistics are calculated separately in Section 19
// using the monotonic log(1 + A) procedure.
//
// This prevents raw-space approximate percentile outputs from being
// accidentally used in the manuscript.
// ============================================================================

var continuousReducer =
  ee.Reducer.mean()

    .combine({
      reducer2: ee.Reducer.stdDev(),
      sharedInputs: true
    })

    .combine({
      reducer2: ee.Reducer.minMax(),
      sharedInputs: true
    });
// ============================================================================
// 7. THRESHOLD-METRIC FUNCTION
// ============================================================================

function thresholdMetricDictionary(
  areaDictionary,
  threshold,
  validAreaHa
) {

  areaDictionary = ee.Dictionary(areaDictionary);

  var tag = 'T' + threshold;

  var vdkArea = numberOrZero(
    areaDictionary.get(tag + '_VDK_ha')
  );

  var linearArea = numberOrZero(
    areaDictionary.get(tag + '_LINEAR_ha')
  );

  var intersectionArea = numberOrZero(
    areaDictionary.get(tag + '_Intersection_ha')
  );

  var unionArea = numberOrZero(
    areaDictionary.get(tag + '_Union_ha')
  );

  var linearOnlyArea = numberOrZero(
    areaDictionary.get(tag + '_LinearOnly_ha')
  );

  var vdkOnlyArea = numberOrZero(
    areaDictionary.get(tag + '_VDKOnly_ha')
  );

  var jaccard = safeDivide(
    intersectionArea,
    unionArea
  );

  var dice = safeDivide(
    intersectionArea.multiply(2),
    vdkArea.add(linearArea)
  );

  return ee.Dictionary()

    // Absolute areas.
    .set(tag + '_VDK_ha', vdkArea)
    .set(tag + '_LINEAR_ha', linearArea)

    .set(
      tag + '_Intersection_ha',
      intersectionArea
    )

    .set(
      tag + '_Union_ha',
      unionArea
    )

    .set(
      tag + '_LinearOnly_ha',
      linearOnlyArea
    )

    .set(
      tag + '_VDKOnly_ha',
      vdkOnlyArea
    )

    // Percentage of valid model area.
    .set(
      tag + '_VDK_pct_valid',
      safePercent(vdkArea, validAreaHa)
    )

    .set(
      tag + '_LINEAR_pct_valid',
      safePercent(linearArea, validAreaHa)
    )

    .set(
      tag + '_Intersection_pct_valid',
      safePercent(intersectionArea, validAreaHa)
    )

    .set(
      tag + '_Union_pct_valid',
      safePercent(unionArea, validAreaHa)
    )

    .set(
      tag + '_LinearOnly_pct_valid',
      safePercent(linearOnlyArea, validAreaHa)
    )

    .set(
      tag + '_VDKOnly_pct_valid',
      safePercent(vdkOnlyArea, validAreaHa)
    )

    // Spatial-overlap indices.
    .set(tag + '_Jaccard', jaccard)
    .set(tag + '_Dice', dice);
}


// ============================================================================
// 8. REGIONAL CONTINUOUS AND AREA STATISTICS
// ============================================================================

var regionContinuousStats =
  comparisonContinuous.reduceRegion({

    reducer: continuousReducer,
    geometry: region,

    scale: ANALYSIS_SCALE,
    crs: ANALYSIS_CRS,

    maxPixels: 1e13,
    tileScale: 4
  });

var regionAreaStats = areaBands.reduceRegion({

  reducer: ee.Reducer.sum(),
  geometry: region,

  scale: ANALYSIS_SCALE,
  crs: ANALYSIS_CRS,

  maxPixels: 1e13,
  tileScale: 4
});

var totalRegionAreaHa = region.area(1)
  .divide(10000);

var regionValidHa = numberOrZero(
  regionAreaStats.get('Valid_RUSLE_Area_ha')
);

var regionRatioValidHa = numberOrZero(
  regionAreaStats.get('Ratio_Valid_Area_ha')
);

var regionSameHa = numberOrZero(
  regionAreaStats.get('Same_Class_Area_ha')
);

var regionLinearHigherHa = numberOrZero(
  regionAreaStats.get(
    'Linear_Higher_Class_Area_ha'
  )
);

var regionVDKHigherHa = numberOrZero(
  regionAreaStats.get(
    'VDK_Higher_Class_Area_ha'
  )
);

var regionLinearUp1Ha = numberOrZero(
  regionAreaStats.get('Linear_Up1_Area_ha')
);

var regionLinearUp2PlusHa = numberOrZero(
  regionAreaStats.get(
    'Linear_Up2Plus_Area_ha'
  )
);

var regionVDKUp1Ha = numberOrZero(
  regionAreaStats.get('VDK_Up1_Area_ha')
);

var regionVDKUp2PlusHa = numberOrZero(
  regionAreaStats.get(
    'VDK_Up2Plus_Area_ha'
  )
);

var regionAbsoluteClassDifferenceSum =
  numberOrZero(
    regionAreaStats.get(
      'Absolute_Class_Difference_x_Area_ha'
    )
  );

var regionSummary = ee.Feature(null)

  .set(regionContinuousStats)
  .set(regionAreaStats)

  .set({
    Spatial_unit: 'Minsk Region',

    Total_region_area_ha:
      totalRegionAreaHa,

    Valid_model_coverage_pct:
      safePercent(
        regionValidHa,
        totalRegionAreaHa
      ),

    Ratio_valid_pct_of_common:
      safePercent(
        regionRatioValidHa,
        regionValidHa
      ),

    Same_class_pct_valid:
      safePercent(
        regionSameHa,
        regionValidHa
      ),

    Linear_higher_class_pct_valid:
      safePercent(
        regionLinearHigherHa,
        regionValidHa
      ),

    VDK_higher_class_pct_valid:
      safePercent(
        regionVDKHigherHa,
        regionValidHa
      ),

    Linear_up1_pct_valid:
      safePercent(
        regionLinearUp1Ha,
        regionValidHa
      ),

    Linear_up2plus_pct_valid:
      safePercent(
        regionLinearUp2PlusHa,
        regionValidHa
      ),

    VDK_up1_pct_valid:
      safePercent(
        regionVDKUp1Ha,
        regionValidHa
      ),

    VDK_up2plus_pct_valid:
      safePercent(
        regionVDKUp2PlusHa,
        regionValidHa
      ),

    Mean_absolute_class_difference:
      safeDivide(
        regionAbsoluteClassDifferenceSum,
        regionValidHa
      )
  })

  .set(
    thresholdMetricDictionary(
      regionAreaStats,
      5,
      regionValidHa
    )
  )

  .set(
    thresholdMetricDictionary(
      regionAreaStats,
      10,
      regionValidHa
    )
  )

  .set(
    thresholdMetricDictionary(
      regionAreaStats,
      20,
      regionValidHa
    )
  )

  .set(
    thresholdMetricDictionary(
      regionAreaStats,
      50,
      regionValidHa
    )
  );

print(
  'Regional VDK–Linear comparison',
  regionSummary
);


// ============================================================================
// 9. DISTRICT-LEVEL COMPARISON
// ============================================================================

var districtComparison = districts.map(function(f) {

  var geom = f.geometry();

  var cont = comparisonContinuous.reduceRegion({

    reducer: continuousReducer,
    geometry: geom,

    scale: ANALYSIS_SCALE,
    crs: ANALYSIS_CRS,

    maxPixels: 1e13,
    tileScale: 4
  });

  var ar = areaBands.reduceRegion({

    reducer: ee.Reducer.sum(),
    geometry: geom,

    scale: ANALYSIS_SCALE,
    crs: ANALYSIS_CRS,

    maxPixels: 1e13,
    tileScale: 4
  });

  var districtAreaHa = geom.area(1)
    .divide(10000);

  var validHa = numberOrZero(
    ar.get('Valid_RUSLE_Area_ha')
  );

  var ratioValidHa = numberOrZero(
    ar.get('Ratio_Valid_Area_ha')
  );

  var sameHa = numberOrZero(
    ar.get('Same_Class_Area_ha')
  );

  var linearHigherHa = numberOrZero(
    ar.get('Linear_Higher_Class_Area_ha')
  );

  var vdkHigherHa = numberOrZero(
    ar.get('VDK_Higher_Class_Area_ha')
  );

  var linearUp1Ha = numberOrZero(
    ar.get('Linear_Up1_Area_ha')
  );

  var linearUp2PlusHa = numberOrZero(
    ar.get('Linear_Up2Plus_Area_ha')
  );

  var vdkUp1Ha = numberOrZero(
    ar.get('VDK_Up1_Area_ha')
  );

  var vdkUp2PlusHa = numberOrZero(
    ar.get('VDK_Up2Plus_Area_ha')
  );

  var absoluteClassDifferenceSum =
    numberOrZero(
      ar.get(
        'Absolute_Class_Difference_x_Area_ha'
      )
    );

  return f

    .set(cont)
    .set(ar)

    .set({
      unit_type: 'district',

      unit_name:
        f.get(DISTRICT_NAME_FIELD),

      District_Area_ha:
        districtAreaHa,

      Valid_model_coverage_pct:
        safePercent(
          validHa,
          districtAreaHa
        ),

      Ratio_valid_pct_of_common:
        safePercent(
          ratioValidHa,
          validHa
        ),

      Same_class_pct_valid:
        safePercent(
          sameHa,
          validHa
        ),

      Linear_higher_class_pct_valid:
        safePercent(
          linearHigherHa,
          validHa
        ),

      VDK_higher_class_pct_valid:
        safePercent(
          vdkHigherHa,
          validHa
        ),

      Linear_up1_pct_valid:
        safePercent(
          linearUp1Ha,
          validHa
        ),

      Linear_up2plus_pct_valid:
        safePercent(
          linearUp2PlusHa,
          validHa
        ),

      VDK_up1_pct_valid:
        safePercent(
          vdkUp1Ha,
          validHa
        ),

      VDK_up2plus_pct_valid:
        safePercent(
          vdkUp2PlusHa,
          validHa
        ),

      Mean_absolute_class_difference:
        safeDivide(
          absoluteClassDifferenceSum,
          validHa
        )
    })

    .set(
      thresholdMetricDictionary(
        ar,
        5,
        validHa
      )
    )

    .set(
      thresholdMetricDictionary(
        ar,
        10,
        validHa
      )
    )

    .set(
      thresholdMetricDictionary(
        ar,
        20,
        validHa
      )
    )

    .set(
      thresholdMetricDictionary(
        ar,
        50,
        validHa
      )
    );
});

print(
  'District VDK–Linear comparison',
  districtComparison.limit(5)
);


// ============================================================================
// 10. FIVE-BY-FIVE EROSION-CLASS TRANSITION MATRIX
// ============================================================================

// Pair codes:
// 11 = VDK class 1, Linear class 1
// 12 = VDK class 1, Linear class 2
// ...
// 55 = VDK class 5, Linear class 5
var pairCode = classVDK
  .multiply(10)
  .add(classLinear)
  .rename('Pair_Code')
  .toInt16();

var classes = ee.List.sequence(1, 5);

var pairAreaTable = ee.FeatureCollection(

  classes.map(function(vClass) {

    vClass = ee.Number(vClass);

    return classes.map(function(lClass) {

      lClass = ee.Number(lClass);

      var pairMask = classVDK.eq(vClass)
        .and(classLinear.eq(lClass));

      var areaDict = pixelAreaKm2
        .updateMask(pairMask)
        .reduceRegion({

          reducer: ee.Reducer.sum(),
          geometry: region,

          scale: ANALYSIS_SCALE,
          crs: ANALYSIS_CRS,

          maxPixels: 1e13,
          tileScale: 4
        });

      var areaKm2 = numberOrZero(
        areaDict.get('area_km2')
      );

      return ee.Feature(null, {

        VDK_class: vClass,
        LINEAR_class: lClass,

        Pair_code:
          vClass.multiply(10).add(lClass),

        Area_km2:
          areaKm2
      });
    });
  }).flatten()
);

var validAreaKm2 = numberOrZero(

  pixelAreaKm2
    .updateMask(commonMask)
    .reduceRegion({

      reducer: ee.Reducer.sum(),
      geometry: region,

      scale: ANALYSIS_SCALE,
      crs: ANALYSIS_CRS,

      maxPixels: 1e13,
      tileScale: 4
    })

    .get('area_km2')
);

pairAreaTable = pairAreaTable.map(function(f) {

  var area = ee.Number(
    f.get('Area_km2')
  );

  return f.set({

    Area_pct_valid:
      safePercent(
        area,
        validAreaKm2
      ),

    Class_difference:
      ee.Number(
        f.get('LINEAR_class')
      ).subtract(
        ee.Number(
          f.get('VDK_class')
        )
      ),

    Absolute_class_difference:
      ee.Number(
        f.get('LINEAR_class')
      ).subtract(
        ee.Number(
          f.get('VDK_class')
        )
      ).abs()
  });
});

print(
  'Five-class transition matrix',
  pairAreaTable
);


// ============================================================================
// 11. CLASS-AREA TABLE
// ============================================================================

var classLabelDictionary = ee.Dictionary({
  '1': '<5',
  '2': '5-<10',
  '3': '10-<20',
  '4': '20-<50',
  '5': '>=50'
});

var classAreaTable = ee.FeatureCollection(

  classes.map(function(erosionClass) {

    erosionClass = ee.Number(erosionClass);

    var vdkAreaKm2 = numberOrZero(

      pairAreaTable
        .filter(
          ee.Filter.eq(
            'VDK_class',
            erosionClass
          )
        )
        .aggregate_sum('Area_km2')
    );

    var linearAreaKm2 = numberOrZero(

      pairAreaTable
        .filter(
          ee.Filter.eq(
            'LINEAR_class',
            erosionClass
          )
        )
        .aggregate_sum('Area_km2')
    );

    return ee.Feature(null, {

      Erosion_class:
        erosionClass,

     Class_interval:
  classLabelDictionary.get(
    erosionClass.format('%.0f')
  ),

      VDK_area_km2:
        vdkAreaKm2,

      VDK_area_pct_valid:
        safePercent(
          vdkAreaKm2,
          validAreaKm2
        ),

      LINEAR_area_km2:
        linearAreaKm2,

      LINEAR_area_pct_valid:
        safePercent(
          linearAreaKm2,
          validAreaKm2
        )
    });
  })
);

print(
  'Erosion-class area table',
  classAreaTable
);


// ============================================================================
// 12. AGREEMENT INDICES AND KAPPA STATISTICS
// ============================================================================

// Add observed and expected agreement components to every matrix cell.
var transitionAgreementComponents =
  pairAreaTable.map(function(f) {

    var vClass = ee.Number(
      f.get('VDK_class')
    );

    var lClass = ee.Number(
      f.get('LINEAR_class')
    );

    var observedProbability = safeDivide(
      ee.Number(f.get('Area_km2')),
      validAreaKm2
    );

    var vdkMarginalFeature = ee.Feature(

      classAreaTable
        .filter(
          ee.Filter.eq(
            'Erosion_class',
            vClass
          )
        )
        .first()
    );

    var linearMarginalFeature = ee.Feature(

      classAreaTable
        .filter(
          ee.Filter.eq(
            'Erosion_class',
            lClass
          )
        )
        .first()
    );

    var vdkMarginalProbability =
      safeDivide(
        vdkMarginalFeature.get(
          'VDK_area_km2'
        ),
        validAreaKm2
      );

    var linearMarginalProbability =
      safeDivide(
        linearMarginalFeature.get(
          'LINEAR_area_km2'
        ),
        validAreaKm2
      );

    var expectedProbability =
      vdkMarginalProbability
        .multiply(
          linearMarginalProbability
        );

    // Exact agreement weight:
    // 1 for the same class, otherwise 0.
    var exactWeight = vClass.eq(lClass);

    // Quadratic agreement weight:
    // 1 for exact agreement;
    // smaller values for increasingly distant classes.
    // Maximum possible distance is 4.
    var quadraticWeight = ee.Number(1)
      .subtract(
        vClass
          .subtract(lClass)
          .pow(2)
          .divide(16)
      );

    return f.set({

      Observed_exact_component:
        observedProbability
          .multiply(exactWeight),

      Expected_exact_component:
        expectedProbability
          .multiply(exactWeight),

      Observed_quadratic_component:
        observedProbability
          .multiply(quadraticWeight),

      Expected_quadratic_component:
        expectedProbability
          .multiply(quadraticWeight)
    });
  });

var observedExactAgreement = ee.Number(

  transitionAgreementComponents
    .aggregate_sum(
      'Observed_exact_component'
    )
);

var expectedExactAgreement = ee.Number(

  transitionAgreementComponents
    .aggregate_sum(
      'Expected_exact_component'
    )
);

var observedQuadraticAgreement = ee.Number(

  transitionAgreementComponents
    .aggregate_sum(
      'Observed_quadratic_component'
    )
);

var expectedQuadraticAgreement = ee.Number(

  transitionAgreementComponents
    .aggregate_sum(
      'Expected_quadratic_component'
    )
);

var cohenKappa = safeDivide(

  observedExactAgreement
    .subtract(expectedExactAgreement),

  ee.Number(1)
    .subtract(expectedExactAgreement)
);

var quadraticWeightedKappa = safeDivide(

  observedQuadraticAgreement
    .subtract(expectedQuadraticAgreement),

  ee.Number(1)
    .subtract(expectedQuadraticAgreement)
);

var transitionSummary = ee.Feature(null, {

  Spatial_unit:
    'Minsk Region',

  Number_of_classes:
    5,

  Valid_area_km2:
    validAreaKm2,

  Overall_exact_agreement:
    observedExactAgreement,

  Overall_exact_agreement_pct:
    observedExactAgreement.multiply(100),

  Expected_exact_agreement:
    expectedExactAgreement,

  Cohens_kappa:
    cohenKappa,

  Observed_quadratic_agreement:
    observedQuadraticAgreement,

  Expected_quadratic_agreement:
    expectedQuadraticAgreement,

  Quadratic_weighted_kappa:
    quadraticWeightedKappa,

  Mean_absolute_class_difference:
    safeDivide(
      regionAbsoluteClassDifferenceSum,
      regionValidHa
    )
});

print(
  'Transition agreement summary',
  transitionSummary
);



// 15. TABLE EXPORTS
// ============================================================================

// 15.1 Regional comparison.
Export.table.toDrive({

  collection:
    ee.FeatureCollection([regionSummary]),

  description:
    'RUSLE_VDK_LINEAR_Region_Comparison_Official5Class_2024',

  fileFormat:
    'CSV'
});


// 15.2 District comparison.
Export.table.toDrive({

  collection:
    districtComparison,

  description:
    'RUSLE_VDK_LINEAR_District_Comparison_Official5Class_2024',

  fileFormat:
    'CSV'
});


// 15.3 Five-by-five transition matrix.
Export.table.toDrive({

  collection:
    pairAreaTable,

  description:
    'RUSLE_VDK_LINEAR_5Class_Transition_Matrix_2024',

  fileFormat:
    'CSV'
});


// 15.4 Class-area table.
Export.table.toDrive({

  collection:
    classAreaTable,

  description:
    'RUSLE_VDK_LINEAR_5Class_Area_Statistics_2024',

  fileFormat:
    'CSV'
});


// 15.5 Transition agreement summary.
Export.table.toDrive({

  collection:
    ee.FeatureCollection([
      transitionSummary
    ]),

  description:
    'RUSLE_VDK_LINEAR_5Class_Agreement_Summary_2024',

  fileFormat:
    'CSV'
});
