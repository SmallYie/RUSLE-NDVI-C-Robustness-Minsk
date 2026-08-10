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


// ============================================================================
// 13. DISTRICT-RANKING ROBUSTNESS
// ============================================================================

// Higher values receive better priority ranks:
// rank 1 = highest value.
function addDescendingRank(
  collection,
  valueField,
  rankField
) {

  var sorted = collection.sort(
    valueField,
    false
  );

  var featureCount = sorted.size();

  var featureList = sorted.toList(
    featureCount
  );

  return ee.FeatureCollection(

    ee.List.sequence(
      0,
      featureCount.subtract(1)
    )

    .map(function(index) {

      index = ee.Number(index);

      return ee.Feature(
        featureList.get(index)
      )

      .set(
        rankField,
        index.add(1)
      );
    })
  );
}


// Build a rank-comparison table for any pair of VDK and Linear indicators.
function buildRankComparison(
  collection,
  vdkField,
  linearField,
  rankBasis
) {

  var rankedVDK = addDescendingRank(
    collection,
    vdkField,
    'VDK_rank'
  );

  var rankedLinear = addDescendingRank(
    collection,
    linearField,
    'LINEAR_rank'
  );

  var joinedRanks = ee.Join.inner().apply({

    primary: rankedVDK,
    secondary: rankedLinear,

    condition: ee.Filter.equals({
      leftField: 'unit_name',
      rightField: 'unit_name'
    })
  });

  var numberOfDistricts =
    collection.size();

  // Manuscript priority group: ranks 1-6.
  // Kept explicit so the code matches the reported top-six analysis.
  var top6CutoffRank = ee.Number(6);

  var rankTable = ee.FeatureCollection(

    joinedRanks.map(function(joinedFeature) {

      joinedFeature =
        ee.Feature(joinedFeature);

      var vdkFeature = ee.Feature(
        joinedFeature.get('primary')
      );

      var linearFeature = ee.Feature(
        joinedFeature.get('secondary')
      );

      var vdkRank = ee.Number(
        vdkFeature.get('VDK_rank')
      );

      var linearRank = ee.Number(
        linearFeature.get('LINEAR_rank')
      );

      var isVDKTop6 =
        vdkRank.lte(top6CutoffRank);

      var isLinearTop6 =
        linearRank.lte(top6CutoffRank);

      var stableTop6 =
        isVDKTop6.and(isLinearTop6);

      var scenarioSensitiveTop6 =
        isVDKTop6.neq(isLinearTop6);

      return vdkFeature.set({

        Rank_basis:
          rankBasis,

        VDK_value:
          vdkFeature.get(vdkField),

        LINEAR_value:
          linearFeature.get(linearField),

        VDK_rank:
          vdkRank,

        LINEAR_rank:
          linearRank,

        Rank_change_LINEAR_minus_VDK:
          linearRank.subtract(vdkRank),

        Absolute_rank_change:
          linearRank
            .subtract(vdkRank)
            .abs(),

        Stable_top6:
          stableTop6,

        Scenario_sensitive_top6:
          scenarioSensitiveTop6
      });
    })
  );

  // Pearson correlation applied to ranks is Spearman rank correlation.
  var rankCorrelation =
    rankTable.reduceColumns({

      reducer:
        ee.Reducer.pearsonsCorrelation(),

      selectors: [
        'VDK_rank',
        'LINEAR_rank'
      ]
    });

  var stableTop6Districts =
    rankTable

      .filter(
        ee.Filter.eq(
          'Stable_top6',
          1
        )
      )

      .aggregate_array('unit_name');

  var sensitiveTop6Districts =
    rankTable

      .filter(
        ee.Filter.eq(
          'Scenario_sensitive_top6',
          1
        )
      )

      .aggregate_array('unit_name');

  var rankSummary = ee.Feature(null, {

    Rank_basis:
      rankBasis,

    VDK_indicator:
      vdkField,

    LINEAR_indicator:
      linearField,

    Number_of_districts:
      numberOfDistricts,

    Top6_cutoff_rank:
      top6CutoffRank,

    Spearman_rank_correlation:
      rankCorrelation.get('correlation'),

    Spearman_p_value:
      rankCorrelation.get('p-value'),

    Mean_absolute_rank_change:
      rankTable.aggregate_mean(
        'Absolute_rank_change'
      ),

    Maximum_absolute_rank_change:
      rankTable.aggregate_max(
        'Absolute_rank_change'
      ),

    Stable_top6_district_count:
      rankTable

        .filter(
          ee.Filter.eq(
            'Stable_top6',
            1
          )
        )

        .size(),

    Scenario_sensitive_top6_district_count:
      rankTable

        .filter(
          ee.Filter.eq(
            'Scenario_sensitive_top6',
            1
          )
        )

        .size(),

    Stable_top6_districts:
      stableTop6Districts,

    Scenario_sensitive_top6_districts:
      sensitiveTop6Districts
  });

  return {
    table: rankTable,
    summary: rankSummary
  };
}


// ----------------------------------------------------------------------------
// 13.1 Ranking based on district mean RUSLE
// ----------------------------------------------------------------------------

var meanRankResult = buildRankComparison(

  districtComparison,

  'A_VDK_mean',
  'A_LINEAR_mean',

  'District mean RUSLE'
);


// ----------------------------------------------------------------------------
// 13.2 Ranking based on area at or above 20
// ----------------------------------------------------------------------------

var t20RankResult = buildRankComparison(

  districtComparison,

  'T20_VDK_pct_valid',
  'T20_LINEAR_pct_valid',

  'Area percentage A >= 20'
);


// ----------------------------------------------------------------------------
// 13.3 Ranking based on area at or above 50
// ----------------------------------------------------------------------------

var t50RankResult = buildRankComparison(

  districtComparison,

  'T50_VDK_pct_valid',
  'T50_LINEAR_pct_valid',

  'Area percentage A >= 50'
);

var rankSummaryCollection =
  ee.FeatureCollection([

    meanRankResult.summary,
    t20RankResult.summary,
    t50RankResult.summary
  ]);

print(
  'District mean ranking',
  meanRankResult.table
);

print(
  'District T20 ranking',
  t20RankResult.table
);

print(
  'District T50 ranking',
  t50RankResult.table
);

print(
  'District ranking summaries',
  rankSummaryCollection
);


// ============================================================================
// 14. COMPARISON RASTER
// ============================================================================
//
// All bands are explicitly converted to Float32 before export.
// This avoids mixed Float32 / Float64 / integer band-type errors.
//
// Integer-valued class bands remain numerically exact:
// 1, 2, 3, 4, 5 and pair codes such as 11, 12, ..., 55.
// ============================================================================

var scenarioComparisonRaster = rusleVDK
  .toFloat()

  .addBands(
    rusleLinear.toFloat()
  )

  .addBands(
    delta.toFloat()
  )

  .addBands(
    absDelta.toFloat()
  )

  .addBands(
    ratio.toFloat()
  )

  .addBands(
    logRatio.toFloat()
  )

  .addBands(
    classVDK.toFloat()
  )

  .addBands(
    classLinear.toFloat()
  )

  .addBands(
    classDifference.toFloat()
  )

  .addBands(
    pairCode.toFloat()
  )

  .toFloat();
// ============================================================================
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


// 15.6 District ranking based on mean RUSLE.
Export.table.toDrive({

  collection:
    meanRankResult.table,

  description:
    'RUSLE_VDK_LINEAR_District_Ranking_Mean_2024',

  fileFormat:
    'CSV'
});


// 15.7 District ranking based on A >= 20.
Export.table.toDrive({

  collection:
    t20RankResult.table,

  description:
    'RUSLE_VDK_LINEAR_District_Ranking_T20_2024',

  fileFormat:
    'CSV'
});


// 15.8 District ranking based on A >= 50.
Export.table.toDrive({

  collection:
    t50RankResult.table,

  description:
    'RUSLE_VDK_LINEAR_District_Ranking_T50_2024',

  fileFormat:
    'CSV'
});


// 15.9 Combined ranking summaries.
Export.table.toDrive({

  collection:
    rankSummaryCollection,

  description:
    'RUSLE_VDK_LINEAR_District_Ranking_Summaries_2024',

  fileFormat:
    'CSV'
});


// ============================================================================
// 16. RASTER EXPORT
// ============================================================================

Export.image.toAsset({

  image:
    scenarioComparisonRaster,

  description:
    'Export_RUSLE_VDK_LINEAR_Comparison_Official5Class_2024',

  assetId:
    OUTPUT_COMPARISON_ASSET,

  region:
    region,

  scale:
    ANALYSIS_SCALE,

  crs:
    ANALYSIS_CRS,

  maxPixels:
    1e13
});


// ============================================================================
// 17. MAP DISPLAY
// ============================================================================

var ruslePalette = [
  '006400',
  '7FFF00',
  'FFFF00',
  'FFA500',
  'FF0000',
  '7F0000'
];

var classPalette = [
  '1A9850',
  '91CF60',
  'FEE08B',
  'FC8D59',
  'D73027'
];

var classDifferencePalette = [
  '081D58', // -4
  '253494', // -3
  '225EA8', // -2
  '41B6C4', // -1
  'F7F7F7', //  0
  'FDAE6B', // +1
  'F16913', // +2
  'D94801', // +3
  '7F0000'  // +4
];

Map.addLayer(
  rusleVDK,
  {
    min: 0,
    max: 50,
    palette: ruslePalette
  },
  'RUSLE VDK',
  false
);

Map.addLayer(
  rusleLinear,
  {
    min: 0,
    max: 50,
    palette: ruslePalette
  },
  'RUSLE Linear',
  false
);

Map.addLayer(
  classVDK,
  {
    min: 1,
    max: 5,
    palette: classPalette
  },
  'VDK five-class erosion',
  true
);

Map.addLayer(
  classLinear,
  {
    min: 1,
    max: 5,
    palette: classPalette
  },
  'Linear five-class erosion',
  true
);

Map.addLayer(
  delta,
  {
    min: -25,
    max: 25,
    palette: [
      '08306B',
      '6BAED6',
      'FFFFFF',
      'FC9272',
      '99000D'
    ]
  },
  'Linear minus VDK',
  false
);

Map.addLayer(
  absDelta,
  {
    min: 0,
    max: 30,
    palette: [
      'FFFFFF',
      'FEE8C8',
      'FDBB84',
      'E34A33',
      '7F0000'
    ]
  },
  'Absolute scenario difference',
  false
);

Map.addLayer(
  classDifference,
  {
    min: -4,
    max: 4,
    palette: classDifferencePalette
  },
  'Erosion-class difference',
  true
);

Map.addLayer(
  pairCode.randomVisualizer(),
  {},
  'VDK–Linear five-class pairs',
  false
);

Map.addLayer(
  districts.style({
    color: '000000',
    fillColor: '00000000',
    width: 1
  }),
  {},
  'District boundaries',
  true
);


// ============================================================================
// 18. FINAL CHECKS
// ============================================================================

print(
  'Class definitions',
  ee.Dictionary({
    Class_1: 'A < 5',
    Class_2: '5 <= A < 10',
    Class_3: '10 <= A < 20',
    Class_4: '20 <= A < 50',
    Class_5: 'A >= 50'
  })
);

print(
  'Analysis CRS',
  ANALYSIS_CRS
);

print(
  'Analysis scale',
  ANALYSIS_SCALE
);

print(
  'Minimum VDK denominator for ratio',
  MIN_VDK_FOR_RATIO
);

// ============================================================================
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
// 20. EXPERIMENT 2
// MODELLING-DOMAIN COMPOSITION AUDIT
//
// Purpose:
// Quantify how the common valid RUSLE domain differs from the excluded
// part of Minsk Region.
//
// Components:
//
// A. Area consistency QA
// B. NDVI composition
// C. Slope composition
// D. Land-cover composition
// E. District coverage versus erosion metrics
//
// IMPORTANT:
// This is a domain-composition / representativeness audit.
// It does NOT claim that the valid domain is an unbiased random sample.
// ============================================================================

// ============================================================================
// SHARED AREA HELPER FOR EXPERIMENT 2
// ============================================================================

var expPixelAreaHa = ee.Image.pixelArea()
  .divide(10000)
  .rename('area_ha');


function expMaskedAreaHa(maskImage, geometry) {

  var areaResult = expPixelAreaHa
    .updateMask(
      ee.Image(maskImage)
        .unmask(0)
        .gt(0)
    )
    .reduceRegion({

      reducer: ee.Reducer.sum(),

      geometry: geometry,

      scale: ANALYSIS_SCALE,
      crs: ANALYSIS_CRS,

      maxPixels: 1e13,
      tileScale: 16
    });

  return numberOrZero(
    areaResult.get('area_ha')
  );
}
// ============================================================================
// 20.1 USER SETTINGS FOR LAND COVER
// ============================================================================

// Existing land-cover raster.
// If your final asset name is different, change ONLY this line.

var EXP2_LANDCOVER_ASSET =
  'projects/YOUR_EE_PROJECT/assets/LANDCOVER_CODE';


var exp2Landcover =
  ee.Image(
    EXP2_LANDCOVER_ASSET
  )
  .select(0)
  .rename('landcover_code')
  .toInt16()
  .clip(region);


// ============================================================================
// 20.2 EXPLICIT VALID / EXCLUDED DOMAIN MASKS
// ============================================================================

// Validity of each RUSLE raster.

var exp2VDKValid =
  rusleVDK0
    .mask()
    .unmask(0)
    .gt(0)

    .and(
      rusleVDK0
        .unmask(-9999)
        .gte(0)
    );


var exp2LinearValid =
  rusleLinear0
    .mask()
    .unmask(0)
    .gt(0)

    .and(
      rusleLinear0
        .unmask(-9999)
        .gte(0)
    );


var exp2CommonValidBinary =
  exp2VDKValid

    .and(
      exp2LinearValid
    )

    .clip(region)

    .rename(
      'valid_domain'
    );


// Region mask.

var exp2RegionMask =
  ee.Image.constant(1)
    .clip(region)
    .selfMask();


// Valid domain.

var exp2ValidMask =
  exp2CommonValidBinary
    .eq(1)
    .updateMask(
      exp2RegionMask
    );


// Excluded domain.

var exp2ExcludedMask =
  exp2CommonValidBinary
    .unmask(0)
    .eq(0)
    .updateMask(
      exp2RegionMask
    );


// ============================================================================
// 20.3 DOMAIN AREA QA
// ============================================================================

// Exact vector area of Minsk Region.

var exp2RegionGeometryAreaHa =
  region.area(1)
    .divide(10000);


// Raster-based valid domain QA.

var exp2ValidRasterAreaHa =
  expMaskedAreaHa(
    exp2ValidMask,
    region
  );


// Raster-based excluded-domain QA.
// This is retained ONLY for rasterisation diagnostics.

var exp2ExcludedRasterAreaHa =
  expMaskedAreaHa(
    exp2ExcludedMask,
    region
  );


// Main valid domain from the principal RUSLE analysis.

var exp2ValidExactAreaHa =
  regionValidHa;


// Exact excluded area used for reporting.

var exp2ExcludedExactAreaHa =
  exp2RegionGeometryAreaHa
    .subtract(
      exp2ValidExactAreaHa
    );


// Raster sum: diagnostic only.

var exp2RasterMaskSumAreaHa =
  exp2ValidRasterAreaHa
    .add(
      exp2ExcludedRasterAreaHa
    );


var exp2AreaQA =
  ee.Feature(null, {

    Region_geometry_area_ha:
      exp2RegionGeometryAreaHa,

    Main_analysis_valid_area_ha:
      regionValidHa,

    EXP2_valid_raster_area_ha:
      exp2ValidRasterAreaHa,

    Exact_excluded_area_ha:
      exp2ExcludedExactAreaHa,

    Raster_excluded_area_QA_ha:
      exp2ExcludedRasterAreaHa,

    Raster_valid_plus_excluded_ha:
      exp2RasterMaskSumAreaHa,

    Raster_minus_geometry_ha:
      exp2RasterMaskSumAreaHa.subtract(
        exp2RegionGeometryAreaHa
      ),

    EXP2_valid_minus_main_valid_ha:
      exp2ValidRasterAreaHa.subtract(
        regionValidHa
      ),

    Valid_pct_region_geometry:
      safePercent(
        regionValidHa,
        exp2RegionGeometryAreaHa
      ),

    Exact_excluded_pct_region_geometry:
      safePercent(
        exp2ExcludedExactAreaHa,
        exp2RegionGeometryAreaHa
      )
  });


// Do NOT force this large calculation interactively.
// Use Export only.

// print(
//   'EXP2 Domain area QA',
//   exp2AreaQA
// );



// ============================================================================
// 20.4 SENTINEL-2 NDVI COLLECTION
//
// Same preprocessing logic as the main RUSLE workflow:
//
// 2024 full year
// CLOUDY_PIXEL_PERCENTAGE < 40
// mask SCL 3, 8, 9, 10, 11
// ============================================================================

function expPrepareS2NDVI(image) {

  var scl =
    image.select('SCL');


  var qualityMask =
    scl.neq(3)

      .and(
        scl.neq(8)
      )

      .and(
        scl.neq(9)
      )

      .and(
        scl.neq(10)
      )

      .and(
        scl.neq(11)
      );


  var red =
    image.select('B4')
      .multiply(0.0001);


  var nir =
    image.select('B8')
      .multiply(0.0001);


  var denominator =
    nir.add(red);


  var ndvi =
    nir.subtract(red)

      .divide(
        denominator
      )

      .rename('NDVI')

      .updateMask(
        denominator.neq(0)
      )

      .updateMask(
        qualityMask
      );


  return ndvi
    .copyProperties(
      image,
      ['system:time_start']
    );
}


var expS2NDVI2024 =
  ee.ImageCollection(
    'COPERNICUS/S2_SR_HARMONIZED'
  )

  .filterBounds(
    region
  )

  .filterDate(
    '2024-01-01',
    '2025-01-01'
  )

  .filter(
    ee.Filter.lt(
      'CLOUDY_PIXEL_PERCENTAGE',
      40
    )
  )

  .map(
    expPrepareS2NDVI
  );


var expAnnualMedianNDVI =
  expS2NDVI2024
    .median()
    .rename('NDVI_2024')
    .clip(region);


// ============================================================================
// 20.5 SRTM SLOPE
// ============================================================================

var exp2DEM =
  ee.Image(
    'USGS/SRTMGL1_003'
  )
  .select('elevation')
  .clip(region);


var exp2Slope =
  ee.Terrain
    .slope(exp2DEM)
    .rename('Slope_deg');


// ============================================================================
// 20.6 CONTINUOUS AUDIT VARIABLES
// ============================================================================

var exp2AuditVariables =
  expAnnualMedianNDVI
    .addBands(
      exp2Slope
    );


var exp2AuditReducer =
  ee.Reducer.mean()

    .combine({

      reducer2:
        ee.Reducer.median(),

      sharedInputs:
        true
    })

    .combine({

      reducer2:
        ee.Reducer.stdDev(),

      sharedInputs:
        true
    })

    .combine({

      reducer2:
        ee.Reducer.percentile([
          25,
          75
        ]),

      sharedInputs:
        true
    });


// ============================================================================
// 20.7 VALID DOMAIN CONTINUOUS STATISTICS
// ============================================================================

var exp2ValidContinuousStats =
  exp2AuditVariables

    .updateMask(
      exp2ValidMask
    )

    .reduceRegion({

      reducer:
        exp2AuditReducer,

      geometry:
        region,

      scale:
        ANALYSIS_SCALE,

      crs:
        ANALYSIS_CRS,

      maxPixels:
        1e13,

      tileScale:
        8
    });


// ============================================================================
// 20.8 EXCLUDED DOMAIN CONTINUOUS STATISTICS
// ============================================================================

var exp2ExcludedContinuousStats =
  exp2AuditVariables

    .updateMask(
      exp2ExcludedMask
    )

    .reduceRegion({

      reducer:
        exp2AuditReducer,

      geometry:
        region,

      scale:
        ANALYSIS_SCALE,

      crs:
        ANALYSIS_CRS,

      maxPixels:
        1e13,

      tileScale:
        8
    });


// ============================================================================
// 20.9 STANDARDISED MEAN DIFFERENCE
//
// SMD =
// (mean_valid - mean_excluded)
// ---------------------------------------------
// sqrt((SD_valid^2 + SD_excluded^2) / 2)
//
// Used as a descriptive effect-size diagnostic only.
// ============================================================================

function exp2SMD(
  validDictionary,
  excludedDictionary,
  bandName
) {

  validDictionary =
    ee.Dictionary(
      validDictionary
    );

  excludedDictionary =
    ee.Dictionary(
      excludedDictionary
    );


  var validMean =
    ee.Number(
      validDictionary.get(
        bandName + '_mean'
      )
    );


  var excludedMean =
    ee.Number(
      excludedDictionary.get(
        bandName + '_mean'
      )
    );


  var validSD =
    ee.Number(
      validDictionary.get(
        bandName + '_stdDev'
      )
    );


  var excludedSD =
    ee.Number(
      excludedDictionary.get(
        bandName + '_stdDev'
      )
    );


  var pooledSD =
    validSD.pow(2)

      .add(
        excludedSD.pow(2)
      )

      .divide(2)

      .sqrt();


  return safeDivide(

    validMean.subtract(
      excludedMean
    ),

    pooledSD
  );
}


var exp2NDVISMD =
  exp2SMD(

    exp2ValidContinuousStats,
    exp2ExcludedContinuousStats,

    'NDVI_2024'
  );


var exp2SlopeSMD =
  exp2SMD(

    exp2ValidContinuousStats,
    exp2ExcludedContinuousStats,

    'Slope_deg'
  );

// ============================================================================
// 20.10 REGIONAL CONTINUOUS AUDIT SUMMARY
//
// Reporting-area convention:
//
// 1. Region_geometry_area_ha
//    = exact vector geometry area.
//
// 2. Valid_area_ha
//    = common valid RUSLE area from the principal analysis.
//
// 3. Excluded_area_ha
//    = exact region geometry area - common valid RUSLE area.
//
// 4. Raster valid/excluded areas are retained only as QA fields.
//
// Continuous NDVI and slope statistics themselves are calculated using
// the raster valid/excluded masks defined in Sections 20.7-20.8.
// ============================================================================

var exp2ContinuousAuditSummary =
  ee.Feature(null, {

    Experiment:
      'Modelling-domain composition audit',

    Spatial_unit:
      'Minsk Region',

    Analysis_year:
      2024,

    Analysis_scale_m:
      ANALYSIS_SCALE,

    Analysis_CRS:
      ANALYSIS_CRS,


    // ========================================================================
    // AREA REPORTING
    // ========================================================================

    Region_geometry_area_ha:
      exp2RegionGeometryAreaHa,

    Valid_area_ha:
      exp2ValidExactAreaHa,

    Excluded_area_ha:
      exp2ExcludedExactAreaHa,

    Valid_area_pct_region:
      safePercent(
        exp2ValidExactAreaHa,
        exp2RegionGeometryAreaHa
      ),

    Excluded_area_pct_region:
      safePercent(
        exp2ExcludedExactAreaHa,
        exp2RegionGeometryAreaHa
      ),


    // ------------------------------------------------------------------------
    // Raster-area QA only
    // ------------------------------------------------------------------------

    Raster_valid_area_QA_ha:
      exp2ValidRasterAreaHa,

    Raster_excluded_area_QA_ha:
      exp2ExcludedRasterAreaHa,

    Raster_valid_minus_main_valid_ha:
      exp2ValidRasterAreaHa
        .subtract(
          exp2ValidExactAreaHa
        ),


    // ========================================================================
    // NDVI
    // ========================================================================

    Valid_NDVI_mean:
      exp2ValidContinuousStats.get(
        'NDVI_2024_mean'
      ),

    Valid_NDVI_median:
      exp2ValidContinuousStats.get(
        'NDVI_2024_median'
      ),

    Valid_NDVI_stdDev:
      exp2ValidContinuousStats.get(
        'NDVI_2024_stdDev'
      ),

    Valid_NDVI_P25:
      exp2ValidContinuousStats.get(
        'NDVI_2024_p25'
      ),

    Valid_NDVI_P75:
      exp2ValidContinuousStats.get(
        'NDVI_2024_p75'
      ),


    Excluded_NDVI_mean:
      exp2ExcludedContinuousStats.get(
        'NDVI_2024_mean'
      ),

    Excluded_NDVI_median:
      exp2ExcludedContinuousStats.get(
        'NDVI_2024_median'
      ),

    Excluded_NDVI_stdDev:
      exp2ExcludedContinuousStats.get(
        'NDVI_2024_stdDev'
      ),

    Excluded_NDVI_P25:
      exp2ExcludedContinuousStats.get(
        'NDVI_2024_p25'
      ),

    Excluded_NDVI_P75:
      exp2ExcludedContinuousStats.get(
        'NDVI_2024_p75'
      ),


    NDVI_standardised_difference:
      exp2NDVISMD,


    // ========================================================================
    // SLOPE
    // ========================================================================

    Valid_slope_mean_deg:
      exp2ValidContinuousStats.get(
        'Slope_deg_mean'
      ),

    Valid_slope_median_deg:
      exp2ValidContinuousStats.get(
        'Slope_deg_median'
      ),

    Valid_slope_stdDev_deg:
      exp2ValidContinuousStats.get(
        'Slope_deg_stdDev'
      ),

    Valid_slope_P25_deg:
      exp2ValidContinuousStats.get(
        'Slope_deg_p25'
      ),

    Valid_slope_P75_deg:
      exp2ValidContinuousStats.get(
        'Slope_deg_p75'
      ),


    Excluded_slope_mean_deg:
      exp2ExcludedContinuousStats.get(
        'Slope_deg_mean'
      ),

    Excluded_slope_median_deg:
      exp2ExcludedContinuousStats.get(
        'Slope_deg_median'
      ),

    Excluded_slope_stdDev_deg:
      exp2ExcludedContinuousStats.get(
        'Slope_deg_stdDev'
      ),

    Excluded_slope_P25_deg:
      exp2ExcludedContinuousStats.get(
        'Slope_deg_p25'
      ),

    Excluded_slope_P75_deg:
      exp2ExcludedContinuousStats.get(
        'Slope_deg_p75'
      ),


    Slope_standardised_difference:
      exp2SlopeSMD
  });


// IMPORTANT:
// Do not force this large reduction in the interactive Console.
// Use batch Export below.

// print(
//   'EXP2 Continuous domain audit',
//   exp2ContinuousAuditSummary
// );




// ============================================================================
// 20.11 LAND-COVER COMPOSITION
//
// Purpose:
//
// Compare land-cover composition between:
//
// 1 = common valid RUSLE domain
// 2 = excluded part of Minsk Region
//
// No semantic class labels are invented here.
// Native landcover_code values are retained.
//
// Percent composition is calculated relative to the land-cover-covered
// raster area of the corresponding domain.
// ============================================================================


// ============================================================================
// 20.11.1 GROUP LAND-COVER AREA BY NATIVE CLASS CODE
// ============================================================================

function exp2LandcoverAreaGroups(
  domainMask
) {

  var groupedResult =
    expPixelAreaHa

      .addBands(
        exp2Landcover
      )

      .updateMask(
        domainMask
      )

      .reduceRegion({

        reducer:
          ee.Reducer.sum()

            .group({

              groupField:
                1,

              groupName:
                'landcover_code'
            }),

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


  return ee.List(
    ee.Dictionary(
      groupedResult
    ).get(
      'groups'
    )
  );
}




// ============================================================================
// 20.11.2 CONVERT GROUPED OUTPUT TO DICTIONARY
//
// Output:
//
// "101" : area_ha
// "..." : area_ha
// ============================================================================

function exp2GroupsToAreaDictionary(
  groups
) {

  groups =
    ee.List(
      groups
    );


  var codes =
    groups.map(
      function(item) {

        item =
          ee.Dictionary(
            item
          );

        return ee.Number(
          item.get(
            'landcover_code'
          )
        )

        .format(
          '%.0f'
        );
      }
    );


  var areas =
    groups.map(
      function(item) {

        item =
          ee.Dictionary(
            item
          );

        return ee.Number(
          item.get(
            'sum'
          )
        );
      }
    );


  return ee.Dictionary.fromLists(
    codes,
    areas
  );
}




// ============================================================================
// 20.11.3 VALID AND EXCLUDED LAND-COVER GROUPS
// ============================================================================

var exp2ValidLCGroups =
  exp2LandcoverAreaGroups(
    exp2ValidMask
  );


var exp2ExcludedLCGroups =
  exp2LandcoverAreaGroups(
    exp2ExcludedMask
  );


var exp2ValidLCDictionary =
  exp2GroupsToAreaDictionary(
    exp2ValidLCGroups
  );


var exp2ExcludedLCDictionary =
  exp2GroupsToAreaDictionary(
    exp2ExcludedLCGroups
  );




// ============================================================================
// 20.11.4 ALL LAND-COVER CODES OCCURRING IN EITHER DOMAIN
// ============================================================================

var exp2AllLCCodes =
  exp2ValidLCDictionary
    .keys()

    .cat(
      exp2ExcludedLCDictionary
        .keys()
    )

    .distinct()

    .sort();




// ============================================================================
// 20.11.5 SUM DICTIONARY VALUES SAFELY
// ============================================================================

function exp2DictionaryValueSum(
  dictionary
) {

  dictionary =
    ee.Dictionary(
      dictionary
    );


  var values =
    ee.List(
      dictionary.values()
    );


  return ee.Number(

    ee.Algorithms.If(

      values.size().gt(0),

      values.reduce(
        ee.Reducer.sum()
      ),

      0
    )
  );
}




// Land-cover-covered raster area.

var exp2ValidLCTotalHa =
  exp2DictionaryValueSum(
    exp2ValidLCDictionary
  );


var exp2ExcludedLCTotalHa =
  exp2DictionaryValueSum(
    exp2ExcludedLCDictionary
  );




// ============================================================================
// 20.11.6 HELPER:
// RETURN CLASS AREA OR ZERO
// ============================================================================

function exp2DictValueOrZero(
  dictionary,
  key
) {

  dictionary =
    ee.Dictionary(
      dictionary
    );


  return ee.Number(

    ee.Algorithms.If(

      dictionary.contains(
        key
      ),

      dictionary.get(
        key
      ),

      0
    )
  );
}




// ============================================================================
// 20.11.7 BUILD LAND-COVER COMPOSITION TABLE
// ============================================================================

var exp2LandcoverCompositionTable =
  ee.FeatureCollection(

    exp2AllLCCodes.map(
      function(code) {

        code =
          ee.String(
            code
          );


        var validClassHa =
          exp2DictValueOrZero(

            exp2ValidLCDictionary,

            code
          );


        var excludedClassHa =
          exp2DictValueOrZero(

            exp2ExcludedLCDictionary,

            code
          );


        var validPct =
          safePercent(

            validClassHa,

            exp2ValidLCTotalHa
          );


        var excludedPct =
          safePercent(

            excludedClassHa,

            exp2ExcludedLCTotalHa
          );


        return ee.Feature(null, {

          Experiment:
            'Modelling-domain land-cover composition',

          Spatial_unit:
            'Minsk Region',

          Analysis_year:
            2024,

          Landcover_code:
            code,


          // ------------------------------------------------------------
          // VALID DOMAIN
          // ------------------------------------------------------------

          Valid_area_ha:
            validClassHa,

          Valid_pct_landcover_domain:
            validPct,


          // ------------------------------------------------------------
          // EXCLUDED DOMAIN
          // ------------------------------------------------------------

          Excluded_area_ha:
            excludedClassHa,

          Excluded_pct_landcover_domain:
            excludedPct,


          // ------------------------------------------------------------
          // COMPOSITION DIFFERENCE
          // ------------------------------------------------------------

          Difference_valid_minus_excluded_pp:
            validPct.subtract(
              excludedPct
            ),

          Absolute_difference_pp:
            validPct
              .subtract(
                excludedPct
              )
              .abs(),

          Valid_to_excluded_share_ratio:
            safeDivide(
              validPct,
              excludedPct
            ),


          // ------------------------------------------------------------
          // LAND-COVER RASTER COVERAGE QA
          // ------------------------------------------------------------

          Valid_LC_total_ha:
            exp2ValidLCTotalHa,

          Excluded_LC_total_ha:
            exp2ExcludedLCTotalHa,


          Valid_LC_coverage_pct_of_valid_raster_domain:
            safePercent(

              exp2ValidLCTotalHa,

              exp2ValidRasterAreaHa
            ),


          Excluded_LC_coverage_pct_of_excluded_raster_domain:
            safePercent(

              exp2ExcludedLCTotalHa,

              exp2ExcludedRasterAreaHa
            )
        });
      }
    )
  );


// Do not expand the full table interactively.
// Export below.

// print(
//   'EXP2 Land-cover composition',
//   exp2LandcoverCompositionTable
// );




// ============================================================================
// 20.12 DISTRICT COVERAGE ASSOCIATIONS
//
// Main question:
//
// Is variation in the percentage of each district included in the common
// RUSLE modelling domain associated with the district erosion metrics?
//
// Principal diagnostic:
// Spearman rank correlation.
//
// GEE p-values are deliberately not reported because previous runs returned
// NaN. Statistical inference can be carried out outside GEE if required.
// ============================================================================

function exp2CoverageAssociation(
  metricField,
  scenario,
  metricLabel
) {

  var result =
    districtComparison

      .reduceColumns({

        reducer:
          ee.Reducer
            .spearmansCorrelation(),

        selectors: [

          'Valid_model_coverage_pct',

          metricField
        ]
      });


  var rho =
    ee.Number(
      result.get(
        'correlation'
      )
    );


  return ee.Feature(null, {

    Experiment:
      'District coverage association',

    Coverage_variable:
      'Valid_model_coverage_pct',

    Metric:
      metricLabel,

    Metric_field:
      metricField,

    Scenario:
      scenario,

    Number_of_districts:
      districtComparison.size(),

    Spearman_rho:
      rho,

    Absolute_Spearman_rho:
      rho.abs()
  });
}




// ============================================================================
// 20.12.1 COVERAGE vs DISTRICT EROSION METRICS
// ============================================================================

var exp2CoverageAssociationTable =
  ee.FeatureCollection([


    // ------------------------------------------------------------------------
    // District mean A
    // ------------------------------------------------------------------------

    exp2CoverageAssociation(

      'A_VDK_mean',

      'VDK',

      'District mean A'
    ),


    exp2CoverageAssociation(

      'A_LINEAR_mean',

      'Linear',

      'District mean A'
    ),


    // ------------------------------------------------------------------------
    // Area A >= 20
    // ------------------------------------------------------------------------

    exp2CoverageAssociation(

      'T20_VDK_pct_valid',

      'VDK',

      'Area percentage A >= 20'
    ),


    exp2CoverageAssociation(

      'T20_LINEAR_pct_valid',

      'Linear',

      'Area percentage A >= 20'
    ),


    // ------------------------------------------------------------------------
    // Area A >= 50
    // ------------------------------------------------------------------------

    exp2CoverageAssociation(

      'T50_VDK_pct_valid',

      'VDK',

      'Area percentage A >= 50'
    ),


    exp2CoverageAssociation(

      'T50_LINEAR_pct_valid',

      'Linear',

      'Area percentage A >= 50'
    )
  ]);


// Avoid interactive full-table evaluation.

// print(
//   'EXP2 Coverage vs district erosion metrics',
//   exp2CoverageAssociationTable
// );




// ============================================================================
// 20.12B DISTRICT COVERAGE vs SCENARIO RANK INSTABILITY
//
// NEW ROBUSTNESS DIAGNOSTIC
//
// Question:
//
// Does uneven valid-domain coverage explain VDK-Linear district-ranking
// instability?
//
// Rank instability = absolute difference between the VDK and Linear ranks.
//
// Three ranking bases:
//
// 1. District mean A
// 2. Percentage of valid domain with A >= 20
// 3. Percentage of valid domain with A >= 50
//
// A weak rho would indicate that uneven coverage is not a major driver of
// scenario-dependent district rank instability.
//
// A stronger rho would indicate that interpretation of ranking disagreement
// should explicitly account for mapped-domain coverage.
// ============================================================================

function exp2CoverageVsRankInstability(
  rankTable,
  rankBasis
) {

  rankTable =
    ee.FeatureCollection(
      rankTable
    );


  var result =
    rankTable

      .reduceColumns({

        reducer:
          ee.Reducer
            .spearmansCorrelation(),

        selectors: [

          'Valid_model_coverage_pct',

          'Absolute_rank_change'
        ]
      });


  var rho =
    ee.Number(
      result.get(
        'correlation'
      )
    );


  return ee.Feature(null, {

    Experiment:
      'Coverage versus scenario rank instability',

    Rank_basis:
      rankBasis,

    Coverage_variable:
      'Valid_model_coverage_pct',

    Instability_variable:
      'Absolute_rank_change',

    Number_of_districts:
      rankTable.size(),

    Spearman_rho:
      rho,

    Absolute_Spearman_rho:
      rho.abs()
  });
}




// ============================================================================
// 20.12B.1 FINAL RANK-INSTABILITY TABLE
// ============================================================================

var exp2CoverageRankInstabilityTable =
  ee.FeatureCollection([


    exp2CoverageVsRankInstability(

      meanRankResult.table,

      'District mean A'
    ),


    exp2CoverageVsRankInstability(

      t20RankResult.table,

      'Area percentage A >= 20'
    ),


    exp2CoverageVsRankInstability(

      t50RankResult.table,

      'Area percentage A >= 50'
    )
  ]);


// Avoid interactive evaluation.

// print(
//   'EXP2 Coverage vs rank instability',
//   exp2CoverageRankInstabilityTable
// );




// ============================================================================
// 20.13 DOMAIN-COMPOSITION RASTER
//
// Raster code:
//
// 1 = common valid modelling domain
// 2 = excluded part of Minsk Region
//
// Intended for:
// QA / supplementary mapping / figure preparation.
// ============================================================================

var exp2DomainRaster =
  ee.Image.constant(0)

    .where(
      exp2ValidMask,
      1
    )

    .where(
      exp2ExcludedMask,
      2
    )

    .updateMask(
      exp2RegionMask
    )

    .rename(
      'Domain_type'
    )

    .toInt8();




// Optional display only.

Map.addLayer(

  exp2DomainRaster,

  {
    min: 1,
    max: 2,

    palette: [
      '2166AC',
      'D9D9D9'
    ]
  },

  'EXP2 Valid vs excluded domain',

  false
);




// ============================================================================
// 20.14 EXPERIMENT 2 EXPORTS
//
// Recommended:
// Run CSV exports as batch tasks.
//
// Do not rely on large interactive print() requests.
// ============================================================================


// ============================================================================
// 20.14.1 AREA QA
// ============================================================================

Export.table.toDrive({

  collection:
    ee.FeatureCollection([
      exp2AreaQA
    ]),

  description:
    'EXP2_RUSLE_Domain_Area_QA_2024',

  fileFormat:
    'CSV'
});




// ============================================================================
// 20.14.2 NDVI / SLOPE CONTINUOUS DOMAIN AUDIT
// ============================================================================

Export.table.toDrive({

  collection:
    ee.FeatureCollection([
      exp2ContinuousAuditSummary
    ]),

  description:
    'EXP2_RUSLE_Domain_NDVI_Slope_Audit_2024',

  fileFormat:
    'CSV'
});




// ============================================================================
// 20.14.3 LAND-COVER COMPOSITION
// ============================================================================

Export.table.toDrive({

  collection:
    exp2LandcoverCompositionTable,

  description:
    'EXP2_RUSLE_Domain_Landcover_Composition_2024',

  fileFormat:
    'CSV'
});




// ============================================================================
// 20.14.4 COVERAGE vs DISTRICT EROSION METRICS
// ============================================================================

Export.table.toDrive({

  collection:
    exp2CoverageAssociationTable,

  description:
    'EXP2_RUSLE_DistrictCoverage_Associations_2024',

  fileFormat:
    'CSV'
});




// ============================================================================
// 20.14.5 COVERAGE vs DISTRICT RANK INSTABILITY
//
// NEW
// ============================================================================

Export.table.toDrive({

  collection:
    exp2CoverageRankInstabilityTable,

  description:
    'EXP2B_RUSLE_Coverage_vs_RankInstability_2024',

  fileFormat:
    'CSV'
});




// ============================================================================
// 20.14.6 OPTIONAL DOMAIN RASTER EXPORT
//
// This raster is NOT required for the numerical manuscript results.
//
// Run only if you want a supplementary map or permanent QA asset.
// ============================================================================

Export.image.toAsset({

  image:
    exp2DomainRaster,

  description:
    'Export_EXP2_RUSLE_Valid_Excluded_Domain_2024',

  assetId:
    'projects/YOUR_EE_PROJECT/assets/' +
    'EXP2_RUSLE_Valid_Excluded_Domain_2024',

  region:
    region,

  scale:
    ANALYSIS_SCALE,

  crs:
    ANALYSIS_CRS,

  maxPixels:
    1e13
});

// ============================================================================
// 20.15 FINAL EXPERIMENT 2 CHECK
// ============================================================================

print(
  'EXP2 prepared:',
  'Use batch exports for domain audit, land-cover composition, ' +
  'coverage associations and coverage-vs-rank-instability.'
);


// ============================================================================
// END OF EXPERIMENT 2
// ============================================================================



// ============================================================================
// 23. FIGURE 3 - OFFICIAL FIVE-CLASS GEOTIFF EXPORTS
// ============================================================================
//
// Purpose:
// Export publication-ready classified rasters for Figure 3.
//
// Official erosion classes:
//   1 = A < 5
//   2 = 5 <= A < 10
//   3 = 10 <= A < 20
//   4 = 20 <= A < 50
//   5 = A >= 50
//
// Class difference:
//   ErosionClass_LINEAR - ErosionClass_VDK
//
//   negative = VDK assigned a higher class
//   0        = same class
//   positive = Linear assigned a higher class
//
// All rasters:
//   CRS   = EPSG:32635
//   Scale = 30 m
//   NoData = -9999
//
// ============================================================================


// ----------------------------------------------------------------------------
// 23.1 SETTINGS
// ----------------------------------------------------------------------------

var FIG3_CRS = 'EPSG:32635';

var FIG3_SCALE = 30;

// GeoTIFF NoData value.
// Int16 is used because -9999 cannot be stored safely in Int8.
var FIG3_NODATA = -9999;

var FIG3_DRIVE_FOLDER = 'RUSLE_Figure3_2024';


// ----------------------------------------------------------------------------
// 23.2 PREPARE FINAL CLASSIFIED RASTERS
// ----------------------------------------------------------------------------
//
// IMPORTANT:
// We use the already calculated classVDK and classLinear images.
// Therefore Figure 3 uses exactly the same classifications as all
// statistics, transition matrices, Jaccard/Dice calculations, etc.
// ----------------------------------------------------------------------------


// VDK official five-class raster.
var fig3VDK = classVDK
  .rename('VDK_RUSLE_Class')
  .toInt16()
  .clip(region);


// Linear official five-class raster.
var fig3Linear = classLinear
  .rename('LINEAR_RUSLE_Class')
  .toInt16()
  .clip(region);


// Difference in ordinal erosion class:
//
// positive: Linear > VDK
// negative: VDK > Linear
// zero:     same class
//
var fig3ClassDifference = classLinear
  .subtract(classVDK)
  .rename('Class_Difference_LINEAR_minus_VDK')
  .toInt16()
  .clip(region);


// ----------------------------------------------------------------------------
// 23.3 SET EXPLICIT NODATA
// ----------------------------------------------------------------------------
//
// Masked pixels are converted to -9999.
// The GeoTIFF metadata also records -9999 as NoData.
//
// This is useful in ArcGIS/QGIS because areas outside the common valid
// modelling domain will automatically be treated as NoData rather than 0.
// ----------------------------------------------------------------------------

var fig3VDKExport = fig3VDK.unmask({
  value: FIG3_NODATA,
  sameFootprint: false
});

var fig3LinearExport = fig3Linear.unmask({
  value: FIG3_NODATA,
  sameFootprint: false
});

var fig3DifferenceExport = fig3ClassDifference.unmask({
  value: FIG3_NODATA,
  sameFootprint: false
});


// ----------------------------------------------------------------------------
// 23.4 EXPORT VDK FIVE-CLASS GEOTIFF
// ----------------------------------------------------------------------------

Export.image.toDrive({

  image: fig3VDKExport,

  description:
    'FIG3A_RUSLE_VDK_Official5Class_2024',

  folder:
    FIG3_DRIVE_FOLDER,

  fileNamePrefix:
    'FIG3A_RUSLE_VDK_Official5Class_2024',

  region:
    region,

  scale:
    FIG3_SCALE,

  crs:
    FIG3_CRS,

  maxPixels:
    1e13,

  fileFormat:
    'GeoTIFF',

  formatOptions: {
    cloudOptimized: true,
    noData: FIG3_NODATA
  }

});


// ----------------------------------------------------------------------------
// 23.5 EXPORT LINEAR FIVE-CLASS GEOTIFF
// ----------------------------------------------------------------------------

Export.image.toDrive({

  image: fig3LinearExport,

  description:
    'FIG3B_RUSLE_LINEAR_Official5Class_2024',

  folder:
    FIG3_DRIVE_FOLDER,

  fileNamePrefix:
    'FIG3B_RUSLE_LINEAR_Official5Class_2024',

  region:
    region,

  scale:
    FIG3_SCALE,

  crs:
    FIG3_CRS,

  maxPixels:
    1e13,

  fileFormat:
    'GeoTIFF',

  formatOptions: {
    cloudOptimized: true,
    noData: FIG3_NODATA
  }

});


// ----------------------------------------------------------------------------
// 23.6 EXPORT CLASS-DIFFERENCE GEOTIFF
// ----------------------------------------------------------------------------

Export.image.toDrive({

  image: fig3DifferenceExport,

  description:
    'FIG3C_RUSLE_ClassDifference_LINEAR_minus_VDK_2024',

  folder:
    FIG3_DRIVE_FOLDER,

  fileNamePrefix:
    'FIG3C_RUSLE_ClassDifference_LINEAR_minus_VDK_2024',

  region:
    region,

  scale:
    FIG3_SCALE,

  crs:
    FIG3_CRS,

  maxPixels:
    1e13,

  fileFormat:
    'GeoTIFF',

  formatOptions: {
    cloudOptimized: true,
    noData: FIG3_NODATA
  }

});


// ----------------------------------------------------------------------------
// 23.7 OPTIONAL: EXPORT A SINGLE THREE-BAND STACK
// ----------------------------------------------------------------------------
//
// This is optional, but useful for archiving and reproducibility.
// For ArcGIS figure preparation I recommend using the three separate TIFFs
// above because symbology is easier to control.
//
// Bands:
//   1. VDK_RUSLE_Class
//   2. LINEAR_RUSLE_Class
//   3. Class_Difference_LINEAR_minus_VDK
// ----------------------------------------------------------------------------

var fig3Stack = fig3VDK
  .addBands(fig3Linear)
  .addBands(fig3ClassDifference)
  .toInt16()
  .unmask({
    value: FIG3_NODATA,
    sameFootprint: false
  });

Export.image.toDrive({

  image:
    fig3Stack,

  description:
    'FIG3_RUSLE_Official5Class_ThreeBandStack_2024',

  folder:
    FIG3_DRIVE_FOLDER,

  fileNamePrefix:
    'FIG3_RUSLE_Official5Class_ThreeBandStack_2024',

  region:
    region,

  scale:
    FIG3_SCALE,

  crs:
    FIG3_CRS,

  maxPixels:
    1e13,

  fileFormat:
    'GeoTIFF',

  formatOptions: {
    cloudOptimized: true,
    noData: FIG3_NODATA
  }

});


// ----------------------------------------------------------------------------
// 23.8 MAP PREVIEW
// ----------------------------------------------------------------------------

// Suggested five-class preview palette.
// This is ONLY for GEE preview.
// GeoTIFF exports contain integer class codes, not embedded symbology.

var FIG3_CLASS_PALETTE = [
  '2C7BB6',   // Class 1
  'ABD9E9',   // Class 2
  'FFFFBF',   // Class 3
  'FDAE61',   // Class 4
  'D7191C'    // Class 5
];


// Diverging class-difference palette:
// -4 ... 0 ... +4
var FIG3_DIFF_PALETTE = [
  '053061',
  '2166AC',
  '4393C3',
  '92C5DE',
  'F7F7F7',
  'F4A582',
  'D6604D',
  'B2182B',
  '67001F'
];


Map.addLayer(
  fig3VDK,
  {
    min: 1,
    max: 5,
    palette: FIG3_CLASS_PALETTE
  },
  'FIG3A - VDK five classes',
  true
);


Map.addLayer(
  fig3Linear,
  {
    min: 1,
    max: 5,
    palette: FIG3_CLASS_PALETTE
  },
  'FIG3B - Linear five classes',
  true
);


Map.addLayer(
  fig3ClassDifference,
  {
    min: -4,
    max: 4,
    palette: FIG3_DIFF_PALETTE
  },
  'FIG3C - Linear minus VDK class',
  true
);


// District boundaries for visual QA.
Map.addLayer(
  districts.style({
    color: '000000',
    fillColor: '00000000',
    width: 1
  }),
  {},
  'District boundaries - Figure 3 QA',
  true
);


// ----------------------------------------------------------------------------
// 23.9 FINAL QA
// ----------------------------------------------------------------------------

print(
  'Figure 3 class definitions',
  ee.Dictionary({
    Class_1: 'A < 5',
    Class_2: '5 <= A < 10',
    Class_3: '10 <= A < 20',
    Class_4: '20 <= A < 50',
    Class_5: 'A >= 50'
  })
);

print(
  'Figure 3 difference interpretation',
  ee.Dictionary({
    Negative: 'VDK class > Linear class',
    Zero: 'Same erosion class',
    Positive: 'Linear class > VDK class'
  })
);

print(
  'Figure 3 export settings',
  ee.Dictionary({
    CRS: FIG3_CRS,
    Scale_m: FIG3_SCALE,
    NoData: FIG3_NODATA
  })
);
