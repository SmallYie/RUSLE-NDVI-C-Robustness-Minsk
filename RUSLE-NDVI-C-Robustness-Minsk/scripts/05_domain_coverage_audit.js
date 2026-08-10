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


// ============================================================================
// DISTRICT PRIORITISATION ROBUSTNESS
// Three ranking criteria: mean A, % area >=20, % area >=50.
// Priority group is explicitly ranks 1-6, matching the manuscript.
// ============================================================================
function districtStats(feature) {
  var geom = feature.geometry();
  var areaHa = geom.area(1).divide(10000);

  var means = rusleVDK.addBands(rusleLinear).reduceRegion({
    reducer: ee.Reducer.mean(), geometry: geom, scale: ANALYSIS_SCALE,
    crs: ANALYSIS_CRS, maxPixels: 1e13, tileScale: 4
  });

  var validHa = numberOrZero(pixelAreaHa.updateMask(commonMask).reduceRegion({
    reducer: ee.Reducer.sum(), geometry: geom, scale: ANALYSIS_SCALE,
    crs: ANALYSIS_CRS, maxPixels: 1e13, tileScale: 4
  }).get('area_ha'));

  function thresholdPct(img, threshold) {
    var ha = numberOrZero(pixelAreaHa.updateMask(img.gte(threshold)).reduceRegion({
      reducer: ee.Reducer.sum(), geometry: geom, scale: ANALYSIS_SCALE,
      crs: ANALYSIS_CRS, maxPixels: 1e13, tileScale: 4
    }).get('area_ha'));
    return safePercent(ha, validHa);
  }

  return feature.set({
    unit_name: feature.get(DISTRICT_NAME_FIELD),
    District_Area_ha: areaHa,
    Valid_model_coverage_pct: safePercent(validHa, areaHa),
    A_VDK_mean: means.get('A_VDK'),
    A_LINEAR_mean: means.get('A_LINEAR'),
    T20_VDK_pct_valid: thresholdPct(rusleVDK, 20),
    T20_LINEAR_pct_valid: thresholdPct(rusleLinear, 20),
    T50_VDK_pct_valid: thresholdPct(rusleVDK, 50),
    T50_LINEAR_pct_valid: thresholdPct(rusleLinear, 50)
  });
}

var districtComparison = districts.map(districtStats);

function addDescendingRank(collection, valueField, rankField) {
  var sorted = collection.sort(valueField, false);
  var n = sorted.size();
  var list = sorted.toList(n);
  return ee.FeatureCollection(ee.List.sequence(0, n.subtract(1)).map(function(i) {
    i = ee.Number(i);
    return ee.Feature(list.get(i)).set(rankField, i.add(1));
  }));
}

function buildRankComparison(collection, vdkField, linearField, rankBasis) {
  var rankedVDK = addDescendingRank(collection, vdkField, 'VDK_rank');
  var rankedLinear = addDescendingRank(collection, linearField, 'LINEAR_rank');
  var joined = ee.Join.inner().apply({
    primary: rankedVDK,
    secondary: rankedLinear,
    condition: ee.Filter.equals({leftField: 'unit_name', rightField: 'unit_name'})
  });

  var top6CutoffRank = ee.Number(6);
  var rankTable = ee.FeatureCollection(joined.map(function(j) {
    j = ee.Feature(j);
    var a = ee.Feature(j.get('primary'));
    var b = ee.Feature(j.get('secondary'));
    var rv = ee.Number(a.get('VDK_rank'));
    var rl = ee.Number(b.get('LINEAR_rank'));
    var vTop = rv.lte(top6CutoffRank);
    var lTop = rl.lte(top6CutoffRank);
    return a.set({
      Rank_basis: rankBasis,
      VDK_value: a.get(vdkField),
      LINEAR_value: b.get(linearField),
      VDK_rank: rv,
      LINEAR_rank: rl,
      Rank_change_LINEAR_minus_VDK: rl.subtract(rv),
      Absolute_rank_change: rl.subtract(rv).abs(),
      Stable_top6: vTop.and(lTop),
      Scenario_sensitive_top6: vTop.neq(lTop)
    });
  }));

  var corr = rankTable.reduceColumns({
    reducer: ee.Reducer.pearsonsCorrelation(),
    selectors: ['VDK_rank','LINEAR_rank']
  });

  var summary = ee.Feature(null, {
    Rank_basis: rankBasis,
    Number_of_districts: rankTable.size(),
    Top6_cutoff_rank: 6,
    Spearman_rank_correlation: corr.get('correlation'),
    Spearman_p_value: corr.get('p-value'),
    Mean_absolute_rank_change: rankTable.aggregate_mean('Absolute_rank_change'),
    Maximum_absolute_rank_change: rankTable.aggregate_max('Absolute_rank_change'),
    Stable_top6_district_count: rankTable.filter(ee.Filter.eq('Stable_top6', 1)).size(),
    Scenario_sensitive_top6_district_count: rankTable.filter(ee.Filter.eq('Scenario_sensitive_top6', 1)).size(),
    Stable_top6_districts: rankTable.filter(ee.Filter.eq('Stable_top6', 1)).aggregate_array('unit_name'),
    Scenario_sensitive_top6_districts: rankTable.filter(ee.Filter.eq('Scenario_sensitive_top6', 1)).aggregate_array('unit_name')
  });
  return {table: rankTable, summary: summary};
}

var meanRankResult = buildRankComparison(districtComparison, 'A_VDK_mean', 'A_LINEAR_mean', 'District mean A');
var t20RankResult = buildRankComparison(districtComparison, 'T20_VDK_pct_valid', 'T20_LINEAR_pct_valid', 'Area percentage A >= 20');
var t50RankResult = buildRankComparison(districtComparison, 'T50_VDK_pct_valid', 'T50_LINEAR_pct_valid', 'Area percentage A >= 50');

var rankSummaryCollection = ee.FeatureCollection([
  meanRankResult.summary, t20RankResult.summary, t50RankResult.summary
]);



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