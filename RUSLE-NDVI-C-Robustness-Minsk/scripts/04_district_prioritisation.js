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

Export.table.toDrive({collection: meanRankResult.table, description: 'RUSLE_District_Ranking_Mean_2024', fileFormat: 'CSV'});
Export.table.toDrive({collection: t20RankResult.table, description: 'RUSLE_District_Ranking_T20_2024', fileFormat: 'CSV'});
Export.table.toDrive({collection: t50RankResult.table, description: 'RUSLE_District_Ranking_T50_2024', fileFormat: 'CSV'});
Export.table.toDrive({collection: rankSummaryCollection, description: 'RUSLE_District_Ranking_Summaries_2024', fileFormat: 'CSV'});
