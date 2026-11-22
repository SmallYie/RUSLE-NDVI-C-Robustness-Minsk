// ============================================================================
// 03_exp3_extra_district_rankings_si2024.js
// ============================================================================
// Standalone modular Google Earth Engine script for EXP3 supplementary
// district-level ranking robustness across 22 Minsk Region districts:
//   - District median RUSLE (A_VDK_median vs. A_LINEAR_median)
//   - Area percentage A >= 5 t ha^-1 yr^-1  (T5)
//   - Area percentage A >= 10 t ha^-1 yr^-1 (T10)
//
// Outputs Spearman's rho, mean/max absolute rank shift, and stable vs.
// scenario-sensitive Top-6 priority districts.
//
// Manuscript: Scale-dependent propagation of NDVI-C parameterisation
//             uncertainty in regional RUSLE soil-erosion modelling:
//             evidence from Minsk Region, Belarus (Earth Science Informatics)
// Authors:    Bochao Zhao, Aliaksandr N. Chervan, Wei Zhang, Dongliang Zhong
// License:    MIT License
// ============================================================================

var DISTRICTS_ASSET = 'projects/YOUR_EE_PROJECT/assets/MINSK_DISTRICTS';
var DISTRICT_NAME_FIELD = 'F09';
var RUSLE_VDK_ASSET = 'projects/YOUR_EE_PROJECT/assets/RUSLE_VDK_SI_2024';
var RUSLE_LIN_ASSET = 'projects/YOUR_EE_PROJECT/assets/RUSLE_LINEAR_SI_2024';
var DRIVE_FOLDER = 'RUSLE_SI_2024_EXP3';

var ANALYSIS_CRS = 'EPSG:32635';
var ANALYSIS_SCALE = 30;
var TOP_K = 6;

var districts = ee.FeatureCollection(DISTRICTS_ASSET);
var region = districts.geometry();

var rusleVDK0 = ee.Image(RUSLE_VDK_ASSET).select(0).toFloat().rename('A_VDK');
var rusleLin0 = ee.Image(RUSLE_LIN_ASSET).select(0).toFloat().rename('A_LINEAR');

var commonMask = rusleVDK0.mask().and(rusleLin0.mask())
  .and(rusleVDK0.gte(0)).and(rusleLin0.gte(0));

var rusleVDK = rusleVDK0.updateMask(commonMask).clip(region);
var rusleLin = rusleLin0.updateMask(commonMask).clip(region);
var areaHa = ee.Image.pixelArea().divide(10000).updateMask(commonMask).clip(region);

function safePct(num, den) {
  num = ee.Number(num);
  den = ee.Number(den);
  return ee.Number(ee.Algorithms.If(den.neq(0), num.divide(den).multiply(100), 0));
}

// Compute per-district median, T5%, and T10%
var districtMetrics = districts.map(function(f) {
  var g = f.geometry();
  var med = rusleVDK.addBands(rusleLin).reduceRegion({
    reducer: ee.Reducer.median(),
    geometry: g,
    crs: ANALYSIS_CRS,
    scale: ANALYSIS_SCALE,
    maxPixels: 1e13,
    tileScale: 4
  });

  var areas = ee.Image.cat([
    areaHa.rename('valid_ha'),
    areaHa.updateMask(rusleVDK.gte(5)).rename('vdk_t5_ha'),
    areaHa.updateMask(rusleLin.gte(5)).rename('lin_t5_ha'),
    areaHa.updateMask(rusleVDK.gte(10)).rename('vdk_t10_ha'),
    areaHa.updateMask(rusleLin.gte(10)).rename('lin_t10_ha')
  ]).reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: g,
    crs: ANALYSIS_CRS,
    scale: ANALYSIS_SCALE,
    maxPixels: 1e13,
    tileScale: 4
  });

  var valid = ee.Number(areas.get('valid_ha', 0));
  return f.set({
    unit_name: f.get(DISTRICT_NAME_FIELD),
    A_VDK_median: ee.Number(med.get('A_VDK', 0)),
    A_LINEAR_median: ee.Number(med.get('A_LINEAR', 0)),
    T5_VDK_pct_valid: safePct(areas.get('vdk_t5_ha', 0), valid),
    T5_LINEAR_pct_valid: safePct(areas.get('lin_t5_ha', 0), valid),
    T10_VDK_pct_valid: safePct(areas.get('vdk_t10_ha', 0), valid),
    T10_LINEAR_pct_valid: safePct(areas.get('lin_t10_ha', 0), valid)
  });
});

function assignDescendingRanks(fc, propName, rankName) {
  var sorted = fc.sort(propName, false).toList(fc.size());
  var n = sorted.size();
  return ee.FeatureCollection(ee.List.sequence(0, n.subtract(1)).map(function(i) {
    i = ee.Number(i);
    return ee.Feature(sorted.get(i)).set(rankName, i.add(1));
  }));
}

function buildRankingTable(fc, vdkField, linField, label) {
  var r1 = assignDescendingRanks(fc, vdkField, 'VDK_rank');
  var r2 = assignDescendingRanks(r1, linField, 'LINEAR_rank');
  return r2.map(function(f) {
    var vr = ee.Number(f.get('VDK_rank'));
    var lr = ee.Number(f.get('LINEAR_rank'));
    var diff = lr.subtract(vr);
    var inV = vr.lte(TOP_K);
    var inL = lr.lte(TOP_K);
    return f.set({
      Rank_basis: label,
      VDK_value: f.get(vdkField),
      LINEAR_value: f.get(linField),
      Rank_change_LINEAR_minus_VDK: diff,
      Absolute_rank_change: diff.abs(),
      Stable_top6: ee.Algorithms.If(inV.and(inL), 1, 0),
      Scenario_sensitive_top6: ee.Algorithms.If(inV.neq(inL), 1, 0)
    });
  });
}

function summarizeRanking(rankedFC, label, vdkField, linField) {
  var corr = rankedFC.reduceColumns({
    reducer: ee.Reducer.spearmansCorrelation(),
    selectors: [vdkField, linField]
  });
  var stableList = rankedFC.filter(ee.Filter.eq('Stable_top6', 1)).aggregate_array('unit_name');
  var sensList = rankedFC.filter(ee.Filter.eq('Scenario_sensitive_top6', 1)).aggregate_array('unit_name');
  return ee.Feature(null, {
    Rank_basis: label,
    VDK_indicator: vdkField,
    LINEAR_indicator: linField,
    Number_of_districts: rankedFC.size(),
    Spearman_rank_correlation: corr.get('correlation'),
    Spearman_p_value: corr.get('p-value'),
    Mean_absolute_rank_change: rankedFC.aggregate_mean('Absolute_rank_change'),
    Maximum_absolute_rank_change: rankedFC.aggregate_max('Absolute_rank_change'),
    Stable_top6_district_count: stableList.size(),
    Stable_top6_districts: stableList,
    Scenario_sensitive_top6_district_count: sensList.size(),
    Scenario_sensitive_top6_districts: sensList
  });
}

var rankMedian = buildRankingTable(districtMetrics, 'A_VDK_median', 'A_LINEAR_median', 'District median RUSLE');
var rankT5 = buildRankingTable(districtMetrics, 'T5_VDK_pct_valid', 'T5_LINEAR_pct_valid', 'Area percentage A >= 5');
var rankT10 = buildRankingTable(districtMetrics, 'T10_VDK_pct_valid', 'T10_LINEAR_pct_valid', 'Area percentage A >= 10');

var extraSummaries = ee.FeatureCollection([
  summarizeRanking(rankMedian, 'District median RUSLE', 'A_VDK_median', 'A_LINEAR_median'),
  summarizeRanking(rankT5, 'Area percentage A >= 5', 'T5_VDK_pct_valid', 'T5_LINEAR_pct_valid'),
  summarizeRanking(rankT10, 'Area percentage A >= 10', 'T10_VDK_pct_valid', 'T10_LINEAR_pct_valid')
]);

Export.table.toDrive({
  collection: rankMedian,
  description: 'EXP3_RUSLE_District_Median_Ranking_Table_SI_2024',
  folder: DRIVE_FOLDER,
  fileFormat: 'CSV'
});

Export.table.toDrive({
  collection: rankT5,
  description: 'EXP3_RUSLE_District_T5_Ranking_Table_SI_2024',
  folder: DRIVE_FOLDER,
  fileFormat: 'CSV'
});

Export.table.toDrive({
  collection: rankT10,
  description: 'EXP3_RUSLE_District_T10_Ranking_Table_SI_2024',
  folder: DRIVE_FOLDER,
  fileFormat: 'CSV'
});

Export.table.toDrive({
  collection: extraSummaries,
  description: 'EXP3_RUSLE_District_Extra_Ranking_Summaries_SI_2024',
  folder: DRIVE_FOLDER,
  fileFormat: 'CSV'
});

print('EXP3 supplementary district ranking tasks queued.');
