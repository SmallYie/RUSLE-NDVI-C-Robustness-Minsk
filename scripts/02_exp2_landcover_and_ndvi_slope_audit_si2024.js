// ============================================================================
// 02_exp2_landcover_and_ndvi_slope_audit_si2024.js
// ============================================================================
// Standalone modular Google Earth Engine script for:
//   1. EXP2  - Modelling-domain land-cover composition (valid vs. excluded)
//   2. EXP2C - Land-cover-stratified VDK vs. Linear disagreement diagnostics
//   3. EXP2  - Continuous NDVI (2024) and SRTM Slope domain audit (with SMD)
//
// Manuscript: Scale-dependent propagation of NDVI-C parameterisation
//             uncertainty in regional RUSLE soil-erosion modelling:
//             evidence from Minsk Region, Belarus (Earth Science Informatics)
// Authors:    Bochao Zhao, Aliaksandr N. Chervan, Wei Zhang, Dongliang Zhong
// License:    MIT License
//
// Usage:
//   Replace 'projects/YOUR_EE_PROJECT/assets/...' with your own GEE assets.
//   Run in Google Earth Engine Code Editor and start the 3 export tasks.
// ============================================================================

// ----------------------------------------------------------------------------
// 0. USER CONFIGURATION
// ----------------------------------------------------------------------------
var DISTRICTS_ASSET = 'projects/YOUR_EE_PROJECT/assets/MINSK_DISTRICTS';
var RUSLE_VDK_ASSET = 'projects/YOUR_EE_PROJECT/assets/RUSLE_VDK_SI_2024';
var RUSLE_LIN_ASSET = 'projects/YOUR_EE_PROJECT/assets/RUSLE_LINEAR_SI_2024';
var LANDCOVER_ASSET = 'projects/YOUR_EE_PROJECT/assets/LANDCOVER_CATEGORIES';

// Numeric land-cover code attribute in the vector layer (Belarusian cadastre: 'LANDTYPE').
var LANDCOVER_CODE_FIELD = 'LANDTYPE';
var DRIVE_FOLDER = 'RUSLE_SI_2024_EXP2';

var ANALYSIS_CRS = 'EPSG:32635';
var ANALYSIS_SCALE = 30;
var YEAR = 2024;

// ----------------------------------------------------------------------------
// 1. LOAD REGION, RUSLE SCENARIOS, AND COMMON VALID DOMAIN
// ----------------------------------------------------------------------------
var districts = ee.FeatureCollection(DISTRICTS_ASSET);
var region = districts.geometry();

var rusleVDK0 = ee.Image(RUSLE_VDK_ASSET).select(0).toFloat().rename('A_VDK');
var rusleLin0 = ee.Image(RUSLE_LIN_ASSET).select(0).toFloat().rename('A_LINEAR');

var validMask = rusleVDK0.mask().unmask(0).gt(0)
  .and(rusleLin0.mask().unmask(0).gt(0))
  .and(rusleVDK0.unmask(-9999).gte(0))
  .and(rusleLin0.unmask(-9999).gte(0))
  .clip(region);

var excludedMask = validMask.not().clip(region);

var rusleVDK = rusleVDK0.updateMask(validMask).clip(region);
var rusleLin = rusleLin0.updateMask(validMask).clip(region);
var delta = rusleLin.subtract(rusleVDK).rename('Delta_LINEAR_minus_VDK');
var pixelAreaHa = ee.Image.pixelArea().divide(10000).rename('area_ha').clip(region);

function classify5(img) {
  return ee.Image(1)
    .where(img.gte(5).and(img.lt(10)), 2)
    .where(img.gte(10).and(img.lt(20)), 3)
    .where(img.gte(20).and(img.lt(50)), 4)
    .where(img.gte(50), 5)
    .updateMask(img.mask())
    .toInt16();
}

var classVDK = classify5(rusleVDK);
var classLin = classify5(rusleLin);
var classDiff = classLin.subtract(classVDK);

// ----------------------------------------------------------------------------
// 2. RASTERISE LAND-COVER VECTOR (NUMERIC FIELD: LANDTYPE)
// ----------------------------------------------------------------------------
var lcFC = ee.FeatureCollection(LANDCOVER_ASSET)
  .filter(ee.Filter.notNull([LANDCOVER_CODE_FIELD]));

print('Land-cover valid feature count:', lcFC.size());

var lcImage = lcFC.reduceToImage({
  properties: [LANDCOVER_CODE_FIELD],
  reducer: ee.Reducer.first()
}).rename('landcover_code').toInt16().clip(region);

function safeDiv(num, den) {
  num = ee.Number(num);
  den = ee.Number(den);
  return ee.Number(ee.Algorithms.If(den.neq(0), num.divide(den), 0));
}

function groupedAreaDict(domainMask) {
  var res = pixelAreaHa.addBands(lcImage).updateMask(domainMask).reduceRegion({
    reducer: ee.Reducer.sum().group({groupField: 1, groupName: 'landcover_code'}),
    geometry: region,
    crs: ANALYSIS_CRS,
    scale: ANALYSIS_SCALE,
    maxPixels: 1e13,
    tileScale: 8
  });
  var groups = ee.List(ee.Dictionary(res).get('groups', ee.List([])));
  return ee.Dictionary(groups.iterate(function(item, acc) {
    item = ee.Dictionary(item);
    return ee.Dictionary(acc).set(ee.Number(item.get('landcover_code')).format('%d'), item.get('sum'));
  }, ee.Dictionary({})));
}

var validLcDict = groupedAreaDict(validMask);
var exclLcDict = groupedAreaDict(excludedMask);
var allCodes = validLcDict.keys().cat(exclLcDict.keys()).distinct().sort();

var totalValidLcHa = ee.Number(validLcDict.values().reduce(ee.Reducer.sum()));
var totalExclLcHa = ee.Number(exclLcDict.values().reduce(ee.Reducer.sum()));

// 2A. EXP2 Land-cover composition table
var exp2LcTable = ee.FeatureCollection(allCodes.map(function(codeStr) {
  codeStr = ee.String(codeStr);
  var vHa = ee.Number(validLcDict.get(codeStr, 0));
  var eHa = ee.Number(exclLcDict.get(codeStr, 0));
  var vPct = safeDiv(vHa, totalValidLcHa).multiply(100);
  var ePct = safeDiv(eHa, totalExclLcHa).multiply(100);
  return ee.Feature(null, {
    Landcover_code: ee.Number.parse(codeStr),
    Valid_area_ha: vHa,
    Valid_pct_of_lc: vPct,
    Excluded_area_ha: eHa,
    Excluded_pct_of_lc: ePct,
    Diff_valid_minus_excluded_pp: vPct.subtract(ePct)
  });
}));

// 2B. EXP2C Land-cover-stratified VDK vs. Linear disagreement table
var areaImg = pixelAreaHa;
var stratifiedStack = ee.Image.cat([
  areaImg.rename('valid_ha'),
  rusleVDK.multiply(areaImg).rename('vdk_x_ha'),
  rusleLin.multiply(areaImg).rename('lin_x_ha'),
  delta.multiply(areaImg).rename('delta_x_ha'),
  classDiff.eq(0).multiply(areaImg).rename('same_class_ha'),
  classDiff.gt(0).multiply(areaImg).rename('lin_higher_ha'),
  classDiff.lt(0).multiply(areaImg).rename('vdk_higher_ha'),
  lcImage
]).updateMask(validMask);

var stratGrouped = stratifiedStack.reduceRegion({
  reducer: ee.Reducer.sum().repeat(7).group({groupField: 7, groupName: 'landcover_code'}),
  geometry: region,
  crs: ANALYSIS_CRS,
  scale: ANALYSIS_SCALE,
  maxPixels: 1e13,
  tileScale: 8
});

var exp2cTable = ee.FeatureCollection(
  ee.List(ee.Dictionary(stratGrouped).get('groups', ee.List([]))).map(function(item) {
    item = ee.Dictionary(item);
    var sums = ee.List(item.get('sum'));
    var validHa = ee.Number(sums.get(0));
    var vdkMean = safeDiv(sums.get(1), validHa);
    var linMean = safeDiv(sums.get(2), validHa);
    var deltaMean = safeDiv(sums.get(3), validHa);
    var samePct = safeDiv(sums.get(4), validHa).multiply(100);
    var linHighPct = safeDiv(sums.get(5), validHa).multiply(100);
    var vdkHighPct = safeDiv(sums.get(6), validHa).multiply(100);
    return ee.Feature(null, {
      Landcover_code: item.get('landcover_code'),
      Valid_area_ha: validHa,
      VDK_mean_SI: vdkMean,
      LINEAR_mean_SI: linMean,
      Delta_LINEAR_minus_VDK_mean: deltaMean,
      Ratio_LINEAR_to_VDK_mean: safeDiv(linMean, vdkMean),
      Same_class_pct: samePct,
      LINEAR_higher_class_pct: linHighPct,
      VDK_higher_class_pct: vdkHighPct
    });
  })
);

// ----------------------------------------------------------------------------
// 3. EXP2 CONTINUOUS NDVI (2024) AND SRTM SLOPE AUDIT (SAFE KEY FALLBACK)
// ----------------------------------------------------------------------------
function maskS2Clouds(image) {
  var qa = image.select('QA60');
  var mask = qa.bitwiseAnd(1 << 10).eq(0).and(qa.bitwiseAnd(1 << 11).eq(0));
  return image.updateMask(mask);
}

var ndvi2024 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(region)
  .filterDate(YEAR + '-01-01', (YEAR + 1) + '-01-01')
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
  .map(maskS2Clouds)
  .map(function(img) {
    return img.normalizedDifference(['B8', 'B4']).rename('NDVI_2024');
  })
  .median()
  .toFloat()
  .clip(region);

var slopeDeg = ee.Terrain.slope(ee.Image('USGS/SRTMGL1_003'))
  .rename('Slope_deg')
  .toFloat()
  .clip(region);

var statReducer = ee.Reducer.mean()
  .combine({reducer2: ee.Reducer.stdDev(), sharedInputs: true})
  .combine({reducer2: ee.Reducer.percentile([25, 50, 75]), sharedInputs: true});

function getStat(dict, bandName, statName) {
  dict = ee.Dictionary(dict);
  var prefixed = bandName + '_' + statName;
  return ee.Number(ee.Algorithms.If(
    dict.contains(prefixed),
    dict.get(prefixed),
    dict.get(statName, 0)
  ));
}

function summarizeVar(img, bandName, mask) {
  return ee.Dictionary(img.updateMask(mask).reduceRegion({
    reducer: statReducer,
    geometry: region,
    crs: ANALYSIS_CRS,
    scale: ANALYSIS_SCALE,
    maxPixels: 1e13,
    tileScale: 8
  }));
}

function calcSMD(vMean, eMean, vSd, eSd) {
  var pooled = vSd.pow(2).add(eSd.pow(2)).divide(2).sqrt();
  return safeDiv(vMean.subtract(eMean), pooled);
}

var vNdvi = summarizeVar(ndvi2024, 'NDVI_2024', validMask);
var eNdvi = summarizeVar(ndvi2024, 'NDVI_2024', excludedMask);
var vSlope = summarizeVar(slopeDeg, 'Slope_deg', validMask);
var eSlope = summarizeVar(slopeDeg, 'Slope_deg', excludedMask);

var exp2NdviSlopeTable = ee.FeatureCollection([ee.Feature(null, {
  Spatial_unit: 'Minsk Region',
  Valid_NDVI_mean: getStat(vNdvi, 'NDVI_2024', 'mean'),
  Valid_NDVI_median: getStat(vNdvi, 'NDVI_2024', 'p50'),
  Valid_NDVI_stdDev: getStat(vNdvi, 'NDVI_2024', 'stdDev'),
  Valid_NDVI_P25: getStat(vNdvi, 'NDVI_2024', 'p25'),
  Valid_NDVI_P75: getStat(vNdvi, 'NDVI_2024', 'p75'),
  Excluded_NDVI_mean: getStat(eNdvi, 'NDVI_2024', 'mean'),
  Excluded_NDVI_median: getStat(eNdvi, 'NDVI_2024', 'p50'),
  Excluded_NDVI_stdDev: getStat(eNdvi, 'NDVI_2024', 'stdDev'),
  Excluded_NDVI_P25: getStat(eNdvi, 'NDVI_2024', 'p25'),
  Excluded_NDVI_P75: getStat(eNdvi, 'NDVI_2024', 'p75'),
  NDVI_SMD: calcSMD(
    getStat(vNdvi, 'NDVI_2024', 'mean'),
    getStat(eNdvi, 'NDVI_2024', 'mean'),
    getStat(vNdvi, 'NDVI_2024', 'stdDev'),
    getStat(eNdvi, 'NDVI_2024', 'stdDev')
  ),
  Valid_Slope_mean: getStat(vSlope, 'Slope_deg', 'mean'),
  Valid_Slope_median: getStat(vSlope, 'Slope_deg', 'p50'),
  Valid_Slope_stdDev: getStat(vSlope, 'Slope_deg', 'stdDev'),
  Valid_Slope_P25: getStat(vSlope, 'Slope_deg', 'p25'),
  Valid_Slope_P75: getStat(vSlope, 'Slope_deg', 'p75'),
  Excluded_Slope_mean: getStat(eSlope, 'Slope_deg', 'mean'),
  Excluded_Slope_median: getStat(eSlope, 'Slope_deg', 'p50'),
  Excluded_Slope_stdDev: getStat(eSlope, 'Slope_deg', 'stdDev'),
  Excluded_Slope_P25: getStat(eSlope, 'Slope_deg', 'p25'),
  Excluded_Slope_P75: getStat(eSlope, 'Slope_deg', 'p75'),
  Slope_SMD: calcSMD(
    getStat(vSlope, 'Slope_deg', 'mean'),
    getStat(eSlope, 'Slope_deg', 'mean'),
    getStat(vSlope, 'Slope_deg', 'stdDev'),
    getStat(eSlope, 'Slope_deg', 'stdDev')
  )
})]);

// ----------------------------------------------------------------------------
// 4. EXPORTS
// ----------------------------------------------------------------------------
Export.table.toDrive({
  collection: exp2LcTable,
  description: 'EXP2_RUSLE_Domain_Landcover_Composition_SI_2024',
  folder: DRIVE_FOLDER,
  fileFormat: 'CSV'
});

Export.table.toDrive({
  collection: exp2cTable,
  description: 'EXP2C_RUSLE_Landcover_Stratified_Disagreement_SI_2024',
  folder: DRIVE_FOLDER,
  fileFormat: 'CSV'
});

Export.table.toDrive({
  collection: exp2NdviSlopeTable,
  description: 'EXP2_RUSLE_Domain_NDVI_Slope_Audit_SI_2024',
  folder: DRIVE_FOLDER,
  fileFormat: 'CSV'
});

print('EXP2 modular audit tasks queued.');
