// ============================================================================
// 01 - NDVI-C AND RUSLE SCENARIO CONSTRUCTION
// ============================================================================
// This script reconstructs the Sentinel-2 NDVI -> C -> RUSLE part of the paper
// from prepared R, K, L and S (or LS) inputs. Hydrological DEM conditioning and
// the final L/S construction were performed outside GEE and are treated as inputs.
// ============================================================================

// ------------------------- USER SETTINGS ------------------------------------
var DISTRICTS_ASSET = 'projects/YOUR_EE_PROJECT/assets/MINSK_DISTRICTS';
var R_ASSET         = 'projects/YOUR_EE_PROJECT/assets/R_FACTOR_2024';
var K_ASSET         = 'projects/YOUR_EE_PROJECT/assets/K_FACTOR';
var L_ASSET         = 'projects/YOUR_EE_PROJECT/assets/L_FACTOR';
var S_ASSET         = 'projects/YOUR_EE_PROJECT/assets/S_FACTOR';

var ANALYSIS_CRS = 'EPSG:32635';
var ANALYSIS_SCALE = 30;
var START_DATE = '2024-01-01';
var END_DATE   = '2025-01-01';
var MAX_SCENE_CLOUD = 40;

var districts = ee.FeatureCollection(DISTRICTS_ASSET);
var region = districts.geometry();

var R = ee.Image(R_ASSET).select(0).toFloat().rename('R');
var K = ee.Image(K_ASSET).select(0).toFloat().rename('K');
var L = ee.Image(L_ASSET).select(0).toFloat().rename('L');
var S = ee.Image(S_ASSET).select(0).toFloat().rename('S');
var LS = L.multiply(S).rename('LS');
var P = ee.Image.constant(1).rename('P').toFloat();

// ------------------------- SENTINEL-2 NDVI ----------------------------------
function prepareNDVI(image) {
  var scl = image.select('SCL');
  var qualityMask = scl.neq(3)
    .and(scl.neq(8))
    .and(scl.neq(9))
    .and(scl.neq(10))
    .and(scl.neq(11));

  var red = image.select('B4').multiply(0.0001);
  var nir = image.select('B8').multiply(0.0001);
  var denominator = nir.add(red);

  return nir.subtract(red)
    .divide(denominator)
    .rename('NDVI')
    .updateMask(denominator.neq(0))
    .updateMask(qualityMask)
    .copyProperties(image, ['system:time_start']);
}

var s2NDVI = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(region)
  .filterDate(START_DATE, END_DATE)
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', MAX_SCENE_CLOUD))
  .map(prepareNDVI);

var ndvi = s2NDVI.median().rename('NDVI_2024').clip(region);
var validNDVI = ndvi.gte(0);
ndvi = ndvi.updateMask(validNDVI);

// ------------------------- C-FACTOR SCENARIOS -------------------------------
// VDK: C = exp[-2*NDVI/(1-NDVI)]
var cVDK = ndvi.expression(
  'exp((-2 * n) / (1 - n))',
  {n: ndvi}
).rename('C_VDK').clamp(0.001, 1);

// Linear: C = (1-NDVI)/2
var cLinear = ee.Image(1)
  .subtract(ndvi)
  .divide(2)
  .rename('C_LINEAR')
  .clamp(0.001, 1);

// ------------------------- RUSLE --------------------------------------------
var base = R.multiply(K).multiply(LS).multiply(P);

var aVDK = base.multiply(cVDK)
  .updateMask(validNDVI)
  .rename('A_VDK')
  .toFloat();

var aLinear = base.multiply(cLinear)
  .updateMask(validNDVI)
  .rename('A_LINEAR')
  .toFloat();

Map.centerObject(region, 7);
Map.addLayer(ndvi, {min: 0, max: 0.9, palette: ['FFFFFF','66BD63','006837']}, 'Annual median NDVI', false);
Map.addLayer(cVDK, {min: 0, max: 0.5, palette: ['FFFFFF','FDAE61','D73027']}, 'C VDK', false);
Map.addLayer(cLinear, {min: 0, max: 0.5, palette: ['FFFFFF','FDAE61','D73027']}, 'C Linear', false);
Map.addLayer(aVDK, {min: 0, max: 50, palette: ['FFFFCC','FED976','FEB24C','FD8D3C','BD0026']}, 'A VDK', false);
Map.addLayer(aLinear, {min: 0, max: 50, palette: ['FFFFCC','FED976','FEB24C','FD8D3C','BD0026']}, 'A Linear', false);

// Optional exports. Uncomment and set destination if needed.
// Export.image.toAsset({image: aVDK, description: 'RUSLE_VDK_2024', assetId: 'projects/YOUR_EE_PROJECT/assets/RUSLE_VDK_2024', region: region, scale: ANALYSIS_SCALE, crs: ANALYSIS_CRS, maxPixels: 1e13});
// Export.image.toAsset({image: aLinear, description: 'RUSLE_LINEAR_2024', assetId: 'projects/YOUR_EE_PROJECT/assets/RUSLE_LINEAR_2024', region: region, scale: ANALYSIS_SCALE, crs: ANALYSIS_CRS, maxPixels: 1e13});
