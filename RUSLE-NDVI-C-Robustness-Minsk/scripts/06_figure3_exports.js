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


function classifyErosionClass5(img) {
  return ee.Image.constant(1)
    .where(img.gte(5).and(img.lt(10)), 2)
    .where(img.gte(10).and(img.lt(20)), 3)
    .where(img.gte(20).and(img.lt(50)), 4)
    .where(img.gte(50), 5)
    .updateMask(commonMask)
    .clip(region)
    .toInt8();
}
var classVDK = classifyErosionClass5(rusleVDK).rename('ErosionClass_VDK');
var classLinear = classifyErosionClass5(rusleLinear).rename('ErosionClass_LINEAR');

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
  'FFFFCC', // <5
  'FED976', // 5-<10
  'FEB24C', // 10-<20
  'FD8D3C', // 20-<50
  'BD0026'  // >=50
];


// Diverging class-difference palette:
// -4 ... 0 ... +4
var FIG3_DIFF_PALETTE = [
  '2166AC', // -4
  '4393C3', // -3
  '92C5DE', // -2
  'D1E5F0', // -1
  'F7F7F7', //  0
  'FDDBC7', // +1
  'F4A582', // +2
  'D6604D', // +3
  'B2182B'  // +4
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