// =============================================================================
// Portfolio project: Soil Erosion Risk with RUSLE in Google Earth Engine
// Upland sub-watershed, East Java (Brantas river basin)
// -----------------------------------------------------------------------------
// Self-initiated portfolio project using public data.
// Author: Utik Tri Wulan Cahya, PhD
//
// A = R x K x LS x C x P   (t ha-1 yr-1)
//   R  rainfall erosivity      CHIRPS annual rainfall, Renard & Freimund (1994)
//   K  soil erodibility        OpenLandMap sand/clay/organic C, Williams (1995) EPIC
//   LS slope length-steepness  Copernicus DEM slope + MERIT Hydro upstream area,
//                              Moore & Burch (1986)
//   C  cover management        Sentinel-2 NDVI, Van der Knijff et al. (1999)
//   P  support practice        contouring on cropland by slope class,
//                              Wischmeier & Smith (1978); 1 elsewhere
//
// Outputs (Google Drive folder RUSLE_erosion):
//   1. RUSLE_layers.tif          bands: A, R, K, LS, C, P, slope_pct, landcover, erosion_class
//   2. RUSLE_area_by_class.csv   area (ha) per land-cover x erosion class
//   3. RUSLE_watershed.geojson   watershed boundary
//
// How to use: paste into the Code Editor, move POINT to the sub-watershed you
// want (NOT your dissertation area), click Run, check the Console, then click
// RUN on the three exports in the Tasks tab.
// =============================================================================


// ----- 1. SETTINGS -----------------------------------------------------------

var POINT        = ee.Geometry.Point([112.78, -8.00]);   // move to your chosen sub-watershed
var BASIN_LEVEL  = 9;          // HydroBASINS level: 8 larger, 10 smaller
var RAIN_YEARS   = [2014, 2023];  // years for mean annual rainfall (inclusive)
var NDVI_YEAR    = 2023;       // year for the C factor
var MAX_SLOPE_LEN = 300;       // m; cap on slope length so river channels are not treated as hillslopes
var SCALE        = 30;         // m; output resolution
var CRS          = 'EPSG:32749';
var FOLDER       = 'RUSLE_erosion';


// ----- 2. STUDY AREA ---------------------------------------------------------

var watershed = ee.FeatureCollection('WWF/HydroSHEDS/v1/Basins/hybas_' + BASIN_LEVEL)
                  .filterBounds(POINT);
var aoi = watershed.geometry();
print('Watershed area (km2):', aoi.area(10).divide(1e6));


// ----- 3. R FACTOR: rainfall erosivity ---------------------------------------
// Mean annual rainfall P (mm) from CHIRPS, then Renard & Freimund (1994):
//   P <= 850 mm : R = 0.0483 * P^1.61
//   P  > 850 mm : R = 587.8 - 1.219 P + 0.004105 P^2      (MJ mm ha-1 h-1 yr-1)

var nYears = RAIN_YEARS[1] - RAIN_YEARS[0] + 1;
var P = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')
          .filterDate(RAIN_YEARS[0] + '-01-01', (RAIN_YEARS[1] + 1) + '-01-01')
          .select('precipitation').sum().divide(nYears)
          .rename('P_mm').clip(aoi.buffer(10000));
var R = P.expression(
  '(P <= 850) ? 0.0483 * pow(P, 1.61) : 587.8 - 1.219 * P + 0.004105 * P * P',
  {P: P}).rename('R');


// ----- 4. K FACTOR: soil erodibility (Williams 1995, EPIC) -------------------
// Topsoil (0 cm) sand and clay in %, organic carbon in g/kg (stored x 5)

var sand = ee.Image('OpenLandMap/SOL/SOL_SAND-WFRACTION_USDA-3A1A1A_M/v02').select('b0');
var clay = ee.Image('OpenLandMap/SOL/SOL_CLAY-WFRACTION_USDA-3A1A1A_M/v02').select('b0');
var ocPct = ee.Image('OpenLandMap/SOL/SOL_ORGANIC-CARBON_USDA-6A1C_M/v02')
              .select('b0').multiply(5).divide(10);   // g/kg -> %
var silt = ee.Image(100).subtract(sand).subtract(clay).max(0);

var K = sand.expression(
  '(0.2 + 0.3 * exp(-0.256 * SAN * (1 - SIL / 100)))' +
  ' * pow(SIL / (CLA + SIL), 0.3)' +
  ' * (1 - 0.25 * C / (C + exp(3.72 - 2.95 * C)))' +
  ' * (1 - 0.7 * (1 - SAN / 100) / ((1 - SAN / 100) + exp(-5.51 + 22.9 * (1 - SAN / 100))))' +
  ' * 0.1317',                                         // US units -> t ha h ha-1 MJ-1 mm-1
  {SAN: sand, SIL: silt, CLA: clay, C: ocPct}).rename('K');


// ----- 5. LS FACTOR (Moore & Burch 1986) -------------------------------------
// LS = (As / 22.13)^0.4 * (sin(slope) / 0.0896)^1.3
// As = specific catchment area = upstream area / cell width, capped at MAX_SLOPE_LEN

var dem = ee.ImageCollection('COPERNICUS/DEM/GLO30_2024_1')
            .filterBounds(aoi).select('DEM').mosaic()
            .setDefaultProjection('EPSG:4326', null, 30);
var slopeDeg = ee.Terrain.slope(dem);
var slopeRad = slopeDeg.multiply(Math.PI / 180);
var slopePct = slopeRad.tan().multiply(100).rename('slope_pct');

var merit = ee.Image('MERIT/Hydro/v1_0_1');
var cellWidth = 92.77;
var As = merit.select('upa').multiply(1e6).divide(cellWidth)   // km2 -> m2, per metre of width
              .min(MAX_SLOPE_LEN);
var LS = As.divide(22.13).pow(0.4)
           .multiply(slopeRad.sin().divide(0.0896).pow(1.3))
           .rename('LS');


// ----- 6. C FACTOR (Van der Knijff et al. 1999) ------------------------------
// C = exp(-2 * NDVI / (1 - NDVI)), from a cloud-masked Sentinel-2 median

var csPlus = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED');
var ndvi = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
             .filterBounds(aoi)
             .filterDate(NDVI_YEAR + '-01-01', (NDVI_YEAR + 1) + '-01-01')
             .linkCollection(csPlus, ['cs_cdf'])
             .map(function (img) {
               return img.normalizedDifference(['B8', 'B4']).toFloat()
                         .updateMask(img.select('cs_cdf').gte(0.6));
             })
             .median().rename('NDVI');
var C = ndvi.expression('exp(-2 * N / (1 - N))', {N: ndvi.clamp(-0.99, 0.99)})
            .clamp(0, 1).rename('C');


// ----- 7. LAND COVER AND P FACTOR --------------------------------------------
// ESA WorldCover 2021: 10 tree, 20 shrub, 30 grass, 40 crop, 50 built-up,
// 60 bare/sparse, 70 snow, 80 water, 90 wetland, 95 mangrove, 100 moss
// P: contouring on cropland by slope (Wischmeier & Smith 1978); 1 elsewhere

var lc = ee.ImageCollection('ESA/WorldCover/v200').first().select('Map').rename('landcover');
var pCrop = slopePct.expression(
  '(s <= 2) ? 0.6 : (s <= 8) ? 0.5 : (s <= 12) ? 0.6 : (s <= 16) ? 0.7 : ' +
  '(s <= 20) ? 0.8 : (s <= 25) ? 0.9 : 1.0', {s: slopePct});
var Pf = ee.Image(1).where(lc.eq(40), pCrop).rename('P');


// ----- 8. SOIL LOSS AND EROSION CLASSES --------------------------------------
// Built-up and water are excluded (no exposed soil).
// Classes (t ha-1 yr-1), as used for erosion hazard in Indonesia:
//   1 very light < 15 | 2 light 15-60 | 3 moderate 60-180 | 4 heavy 180-480 | 5 very heavy > 480

var soilMask = lc.neq(50).and(lc.neq(80));
var A = R.multiply(K).multiply(LS).multiply(C).multiply(Pf)
         .updateMask(soilMask).rename('A');
var erosionClass = ee.Image(1)
  .where(A.gte(15), 2).where(A.gte(60), 3).where(A.gte(180), 4).where(A.gte(480), 5)
  .updateMask(A.mask()).rename('erosion_class');

var layers = A.addBands([R, K, LS, C, Pf, slopePct, lc, erosionClass]).toFloat().clip(aoi);


// ----- 9. QUICK CHECKS IN THE CONSOLE ----------------------------------------

var stats = layers.select(['A', 'R', 'K', 'LS', 'C', 'P']).reduceRegion({
  reducer: ee.Reducer.mean().combine(ee.Reducer.percentile([50, 95]), '', true),
  geometry: aoi, scale: 90, maxPixels: 1e10, bestEffort: true
});
print('Factor summary (mean, median p50, p95):', stats);


// ----- 10. AREA BY LAND COVER x EROSION CLASS --------------------------------

var combo = lc.multiply(10).add(erosionClass).rename('combo');
var areaHa = ee.Image.pixelArea().divide(1e4).addBands(combo).updateMask(A.mask());
var grouped = areaHa.reduceRegion({
  reducer: ee.Reducer.sum().group({groupField: 1, groupName: 'combo'}),
  geometry: aoi, scale: SCALE, crs: CRS, maxPixels: 1e10, tileScale: 4
});
var areaTable = ee.FeatureCollection(ee.List(grouped.get('groups')).map(function (g) {
  g = ee.Dictionary(g);
  var code = ee.Number(g.get('combo'));
  return ee.Feature(null, {
    landcover_code: code.divide(10).floor(),
    erosion_class: code.mod(10),
    area_ha: g.get('sum')
  });
}));
print('Area by land cover x erosion class (first rows):', areaTable.limit(10));


// ----- 11. PREVIEW -----------------------------------------------------------

Map.centerObject(aoi, 11);
Map.addLayer(C, {min: 0, max: 0.5, palette: ['#1f6b33', '#f3efe0']}, 'C factor', false);
Map.addLayer(LS.clip(aoi), {min: 0, max: 20, palette: ['#f3efe0', '#7a7a7a', '#0b0b0b']}, 'LS factor', false);
Map.addLayer(erosionClass.clip(aoi),
  {min: 1, max: 5, palette: ['#1a9850', '#a6d96a', '#fee08b', '#f46d43', '#a50026']}, 'Erosion class');
Map.addLayer(watershed.style({color: '#000000', fillColor: '00000000', width: 2}), {}, 'Watershed');


// ----- 12. EXPORTS -----------------------------------------------------------

Export.image.toDrive({
  image: layers, description: 'RUSLE_layers', folder: FOLDER,
  region: aoi.bounds(), scale: SCALE, crs: CRS, maxPixels: 1e10, fileFormat: 'GeoTIFF'
});
Export.table.toDrive({
  collection: areaTable, description: 'RUSLE_area_by_class', folder: FOLDER, fileFormat: 'CSV',
  selectors: ['landcover_code', 'erosion_class', 'area_ha']
});
Export.table.toDrive({
  collection: watershed.select(['HYBAS_ID', 'SUB_AREA']), description: 'RUSLE_watershed',
  folder: FOLDER, fileFormat: 'GeoJSON'
});

// Data credits: CHIRPS (Funk et al. 2015); OpenLandMap soil (Hengl 2018, CC BY-SA 4.0);
// Copernicus WorldDEM-30 (c) DLR e.V. 2010-2014 and (c) Airbus Defence and Space GmbH
// 2014-2018, provided under COPERNICUS by the EU and ESA; MERIT Hydro (Yamazaki et al.
// 2019; CC BY-NC 4.0 or ODbL 1.0); Sentinel-2 (modified Copernicus Sentinel data 2023);
// Cloud Score+; ESA WorldCover 2021 (CC BY 4.0); HydroSHEDS (Lehner & Grill 2013).
