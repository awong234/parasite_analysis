# Setup ------------------------
library(ggplot2)
library(unmarked)
library(dplyr)
library(sp)
library(lubridate)

source('functions.R')

# Covariates for intensity of parasite ----------------------------------------

# Prep spatial data ------------------------------------------------------------

metadata = read.csv(file = 'metadata_adjusted.csv', stringsAsFactors = F)

metadata$Date = as.Date(metadata$Date)

ff = read.csv(file = 'fecal_flukefinder.csv', stringsAsFactors = F)
mb = read.csv(file = 'fecal_MB.csv', stringsAsFactors = F)
quant = read.csv(file = 'fecal_quant.csv', stringsAsFactors = F)

metadata_join = metadata %>% left_join(ff %>% select(PK, Total_eggs)) %>% rename(fmagna_ff = Total_eggs) %>% 
  left_join(mb %>% select(PK, Total_larvae_dsl, total_eggs_fmagna)) %>% rename(dsl_mb = Total_larvae_dsl, fmagna_mb = total_eggs_fmagna) %>% 
  left_join(quant %>% select(PK, Fascioloides_magna, Protostrongylid_DSL)) %>% rename(fmagna_quant = Fascioloides_magna, dsl_quant = Protostrongylid_DSL)

# Which na

na_location_index = is.na(metadata$Easting) | is.na(metadata$Northing)

# Easting and Northing ----------------------------------------------------------------

metadata_sp = metadata[!na_location_index,]

coordinates(metadata_sp) = ~Easting + Northing

covariate = data.frame(Easting = metadata$Easting, Northing = metadata$Northing)

# Prep precipitation data ------------------------------------------------------------

# NOTICE: The precipitation variable is a derived variable. It references the spatial estimate of precipitation for each data point for the month *preceding* the data collection month. 

precip_names = dir(path = '../../GIS/PRISM_precip_data/', pattern = 'PRISM.+tif$', full.names = F)
precip_files = dir(path = '../../GIS/PRISM_precip_data/', pattern = 'PRISM.+tif$', full.names = T)


precip_dates = data.frame(file = precip_names, year = precip_names %>% {regmatches(x = ., m = regexec(pattern = '\\d{4}', text = ., perl = T))} %>% as.integer)

precip_data = matrix(NA, nrow = nrow(metadata), ncol = length(precip_files))

for(f in 1:length(precip_files)){
  
  precip = raster::raster(precip_files[f])
  
  precip_data[!na_location_index,f] = raster::extract(precip, metadata_sp)
  
}

colnames(precip_data) = precip_names

precip_data = as.data.frame(precip_data)

# Using existing function for snowfall from previous year.
precip_covariate = getSnowFromDate(dates = metadata$Date, snow_data = precip_data, snow_dates = precip_dates)

metadata$Precipitation = precip_covariate

covariate$Precipitation = precip_covariate

# Check a few rows

rand = sample(x = seq(1,nrow(metadata)), size = 1)

test = metadata[rand,]

date_test = test$Date

test$Precip == precip_data[rand,which(year(date_test) - 1 == precip_dates$year)]

# Prep snowcover data     ------------------------------------------------------------

snowcover_names = dir(path = '../../GIS/NWS_snowfall_data/', pattern = 'snowfall.*.tif$')
snowcover_paths = dir(path = '../../GIS/NWS_snowfall_data/', pattern = 'snowfall.*.tif$', full.names = T)

snowcover_dates = data.frame(file = snowcover_names, year = snowcover_names %>% {regmatches(x = ., m = regexec(pattern = '\\d{4}', text = ., perl = T))} %>% as.integer)

# Extract snowcover values for year prior from location

snowcover_data = matrix(NA, nrow = nrow(metadata), ncol = length(snowcover_paths))

for(f in 1:length(snowcover_paths)){
  
  snowcover = raster::raster(snowcover_paths[f])
  
  snowcover_data[!na_location_index,f] = raster::extract(snowcover, metadata_sp)
  
}

snow_covariate = getSnowFromDate(dates = metadata$Date, snow_data = snowcover_data, snow_dates = snowcover_dates)

metadata$Snow = snow_covariate

covariate$Snow = snow_covariate

# Check a few rows

rand = sample(x = seq(1,nrow(metadata)), size = 1)

test = metadata[rand,]

date_test = test$Date

test$Snow == snowcover_data[rand,which(year(date_test) - 1 == snowcover_dates$year)]

# Prep elevation data --------------------------------------------------------

elev = raster::raster(x = '../../GIS/Extract_Elev1/Extract_Elev1.tif')

elev_data = rep(NA, times = nrow(metadata)) 

elev_data[!na_location_index] = raster::extract(elev, metadata_sp)

metadata$Elevation = elev_data
covariate$Elevation = elev_data

# Prep dist to wetland data --------------------------------------------------
# 
# # IF it exists already, load it.
# 
# if(file.exists('../../GIS/NY_shapefile_wetlands/Dist_to_wetlands.Rdata')){
#   
#   load('../../GIS/NY_shapefile_wetlands/Dist_to_wetlands.Rdata')
#   
# } else {
#   
#   # Make a raster upon which to calculate distance to nearest wetland boundary
#   
#   cellSize = 1000
#   
#   # Load adk boundary
#   
#   rgdal::ogrListLayers('../../GIS/adkParkPoly/adkParkPoly.shp')
#   
#   adkbound = rgdal::readOGR(dsn = '../../GIS/adkParkPoly/adkParkPoly.shp')
#   
#   # Load uninhabitable areas
#   
#   uhm = rgdal::readOGR(dsn = '../../GIS/uninhabitable_mask/uninhabitable_mask.shp')
#   
#   predict_grid = makegrid(x = adkbound, cellsize = cellSize) %>% rename(x = x1, y = x2)
#   
#   coordinates(predict_grid) = ~x + y
#   proj4string(predict_grid) = proj4string(adkbound)
#   # Promote to spatialgriddataframe
#   gridded(predict_grid) = TRUE
#   predict_grid = as(predict_grid, "SpatialGrid")
#   
#   # Convert to raster to mask by polygons
#   predict_grid = raster::raster(predict_grid)
#   predict_grid[] = 1
#   predict_grid = raster::mask(x = predict_grid, mask = adkbound)
#   
#   predict_grid_raster = raster::mask(x = predict_grid, mask = uhm, inverse = T)
#   
#   # Convert back to grid
#   
#   predict_grid = as(predict_grid_raster, 'SpatialPointsDataFrame')
#   predict_grid = predict_grid[!is.na(predict_grid@data$layer),]
#   gridded(predict_grid) = T
#   summary(predict_grid)
#   plot(predict_grid)
#   
#   # Obtain wetland features
#   
#   wetlands = rgdal::readOGR(dsn = '../../GIS/NY_shapefile_wetlands/NY_shapefile_wetlands_proj_clip.shp')
#   
#   # Condense into a single feature
#   
#   wetlands_union = rgeos::gUnaryUnion(wetlands)
#   
#   # Calculate distance to nearest wetland feature, upon grid
#   
#   dists = rgeos::gDistance(spgeom1 = wetlands_union, spgeom2 = as(predict_grid, "SpatialPoints"), byid = T)
#   
#   save('dists', file = '../../GIS/NY_shapefile_wetlands/Dist_to_wetlands.Rdata')
#   
#   
# }
# 
# if(file.exists('../../GIS/NY_shapefile_wetlands/dist_to_wetlands.tif')){
#   
#   dist_to_wetland_rast = raster::raster('../../GIS/NY_shapefile_wetlands/dist_to_wetlands.tif')
#   
# } else {
#   
#   dist_to_wetland_rast = predict_grid_raster
#   
#   dist_to_wetland_rast[dist_to_wetland_rast > 0] = dists
#   
#   raster::writeRaster(dist_to_wetland_rast, filename = '../../GIS/NY_shapefile_wetlands/dist_to_wetlands.tif', format = 'GTiff')
#   
# }
# 
# dist_to_wetland_data = rep(NA, nrow(metadata))
# extraction = raster::extract(dist_to_wetland_rast, metadata_sp)
# extraction[is.na(extraction)] = 0
# 
# dist_to_wetland_data[!na_location_index] = extraction
# 
# metadata$Distance_to_wetland = dist_to_wetland_data
# covariate$Distance_to_wetland = dist_to_wetland_data


# Prep geomorphon data ----------------------------------------------------------------

files = c('../../GIS/geomorphon/geomorphon_40.tif', '../../GIS/geomorphon/geomorphon_20.tif')
names = regmatches(x = files, m = regexpr(pattern = '\\w+(?=\\.)', text = files, ignore.case = T, perl = T))


for(f in seq_along(files)){
  geomorphon = raster::raster(files[f])
  geomorph_data = as.integer(raster::extract(geomorphon, metadata_sp, method = 'simple'))
  metadata[!na_location_index, names[f]] = geomorph_data
  covariate[!na_location_index, names[f]] = geomorph_data
}

var(metadata$geomorphon_20, na.rm = T)
var(metadata$geomorphon_40, na.rm = T)
hist(metadata$geomorphon_20, breaks = seq(1,9))
hist(metadata$geomorphon_40, breaks = seq(1,9))

# Save covariates

cor(covariate, use = 'complete.obs')

save(covariate, file = 'raw_covariates.Rdata')

# Scale covariate data -------------------------------------------------------------------

# Only scale numeric values; 

integer_test = lapply(covariate, is.integer)

covariate_scaled = scale(covariate[, integer_test == F])
covariate_scaled_df = cbind(covariate_scaled, covariate[, integer_test == T])

covariate_scaled_attr = data.frame(Center = attr(covariate_scaled, which = 'scaled:center'),
                                  Scale  = attr(covariate_scaled, which = 'scaled:scale')
                                  )

save(covariate_scaled_df, file = 'scaled_covariates.Rdata')
save(covariate_scaled_attr, file = 'scaled_covariates_attr.Rdata')
