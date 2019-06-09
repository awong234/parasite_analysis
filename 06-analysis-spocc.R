# This file was made after 06-analysis_HDS.R

# The overdispersed zero-inflated-type data were not able to be modeled well with INLA tools. This
# document is taking the observations and performing a simpler binary (spatial, if possible)
# regression.



# Setup ----------------------------------------

library(INLA)
library(inlabru)
library(magrittr)
library(ggplot2)
library(viridis)
library(reshape2)
library(dplyr)

source('functions.R')
load(file = 'scaled_covariates.Rdata')
load(file = 'scaled_covariates_attr.Rdata')
# load(file = '.RData')

# Load data ----------------------------------------

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

data = bind_cols(metadata_join %>% rename(Easting_real = Easting, Northing_real = Northing), covariate_scaled_df)

# missing_data_rows = is.na(data$Easting) | is.na(data$Northing)

data = data[complete.cases(data %>% select(fmagna_ff, dsl_mb, Easting:Distance_to_wetland)),]

# Change data to binary

data %<>% mutate(fmagna_ff_occ = fmagna_ff > 0,
                 dsl_mb_occ    = dsl_mb > 0)



# Data points as sp

data_sp = data
coordinates(data_sp) = ~Easting_real + Northing_real

# Create mesh

adkbound = rgdal::readOGR(dsn = '../../GIS/adkParkPoly/adkParkPoly.shp')

boundary <- list(
  inla.nonconvex.hull(coordinates(data_sp), 10000),
  inla.nonconvex.hull(coordinates(data_sp), 30000))


adk_mesh = inla.mesh.2d(boundary=boundary,
                        max.edge=c(2000, 5400),
                        min.angle=c(30, 21),
                        max.n=c(48000, 16000), ## Safeguard against large meshes.
                        max.n.strict=c(128000, 128000), ## Don't build a huge mesh!
                        cutoff=1050, ## Filter away adjacent points.
                        offset=c(5200, 15000)) ## Offset for extra boundaries, if needed.

# ggplot() +
#   gg(data_sp) +
#   gg(adk_mesh) +
#   gg(adkbound) +
#   coord_fixed(ratio = 1)



# Load prediction grid --------------------------------------------------------------

load('predict_grid_1000.Rdata')

# Add new predictors

# Add Precipitation --------------------------------------------------------

# Precipitation will be the mean of the three years.

precip_names = dir(path = '../../GIS/PRISM_precip_data/', pattern = 'PRISM.+tif$', full.names = F)
precip_files = dir(path = '../../GIS/PRISM_precip_data/', pattern = 'PRISM.+tif$', full.names = T)


precip_dates = data.frame(file = precip_names, year = precip_names %>% {regmatches(x = ., m = regexec(pattern = '\\d{4}', text = ., perl = T))} %>% as.integer)

precip_data = matrix(NA, nrow = NROW(predict_grid), ncol = length(precip_files))

for(f in 1:length(precip_files)){
  
  precip = raster::raster(precip_files[f])
  
  precip_data[,f] = raster::extract(precip, predict_grid)
  
}

precip_means = precip_data %>% rowMeans()

predict_grid@data$Precipitation = precip_means

# Add snow ----------------------------------------------------------------------


snowcover_names = dir(path = '../../GIS/NWS_snowfall_data/', pattern = 'snowfall.*.tif$')
snowcover_paths = dir(path = '../../GIS/NWS_snowfall_data/', pattern = 'snowfall.*.tif$', full.names = T)

snowcover_dates = data.frame(file = snowcover_names, year = snowcover_names %>% {regmatches(x = ., m = regexec(pattern = '\\d{4}', text = ., perl = T))} %>% as.integer)

# Extract snowcover values for year prior from location

snowcover_data = matrix(NA, nrow = NROW(predict_grid), ncol = length(snowcover_paths))

for(f in 1:length(snowcover_paths)){
  
  snowcover = raster::raster(snowcover_paths[f])
  
  snowcover_data[,f] = raster::extract(snowcover, predict_grid)
  
}

# Mean snowcover
snowcover_means = snowcover_data %>% rowMeans()

predict_grid@data$Snow = snowcover_means

# Add wetland distance

wetland_dist = raster::raster('../../GIS/NY_shapefile_wetlands/dist_to_wetlands.tif')

predict_grid@data$Distance_to_wetland = raster::extract(wetland_dist, predict_grid)

# Scale prediction

predict_grid_complete = predict_grid[complete.cases(predict_grid@data),]

covariate_scaled_attr$covar = row.names(covariate_scaled_attr)
rownames(covariate_scaled_attr) = NULL

predict_grid_scaled = predict_grid_complete

rel_cols = covariate_scaled_attr$covar

# Center and scale 
for(c in rel_cols){
  print(c)
  print(head(predict_grid_scaled@data[[c]]))
  print(covariate_scaled_attr %>% filter(covar == c) %>% pull(Center))
  print(covariate_scaled_attr %>% filter(covar == c) %>% pull(Scale))
  # Center
  predict_grid_scaled@data[[c]] = predict_grid_scaled@data[[c]] - covariate_scaled_attr %>% filter(covar == c) %>% pull(Center)
  # Scale
  predict_grid_scaled@data[[c]] = predict_grid_scaled@data[[c]] / covariate_scaled_attr %>% filter(covar == c) %>% pull(Scale)
  
}
