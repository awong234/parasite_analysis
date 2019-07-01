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

data = data[complete.cases(data %>% select(fmagna_ff, dsl_mb, Easting:geomorphon_20)),]

# Change data to binary

data %<>% mutate(fmagna_ff = as.integer(fmagna_ff > 0),
                 dsl_mb    = as.integer(dsl_mb > 0))

# Select which geomorphon to use
# Select 40-cell scale

s_geomorph = 'geomorphon_40'

data$geomorphon = data[, s_geomorph]

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

# # No longer interested in distance to wetland now that we have geomorphon
# # Add wetland distance
# 
# wetland_dist = raster::raster('../../GIS/NY_shapefile_wetlands/dist_to_wetlands.tif')
# 
# predict_grid@data$Distance_to_wetland = raster::extract(wetland_dist, predict_grid)

# Add geomorphon
# use selection from before to parameterize

geomorphon = raster::raster(paste0('../../GIS/geomorphon/', s_geomorph, '.tif'))
predict_grid@data$geomorphon = raster::extract(geomorphon, predict_grid, method = 'simple', fun =)

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



# List all models to run and run them! -----------------------------------------

f_magna_models = list()

# Null model

f_magna_models[['null_model']] = as.formula("fmagna_ff ~ 1")

# Full model with spatial effect
f_magna_models[['full_model']] = as.formula("fmagna_ff ~ Easting + Northing + Precipitation + Snow + geomorphon + Elevation + f(i, model = spde)")

# Full model with elevation random walk
f_magna_models[['full_model_elev_rw']] = as.formula("fmagna_ff ~ Easting + Northing + Precipitation + Snow + geomorphon + f(Elevation, model = 'rw1') + f(i, model = spde)")

# Reduced models

# Remove spatial effect, all covariates and rw elevation
f_magna_models[['red_mod_minus_speff']] = as.formula("fmagna_ff ~ Easting + Northing + Precipitation + Snow + geomorphon + f(Elevation, model = 'rw1')")

# Remove elevation rw; linear in all covariates
f_magna_models[['red_mod_linear_all_cov']] = as.formula("fmagna_ff ~ Easting + Northing + Precipitation + Snow + geomorphon + Elevation")

# Survival hypothesis -- precip, snow, wetland
f_magna_models[['red_mod_survival']] = as.formula("fmagna_ff ~ Precipitation + Snow + geomorphon")

# No survival -- large-scale and small-scale spatial variation

f_magna_models[['red_mod_spatial']] = as.formula("fmagna_ff ~ Easting + Northing + f(i, model = spde)")

# No survival -- with linear elevation also

f_magna_models[['red_mod_spatial_elev']] = as.formula("fmagna_ff ~ Easting + Northing + f(i, model = spde) + Elevation")

#'###############################################
# MODEL RUNS ###################################
#'###############################################
#
#

# Make the A matrices and spde for all subsequent stacks
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

projector_A = inla.spde.make.A(adk_mesh, loc=data %>% select(Easting_real, Northing_real) %>% as.matrix())
predictor_A = inla.spde.make.A(adk_mesh, loc = predict_grid_scaled@coords)


spde <- inla.spde2.pcmatern(
  mesh = adk_mesh,
  alpha = 2,
  ### mesh and smoothness parameter
  prior.range = c(1000, 0.01),
  ### P(practic.range<0.3)=0.5
  prior.sigma = c(1, 0.01)
  ### P(sigma>1)=0.01
)

# Null model -------------------------------

null_model = inla(formula = f_magna_models$null_model, family = 'binomial', data = data, control.compute = list(waic = TRUE, cpo  = TRUE))

null_model_ls = list(
  "data_stack" = NA,
  "join_stack" = NA,
  "model"      = null_model
)

save(null_model_ls, file = "model_outputs/logreg/fmagna/null_model.Rdata")

# full_model -------------------------------

# Set up stacks
{
  
  
  # ggplot() +
  #   gg(data_sp) +
  #   gg(adk_mesh) +
  #   gg(adkbound) +
  #   coord_fixed(ratio = 1)
  
  
  stk <- inla.stack(
    data = list(fmagna_ff = data$fmagna_ff),
    A = list(projector_A, 1),
    effects = list(
      list(i = 1:spde$n.spde),
      data.frame(
        Intercept = 1,
        Easting = data$Easting,
        Elevation = data$Elevation,
        Northing = data$Northing,
        Precipitation = data$Precipitation,
        Snow = data$Snow,
        geomorphon = data$geomorphon
      )
    ),
    tag = 'dat'
  )
  
  stk_predictor = inla.stack(
    data = list(fmagna_ff = NA),
    A = list(predictor_A, 1),
    effects = list(
      list(i = 1:spde$n.spde),
      data.frame(
        Intercept = 1,
        Easting = predict_grid_scaled@data$Easting,
        Elevation = predict_grid_scaled@data$Elevation,
        Northing = predict_grid_scaled@data$Northing,
        Precipitation = predict_grid_scaled@data$Precipitation,
        Snow = predict_grid_scaled@data$Snow,
        geomorphon = predict_grid_scaled$geomorphon
      )
    ),
    tag = 'pred'
  )
  
  stk_jp = inla.stack(stk, stk_predictor)
  
}

full_model = inla(formula = f_magna_models$full_model,
                  family = 'binomial',
                  data = inla.stack.data(stk),
                  control.compute = list(waic = TRUE,
                                         cpo  = TRUE),
                  control.predictor = list(A = inla.stack.A(stk))
)

full_model_ls = list(
  "data_stack" = stk,
  "joint_stack" = stk_jp,
  "model" = full_model
)

save(full_model_ls, file = "model_outputs/logreg/full_model.Rdata")

Cairo::Cairo(width = 1024*3, height = 768*3, file = 'images/logreg_full_cpo.png', dpi = 100*4)
ggplot() + 
  geom_point(aes(x = seq_along(full_model$cpo$cpo), y = full_model$cpo$cpo)) +
  theme_bw() + xlab("Index") + ylab("Conditional Predictive Ordinate") + ggtitle("CPO values for the full model fit to binary data")
dev.off()

Cairo::Cairo(width = 1024*3, height = 768*3, file = 'images/logreg_full_pit.png', dpi = 100*4)
ggplot() + 
  geom_point(aes(x = (1:n) / (n+1), y = sort(full_model$cpo$pit))) +
  geom_abline(slope = 1, intercept = 0) +
  xlim(c(0,1)) + ylim(c(0,1)) + coord_equal() + 
  theme_bw() + xlab("Theoretical quantiles") + ylab("PIT probabilities") + ggtitle("PIT values for the full model fit to binary data")
dev.off()

