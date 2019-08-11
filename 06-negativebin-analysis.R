
# Setup ----------------------------------------

library(INLA)
library(inlabru)
library(doParallel)
library(magrittr)
library(ggplot2)
library(MASS)
library(GGally)
library(MuMIn)
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

data %<>% mutate(JulianDay = lubridate::yday(Date) %>% scale())

# missing_data_rows = is.na(data$Easting) | is.na(data$Northing)

data = data[complete.cases(data %>% select(fmagna_ff, dsl_mb, Easting:JulianDay)),]

# Select geomorphon 
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

# Add wetland distance

wetland_dist = raster::raster('../../GIS/NY_shapefile_wetlands/dist_to_wetlands.tif')

predict_grid@data$Distance_to_wetland = raster::extract(wetland_dist, predict_grid)


# Add geomorphon values
# Select 40-cell scale

s_geomorph = 'geomorphon_40'
geomorphon = raster::raster(paste0('../../GIS/geomorphon/', s_geomorph, '.tif'))
temp = raster::extract(geomorphon, predict_grid, method = 'simple')
temp = factor(temp, levels = as.integer(seq(1,10)), labels = c('flat', 'summit', 'ridge', 'shoulder', 
                                                               'spur', 'slope', 'hollow', 'footslope',
                                                               'valley', 'depression'))

predict_grid@data$geomorphon = temp

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

save(adk_mesh, file = 'model_outputs/adk_mesh.Rdata')

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

null_model = inla(formula = f_magna_models$null_model, family = 'nbinomial', data = data, control.compute = list(waic = TRUE, cpo  = TRUE))

null_model_ls = list(
  "data_stack" = NA,
  "join_stack" = NA,
  "model"      = null_model
)

save(null_model_ls, file = "model_outputs/fmagna/null_model.Rdata")

# Full model with spatial effect -------------------------------

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
                  family = 'nbinomial',
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

save(full_model_ls, file = "model_outputs/fmagna/full_model.Rdata")

# Full model with elevation random walk -----------------------------------------------

# Set up stack --

{
  stk <- inla.stack(
    data = list(fmagna_ff = data$fmagna_ff),
    A = list(projector_A, 1),
    effects = list(
      list(i = 1:spde$n.spde),
      data.frame(
        Intercept = 1,
        Easting = data$Easting,
        Elevation = inla.group(data$Elevation, method = 'quantile', n = 30),
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
        Elevation = inla.group(predict_grid_scaled$Elevation, method = 'quantile', n = 30),
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

full_model_elev_rw = inla(formula = f_magna_models$full_model_elev_rw,
                          family = 'nbinomial',
                          data = inla.stack.data(stk),
                          control.compute = list(waic = TRUE, cpo  = TRUE),
                          control.predictor = list(A = inla.stack.A(stk))
)

full_model_elev_rw_ls = list(
  "data_stack" = stk,
  "joint_stack" = stk_jp,
  "model"      = full_model_elev_rw
)

save(full_model_elev_rw_ls, file = "model_outputs/fmagna/full_model_elev_rw_ls.Rdata")

# Reduced: Remove spatial effect, all covariates and rw elevation --------------

# set up stacks
{
  stk <- inla.stack(
    data = list(fmagna_ff = data$fmagna_ff),
    A = list(1),
    effects = list(
      data.frame(
        Intercept = 1,
        Easting = data$Easting,
        Elevation = inla.group(data$Elevation, method = 'quantile', n = 30),
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
    A = list(1),
    effects = list(
      data.frame(
        Intercept = 1,
        Easting = predict_grid_scaled@data$Easting,
        Elevation = inla.group(predict_grid_scaled$Elevation, method = 'quantile', n = 30),
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

red_mod_minus_speff = inla(formula = f_magna_models$red_mod_minus_speff,
                          family = 'nbinomial',
                          data = inla.stack.data(stk),
                          control.compute = list(waic = TRUE, cpo  = TRUE),
                          control.predictor = list(A = inla.stack.A(stk))
)

red_mod_minus_speff_ls = list(
  "data_stack" = stk,
  "joint_stack" = stk_jp,
  "model"      = red_mod_minus_speff
)

save(red_mod_minus_speff_ls, file = "model_outputs/fmagna/red_mod_minus_speff_ls.Rdata")


# Reduced: Remove elevation rw; linear in all covariates ------------------

# No stacks necessary

red_mod_linear_all_cov = inla(formula = f_magna_models$red_mod_linear_all_cov,
                              family = 'nbinomial',
                              data = data,
                              control.compute = list(waic = TRUE, cpo  = TRUE)
                              )

red_mod_linear_all_cov_ls = list(
  "data_stack" = NA,
  "join_stack" = NA,
  "model"      = red_mod_linear_all_cov
)

save(red_mod_linear_all_cov_ls, file = 'model_outputs/fmagna/red_mod_linear_all_cov_ls.Rdata')


# Reduced: Survival hypothesis -- precip, snow, wetland -------------------

red_mod_survival = inla(formula = f_magna_models$red_mod_survival,
                              family = 'nbinomial',
                              data = data,
                              control.compute = list(waic = TRUE, cpo  = TRUE)
)

red_mod_survival_ls = list(
  "data_stack" = NA,
  "join_stack" = NA,
  "model"      = red_mod_survival
)

save(red_mod_survival_ls, file = 'model_outputs/fmagna/red_mod_survival_ls.Rdata')


# Reduced: no survival -- large-scale and small-scale spatial vari --------

# Set up stacks
{

  stk <- inla.stack(
    data = list(fmagna_ff = data$fmagna_ff),
    A = list(projector_A, 1),
    effects = list(
      list(i = 1:spde$n.spde),
      data.frame(
        Intercept = 1,
        Easting = data$Easting,
        Northing = data$Northing
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
        Northing = predict_grid_scaled@data$Northing
      )
    ),
    tag = 'pred'
  )

  stk_jp = inla.stack(stk, stk_predictor)

}

red_mod_spatial = inla(formula = f_magna_models$red_mod_spatial,
                  family = 'nbinomial',
                  data = inla.stack.data(stk),
                  control.compute = list(waic = TRUE, cpo  = TRUE),
                  control.predictor = list(A = inla.stack.A(stk))
)

red_mod_spatial_ls = list(
  "data_stack" = stk,
  "joint_stack" = stk_jp,
  "model" = red_mod_spatial
)

save(red_mod_spatial_ls, file = "model_outputs/fmagna/red_mod_spatial.Rdata")


# Reduced: No survival -- with linear elevation also ----------------------

# Set up stacks
{

  stk <- inla.stack(
    data = list(fmagna_ff = data$fmagna_ff),
    A = list(projector_A, 1),
    effects = list(
      list(i = 1:spde$n.spde),
      data.frame(
        Intercept = 1,
        Easting = data$Easting,
        Northing = data$Northing,
        Elevation = data$Elevation
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
        Northing = predict_grid_scaled@data$Northing,
        Elevation = predict_grid_scaled@data$Elevation
      )
    ),
    tag = 'pred'
  )

  stk_jp = inla.stack(stk, stk_predictor)

}

red_mod_spatial_elev = inla(formula = f_magna_models$red_mod_spatial_elev,
                       family = 'nbinomial',
                       data = inla.stack.data(stk),
                       control.compute = list(waic = TRUE, cpo  = TRUE),
                       control.predictor = list(A = inla.stack.A(stk))
)

red_mod_spatial_elev_ls = list(
  "data_stack" = stk,
  "joint_stack" = stk_jp,
  "model" = red_mod_spatial_elev
)

save(red_mod_spatial_elev_ls, file = "model_outputs/fmagna/red_mod_spatial_elev.Rdata")

# P tenuis models -----------------------------------------------

p_tenuis_models = list()

# Null model

p_tenuis_models[['null_model']] = as.formula("dsl_mb ~ 1")

# Full model with spatial effect
p_tenuis_models[['full_model']] = as.formula("dsl_mb ~ Easting + Northing + Precipitation + Snow + geomorphon + Elevation + f(i, model = spde)")

# Full model with elevation random walk
p_tenuis_models[['full_model_elev_rw']] = as.formula("dsl_mb ~ Easting + Northing + Precipitation + Snow + geomorphon + f(Elevation, model = 'rw1') + f(i, model = spde)")

# Reduced models

# Remove spatial effect, all covariates and rw elevation
p_tenuis_models[['red_mod_minus_speff']] = as.formula("dsl_mb ~ Easting + Northing + Precipitation + Snow + geomorphon + f(Elevation, model = 'rw1')")

# Remove elevation rw; linear in all covariates
p_tenuis_models[['red_mod_linear_all_cov']] = as.formula("dsl_mb ~ Easting + Northing + Precipitation + Snow + geomorphon + Elevation")

# Survival hypothesis -- precip, snow, wetland
p_tenuis_models[['red_mod_survival']] = as.formula("dsl_mb ~ Precipitation + Snow + geomorphon")

# No survival -- large-scale and small-scale spatial variation

p_tenuis_models[['red_mod_spatial']] = as.formula("dsl_mb ~ Easting + Northing + f(i, model = spde)")

# No survival -- with linear elevation also

p_tenuis_models[['red_mod_spatial_elev']] = as.formula("dsl_mb ~ Easting + Northing + f(i, model = spde) + Elevation")

# P tenuis Null model -------------------------------

null_model = inla(formula = p_tenuis_models$null_model, family = 'nbinomial', data = data, control.compute = list(waic = TRUE, cpo  = TRUE))

null_model_ls = list(
  "data_stack" = NA,
  "joint_stack" = NA,
  "model"      = null_model
)

save(null_model_ls, file = "model_outputs/ptenuis/null_model.Rdata")

# P tenuis Full model with spatial effect -------------------------------

# Set up stacks
{


  # ggplot() +
  #   gg(data_sp) +
  #   gg(adk_mesh) +
  #   gg(adkbound) +
  #   coord_fixed(ratio = 1)


  stk <- inla.stack(
    data = list(dsl_mb = data$dsl_mb),
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
    data = list(dsl_mb = NA),
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

full_model = inla(formula = p_tenuis_models$full_model,
                  family = 'nbinomial',
                  data = inla.stack.data(stk),
                  control.compute = list(waic = TRUE, cpo  = TRUE),
                  control.predictor = list(A = inla.stack.A(stk))
)

full_model_ls = list(
  "data_stack" = stk,
  "joint_stack" = stk_jp,
  "model" = full_model
)

save(full_model_ls, file = "model_outputs/ptenuis/full_model.Rdata")

# P tenuis Full model with elevation random walk -----------------------------------------------

# Set up stack --

{
  stk <- inla.stack(
    data = list(dsl_mb = data$dsl_mb),
    A = list(projector_A, 1),
    effects = list(
      list(i = 1:spde$n.spde),
      data.frame(
        Intercept = 1,
        Easting = data$Easting,
        Elevation = inla.group(data$Elevation, method = 'quantile', n = 30),
        Northing = data$Northing,
        Precipitation = data$Precipitation,
        Snow = data$Snow,
        geomorphon = data$geomorphon
      )
    ),
    tag = 'dat'
  )

  stk_predictor = inla.stack(
    data = list(dsl_mb = NA),
    A = list(predictor_A, 1),
    effects = list(
      list(i = 1:spde$n.spde),
      data.frame(
        Intercept = 1,
        Easting = predict_grid_scaled@data$Easting,
        Elevation = inla.group(predict_grid_scaled$Elevation, method = 'quantile', n = 30),
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

full_model_elev_rw = inla(formula = p_tenuis_models$full_model_elev_rw,
                          family = 'nbinomial',
                          data = inla.stack.data(stk),
                          control.compute = list(waic = TRUE, cpo  = TRUE),
                          control.predictor = list(A = inla.stack.A(stk))
)

full_model_elev_rw_ls = list(
  "data_stack" = stk,
  "joint_stack" = stk_jp,
  "model"      = full_model_elev_rw
)

save(full_model_elev_rw_ls, file = "model_outputs/ptenuis/full_model_elev_rw_ls.Rdata")

# P tenuis Reduced: Remove spatial effect, all covariates and rw elevation --------------

# set up stacks
{
  stk <- inla.stack(
    data = list(dsl_mb = data$dsl_mb),
    A = list(1),
    effects = list(
      data.frame(
        Intercept = 1,
        Easting = data$Easting,
        Elevation = inla.group(data$Elevation, method = 'quantile', n = 30),
        Northing = data$Northing,
        Precipitation = data$Precipitation,
        Snow = data$Snow,
        geomorphon = data$geomorphon
      )
    ),
    tag = 'dat'
  )

  stk_predictor = inla.stack(
    data = list(dsl_mb = NA),
    A = list(1),
    effects = list(
      data.frame(
        Intercept = 1,
        Easting = predict_grid_scaled@data$Easting,
        Elevation = inla.group(predict_grid_scaled$Elevation, method = 'quantile', n = 30),
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

red_mod_minus_speff = inla(formula = p_tenuis_models$red_mod_minus_speff,
                           family = 'nbinomial',
                           data = inla.stack.data(stk),
                           control.compute = list(waic = TRUE, cpo  = TRUE),
                           control.predictor = list(A = inla.stack.A(stk))
)

red_mod_minus_speff_ls = list(
  "data_stack" = stk,
  "joint_stack" = stk_jp,
  "model"      = red_mod_minus_speff
)

save(red_mod_minus_speff_ls, file = "model_outputs/ptenuis/red_mod_minus_speff_ls.Rdata")

# P tenuis Reduced: Remove elevation rw; linear in all covariates ------------------

# No stacks necessary

red_mod_linear_all_cov = inla(formula = p_tenuis_models$red_mod_linear_all_cov,
                              family = 'nbinomial',
                              data = data,
                              control.compute = list(waic = TRUE, cpo  = TRUE)
)

red_mod_linear_all_cov_ls = list(
  "data_stack" = NA,
  "join_stack" = NA,
  "model"      = red_mod_linear_all_cov
)

save(red_mod_linear_all_cov_ls, file = 'model_outputs/ptenuis/red_mod_linear_all_cov_ls.Rdata')


# P tenuis Reduced: Survival hypothesis -- precip, snow, wetland -------------------

red_mod_survival = inla(formula = p_tenuis_models$red_mod_survival,
                        family = 'nbinomial',
                        data = data,
                        control.compute = list(waic = TRUE, cpo  = TRUE)
)

red_mod_survival_ls = list(
  "data_stack" = NA,
  "join_stack" = NA,
  "model"      = red_mod_survival
)

save(red_mod_survival_ls, file = 'model_outputs/ptenuis/red_mod_survival_ls.Rdata')

# P tenuis Reduced: no survival -- large-scale and small-scale spatial vari --------

# Set up stacks
{

  stk <- inla.stack(
    data = list(dsl_mb = data$dsl_mb),
    A = list(projector_A, 1),
    effects = list(
      list(i = 1:spde$n.spde),
      data.frame(
        Intercept = 1,
        Easting = data$Easting,
        Northing = data$Northing
      )
    ),
    tag = 'dat'
  )

  stk_predictor = inla.stack(
    data = list(dsl_mb = NA),
    A = list(predictor_A, 1),
    effects = list(
      list(i = 1:spde$n.spde),
      data.frame(
        Intercept = 1,
        Easting = predict_grid_scaled@data$Easting,
        Northing = predict_grid_scaled@data$Northing
      )
    ),
    tag = 'pred'
  )

  stk_jp = inla.stack(stk, stk_predictor)

}

red_mod_spatial = inla(formula = p_tenuis_models$red_mod_spatial,
                       family = 'nbinomial',
                       data = inla.stack.data(stk),
                       control.compute = list(waic = TRUE, cpo  = TRUE),
                       control.predictor = list(A = inla.stack.A(stk))
)

red_mod_spatial_ls = list(
  "data_stack" = stk,
  "joint_stack" = stk_jp,
  "model" = red_mod_spatial
)

save(red_mod_spatial_ls, file = "model_outputs/ptenuis/red_mod_spatial.Rdata")

# P tenuis Reduced: No survival -- with linear elevation also ----------------------


# Set up stacks
{

  stk <- inla.stack(
    data = list(dsl_mb = data$dsl_mb),
    A = list(projector_A, 1),
    effects = list(
      list(i = 1:spde$n.spde),
      data.frame(
        Intercept = 1,
        Easting = data$Easting,
        Northing = data$Northing,
        Elevation = data$Elevation
      )
    ),
    tag = 'dat'
  )

  stk_predictor = inla.stack(
    data = list(dsl_mb = NA),
    A = list(predictor_A, 1),
    effects = list(
      list(i = 1:spde$n.spde),
      data.frame(
        Intercept = 1,
        Easting = predict_grid_scaled@data$Easting,
        Northing = predict_grid_scaled@data$Northing,
        Elevation = predict_grid_scaled@data$Elevation
      )
    ),
    tag = 'pred'
  )

  stk_jp = inla.stack(stk, stk_predictor)

}

red_mod_spatial_elev = inla(formula = p_tenuis_models$red_mod_spatial_elev,
                                 family = 'nbinomial',
                                 data = inla.stack.data(stk),
                                 control.compute = list(waic = TRUE, cpo  = TRUE),
                                 control.predictor = list(A = inla.stack.A(stk))
)

red_mod_spatial_elev_ls = list(
  "data_stack" = stk,
  "joint_stack" = stk_jp,
  "model" = red_mod_spatial_elev
)

save(red_mod_spatial_elev_ls, file = "model_outputs/ptenuis/red_mod_spatial_elev.Rdata")

# Test HDS Model --------------------------------------------------------------------

# Keep mesh model from previous.

# Need transect lines, ds data

deer_transects_2018 = rgdal::readOGR(dsn = '../../GIS/DistanceSampling2018/ds_transects_2018.shp')

ds_data = read.csv('ds_data_adjusted.csv')

ds_data_sp = ds_data

coordinates(ds_data_sp) = ~Easting + Northing

proj4string(ds_data_sp) = proj4string(predict_grid)

ds_data_sp$distance = ds_data_sp$PERP_DIST_M

# What is the strip half-width? Set to a little larger than max distance

W = ceiling(max(ds_data$PERP_DIST_M)) # 3 meters

# Define half-normal detection function

hn = function(distance, lsig){
  exp(-0.5*(distance/exp(lsig))^2)}

# Define matern SPDE function for deer scats

matern <- inla.spde2.pcmatern(adk_mesh,
                              prior.sigma = c(2, 0.01),
                              prior.range = c(1000, 0.5))

# Define components of SPDE model

cmp = ~ mySPDE(map = coordinates, model = matern) +
  lsig + Intercept

formula = coordinates + distance ~ mySPDE +
  log(hn(distance, lsig)) +
  log(1/W) +
  Intercept

fit = lgcp(components = cmp,
           data = ds_data_sp,
           samplers = deer_transects_2018,
           formula = formula)


spde.range <- spde.posterior(fit, "mySPDE", what = "range"); plot(spde.range)
spde.logvar <- spde.posterior(fit, "mySPDE", what = "log.variance"); plot(spde.logvar)


# Predict over range

pxl = pixels(adk_mesh, nx = 200, ny = 200)
pr.int <- predict(fit, pxl, ~ exp(mySPDE))

ggplot() + gg(pr.int) + gg(adkbound) +
  gg(deer_transects_2018, color = "red") +
  gg(ds_data_sp, size = 0.2, alpha = 1) +
  noyticks + noxticks  +
  theme(legend.key.width = unit(x = 0.2,"cm"), legend.key.height = unit(x = 0.3,"cm")) +
  theme(legend.text=element_text(size=6)) +
  # guides(fill=FALSE) +
  coord_equal()


# View distance function
distdf <- data.frame(distance = seq(0,8,length=100))
dfun <- predict(fit, distdf, ~ hn(distance,lsig))
plot(dfun)


# Setup to fit covariate models -------------------------------------------

# Null model is meaningless here, since it will be multiplied against parasite models. Will be using
# covariate models only.

# Combine in with ds_data

ds_data_scaled = ds_data_sp

ds_data_scaled@data = cbind.data.frame(ds_data_sp@data,
                                       over(x = ds_data_sp, y = predict_grid_scaled))

ds_data_sp@data = cbind.data.frame(ds_data_sp@data,
                                   over(x = ds_data_sp, y = predict_grid))


