
# Setup ----------------------------------------

library(INLA)
library(inlabru)
library(magrittr)
library(ggplot2)
library(MASS)
library(GGally)
library(MuMIn)
library(viridis)
library(dplyr)

source('functions.R')
load(file = 'scaled_covariates.Rdata')
load(file = 'scaled_covariates_attr.Rdata')

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

# Load maxlik outputs ------------------------------------------------------------------------------------------------------------

load('model_combos_dsl.Rdata')
load('model_combos_fmagna.Rdata')

# Fit INLA model same effects as last ---------------------------------------

# Best model was Elevation^2, Northing, Precipitation^2, Snow.

inla_mod = fmagna_ff ~ Easting + Elevation + I(Elevation*Elevation) + Northing + Precipitation + I(Precipitation*Precipitation) + Snow + JulianDay

inla_out = inla(formula = inla_mod, family = "nbinomial", data = data, control.compute = list(waic=TRUE))

summary(inla_out) # All of the estimates are identical.

# Try fitting same model with bru() -- yes!

cmp_test = fmagna_ff ~ Easting + Elevation + ElevQuad(map = Elevation^2, model = 'linear') + Northing + Precipitation + PrecipQuad(map = Precipitation^2, model = 'linear') + Snow + JulianDay
bru_out = bru(components = cmp_test, family = 'nbinomial', data = data)

summary(bru_out)

bru_out$summary.fixed$mean

# Fit model with spatial random effect

# Create mesh

adkbound = rgdal::readOGR(dsn = '../../GIS/adkParkPoly/adkParkPoly.shp')

# Data points as sp

data_sp = data
coordinates(data_sp) = ~Easting_real + Northing_real

# Look at variogram
# brk = seq(0, 10000, by = 1)
# vgram_out = fields::vgram(loc = data_sp@coords, y = data_sp@data$fmagna_ff, breaks = brk)

# No evidence of spatial correlation.
# plot(vgram_out, breaks = seq(0,2000, length.out = 50))

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

ggplot() +
  gg(data_sp) +
  gg(adk_mesh) +
  gg(adkbound) +
  coord_fixed(ratio = 1)

spde <- inla.spde2.pcmatern(
  mesh = adk_mesh,
  alpha = 2,
  ### mesh and smoothness parameter
  prior.range = c(1000, 0.01),
  ### P(practic.range<0.3)=0.5
  prior.sigma = c(1, 0.01)
  ### P(sigma>1)=0.01
) 

projector_A <- inla.spde.make.A(adk_mesh, loc=data %>% select(Easting_real, Northing_real) %>% as.matrix())

stk <- inla.stack(
  data = list(y = data$fmagna_ff),
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
      JulianDay = data$JulianDay
    )
  ),
  tag = 'dat'
)

# Load prediction grid --------------------------------------------------------------

load('predict_grid_1000.Rdata')

predict_grid@data = predict_grid@data %>% select(Northing, Easting, Elevation)

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

center = covariate_scaled_attr$Center[c(2,1,6,3,4,5)]
scale  = covariate_scaled_attr$Scale[c(2,1,6,3,4,5)]

temp = apply(X = predict_grid_complete@data, MARGIN = 1, FUN = function(x){x - center}) %>% t
temp = apply(X = temp, MARGIN = 1, FUN = function(x){x / scale}) %>% t %>% as.data.frame

predict_grid_scaled = predict_grid_complete
predict_grid_scaled@data = temp

# Stack predict grid projection into inla stack --------------------

predictor_A = inla.spde.make.A(adk_mesh, loc = predict_grid_scaled@coords)

stk_predictor = inla.stack(
  data = list(y = NA),
  A = list(predictor_A, 1),
  effects = list(
    list(i = 1:spde$n.spde),
    data.frame(
      Intercept = 1,
      Easting = predict_grid_scaled@data$Easting,
      Elevation = predict_grid_scaled@data$Elevation,
      Northing = predict_grid_scaled@data$Northing,
      Precipitation = predict_grid_scaled@data$Precipitation,
      Snow = predict_grid_scaled@data$Snow
    )
  ),
  tag = 'pred'
)

stk_jp = inla.stack(stk, stk_predictor)

# Run models -------------------

# Examples

fmagna_spatial_mod = y ~ 0 + Intercept + Easting + Elevation + I(Elevation^2) + Northing + Precipitation + I(Precipitation^2) + Snow + JulianDay + f(i, model = spde)

fmagna_spatial_out = inla(fmagna_spatial_mod, family = 'nbinomial', 
                          control.compute = list(waic = TRUE),
                          data = inla.stack.data(stk),
                          control.predictor = list(A = inla.stack.A(stk))
                          )

summary(fmagna_spatial_out)

# Predict on grid using posterior mode

fmagna_spatial_out_md = inla(fmagna_spatial_mod, family = 'nbinomial',
                             data = inla.stack.data(stk_jp),
                             control.predictor = list(A = inla.stack.A(stk_jp),
                                                      compute = TRUE,
                                                      link = 1),
                             quantiles = NULL,
                             control.inla = list(strategy = 'gaussian'),
                             control.results = list(return.marginals.random = F,
                                                    return.marginals.predictor=F),
                             control.mode = list(theta = fmagna_spatial_out$mode$theta, restart = F)
)

summary(fmagna_spatial_out_md)

pred_index = inla.stack.index(stk_jp, tag = 'pred')$data
pred_response = fmagna_spatial_out_md$summary.fitted.values$mean[pred_index]

predict_grid_complete@data$INLA_pred_fmagna = pred_response

summary(pred_response)

outlier_index = pred_response > 1e20


ggplot(predict_grid_complete@data) + 
  geom_raster(aes(x = Easting, y = Northing, fill = INLA_pred_fmagna)) +
  scale_fill_viridis(option = 'C', limits = c(0,1000)) + 
  gg(adkbound) + coord_fixed(ratio = 1)
  # geom_raster(aes(x = Easting, y = Northing, fill = Elevation))


# Example with random walk elevation effect -------------------------------

stk <- inla.stack(
  data = list(y = data$fmagna_ff),
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
      JulianDay = data$JulianDay
    )
  ),
  tag = 'dat'
)


fmagna_spatial_mod = y ~ 0 + Intercept + Easting + Northing + Precipitation + I(Precipitation^2) + Snow + JulianDay + f(i, model = spde) + f(Elevation, model = "rw1")

fmagna_spatial_out = inla(fmagna_spatial_mod, family = 'nbinomial', 
                          control.compute = list(waic = TRUE),
                          data = inla.stack.data(stk),
                          control.predictor = list(A = inla.stack.A(stk))
)

summary(fmagna_spatial_out)


# Predict on grid using posterior mode

stk_predictor = inla.stack(
  data = list(y = NA),
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
      JulianDay = predict_grid_scaled@data$JulianDay
    )
  ),
  tag = 'pred'
)

stk_jp = inla.stack(stk, stk_predictor)

fmagna_spatial_out_md = inla(fmagna_spatial_mod, family = 'nbinomial',
                             data = inla.stack.data(stk_jp),
                             control.predictor = list(A = inla.stack.A(stk_jp),
                                                      compute = TRUE,
                                                      link = 1),
                             quantiles = NULL,
                             control.inla = list(strategy = 'gaussian'),
                             control.results = list(return.marginals.random = F,
                                                    return.marginals.predictor=F),
                             control.mode = list(theta = fmagna_spatial_out$mode$theta, restart = F)
)

summary(fmagna_spatial_out_md)

pred_index = inla.stack.index(stk_jp, tag = 'pred')$data
pred_response = fmagna_spatial_out_md$summary.fitted.values$mean[pred_index]
pred_sd = fmagna_spatial_out_md$summary.fitted.values$sd[pred_index]

predict_grid_complete@data$INLA_pred_fmagna = pred_response
predict_grid_complete@data$INLA_sd = pred_sd

summary(pred_response)

outlier_index = pred_response > 1e20


ggplot(predict_grid_complete@data) + 
  geom_raster(aes(x = Easting, y = Northing, fill = INLA_pred_fmagna)) +
  scale_fill_viridis(option = 'C') + 
  gg(adkbound) + coord_fixed(ratio = 1)
# geom_raster(aes(x = Easting, y = Northing, fill = Elevation))

# List all models to run and run them! -----------------------------------------



# HDS Model --------------------------------------------------------------------

# Keep mesh model from previous. 

# Need transect lines, ds data

deer_transects_2018 = rgdal::readOGR(dsn = '../../GIS/DistanceSampling2018/ds_transects_2018.shp')

ds_data = read.csv('ds_data_adjusted.csv')

ds_data_sp = ds_data

coordinates(ds_data_sp) = ~Easting + Northing

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

# Very cool ---- now combine the two simply -----------------------

lik_parasite = like(family = 'nbinomial', 
     formula = fmagna_ff ~ Easting + Elevation + I(Elevation^2) + Northing + Precipitation + I(Precipitation^2) + Snow + JulianDay, 
     data = data)

cmp_deer = scat_dens ~ mySPDE(map = coordinates, model = matern) + 
  lsig + Intercept

lik_deer = like(family = 'cp',
                formula = coordinates + distance ~ mySPDE +
                  log(hn(distance, lsig)) + 
                  log(1/W) +
                  Intercept, data = ds_data_sp,
                components = cmp_deer,
                samplers = deer_transects_2018)

bru(cmp_deer, lik_deer)

jcmp = fmagna_ff ~ Easting + Elevation + Northing + Precipitation + Snow + JulianDay + Intercept + coordinates + distance + mySPDE(map = coordinates, model = matern)

test = bru(components = jcmp, lik_parasite, lik_deer)

# Practice ---------------

data(Seeds)

df = data.frame(y = Seeds$r, Ntrials = Seeds$n, Seeds[, 3:5])

family1 = "binomial"
control.family1 = list(control.link=list(model="logit"))
# number of trials is df$Ntrials

hyper1 = list(theta = list(prior="pc.prec", param=c(1,0.01)))
formula1 = y ~ x1 + x2 + f(plate, model="iid", hyper=hyper1)

q1 = INLA:::inla.rw1(n = 5)

crossprod(q1)

INLA:::inla.rw2(n = 5)

# Mesh example

s <- 3 ### this factor will only changes C, not G
pts <- rbind(c(1,1), c(2,1),
             c(2.6, 1), c(0.7,1.7), 4:5/3, c(2,1.7))*s
n <- nrow(pts)
mesh0 <- inla.mesh.2d(pts[1:2,], max.edge=3*s,
                      offset=1*s, n=6, cutoff=s*1/2)
mesh <- inla.mesh.2d(rbind(c(3.3,1)*s, c(2.4,2)*s,
                           mesh0$loc[-c(3:4),1:2]),
                     max.edge=3*s, offset=1e-5, cutoff=s*1/2, n=100)
(m <- mesh$n)
dmesh <- inla.mesh.dual(mesh)
fem <- inla.mesh.fem(mesh, order=1)
A <- inla.spde.make.A(mesh, pts)

# Geostats example

data("SPDEtoy")

str(SPDEtoy)

pl.dom <- cbind(c(0, 1, 1, 0.7, 0), c(0, 0, 0.7, 1, 1))

mesh5 <- inla.mesh.2d(pl.dom, max.e = c(0.092, 0.2))


spde5 <- inla.spde2.pcmatern(
  mesh = mesh5,
  alpha = 2,
  ### mesh and smoothness parameter
  prior.range = c(0.3, 0.5),
  ### P(practic.range<0.3)=0.5
  prior.sigma = c(1, 0.01)
) ### P(sigma>1)=0.01

coords <- as.matrix(SPDEtoy[,1:2])
A5 <- inla.spde.make.A(mesh5, loc=coords)

stk5 <- inla.stack(
  data = list(resp = SPDEtoy$y),
  A = list(A5, 1),
  effects = list(i = 1:spde5$n.spde,
                 m = rep(1, nrow(SPDEtoy))),
  tag = 'est'
)

res5 <- inla(
  resp ~ 0 + m + f(i, model = spde5),
  data = inla.stack.data(stk5),
  control.predictor = list(A = inla.stack.A(stk5))
)

res5.field <- inla.spde2.result(res5, 'i', spde5, do.transf=TRUE)


# Prediction during estimation
pts3 <- rbind(c(.1,.1), c(.5,.55), c(.7,.9))

dim(A5pts3 <- inla.spde.make.A(mesh5, loc=pts3))

(jj3 <- which(colSums(A5pts3)>0))

round(A5pts3[, jj3],3)

stk5p.rf <- inla.stack(
  data = list(resp = NA),
  A = list(A5pts3),
  effects = list(i = 1:spde5$n.spde),
  tag = 'prd5r'
)

stk5.jp <- inla.stack(stk5, stk5p.rf)

res5p <- inla(
  resp ~ 0 + m + f(i, model = spde5),
  data = inla.stack.data(stk5.jp),
  control.predictor = list(A = inla.stack.A(stk5.jp), compute = TRUE)
)

# Prediction after estimation

drop(A5pts3%*%res5$summary.random$i$mean)

inla.mesh.project(inla.mesh.projector(mesh5, loc = pts3),
                  res5$summary.random$i$mean)

# Project onto mesh

# project on grid with domain [0,1] x [0,1]
pgrid0 <- inla.mesh.projector(mesh5, xlim=0:1, ylim=0:1, dims=c(101,101))

# project posterior mean and standard deviation
prd0.m <- inla.mesh.project(pgrid0, res5$summary.ran$i$mean)
prd0.s <- inla.mesh.project(pgrid0, res5$summary.ran$i$s)


# Predict response
stk5.presp <- inla.stack(
  data = list(resp = NA),
  A = list(A5pts3, 1),
  effects = list(i = 1:spde5$n.spde, m = rep(1, 3)),
  tag = 'prd5.resp'
)

stk5.full <- inla.stack(stk5, stk5.presp)

r5presp <- inla(
  resp ~ 0 + m + f(i, model = spde5),
  data = inla.stack.data(stk5.full),
  control.predictor = list(A = inla.stack.A(stk5.full), compute = TRUE)
)

(indd3r <- inla.stack.index(stk5.full, 'prd5.resp')$data)

round(r5presp$summary.fitted.values[indd3r,], 3)

# Marginal distribution of point response
marg3r <- r5presp$marginals.fitted.values[indd3r]

# Response on a grid, disabling marginal distributions.
stkgrid <- inla.stack(
  data = list(resp = NA),
  A = list(pgrid0$proj$A, 1),
  effects = list(i = 1:spde5$n.spde,
                 m = rep(1, 101 * 101)),
  tag = 'prd.gr'
)

stk.all <- inla.stack(stk5, stkgrid)
res5g <- inla(
  resp ~ 0 + m + f(i, model = spde5),
  data = inla.stack.data(stk.all),
  control.predictor = list(A = inla.stack.A(stk.all),
                           compute = TRUE),
  quantiles = NULL,
  control.results = list(
    return.marginals.random = FALSE,
    return.marginals.predictor = FALSE
  )
)
res5g$cpu

# Indexes of prediction grid
igr <- inla.stack.index(stk.all, 'prd.gr')$data


