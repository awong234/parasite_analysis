# Analysis functions ------------


inla_model_fn = function(data, 
                         formulas, 
                         prior.range = c(1000, 0.01), prior.sigma = c(1, 0.01), 
                         par = FALSE, par_cores = 4,
                         prediction){
  
  
  # Classify models -- spde? rw?
  # browser()
  rw_present = sapply(X = formulas, FUN = function(x){grepl(pattern = 'rw', x = x[3], perl = T)})
  spde_present = sapply(X = formulas, FUN = function(x){grepl(pattern = 'spde', x = x[3], perl = T)})
  
  spde_rw =  rw_present &   spde_present
  spde    = !rw_present &   spde_present
  rw      =  rw_present &  !spde_present
  clear   = !rw_present &  !spde_present
  
  model_class = bind_rows(spde_rw, spde, rw, clear)
  model_class = apply(X = model_class, MARGIN = 2, FUN = function(x){which(x == 1)})
  # browser()
  # Set up mesh, spde
  
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
  
  spde <- inla.spde2.pcmatern(
    mesh = adk_mesh,
    alpha = 2,
    ### mesh and smoothness parameter
    prior.range = prior.range,
    ### P(practic.range<1000)=0.01 -- conservative, states that there is a very small chance that there is autocorrelation at ranges below 1000m. 
    prior.sigma = prior.sigma
    ### P(sigma>1)=0.01 -- conservative, states that there is a very small chance that the standard deviation is larger than 1.
  ) 
  
  projector_A <- inla.spde.make.A(adk_mesh, loc=data %>% select(Easting_real, Northing_real) %>% as.matrix())

  browser()
  # Classify stacks for each situation ------------------------------------------------------------
  {
  # SPDE_RW ---------------------------------------------------
  
  stk_spde_rw = inla.stack(
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
        Distance_to_wetland = data$Distance_to_wetland
      )
    ),
    tag = 'dat'
  )
    
    
  # Predictor stack
  
  predictor_A = inla.spde.make.A(adk_mesh, loc = prediction@coords) # Use predict_grid_scaled generated near the top of the HDS.R file
  
  stk_predictor_spde_rw = inla.stack(
    data = list(y = NA),
    A = list(predictor_A, 1),
    effects = list(
      list(i = 1:spde$n.spde),
      data.frame(
        Intercept = 1,
        Easting = prediction@data$Easting,
        Elevation = inla.group(prediction$Elevation, method = 'quantile', n = 30),
        Northing = prediction@data$Northing,
        Precipitation = prediction@data$Precipitation,
        Snow = prediction@data$Snow,
        Distance_to_wetland = prediction@data$Distance_to_wetland
      )
    ),
    tag = 'pred'
  )
  
  stk_jp_spde_rw = inla.stack(stk_spde_rw, stk_predictor_spde_rw)
  
  # SPDE ------------------------------------------------
  
  stk_spde = inla.stack(
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
        Distance_to_wetland = data$Distance_to_wetland
      )
    ),
    tag = 'dat'
  )
  
  
  predictor_A = inla.spde.make.A(adk_mesh, loc = prediction@coords) # Use predict_grid_scaled generated near the top of the HDS.R file
  
  stk_predictor_spde = inla.stack(
    data = list(y = NA),
    A = list(predictor_A, 1),
    effects = list(
      list(i = 1:spde$n.spde),
      data.frame(
        Intercept = 1,
        Easting = prediction@data$Easting,
        Elevation = prediction$Elevation,
        Northing = prediction@data$Northing,
        Precipitation = prediction@data$Precipitation,
        Snow = prediction@data$Snow,
        Distance_to_wetland = prediction@data$Distance_to_wetland
        
      )
    ),
    tag = 'pred'
  )
  
  stk_jp_spde = inla.stack(stk_spde, stk_predictor_spde)
  
  # RW ----------------------------------------------------
  
  stk_rw = inla.stack(
    data = list(y = data$fmagna_ff),
    A = list(1),
    effects = list(
      list(data.frame(
        Intercept = 1,
        Easting = data$Easting,
        Elevation = inla.group(data$Elevation, method = 'quantile', n = 30),
        Northing = data$Northing,
        Precipitation = data$Precipitation,
        Snow = data$Snow,
        Distance_to_wetland = data$Distance_to_wetland
      )
      ),
      tag = 'dat'
    )
  )
  
  
  predictor_A = inla.spde.make.A(adk_mesh, loc = prediction@coords) # Use predict_grid_scaled generated near the top of the HDS.R file
  
  stk_predictor_rw = inla.stack(
    data = list(y = NA),
    A = list(1),
    effects = list(
      data.frame(
        Intercept = 1,
        Easting = prediction@data$Easting,
        Elevation = inla.group(prediction$Elevation, method = 'quantile', n = 30),
        Northing = prediction@data$Northing,
        Precipitation = prediction@data$Precipitation,
        Snow = prediction@data$Snow,
        Distance_to_wetland = prediction@data$Distance_to_wetland
        
      )
    ),
    tag = 'pred'
  )
  
  stk_jp_rw = inla.stack(stk_rw, stk_predictor_rw)
  
  }
  
  # Run through each model. 
  # Recall model_class construction   
  # 1 spde_rw =  rw_present &   spde_present
  # 2 spde    = !rw_present &   spde_present
  # 3 rw      =  rw_present &  !spde_present
  # 4 clear   = !rw_present &  !spde_present
  
  # Stacks available are
  # stk_jp_spde_rw
  # stk_jp_spde
  # stk_jp_rw
  # none
  
  for(i in 1:length(formulas)){
    browser()
    
    type = model_class[i] %>% as.character()
    
    stack_object = switch(type,
                          "1" = stk_jp_spde_rw,
                          "2" = stk_jp_spde,
                          "3" = stk_jp_rw,
                          "4" = NULL)
    
    str(stack_object)
    
  }
#   
#   if(par){doParallel::registerDoParallel(cores = par_cores)}
#   
#   foreach(i = 1:length(formulas)) %dopar% {
#     
#     
#     
#     mod_out = inla(formula = formula[[i]], 
#                    family  = 'nbinomial',
#                    data    = data,
#                    control.compute = list(waic=TRUE))
#     
#   }
#     
}



getAngle = function(start_pt, end_pt){

  # Obtain an angle from two points

  # Translation

  translation = end_pt - start_pt
  
  dx = translation[1]
  
  # Figure out quadrant based off translation

  signs = sign(translation)
  
  # How to calculate depends on signs.
  
  if(signs[2] >= 0){
    
    dist_to_end = fields::rdist(matrix(start_pt, nrow = 1), matrix(end_pt, nrow = 1)) %>% as.numeric()
    
    fake_east_pt = start_pt + c(1, 0)
    
    dist_to_fake = fields::rdist(matrix(end_pt, nrow = 1), matrix(fake_east_pt, nrow = 1)) %>% as.numeric()
    
    # Use cosign law
    
    angle = acos((1 + dist_to_end^2 - dist_to_fake^2) / (2 * 1 * dist_to_end))
    
    return(angle * (180 / pi))
    
  } else {
    
    dist_to_end = fields::rdist(matrix(start_pt, nrow = 1), matrix(end_pt, nrow = 1)) %>% as.numeric()
    
    fake_east_pt = start_pt + c(1, 0)
    
    dist_to_fake = fields::rdist(matrix(end_pt, nrow = 1), matrix(fake_east_pt, nrow = 1)) %>% as.numeric()
    
    # Use cosign law
    
    angle = acos((1 + dist_to_end^2 - dist_to_fake^2) / (2 * 1 * dist_to_end))
    
    return(180 - angle * (180 / pi) + 180)
    
    
  }

}

adjustPoint = function(tsct_start, angle, parallel_dist, perp_dist){
  
  # NOTE: IT IS NOT KNOWN WHETHER THE SCATS WERE ENCOUNTERED TO THE LEFT OR
  # RIGHT OF THE TRANSECTS. THEREFORE, WE MUST RANDOMIZE THE PERPENDICULAR
  # DISTANCE TO BE NEGATIVE OR POSITIVE.
  
  rad_angle = angle * pi / 180
  
  perp_dist = perp_dist * sample(x = c(-1, 1), size = 1)
  
  untransformed_new_pt = c(parallel_dist, perp_dist)
  
  transformed_new_pt = matrix(c(cos(rad_angle), sin(rad_angle), -sin(rad_angle), cos(rad_angle)), nrow = 2, ncol = 2) %*% untransformed_new_pt
  
  # Now translate to tsct_start 
  
  transformed_new_pt = transformed_new_pt + tsct_start
  
  return(transformed_new_pt)
  
}

chisq <- function(fm) {
  umf <- getData(fm)
  y <- getY(umf)
  y[y>1] <- 1
  sr <- fm@sitesRemoved
  if(length(sr)>0)
    y <- y[-sr,,drop=FALSE]
  fv <- fitted(fm, na.rm=TRUE)
  y[is.na(fv)] <- NA
  sum((y-fv)^2/(fv*(1-fv)), na.rm=TRUE)
}

# Data formatting ---------------

getSnowFromDate = function(dates, snow_data, snow_dates){
  
  date_year = lubridate::year(dates)
  
  n = length(date_year)
  
  data_out = rep(NA_real_, n)
  
  for(i in 1:n){
    
    index = which(date_year[i] - 1 == snow_dates$year)
    tryCatch(expr = {data_out[i] = snow_data[i,index]},
             error = function(m){data_out[i] = NA; return(data_out[i])})
    
  }
  
  return(data_out)
  
}

# Deprecated ------------------

# getPrecipFromDate = function(dates, precip_data, precip_dates){
#   
#   data_month = lubridate::month(dates)
#   data_year  = lubridate::year(dates)
#   data_dates = data.frame(year = data_year, month = data_month)
#   
#   n = nrow(data_dates)
#   
#   data_out = rep(NA_real_, n)
#   
#   for(i in 1:n){
#     index = which(data_dates[i,'year'] == precip_dates$year & data_dates[i,'month'] - 1 == precip_dates$month)
#     tryCatch(expr = {data_out[i]  = precip_data[i,index]},
#              error = function(m){data_out[i] = NA; return(data_out[i])}
#     )
#   }
#   
#   return(data_out)
#   
# }
