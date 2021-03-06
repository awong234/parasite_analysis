# Analysis util functions ----------

effects_extr = function(model_list, model_names, var_order) {
  
  if(var_order[1] != "(Intercept)") var_order = c("(Intercept)", var_order)
  
  fixed_eff = Map(f = function(x, nam){
    # browser()
    out = x$model$summary.fixed
    out$name = nam
    out$variable = rownames(out) %>% gsub(pattern = "geomorphon", replacement = "")
    out$variable = factor(out$variable, levels = rev(var_order) %>% gsub(pattern = "geomorphon", replacement = ""))
    rownames(out) = NULL
    return(out)
    },
    x = model_list,
    nam = model_names
    )
  # browser()
  fixed_eff = do.call(what = rbind, args = fixed_eff)
  rownames(fixed_eff) = NULL
  
  ranef = Map(f = function(x, nam){
    out = x$model$summary.hyperpar
    out$name = nam
    out$variable = rownames(out)
    rownames(out) = NULL
    return(out)
    },
    x = model_list,
    nam = model_names
    )
  
  ranef = do.call(what = rbind, args = ranef)
  rownames(ranef) = NULL
  
  out = list("Fixed" = fixed_eff,
             "Hyperpars" = ranef)
  
  return(out)
  
}

# Aids in extracting cpo from models

models_cpo = function(model_names, model_list){
  cpo_out = matrix(data = NA, nrow = nrow(data), ncol = nrow(model_list))
  names(cpo_out) = model_names
  for( name in model_names ) {
    mod_to_load = models_list[grep(models_list, pattern = paste0("\\/",name, "\\."), ignore.case = T)]
    if(length(mod_to_load) == 0){next}
    mod = readRDS(mod_to_load)
    mod_cpo = mod$cpo$cpo
    cpo_out[[name]] = mod_cpo
  }
  
  return(cpo_out)
}

# Aids in extracting waic from models located on disk.

hds_model_waic = function(hds_models_name, models_list) {
  aic_list = vector(mode = 'list', length = length(hds_models_name))
  names(aic_list) = hds_models_name
  for(name in hds_models_name){
    mod_to_load = models_list[grep(models_list, pattern = paste0("\\/",name, "\\."), ignore.case = T)]
    if(length(mod_to_load) == 0){next}
    mod = readRDS(mod_to_load)
    mod_waic = mod$waic$waic
    aic_list[[name]] = mod_waic
  }
  return(aic_list)
}

# aids in extracting waic from loaded models

aicFunc = function(model_ls){
  
  item = model_ls[[3]]
  
  waic = item$waic$`waic`
  
  return(waic)
  
}

# aids in extracting cpo from loaded models

cpoFunc = function(model_ls, name){
  
  item = model_ls[[3]]
  
  out1 = item$cpo$cpo
  out2 = item$cpo$pit
  
  out = data.frame("CPO" = out1, "PIT" = out2, "model" = name, stringsAsFactors = F)
  
  return(out)
  
}




# Analysis functions ------------



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

# Utility -------------------------------

ggmeshcust = function (data, color = NULL, alpha = NULL, edge.color = "grey", 
                       interior = TRUE, exterior = TRUE, ext.color = "black", 
                       crs = NULL, mask = NULL, nx = 500, ny = 500, ...) 
{
  if (!is.null(color)) {
    px = pixels(data, nx = nx, ny = ny)
    A = INLA::inla.spde.make.A(data, px)
    px$color = as.vector(A %*% color)
    if (!is.null(alpha)) {
      px$alpha = as.vector(A %*% alpha)
      gg = gg(px, mask = mask, alpha = "alpha")
    }
    else {
      gg = gg(px, mask = mask)
    }
    return(gg)
  }
  else {
    if (data$manifold == "S2") {
      stop("Geom not implemented for spherical meshes (manifold = S2)")
    }
    if (!is.null(crs)) {
      data = INLA::inla.spTransform(data, CRSobj = crs)
    }
    df = rbind(data.frame(a = data$loc[data$graph$tv[, 1], 
                                       c(1, 2)], b = data$loc[data$graph$tv[, 2], c(1, 
                                                                                    2)]), data.frame(a = data$loc[data$graph$tv[, 2], 
                                                                                                                  c(1, 2)], b = data$loc[data$graph$tv[, 3], c(1, 
                                                                                                                                                               2)]), data.frame(a = data$loc[data$graph$tv[, 1], 
                                                                                                                                                                                             c(1, 2)], b = data$loc[data$graph$tv[, 3], c(1, 
                                                                                                                                                                                                                                          2)]))
    colnames(df) = c("x", "y", "xend", "yend")
    mp = aes_string(x = "x", y = "y", xend = "xend", yend = "yend")
    msh = geom_segment(data = df, mapping = mp, color = edge.color)
    if (exterior) {
      df = data.frame(data$loc[data$segm$bnd$idx[, 1], 
                               1:2], data$loc[data$segm$bnd$idx[, 2], 1:2])
      colnames(df) = c("x", "y", "xend", "yend")
      bnd = geom_segment(data = df, mapping = mp, color = ext.color)
    }
    else {
      bnd = NULL
    }
    if (interior) {
      df = data.frame(data$loc[data$segm$int$idx[, 1], 
                               1:2], data$loc[data$segm$int$idx[, 2], 1:2])
      colnames(df) = c("x", "y", "xend", "yend")
      if (nrow(df) == 0) {
        int = NULL
      }
      else {
        int = geom_segment(data = df, mapping = aes_string(x = "x", y = "y", xend = "xend", yend = "yend", color = '"interior"'))
      }
    }
    else {
      int = NULL
    }
    c(msh, bnd, int)
  }
}
