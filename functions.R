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

getPrecipFromDate = function(dates, precip_data, precip_dates){
  
  data_month = lubridate::month(dates)
  data_year  = lubridate::year(dates)
  data_dates = data.frame(year = data_year, month = data_month)
  
  n = nrow(data_dates)
  
  data_out = rep(NA_real_, n)
  
  for(i in 1:n){
    index = which(data_dates[i,'year'] == precip_dates$year & data_dates[i,'month'] - 1 == precip_dates$month)
    tryCatch(expr = {data_out[i]  = precip_data[i,index]},
             error = function(m){data_out[i] = NA; return(data_out[i])}
    )
  }
  
  return(data_out)
  
}

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