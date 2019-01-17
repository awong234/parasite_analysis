# Setup -------------------------------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(rgdal)

source('functions.R')

# Clean up metadata  ------------------------------------------------------------

metadata = read.csv(file = 'metadata.csv', stringsAsFactors = F)

metadata %>% head
metadata %>% str

metadata$Date = as.Date(metadata$Date, format = '%m/%d/%Y')

# Identify missing location values -- not removing just yet because they may yet be useful for site-specific occupancy.

missingLocRows = with(data = metadata, expr = {(is.na(Easting) | is.na(Northing))})

# Standardize fecal sample quality to integer. 

metadata$Condition[metadata$Condition == "Old" & !is.na(metadata$Condition)] = 3
metadata$Condition[metadata$Condition == "Fresh" & !is.na(metadata$Condition)] = 2
metadata$Condition[grepl(metadata$Condition, pattern = 'Very \\wresh')] = 1
metadata$Condition[metadata$Condition == '2.5' & !is.na(metadata$Condition)] = 2
metadata$Condition[metadata$Condition == 'Semi-fresh' & !is.na(metadata$Condition)] = 2
metadata$Condition[metadata$Condition == '' & !is.na(metadata$Condition)] = NA

metadata$Condition = as.integer(metadata$Condition)

# Verify locations - no large deviations

transects = readOGR(dsn = 'DeerTransects2018/DeerTransects2018.shp', layer = 'DeerTransects2018', stringsAsFactors = F)
transects@coords = transects@coords[,c(1,2)]
attr(x = transects@coords, which = 'dimnames') = list(NULL, c("Easting", "Northing"))

# Convert transect points to lines

transects_lines = Map(f = function(x,y){sp::Lines(slinelist = sp::Line(x),ID = y)},
                      x = split(x = transects, f = transects@data$ident),
                      y = {transects@data$ident %>% unique}
)

transects_lines = SpatialLines(transects_lines, proj4string = CRS(proj4string(transects)))

# Export

#rgdal::writeOGR(obj = transects_lines, )

# Get distances

metadata_sp = metadata[!missingLocRows,]

coordinates(metadata_sp) = ~Easting + Northing
proj4string(metadata_sp) = proj4string(transects)

dist_to_scat = rgeos::gDistance(
  metadata_sp[metadata_sp@data$Year == 2018,],
  transects_lines,
  byid = T
  )


dist_to_scat = dist_to_scat %>% apply(X = ., MARGIN = 2, FUN = sort)

# Adjusting 2018 fecal locations --------------------------------------------------------------------------------------------------------------------------------------------

# Deer fecal locations in 2016 and 2018 can potentially be improved by using the transect start location, and displacement along the transect and off of the transect.
# This would require gathering the transect angles (should be easy from access db), and calculating a vector by which to transform the original starting point.
# Assumes *better* knowledge of starting points than collection points. This is likely reasonable in 2018 because of averaging performed during site survey.

# Actually this MUST be done for the 2018 data, because INLABRU uses the location of the points supplied, and predicts intensity off of that, can't just feed it distances I think.

# Load in transect points

ds_transects_2018 = read.csv(file = 'distance_sampling_2018_transects.csv', stringsAsFactors = F) %>% select(-X)

start_points = ds_transects_2018 %>% select(StartEasting, StartNorthing) %>% rename(Easting = StartEasting, Northing = StartNorthing)
end_points = ds_transects_2018 %>% select(EndEasting, EndNorthing) %>% rename(Easting = EndEasting, Northing = EndNorthing)

# Generate angles

angles = vector(mode = 'list', length = nrow(start_points))

for(i in 1:nrow(start_points)){
  
  angles[[i]] = getAngle(start_pt = start_points[i,], end_pt = end_points[i,])
  
}

angles = do.call(what = c, args = angles)

ds_transects_2018 = ds_transects_2018 %>% mutate(Angles = angles)

# Load in distance sampling encounter data

ds_obs_2018 = read.csv(file = 'distance_sampling_2018_observations.csv', stringsAsFactors = F)

# Compare to metadata

ds_compare = ds_obs_2018 %>% left_join(metadata %>% select(PK, Specimen_ID, Grid, Date), by = c("SCAT_ID" = "Specimen_ID")) %>% select(PK, TSCT, Date, SCAT_ID, Grid)

ds_compare[ds_compare$TSCT != ds_compare$Grid,] # Now there are no discrepancies.

# Re-specify encounter coordinates by displacing from start using parallel length, perpendicular length, and angle. Adjust metadata.

ds_data_join = ds_obs_2018 %>% left_join(ds_transects_2018)                 # adds start coordinates recorded for transect, and measured angles
ds_data_join = ds_data_join %>% left_join(metadata %>% select(Specimen_ID, Easting, Northing), by = c("SCAT_ID" = "Specimen_ID")) # adds RECORDED Easting and Northing for samples. This is only used inthe next part

# Ensure first that all collection locations are labeled correctly -- scat should not be too far from labeled transect.

ds_data_sp = ds_data_join
coordinates(ds_data_sp) = ~Easting + Northing
proj4string(ds_data_sp) = proj4string(transects)

dist_to_scat = rgeos::gDistance(
  ds_data_sp,
  transects_lines,
  byid = T
)

attr(dist_to_scat, which = 'dimnames') = NULL

dist_to_scat = dist_to_scat %>% data.frame %>% mutate(Grid = transects_lines@lines %>% names)

mins = apply(X = dist_to_scat[,-136], MARGIN = 2, FUN = which.min)

dist_to_scat$Grid[mins]

test = cbind(ds_data_join, "DerGrid" = dist_to_scat$Grid[mins])

test[which(test$TSCT != test$DerGrid),] # Discrepancy in grid location is anomalous -- the recorded data is 7B5-2. Will fix in readjustment of locations.

# DO THE ADJUSTMENT!!!!

ds_data_join = ds_data_join %>% mutate(PERP_DIST_M = PERP_DIST / 100)

adjustedLoc = vector(mode = 'list', length = nrow(ds_data_join))

for(i in 1:nrow(ds_data_join)){
  
  adjustedLoc[[i]] = adjustPoint(tsct_start = ds_data_join %>% select(StartEasting, StartNorthing) %>% slice(i), 
                                 angle = ds_data_join$Angles[i], 
                                 parallel_dist = ds_data_join$PARALLEL_DIST[i],
                                 perp_dist = ds_data_join$PERP_DIST_M[i]
  )
  
}

adjustedLoc = do.call(what = rbind, args = adjustedLoc)

ds_data_join = ds_data_join %>% mutate(NewEasting = adjustedLoc[,1], NewNorthing = adjustedLoc[,2])

ds_data_join %>% head

ds_data_out = ds_data_join %>% select(TSCT:PARALLEL_DIST, PERP_DIST_M, APPROX_N:EndTime, Obs1:Angles, NewEasting, NewNorthing)

ds_data_out = ds_data_out %>% rename(Easting = NewEasting, Northing = NewNorthing)

write.csv(ds_data_out, file = 'ds_data_adjusted.csv')

# Should move the end coordinates 100m from start using the same angles -- but only if INLAbru doesn't like the mismatch in locations. 

# Update metadata with new coordinates

metadata_new = metadata

metadata_new[metadata_new$Specimen_ID %in% ds_data_out$SCAT_ID,]$Easting = ds_data_out$Easting
metadata_new[metadata_new$Specimen_ID %in% ds_data_out$SCAT_ID,]$Northing = ds_data_out$Northing

write.csv(metadata_new, file = 'metadata_adjusted.csv', row.names = F)
