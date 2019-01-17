# Setup ----------------------------------------

library(dplyr)
library(sp)

select = dplyr::select

metadata = read.csv(file = 'metadata_adjusted.csv', stringsAsFactors = F)
metadata$Date = as.Date(metadata$Date)


metadata_sp = metadata
na_location_index = is.na(metadata$Easting) | is.na(metadata$Northing)

# Easting and Northing ----------------------------------------------------------------

metadata_sp = metadata[!na_location_index,]

coordinates(metadata_sp) = ~Easting + Northing

ff = read.csv(file = 'fecal_flukefinder.csv', stringsAsFactors = F)
mb = read.csv(file = 'fecal_MB.csv', stringsAsFactors = F)
quant = read.csv(file = 'fecal_quant.csv', stringsAsFactors = F)

metadata_join = metadata %>% left_join(ff %>% select(PK, Total_eggs)) %>% rename(fmagna_ff = Total_eggs) %>% 
  left_join(mb %>% select(PK, Total_larvae_dsl, total_eggs_fmagna)) %>% rename(dsl_mb = Total_larvae_dsl, fmagna_mb = total_eggs_fmagna) %>% 
  left_join(quant %>% select(PK, Fascioloides_magna, Protostrongylid_DSL)) %>% rename(fmagna_quant = Fascioloides_magna, dsl_quant = Protostrongylid_DSL)

load('data_no_na.Rdata')

load('model_combos_fmagna.Rdata')
load('model_combos_poisson_fmagna.Rdata')
load('model_combos_dsl.Rdata')
load('model_combos_poisson_dsl.Rdata')

load('predict_grid_1000.Rdata')

# Prediction ------------------------------------

# Remove old predictors

predict_grid@data = predict_grid@data %>% select(Northing, Easting, Elevation)

# Add new predictors

# Add Precipitation --------------------------------------------------------

# Precipitation will be the mean of the three years.

precip_names = dir(path = '../../GIS/PRISM_precip_data/', pattern = 'PRISM.+tif$', full.names = F)
precip_files = dir(path = '../../GIS/PRISM_precip_data/', pattern = 'PRISM.+tif$', full.names = T)


precip_dates = data.frame(file = precip_names, year = precip_names %>% {regmatches(x = ., m = regexec(pattern = '\\d{4}', text = ., perl = T))} %>% as.integer,
                          month = precip_names %>% {regmatches(x = ., m = regexec(text = ., pattern = '(?!\\d{3})\\d{2}', perl = T))} %>% as.integer)

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

ggplot(predict_grid@data) +
  geom_raster(aes(x = Easting, y = Northing, fill = Snow))

# Perform prediction ------------------------------------------------------------------------------------------------

load('scaled_covariates_attr.Rdata')

predict_grid_complete = predict_grid_scaled = predict_grid[complete.cases(predict_grid@data),]

center = covariate_scaled_attr$Center[c(2,1,6,3,4,5)]
scale  = covariate_scaled_attr$Scale[c(2,1,6,3,4,5)]

temp = apply(X = predict_grid_complete@data, MARGIN = 1, FUN = function(x){x - center}) %>% t
predict_grid_scaled@data = apply(X = temp, MARGIN = 1, FUN = function(x){x / scale}) %>% t %>% as.data.frame

top_model = MuMIn::get.models(out_nb_fmagna, subset = delta == 0)
top_model_dsl = MuMIn::get.models(out_nb_dsl, subset = delta == 0)

# Predict fmagna
pr_out = predict(top_model$`255`, newdata = predict_grid_scaled@data, type = 'response')
# Predict dsl
pr_out_dsl = predict(top_model_dsl$`111`, newdata = predict_grid_scaled@data, type = 'response')

# Fmagna assign
predict_grid_scaled@data$PredictedIntensity = pr_out
predict_grid_complete@data$PredictedIntensity = pr_out

# DSL assign
predict_grid_scaled@data$PredictedIntensity_dsl = pr_out_dsl
predict_grid_complete@data$PredictedIntensity_dsl = pr_out_dsl

# Fmagna
Cairo::Cairo(width = 768*3, height = 1024*3, dpi = 150*2, file = 'images/intensity_fmagna_top_nb.png')
ggplot(predict_grid_complete@data) +
  geom_raster(aes(x = Easting, y = Northing, fill = PredictedIntensity)) + 
  viridis::scale_fill_viridis(option = "B") + 
  theme_bw() + 
  theme(panel.grid = element_blank()) + 
  coord_equal()
dev.off()

# DSL
Cairo::Cairo(width = 768*3, height = 1024*3, dpi = 150*2, file = 'images/intensity_dsl_top_nb.png')
ggplot(predict_grid_complete@data) +
  geom_raster(aes(x = Easting, y = Northing, fill = PredictedIntensity_dsl)) + 
  viridis::scale_fill_viridis(option = "B") + 
  theme_bw() + 
  theme(panel.grid = element_blank()) + 
  coord_equal()
dev.off()


# Relationships -----------------------------

# What's the relationship with elevation

sampled_elevation = data$Elevation * covariate_scaled_attr$Scale[6] + covariate_scaled_attr$Center[6]

new_dat = with(data, data.frame(Northing = mean(Northing), Easting = mean(Easting), 
                                                    Elevation = seq(min(Elevation), max(Elevation), length = 1000),
                                                    Precipitation = mean(Precipitation), Snow = mean(Snow),
                                                    Distance_to_wetland = mean(Distance_to_wetland)
                                                    )
)

predicted_intensity = predict(top_model$`255`, newdata = new_dat, type = 'response')


plot_data = data.frame(Elevation = with(data, seq(min(sampled_elevation), max(sampled_elevation), length = 1000)), Intensity = predicted_intensity)

plot_points = data.frame(Elevation = sampled_elevation,
                         Intensity = data$fmagna_ff)

Cairo::Cairo(file = 'images/elevation_rel_fmagna.png', dpi = 150)
ggplot(plot_data) + 
  geom_line(aes(x = Elevation, y = predicted_intensity)) + 
  geom_point(data = plot_points, aes(x = Elevation, y = Intensity), alpha = 0.05) + 
  theme_bw() + ylab("Intensity per Fecal Sample")
dev.off()

# What's the relationship with precipitation

new_dat = with(data, data.frame(Northing = mean(Northing), Easting = mean(Easting), 
                                                    Elevation = mean(Elevation),
                                                    Precipitation = seq(min(Precipitation), max(Precipitation), length.out = 1000), 
                                                    Snow = mean(Snow),
                                                    Distance_to_wetland = mean(Distance_to_wetland)
                                                    )
)

predicted_intensity = predict(top_model$`255`, newdata = new_dat, type = 'response')

sampled_precip = data$Precipitation * covariate_scaled_attr$Scale[3] + covariate_scaled_attr$Center[3]

plot_data = data.frame(Precipitation = with(data, seq(min(sampled_precip), max(sampled_precip), length = 1000)), `Predicted Intensity` = predicted_intensity)

plot_points = data.frame(Precipitation = sampled_precip,
                         Intensity = data$fmagna_ff)

Cairo::Cairo(file = 'images/precipitation_rel_fmagna.png', dpi = 150)
ggplot(plot_data) + 
  geom_line(aes(x = Precipitation, y = predicted_intensity)) + 
  geom_point(data = plot_points, aes(x = Precipitation, y = Intensity), alpha = 0.05) + 
  theme_bw() + ylab("Intensity per Fecal Sample")
dev.off()

# What's the relationship with snow

sampled_snow = data$Snow * covariate_scaled_attr$Scale[4] + covariate_scaled_attr$Center[4]

new_dat = with(data, data.frame(Northing = mean(Northing), Easting = mean(Easting), 
                                                    Elevation = mean(Elevation),
                                                    Precipitation = mean(Precipitation), 
                                                    Snow = seq(min(Snow), max(Snow), length = 1000),
                                                    Distance_to_wetland = mean(Distance_to_wetland)
                                                    ))

predicted_intensity = predict(top_model$`255`, newdata = new_dat, type = 'response')

plot_data = data.frame(Snow = with(data, seq(min(sampled_snow), max(sampled_snow), length = 1000)), `Predicted Intensity` = predicted_intensity)

plot_points = data.frame(Snow = sampled_snow,
                         Intensity = data$fmagna_ff)

Cairo::Cairo(file = 'images/snow_rel_fmagna.png', dpi = 150)
ggplot(plot_data) + 
  geom_line(aes(x = Snow, y = predicted_intensity)) + 
  geom_point(data = plot_points, aes(x = Snow, y = Intensity), alpha = 0.05) + 
  theme_bw() + ylab("Intensity per Fecal Sample")
dev.off()

# What's the relationship with Easting?

sampled_easting = data$Easting * covariate_scaled_attr$Scale[1] + covariate_scaled_attr$Center[1]

new_dat = with(data, data.frame(Northing = mean(Northing), Easting = seq(min(Easting), max(Easting), length = 1000), 
                                                    Elevation = mean(Elevation),
                                                    Precipitation = mean(Precipitation), 
                                                    Snow = mean(Snow),
                                                    Distance_to_wetland = mean(Distance_to_wetland)
                                                    ))

predicted_intensity = predict(top_model$`255`, newdata = new_dat, type = 'response')

plot_data = data.frame(Easting = with(data, seq(min(sampled_easting), max(sampled_easting), length = 1000)), `Predicted Intensity` = predicted_intensity)

plot_points = data.frame(Easting = sampled_easting,
                         Intensity = data$fmagna_ff)

Cairo::Cairo(file = 'images/easting_rel_fmagna.png', dpi = 150)
ggplot(plot_data) + 
  geom_line(aes(x = Easting, y = predicted_intensity)) + 
  geom_point(data = plot_points, aes(x = Easting, y = Intensity), alpha = 0.05) + 
  theme_bw() + ylab("Intensity per Fecal Sample")
dev.off()

# Northing?

sampled_northing = data$Northing * covariate_scaled_attr$Scale[2] + covariate_scaled_attr$Center[2]

new_dat = with(data, data.frame(Northing = seq(min(Northing), max(Northing), length = 1000), Easting = mean(Easting), 
                                Elevation = mean(Elevation),
                                Precipitation = mean(Precipitation), 
                                Snow = mean(Snow),
                                Distance_to_wetland = mean(Distance_to_wetland)
))

predicted_intensity = predict(top_model$`255`, newdata = new_dat, type = 'response')

plot_data = data.frame(Northing = with(data, seq(min(sampled_northing), max(sampled_northing), length = 1000)), `Predicted Intensity` = predicted_intensity)

plot_points = data.frame(Northing = sampled_northing,
                         Intensity = data$fmagna_ff)

Cairo::Cairo(file = 'images/northing_rel_fmagna.png', dpi = 150)
ggplot(plot_data) + 
  geom_line(aes(x = Northing, y = predicted_intensity)) + 
  geom_point(data = plot_points, aes(x = Northing, y = Intensity), alpha = 0.05) + 
  theme_bw() + ylab("Intensity per Fecal Sample")
dev.off()


# Mean prediction

new_dat = with(data, data.frame(Northing = mean(Northing), Easting = mean(Easting), 
                                                    Elevation = mean(Elevation),
                                                    Precipitation = mean(Precipitation), 
                                                    Snow = mean(Snow),
                                                    Distance_to_wetland = mean(Distance_to_wetland)
                                                    ))

predicted_intensity = predict(top_model$`255`, newdata = new_dat, type = 'response', se.fit = T)

