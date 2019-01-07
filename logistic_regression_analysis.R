# Setup

library(magrittr)
library(dplyr)
library(ggplot2)
library(MASS)
library(GGally)

select = dplyr::select

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

# Exploration -----------------------------------------------------------------------------

data = bind_cols(metadata_join %>% select(-Easting, -Northing), covariate_scaled_df)

data %<>% mutate(JulianDay = lubridate::yday(Date))

missing_data_rows = is.na(data$Easting) | is.na(data$Northing)

my_bin <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_bin2d(...) +
    viridis::scale_fill_viridis(option = "D") + 
    theme_bw() + 
    theme(axis.text = element_blank(),
          panel.background = element_rect(fill = 'gray20'),
          panel.grid = element_blank())
}

center = function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) + 
    geom_density(...) + 
    theme_bw() + 
    theme(axis.text = element_blank())
}

pairs_plot = data %>% filter(!missing_data_rows) %>% dplyr::select(JulianDay, Year, Easting, Northing, Precipitation:Elevation, fmagna_ff, dsl_mb) %>% 
  ggpairs(lower = list(continuous = my_bin), diag = list(continuous = center)) + theme(strip.text = element_text(size = 5))

Cairo::Cairo(width = 1000, height = 1000, file = 'images/pairs_plot.png', dpi = 150)
pairs_plot
dev.off()


# Poisson regression ----------------------------------------------------------------------

# F magna

full_model_poisson_fmagna = glm(formula = fmagna_ff ~ Easting + Northing + Year + Precipitation + Snow + Distance_to_wetland + Elevation, family = "poisson", data = data)

full_model_nb_fmagna = glm.nb(formula = fmagna_ff ~ Easting + Northing + Year + Precipitation + Snow + Distance_to_wetland + Elevation, data = data)

# P tenuis

full_model_poisson_dsl = glm(formula = dsl_mb ~ Easting + Northing + Year + Precipitation + Snow + Distance_to_wetland + Elevation, family = "poisson", data = data)

full_model_nb_dsl
