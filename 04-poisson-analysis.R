# Setup

library(magrittr)
library(ggplot2)
library(MASS)
library(GGally)
library(MuMIn)
library(dplyr)

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

data = bind_cols(metadata_join %>% rename(Easting_real = Easting, Northing_real = Northing), covariate_scaled_df)

data %<>% mutate(JulianDay = lubridate::yday(Date) %>% scale() %>% c)

data$Year = factor(data$Year)

# missing_data_rows = is.na(data$Easting) | is.na(data$Northing)

data = data[complete.cases(data %>% select(fmagna_ff, dsl_mb, Easting:JulianDay)),]

save(data, file = 'data_no_na.Rdata')

my_bin <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_bin2d(...) +
    geom_smooth(..., color = 'red') +
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

pairs_plot = data %>% dplyr::select(Easting, Northing, Precipitation:Elevation, fmagna_ff, dsl_mb) %>% 
  ggpairs(lower = list(continuous = my_bin), diag = list(continuous = center)) + theme(strip.text = element_text(size = 5))

Cairo::Cairo(width = 1000, height = 1000, file = 'images/pairs_plot.png', dpi = 150)
pairs_plot
dev.off()

# What's that relationship between fmagna and dsl ?

data %>% select(fmagna_ff, dsl_mb) %>% plot
({data %>% select(fmagna_ff, dsl_mb)} + 0.001) %>% log %>% plot



# Poisson regression ----------------------------------------------------------------------

# F magna

full_model_poisson_fmagna = glm(formula = fmagna_ff ~ Easting + Northing + Year + Precipitation + I(Precipitation^2) + Snow + I(Snow^2) + geomorphon_40 + Elevation+ I(Elevation^2) + JulianDay + I(JulianDay^2), family = "poisson", data = data, na.action = 'na.fail')

full_model_nb_fmagna = glm.nb(formula = fmagna_ff ~ Easting + Northing + Year + Precipitation + I(Precipitation^2) + Snow + I(Snow^2) + geomorphon_40 + Elevation + I(Elevation^2) + JulianDay + I(JulianDay^2), data = data, na.action = 'na.fail')

# Absolutely no support for poisson model
AIC(full_model_poisson_fmagna, full_model_nb_fmagna)
pchisq(-2 * (as.numeric(stats::logLik(full_model_poisson_fmagna)) - as.numeric(stats::logLik(full_model_nb_fmagna))), df = 1, lower.tail = F)

summary(full_model_nb_fmagna)

out_nb_fmagna = dredge(global.model = full_model_nb_fmagna, 
                       subset = dc(Elevation, I(Elevation^2), Precipitation, I(Precipitation^2), Snow, I(Snow^2), JulianDay, I(JulianDay^2)))

out_poisson_fmagna = dredge(global.model = full_model_poisson_fmagna, 
                            subset = dc(Elevation, I(Elevation^2), Precipitation, I(Precipitation^2), Snow, I(Snow^2), JulianDay, I(JulianDay^2))) 

save(out_nb_fmagna, file = 'model_combos_fmagna.Rdata')
save(out_poisson_fmagna, file = 'model_combos_poisson_fmagna.Rdata')

# P tenuis

full_model_poisson_dsl = glm(formula = dsl_mb ~ Easting + Northing + Year + Precipitation + I(Precipitation^2) + Snow + I(Snow^2) + geomorphon_40 + Elevation+ I(Elevation^2) + JulianDay + I(JulianDay^2), family = "poisson", data = data, na.action = 'na.fail')

full_model_nb_dsl = glm.nb(formula = dsl_mb ~ Easting + Northing + Year + Precipitation + I(Precipitation^2) + Snow + I(Snow^2) + geomorphon_40 + Elevation + I(Elevation^2) + JulianDay + I(JulianDay^2), data = data, na.action = 'na.fail', init.theta = 0.35)

out_nb_dsl = dredge(global.model = full_model_nb_dsl, subset = dc(Elevation, I(Elevation^2), Precipitation, I(Precipitation^2), Snow, I(Snow^2), JulianDay, I(JulianDay^2)))
out_poisson_dsl = dredge(global.model = full_model_poisson_dsl, subset = dc(Elevation, I(Elevation^2), Precipitation, I(Precipitation^2), Snow, I(Snow^2), JulianDay, I(JulianDay^2)))

save(out_nb_dsl, file = 'model_combos_dsl.Rdata')
save(out_poisson_dsl, file = 'model_combos_poisson_dsl.Rdata')
