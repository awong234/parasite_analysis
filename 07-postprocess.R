# Setup -----------
source('functions.R')

library(ggplot2)
library(dplyr)


# Summaries to refresh waic -----------

hds_models = readr::read_csv('hds_models.csv')

models_list = list.files(path = 'model_outputs/hds/', full.names = T)

hds_models_name = hds_models$Name

waic_ls = hds_model_waic(hds_models_name = hds_models_name, models_list = models_list)

waic_vec = do.call(what = c, args = waic_ls)

waic_df = data.frame(model = names(waic_vec), waic = waic_vec, row.names = NULL) %>% arrange(waic)

waic_df

# Obtain parasite distribution models ----------------------------------


# Summaries f magna -------------------------------------

m_out = list.files('model_outputs/fmagna/', full.names = T, pattern = 'Rdata')

for(m in m_out){
  load(m)
}

mod_list = list(
  'null_model_ls'             = null_model_ls,
  'full_model_ls'             = full_model_ls,
  'full_model_elev_rw_ls'     = full_model_elev_rw_ls,
  'red_mod_spatial_ls'        = red_mod_spatial_ls,
  # 'red_mod_minus_speff_ls'    = red_mod_minus_speff_ls,
  'red_mod_linear_all_cov_ls' = red_mod_linear_all_cov_ls,
  'red_mod_survival_ls'       = red_mod_survival_ls
)

# WAIC

aic_vals_fmagna = sapply(mod_list, aicFunc)

aic_vals_fmagna = data.frame("model" = names(mod_list), "waic" = aic_vals_fmagna) %>% arrange(waic)

# Get all the cpo and pit values for charts ------------

cpo_vals_fmagna = Map(f = function(mod, nam) {
  cpoFunc(model_ls = mod, name = nam)
  }, mod = mod_list, nam = names(mod_list)
)

cpo_vals_fmagna = do.call(what = rbind, args = cpo_vals_fmagna)
rownames(cpo_vals_fmagna) = NULL

cpo_vals_fmagna = reshape2::melt(cpo_vals_fmagna, id.vars = "model")

ggplot(cpo_vals_fmagna) + 
  geom_density(aes(x = value, fill = model), alpha = 0.3) + 
  facet_wrap(~variable) + 
  scale_fill_viridis_d(direction = -1) + 
  theme_bw()

rm(list = mod_list %>% names) ; gc()  

# Compare against fmagna values

combos = cpo_vals_fmagna %>% select(model, variable) %>% unique %>% nrow

cpo_vals_fmagna %>% mutate(response = rep(data$fmagna_ff, times = combos)) %>% 
  ggplot() + 
  geom_point(aes(x = response, y = value), shape = 1, alpha = 0.5) + 
  facet_grid(model ~ variable) + 
  theme_bw()

# Summaries p tenuis models' waic -------------------------------------------------

m_out = list.files('model_outputs/ptenuis/', full.names = T, pattern = 'Rdata')

for(m in m_out){
  load(m)
}

mod_list = list(
  'null_model_ls'             = null_model_ls,
  'full_model_ls'             = full_model_ls,
  # 'full_model_elev_rw_ls'     = full_model_elev_rw_ls,
  'red_mod_spatial_ls'        = red_mod_spatial_ls,
  'red_mod_minus_speff_ls'    = red_mod_minus_speff_ls,
  'red_mod_linear_all_cov_ls' = red_mod_linear_all_cov_ls,
  'red_mod_survival_ls'       = red_mod_survival_ls
)

# WAIC

aic_vals_ptenuis = sapply(mod_list, aicFunc)

aic_vals_ptenuis = data.frame("model" = names(mod_list), "waic" = aic_vals_ptenuis) %>% arrange(waic)

# Get all the cpo and pit values for charts ------------

cpo_vals_ptenuis = Map(f = function(mod, nam) {
  cpoFunc(model_ls = mod, name = nam)
}, mod = mod_list, nam = names(mod_list)
)

cpo_vals_ptenuis = do.call(what = rbind, args = cpo_vals_ptenuis)
rownames(cpo_vals_ptenuis) = NULL

rm(list = mod_list %>% names) ; gc()

cpo_vals_ptenuis = reshape2::melt(cpo_vals_ptenuis, id.vars = "model")

ggplot(cpo_vals_ptenuis) + 
  geom_density(aes(x = value, fill = model), alpha = 0.3) + 
  facet_wrap(~variable) + 
  scale_fill_viridis_d(direction = -1) + 
  theme_bw()
