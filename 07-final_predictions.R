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

m_out = list.files('model_outputs/fmagna/', full.names = T)

for(m in m_out){
  load(m)
}

mod_list = list(
  'null_model_ls'             = null_model_ls,
  'full_model_ls'             = full_model_ls,
  'full_model_elev_rw_ls'     = full_model_elev_rw_ls,
  'red_mod_spatial_ls'        = red_mod_spatial_ls,
  'red_mod_minus_speff_ls'    = red_mod_minus_speff_ls,
  'red_mod_linear_all_cov_ls' = red_mod_linear_all_cov_ls,
  'red_mod_survival_ls'       = red_mod_survival_ls
)


aic_vals_fmagna = sapply(mod_list, aicFunc)

aic_vals_fmagna = data.frame("model" = names(mod_list), "waic" = aic_vals_fmagna) %>% arrange(waic)

rm(list = mod_list %>% names) ; gc()

aic_vals_fmagna

# Summaries p tenuis models' waic -------------------------------------------------

m_out = list.files('model_outputs/ptenuis/', full.names = T)

for(m in m_out){
  load(m)
}

mod_list = list(
  'null_model_ls'             = null_model_ls,
  'full_model_ls'             = full_model_ls,
  'full_model_elev_rw_ls'     = full_model_elev_rw_ls,
  'red_mod_spatial_ls'        = red_mod_spatial_ls,
  'red_mod_minus_speff_ls'    = red_mod_minus_speff_ls,
  'red_mod_linear_all_cov_ls' = red_mod_linear_all_cov_ls,
  'red_mod_survival_ls'       = red_mod_survival_ls
)


aic_vals_ptenuis = sapply(mod_list, aicFunc)

aic_vals_ptenuis = data.frame("model" = names(mod_list), "waic" = aic_vals_ptenuis) %>% arrange(waic)

rm(list = mod_list %>% names) ; gc()

aic_vals_ptenuis

# Select top models and multiply posterior mean together ---------------

# Note that for the parasite models prediction has already been done. In the stack, the 'pred' tag
# should extract the mean predictions for the model.