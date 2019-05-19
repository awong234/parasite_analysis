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

# Select top models and multiply posterior dist together ---------------