# Setup -----------
source('functions.R')

library(ggplot2)
library(dplyr) ; `%<>%` = magrittr::`%<>%`
library(reshape2)



load('data_no_na.Rdata')

model_matches = data.frame(fitnames = c("null_model", 
                                        "full_model", 
                                        "red_mod_linear_all_cov", 
                                        "red_mod_survival",
                                        "red_mod_spatial", 
                                        "red_mod_spatial_elev"),
                           pubnames = c("Null Model", 
                                        "Full model", 
                                        "Fixed effects model", 
                                        "Fixed Effects Parasite Survival Model", 
                                        "Two-scale spatial model", 
                                        "Two-scale spatial model and elevation"),
                           stringsAsFactors = F)

# Summaries f magna -------------------------------------

m_out = list.files('model_outputs/fmagna/', full.names = T, pattern = 'Rdata')

for(m in m_out){
  load(m)
}

mod_list = list(
  'null_model'             = null_model_ls,
  'full_model'             = full_model_ls,
  # 'full_model_elev_rw'     = full_model_elev_rw_ls,
  'red_mod_spatial'        = red_mod_spatial_ls,
  'red_mod_spatial_elev'   = red_mod_spatial_elev_ls,
  # 'red_mod_minus_speff'    = red_mod_minus_speff_ls,
  'red_mod_linear_all_cov' = red_mod_linear_all_cov_ls,
  'red_mod_survival'       = red_mod_survival_ls
)

# WAIC

aic_vals_fmagna = sapply(mod_list, aicFunc)

aic_vals_fmagna = data.frame("model" = names(mod_list), "waic" = aic_vals_fmagna, stringsAsFactors = F) %>% arrange(waic)

save(aic_vals_fmagna, file = 'model_outputs/aicvals_fmagna.Rdata')

# Plot fixed effects -----------------------------------

varorder =  rownames(full_model_ls$model$summary.fixed)

fmagna_effects = effects_extr(mod_list, model_names = names(mod_list), var_order = varorder)

fmagna_effects$Fixed %<>% left_join(model_matches, by = c("name" = "fitnames"))

save(fmagna_effects, file = 'model_outputs/fmagna_effects.Rdata')

Cairo::Cairo(width = 650, height = 1600, file = 'images/parplot_fmagna.pdf', dpi = 110, type = 'pdf')
ggplot(fmagna_effects$Fixed) + 
  geom_hline(yintercept = 0) + 
  geom_point(aes(x = variable, y = mode, color = pubnames, shape = pubnames), position = position_dodge(0.3)) + 
  geom_errorbar(aes(x = variable, ymin = `0.025quant`, ymax = `0.975quant`, color = pubnames, linetype = pubnames), width = 0, position = position_dodge(0.3)) + 
  # facet_wrap(facets = ~name, nrow = length(names(mod_list))) + 
  scale_color_viridis_d(option = 'D') + 
  # guides(linetype = guide_legend(title = "Model"), shape = guide_legend(title = "Model"),  color = guide_legend(title = "Model")) +
  # coord_cartesian(ylim = c(-5, 5)) + 
  coord_flip(ylim = c(-5, 5)) +
  theme_bw() + theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  )
dev.off()


# Get all the cpo and pit values for charts ------------

cpo_vals_fmagna = Map(f = function(mod, nam) {
  cpoFunc(model_ls = mod, name = nam)
  }, mod = mod_list, nam = names(mod_list)
)

cpo_vals_fmagna = do.call(what = rbind, args = cpo_vals_fmagna) %>% left_join(model_matches, by = c("model" = "fitnames"))

rownames(cpo_vals_fmagna) = NULL

cpo_vals_fmagna = reshape2::melt(cpo_vals_fmagna, id.vars = "pubnames", measure.vars = c("cpo", "pit"))

Cairo::Cairo(width = 1024*2, height = 768*2, file = "images/cpo_fmagna.pdf", type = "pdf", dpi = 120)
ggplot(cpo_vals_fmagna) + 
  geom_density(aes(x = value, fill = pubnames), alpha = 0.3) + 
  facet_wrap(~variable) + 
  scale_fill_viridis_d(direction = -1) + 
  theme_bw()
dev.off()

# Compare against fmagna values

combos = cpo_vals_fmagna %>% select(pubnames, variable) %>% unique %>% nrow

Cairo::Cairo(width = 1024, height = 768*2, file = "images/cpo_by_response_fmagna.pdf", type = "pdf", dpi = 80)
cpo_vals_fmagna %>% mutate(response = rep(data$fmagna_ff, times = combos)) %>% 
  ggplot() + 
  geom_point(aes(x = response, y = value), shape = 1, alpha = 0.5) + 
  facet_grid(pubnames ~ variable) + 
  theme_bw()
dev.off()

rm(list = ls(pattern = '_ls'))

# Summaries p tenuis models' waic -------------------------------------------------

m_out = list.files('model_outputs/ptenuis/', full.names = T, pattern = 'Rdata')

for(m in m_out){
  load(m)
}

mod_list = list(
  'null_model'             = null_model_ls,
  'full_model'             = full_model_ls,
  # 'full_model_elev_rw'     = full_model_elev_rw_ls,
  'red_mod_spatial'        = red_mod_spatial_ls,
  'red_mod_spatial_elev'   = red_mod_spatial_elev_ls,
  # 'red_mod_minus_speff'    = red_mod_minus_speff_ls,
  'red_mod_linear_all_cov' = red_mod_linear_all_cov_ls,
  'red_mod_survival'       = red_mod_survival_ls
)

# WAIC

aic_vals_ptenuis = sapply(mod_list, aicFunc)

aic_vals_ptenuis = data.frame("model" = names(mod_list), "waic" = aic_vals_ptenuis) %>% arrange(waic)

save(aic_vals_ptenuis, file = 'model_outputs/aicvals_ptenuis.Rdata')

# Plot fixed effects -----------------------------------

ptenuis_effects = effects_extr(mod_list, model_names = names(mod_list), var_order = varorder)

ptenuis_effects$Fixed %<>% left_join(model_matches, by = c("name" = "fitnames"))

save(ptenuis_effects, file = 'model_outputs/ptenuis_effects.Rdata')

Cairo::Cairo(width = 1100, height = 1600, file = 'images/parplot_ptenuis.pdf', dpi = 110, type = 'pdf')
ggplot(ptenuis_effects$Fixed) + 
  geom_hline(yintercept = 0) + 
  geom_point(aes(x = variable, y = mode, color = pubnames, shape = pubnames), position = position_dodge(0.3)) + 
  geom_errorbar(aes(x = variable, ymin = `0.025quant`, ymax = `0.975quant`, color = pubnames, linetype = pubnames), width = 0, position = position_dodge(0.3)) + 
  # facet_wrap(facets = ~name, nrow = length(names(mod_list))) + 
  scale_color_viridis_d(option = 'D') + 
  guides(linetype = guide_legend(title = "Model"), shape = guide_legend(title = "Model"),  color = guide_legend(title = "Model")) +
  # coord_cartesian(ylim = c(-5, 5)) + 
  coord_flip(ylim = c(-5, 5)) +
  theme_bw() + theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )
dev.off()


# Get all the cpo and pit values for charts ------------

cpo_vals_ptenuis = Map(f = function(mod, nam) {
  cpoFunc(model_ls = mod, name = nam)
}, mod = mod_list, nam = names(mod_list)
)

cpo_vals_ptenuis = do.call(what = rbind, args = cpo_vals_ptenuis) %>% left_join(model_matches, by = c("model" = "fitnames"))
rownames(cpo_vals_ptenuis) = NULL

cpo_vals_ptenuis = melt(cpo_vals_ptenuis, id.vars = "pubnames", measure.vars = c("cpo", "pit"))

Cairo::Cairo(width = 1024*2, height = 768*2, file = "images/cpo_ptenuis.pdf", type = "pdf", dpi = 120)
ggplot(cpo_vals_ptenuis) + 
  geom_density(aes(x = value, fill = pubnames), alpha = 0.3) + 
  facet_wrap(~variable) + 
  scale_fill_viridis_d(direction = -1) + 
  theme_bw()
dev.off()

# Compare against ptenuis values

combos = cpo_vals_ptenuis %>% select(pubnames, variable) %>% unique %>% nrow

Cairo::Cairo(width = 1024, height = 768*2, file = "images/cpo_by_response_ptenuis.pdf", type = "pdf", dpi = 80)
cpo_vals_ptenuis %>% mutate(response = rep(data$dsl_mb, times = combos)) %>% 
  ggplot() + 
  geom_point(aes(x = response, y = value), shape = 1, alpha = 0.5) + 
  facet_grid(pubnames ~ variable) + 
  theme_bw()
dev.off()