library(ggplot2)
library(unmarked)
library(dplyr)

metadata = read.csv(file = 'metadata_adjusted.csv', stringsAsFactors = F)

ff = read.csv(file = 'fecal_flukefinder.csv', stringsAsFactors = F)
mb = read.csv(file = 'fecal_MB.csv', stringsAsFactors = F)
quant = read.csv(file = 'fecal_quant.csv', stringsAsFactors = F)

metadata_join = metadata %>% left_join(ff %>% select(PK, Total_eggs)) %>% rename(fmagna_ff = Total_eggs) %>% 
  left_join(mb %>% select(PK, Total_larvae_dsl, total_eggs_fmagna)) %>% rename(dsl_mb = Total_larvae_dsl, fmagna_mb = total_eggs_fmagna) %>% 
  left_join(quant %>% select(PK, Fascioloides_magna, Protostrongylid_DSL)) %>% rename(fmagna_quant = Fascioloides_magna, dsl_quant = Protostrongylid_DSL)

# Comparison of methods ---------------------------------------------------

# What is the sensitivity of the tests -- i.e. what is the probability of detection for the samples, given that the parasite occurs within?

metadata_join$fmagna_ff_y = metadata_join$fmagna_ff > 0
metadata_join$fmagna_mb_y = metadata_join$fmagna_mb > 0
metadata_join$fmagna_quant_y = metadata_join$fmagna_quant > 0

metadata_join$dsl_mb_y = metadata_join$dsl_mb > 0
metadata_join$dsl_quant_y = metadata_join$dsl_quant > 0

## Occupancy analysis

## Assess f magna

### Covariates

conditions = factor(metadata_join$Condition %>% na.omit)

years = factor(metadata_join$Year[!is.na(metadata_join$Condition)])

site_covs = data.frame(Condition = conditions, Year = years)

obs_covs = list(Method = matrix(data = c("FF", "MB", "QUANT"), nrow = length(metadata_join$Condition %>% na.omit), ncol = 3, byrow = T))

### Observations

yy_fmagna = (metadata_join %>% select(fmagna_ff_y : fmagna_quant_y))[!is.na(metadata_join$Condition),]

## Put data together

data = unmarkedFrameOccu(y = yy_fmagna, siteCovs = site_covs, obsCovs = obs_covs)

# Assess

out = list()

out[['null']]                            = occu(data = data, formula = ~1 ~1)
out[["year_null"]]                       = occu(data = data, formula = ~1 ~Year)
out[['null_year']]                       = occu(data = data, formula = ~Year ~1)
out[['year_condition']]                  = occu(data = data, formula = ~Condition ~Year)
out[['null_condition']]                  = occu(data = data, formula = ~Condition ~1)
out[['null_condition+method']]           = occu(data = data, formula = ~Condition + Method ~1)
out[['year_condition+method']]           = occu(data = data, formula = ~Condition + Method ~Year)
out[['null_method']]                     = occu(data = data, formula = ~Method ~1)
out[['year_method']]                     = occu(data = data, formula = ~Method ~Year)
out[['condition_null']]                  = occu(data = data, formula = ~1 ~Condition)
out[['condition_condition']]             = occu(data = data, formula = ~Condition ~Condition)
out[['condition_method']]                = occu(data = data, formula = ~Method ~Condition)
out[['condition_method+condition']]      = occu(data = data, formula = ~Method + Condition ~Condition)
out[['condition_method+condition+year']]      = occu(data = data, formula = ~Method + Condition + Year ~Condition)

models = names(out)

AIC_vals = lapply(X = out, FUN = function(x){x@AIC}) %>% do.call(rbind, args = .) %>% as.data.frame %>% cbind.data.frame(., models) %>% rename(AIC = V1)

rownames(AIC_vals) = NULL

(AIC_vals = AIC_vals %>% arrange(AIC))

fits = fitList(fits = out)
selection = modSel(fits)

save(out, file = 'occ_fmagna_1127.Rdata')

with(yy_fmagna, expr = fmagna_ff_y == T &
       (fmagna_mb_y == F | !is.na(fmagna_mb_y)) &
       (fmagna_quant_y == F | !is.na(fmagna_quant_y))) %>% sum(na.rm = T)


# Assess model fit

parboot(out$null, statistic=chisq, nsim=100, parallel=FALSE)

parboots = lapply(out, FUN = function(x){parboot(x, statistic = chisq, nsim = 10000, parallel = T)})

# Models with null components appear to be poor fits to the data.

## Assess dsl

### Observations

yy_dsl = (metadata_join %>% select(dsl_mb_y, dsl_quant_y))[!is.na(metadata_join$Condition),]

### Covariates

conditions = factor(metadata_join$Condition %>% na.omit)

years = factor(metadata_join$Year[!is.na(metadata_join$Condition)])

site_covs = data.frame(Condition = conditions, Year = years)

obs_covs = list(Method = matrix(data = c("MB", "QUANT"), nrow = length(metadata_join$Condition %>% na.omit), ncol = 2, byrow = T))

data = unmarkedFrameOccu(y = yy_dsl, siteCovs = site_covs, obsCovs = obs_covs)

# Assess

out = list()

out[['null']]                            = occu(data = data, formula = ~1 ~1)
out[["year_null"]]                       = occu(data = data, formula = ~1 ~Year)
out[['null_year']]                       = occu(data = data, formula = ~Year ~1)
out[['year_condition']]                  = occu(data = data, formula = ~Condition ~Year)
out[['null_condition']]                  = occu(data = data, formula = ~Condition ~1)
out[['null_condition+method']]           = occu(data = data, formula = ~Condition + Method ~1)
out[['year_condition+method']]           = occu(data = data, formula = ~Condition + Method ~Year)
out[['null_method']]                     = occu(data = data, formula = ~Method ~1)
out[['year_method']]                     = occu(data = data, formula = ~Method ~Year)
out[['condition_null']]                  = occu(data = data, formula = ~1 ~Condition)
out[['condition_condition']]             = occu(data = data, formula = ~Condition ~Condition)
out[['condition_method']]                = occu(data = data, formula = ~Method ~Condition)
out[['condition_method+condition']]      = occu(data = data, formula = ~Method + Condition ~Condition)
out[['condition_method+condition+year']] = occu(data = data, formula = ~Method + Condition + Year ~Condition)


models = names(out)

AIC_vals = lapply(X = out, FUN = function(x){x@AIC}) %>% do.call(rbind, args = .) %>% as.data.frame %>% cbind.data.frame(., models) %>% rename(AIC = V1)

rownames(AIC_vals) = NULL

(AIC_vals = AIC_vals %>% arrange(AIC))

fits = fitList(fits = out)
selection = modSel(fits)

save(out, file = 'occ_dsl_1127.Rdata')

backTransform(linearComb(out$`year_condition+method`, type = 'state', coefficients = c(1,0,0)))
backTransform(linearComb(out$`year_condition+method`, type = 'state', coefficients = c(1,1,0)))
backTransform(linearComb(out$`year_condition+method`, type = 'state', coefficients = c(1,0,1)))

metadata_join %>% filter(Year == 2016) %>% select(dsl_mb_y, dsl_quant_y) %>% rowSums() %>% sum /
metadata_join %>% filter(Year == 2016) %>% nrow

