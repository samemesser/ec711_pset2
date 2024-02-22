# Samuel Messer
# EC711 Problem Set 2
# Worked closely on code with Erin Eidschun
library(tidyverse) # Data manipulation
library(readstata13) # To read in stata data file
library(xtable) # To make output convenient

# Clear environment
rm(list = ls())

# Read in and clean data
sipp_raw <- read.dta13("data/sipp1991.dta")

# Split income into 7 quantiles
break_points = as.numeric(quantile(sipp_raw$inc, probs = ((1:6) / 7) ))
sipp <- sipp_raw %>% 
  mutate(inc_bin = as.factor(findInterval(inc, break_points)),
         inc2 = inc^2, # Quadratic in income
         age2 = age^2, age3 = age^3, # Cubic in age
         fsize2 = fsize^2, # Quadratic in fsize
         educ2 = educ^2) # Quadratic in educ

# Generate low-dim control matrix
# Need this to calculate fitted values later
ld_control <- model.matrix(~ marr + twoearn + db + pira + hown +
                             fsize + fsize2 + educ + educ2 + age + age2 + age3 +
                             inc_bin + inc_bin:inc + inc_bin:inc2, data = sipp)

# Remove raw data
rm(sipp_raw)

### Eligibility as treatment

## Non-parametric 
# Covariate control formula
ctrl_form = net_tfa ~ marr + twoearn + db + pira + hown +
                      fsize + fsize2 + educ + educ2 + age + age2 + age3 +
                      inc_bin + inc_bin:inc + inc_bin:inc2

# Fit OLS model only on treated
elg_tr <- sipp[sipp$e401 == 1,]
mod_elg_tr <- lm(ctrl_form, data = elg_tr)

# Fit OLS model only on untreated
elg_untr <- sipp[sipp$e401 == 0,]
mod_elg_untr <- lm(ctrl_form, data = elg_untr)

# Calculate fitted Y vals
elg_coef = mod_elg_tr$coefficients - mod_elg_untr$coefficients
y_hat_np_elg = as.matrix(ld_control) %*% as.matrix(elg_coef)

ATE_np_elg = mean(y_hat_np_elg)

## Propensity score reweighting
mod_elg_ps = glm(e401 ~ marr + twoearn + db + pira + hown +
                 fsize + fsize2 + educ + educ2 + age + age2 + age3 +
                 inc_bin + inc_bin:inc + inc_bin:inc2, data = sipp, 
                 family = binomial(link = "logit"))

# Fitted propensity scores
ps = fitted(mod_elg_ps) 

# Check if common support holds
min(ps) > 0
max(ps) < 1

# Calculate ATE according to HIR
ATE_ps_elg = mean((sipp$net_tfa * sipp$e401 / ps) - (sipp$net_tfa * (1 - sipp$e401) / (1 - ps)))


## Doubly robust estimator
# We calculated the expectations for the adjustment already in the np estimator
# Treated
tr_fit = as.matrix(ld_control) %*% as.matrix(mod_elg_tr$coefficients)

# Untreated
untr_fit = as.matrix(ld_control) %*% as.matrix(mod_elg_untr$coefficients)

# Adjustment term in dr estimator
dr_elg_adj = mean((sipp$e401 * tr_fit) / ps - (1 - sipp$e401) * untr_fit / (1 - ps))
ATE_dr_elg = ATE_np_elg + ATE_ps_elg  - dr_elg_adj


print(xtable(data.frame(Estimator = c("Non-parametric", "Propensity Score", "Doubly Robust"), 
                        ATE = c(ATE_np_elg, ATE_ps_elg, ATE_dr_elg))), 
      include.rownames = F)

### Participation as treatment

## Non-parametric
# Fit OLS model only on treated
par_tr <- sipp[sipp$p401 == 1,]
mod_par_tr <- lm(ctrl_form, data = par_tr)

# Fit OLS model only on untreated
par_untr <- sipp[sipp$p401 == 0,]
mod_par_untr <- lm(ctrl_form, data = par_untr)

# Plug into framework to predict Y
par_coef = mod_par_tr$coefficients - mod_par_untr$coefficients
y_hat_np_par = as.matrix(ld_control) %*% as.matrix(par_coef)

ATE_np_par = mean(y_hat_np_par)

## Propensity score reweighting
mod_par_ps = glm(p401 ~ marr + twoearn + db + pira + hown +
                   fsize + fsize2 + educ + educ2 + age + age2 + age3 +
                   inc_bin + inc_bin:inc + inc_bin:inc2, data = sipp, 
                 family = binomial(link = "logit"))

# Fitted propensity scores
ps = fitted(mod_par_ps) 

# Check if common support holds
min(ps) > 0
max(ps) < 1

# Calculate ATE according to HIR
ATE_ps_par = mean((sipp$net_tfa * sipp$p401 / ps) - (sipp$net_tfa * (1 - sipp$p401) / (1 - ps)))


## Doubly robust estimator
# We calculated the expectations for the adjustment already in the np estimator
# Treated
tr_fit = as.matrix(ld_control) %*% as.matrix(mod_par_tr$coefficients)

# Untreated
untr_fit = as.matrix(ld_control) %*% as.matrix(mod_par_untr$coefficients)

# Adjustment term in dr estimator
dr_par_adj = mean((sipp$p401 * tr_fit) / ps - (1 - sipp$p401) * untr_fit / (1 - ps))
ATE_dr_par = ATE_np_par + ATE_ps_par  - dr_par_adj


print(xtable(data.frame(Estimator = c("Non-parametric", "Propensity Score", "Doubly Robust"), 
                        ATE = c(ATE_np_par, ATE_ps_par, ATE_dr_par))), 
      include.rownames = F)

