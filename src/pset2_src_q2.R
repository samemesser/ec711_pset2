# Samuel Messer
# EC711 Problem Set 2
# This code draws heavily on the replication package from Cattaneo, Titiunik, 
# and Vazquez-Bare (2017)
rm(list = ls())

library(readstata13)
library(rdlocrand)
library(rdrobust)
library(rddensity)
library(stargazer)

# Read data
headstart <- read.dta13("data/headstart.dta")

# Set cutoff poverty rate
cutoff = 59.1984
# Outcome variable of interest
Y = headstart$mort_age59_related_postHS

# Generate running variable
R = headstart$povrate60 - cutoff

# Generate treatment indicator
D = as.numeric(R >= 0)

# Placebo outcomes
# These people should not be affected by head start, so we would expect to see no effect here
Plac <- cbind(headstart$mort_age25plus_related_postHS, headstart$mort_age25plus_injuries_postHS)

#Initialize table to store results 
local_linear <- array(NA, dim = c(4, 3))

# Local Linear Regression
# masspoints and stdvars are applied according to instructions in replication package
# Outcome of interest
tmp <- rdrobust(Y, R, p = 1, masspoints = "off", stdvars = "on")
local_linear[1, 1] = tmp$coef[1]
local_linear[2, 1] = tmp$se[1]
local_linear[3, 1] = tmp$pv[3]
local_linear[4, 1] = tmp$bws[1,1]

# Placebo outcomes
tmp <- rdrobust(Plac[, 1], R, p = 1, masspoints = "off", stdvars = "on")
local_linear[1, 2] = tmp$coef[1]
local_linear[2, 2] = tmp$se[1]
local_linear[3, 2] = tmp$pv[3]
local_linear[4, 2] = tmp$bws[1,1]

# Placebo outcome 2
tmp <- rdrobust(Plac[, 2], R, p = 1, masspoints = "off", stdvars = "on")
local_linear[1, 3] = tmp$coef[1]
local_linear[2, 3] = tmp$se[1]
local_linear[3, 3] = tmp$pv[3]
local_linear[4, 3] = tmp$bws[1,1]

round(local_linear,3)

stargazer(local_linear)

#Initialize table to store results 
local_quadratic <- array(NA, dim = c(4, 3))

# Local quadratic Regression
# masspoints and stdvars are applied according to instructions in replication package
# Outcome of interest
tmp <- rdrobust(Y, R, p = 2, masspoints = "off", stdvars = "on")
local_quadratic[1, 1] = tmp$coef[1]
local_quadratic[2, 1] = tmp$se[1]
local_quadratic[3, 1] = tmp$pv[3]
local_quadratic[4, 1] = tmp$bws[1,1]

# Placebo outcomes
tmp <- rdrobust(Plac[, 1], R, p = 2, masspoints = "off", stdvars = "on")
local_quadratic[1, 2] = tmp$coef[1]
local_quadratic[2, 2] = tmp$se[1]
local_quadratic[3, 2] = tmp$pv[3]
local_quadratic[4, 2] = tmp$bws[1,1]

# Placebo outcome 2
tmp <- rdrobust(Plac[, 2], R, p = 2, masspoints = "off", stdvars = "on")
local_quadratic[1, 3] = tmp$coef[1]
local_quadratic[2, 3] = tmp$se[1]
local_quadratic[3, 3] = tmp$pv[3]
local_quadratic[4, 3] = tmp$bws[1,1]

round(local_quadratic,3)

stargazer(local_quadratic)

# Specify covariates for randomization window selection
cov_1960 <- cbind(headstart$census1960_pop,
                  headstart$census1960_pctsch1417,
                  headstart$census1960_pctsch534, 
                  headstart$census1960_pctsch25plus,
                  headstart$census1960_pop1417,
                  headstart$census1960_pop534,
                  headstart$census1960_pop25plus,
                  headstart$census1960_pcturban,
                  headstart$census1960_pctblack)


cov_test <- cbind(headstart$mort_age59_related_preHS, cov_1960)
# Plot p-values for winselect test
tmp <- rdwinselect(R, cov_test, reps = 1000, statistic = "ksmirnov", wmin = .3 , wstep = .05, level = .15, nwin = 34, plot = TRUE, quietly = F)

# Set optimal window
w = 1.1

# Initialize otput array
rand_inf <- array(NA, dim = c(4,2))

# Outcome

tmp <- rdrandinf(Y, R, wl = -w, wr = w, reps = 1000, quietly = TRUE)
rand_inf[1:4,1] <- c(0, tmp$window[2], tmp$obs.stat, tmp$p.value)


tmp <- rdrandinf(Y, R, wl = -w, wr = w, reps = 1000, p = 1, quietly = TRUE)
rand_inf[1:4,2] <- c(1, tmp$window[2], tmp$obs.stat, tmp$p.value)

round(rand_inf, 3)

stargazer(rand_inf)
