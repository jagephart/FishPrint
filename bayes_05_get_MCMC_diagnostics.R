# Author: Kelvin Gorospe
# Get convergence diagnostics

rm(list=ls())
library(tidyverse)
library(bayesplot)
datadir <- "/Volumes/jgephart/BFA Environment 2/Data"
outdir <- "/Volumes/jgephart/BFA Environment 2/Outputs"

# Set filenames:
# Mass allocation
c_results <- "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/GHG/2021-04-27_full-model-posterior_Global warming potential_Mass-allocation.RData"
n_results <- "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Nitrogen/2021-04-28_full-model-posterior_Marine eutrophication_Mass-allocation.RData"
p_results <- "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Phosphorus/2021-04-28_full-model-posterior_Freshwater eutrophication_Mass-allocation.RData"
land_results <- "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Land-FCR-priors-only/2021-04-28_full-model-posterior_Land Use_Mass-allocation.RData"
water_results <- "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Water/2021-04-28_full-model-posterior_Water Consumption_Mass-allocation.RData"
wild_results <- "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Wild/2021-04-27_full-model-posterior_Wild-capture-ghg.RData"

# Update filepaths as needed for various model outputs (i.e., allocation method, priors vs no priors, edible weight vs live weight)

####################################################
# Carbon
load(file.path(outdir, c_results))

# Get neff ratios:
mod_neff <- neff_ratio(fit_no_na)
min(mod_neff)

# Get Rhat:
mod_rhat <- rhat(fit_no_na)
min(mod_rhat)
max(mod_rhat)

rm(fit_no_na)

####################################################
# Nitrogen
load(file.path(outdir, n_results))

# Get neff ratios:
mod_neff <- neff_ratio(fit_no_na)
min(mod_neff)

# Get Rhat:
mod_rhat <- rhat(fit_no_na)
min(mod_rhat)
max(mod_rhat)

rm(fit_no_na)

####################################################
# Phosphorus
load(file.path(outdir, p_results))

# Get neff ratios:
mod_neff <- neff_ratio(fit_no_na)
min(mod_neff)

# Get Rhat:
mod_rhat <- rhat(fit_no_na)
min(mod_rhat)
max(mod_rhat)

rm(fit_no_na)

####################################################
# Water
load(file.path(outdir, water_results))

# Get neff ratios:
mod_neff <- neff_ratio(fit_no_na)
min(mod_neff)

# Get Rhat:
mod_rhat <- rhat(fit_no_na)
min(mod_rhat)
max(mod_rhat)

rm(fit_no_na)

####################################################
# Land
load(file.path(outdir, land_results))

# Get neff ratios:
mod_neff <- neff_ratio(fit_no_na)
min(mod_neff)

# Get Rhat:
mod_rhat <- rhat(fit_no_na)
min(mod_rhat)
max(mod_rhat)

rm(fit_no_na)

####################################################
# Wild capture
load(file.path(outdir, wild_results))

# Get neff ratios:
mod_neff <- neff_ratio(fit_no_na)
min(mod_neff)

# Get Rhat:
mod_rhat <- rhat(fit_no_na)
min(mod_rhat)
max(mod_rhat)

rm(fit_no_na)