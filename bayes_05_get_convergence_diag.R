# Get convergence diagnostics

rm(list=ls())
library(tidyverse)
library(bayesplot)
datadir <- "/Volumes/jgephart/BFA Environment 2/Data"
outdir <- "/Volumes/jgephart/BFA Environment 2/Outputs"

# Set filenames:
# Mass allocation
c_results <- "PRIORS/Carbon/2021-01-11_full-model-posterior_Global warming potential_Mass-allocation.RData"
n_results <- "PRIORS/Nitrogen/2021-01-11_full-model-posterior_Marine eutrophication_Mass-allocation.RData"
p_results <- "PRIORS/Phosphorus/2021-01-11_full-model-posterior_Freshwater eutrophication_Mass-allocation.RData"
land_results <- "PRIORS/Land/2021-01-11_full-model-posterior_Land Use_Mass-allocation.RData"
water_results <- "PRIORS/Water/2021-01-11_full-model-posterior_Water Consumption_Mass-allocation.RData"
wild_results <- "PRIORS/Wild/2021-01-12_full-model-posterior_Wild-Capture-ghg.RData"

# Gross energy allocation (no wild capture results)
# c_results <- "PRIORS/Carbon/2021-01-11_full-model-posterior_Global warming potential_Gross energy content-allocation.RData"
# n_results <- "PRIORS/Nitrogen/2021-01-11_full-model-posterior_Marine eutrophication_Gross energy content-allocation.RData"
# p_results <- "PRIORS/Phosphorus/2021-01-11_full-model-posterior_Freshwater eutrophication_Gross energy content-allocation.RData"
# land_results <- "PRIORS/Land/2021-01-12_full-model-posterior_Land Use_Gross energy content-allocation.RData"
# water_results <- "PRIORS/Water/2021-01-11_full-model-posterior_Water Consumption_Gross energy content-allocation.RData"

# Economic allocation (no wild capture results)
# c_results <- "PRIORS/Carbon/2021-01-11_full-model-posterior_Global warming potential_Economic-allocation.RData"
# n_results <- "PRIORS/Nitrogen/2021-01-11_full-model-posterior_Marine eutrophication_Economic-allocation.RData"
# p_results <- "PRIORS/Phosphorus/2021-01-11_full-model-posterior_Freshwater eutrophication_Economic-allocation.RData"
# land_results <- "PRIORS/Land/2021-01-13_full-model-posterior_Land Use_Economic-allocation.RData"
# water_results <- "PRIORS/Water/2021-01-11_full-model-posterior_Water Consumption_Economic-allocation.RData"

# Mass allocation NO PRIORS
# c_results <- "NO PRIORS/Carbon/2021-01-11_full-model-posterior_Global warming potential_Mass-allocation.RData"
# n_results <- "NO PRIORS/Nitrogen/2021-01-11_full-model-posterior_Marine eutrophication_Mass-allocation.RData"
# p_results <- "NO PRIORS/Phosphorus/2021-01-11_full-model-posterior_Freshwater eutrophication_Mass-allocation.RData"
# land_results <- "NO PRIORS/Land/2021-01-12_full-model-posterior_Land Use_Mass-allocation.RData"
# water_results <- "NO PRIORS/Water/2021-01-11_full-model-posterior_Water Consumption_Mass-allocation.RData"
# wild_results <- "NO PRIORS/Wild/2021-01-13_full-model-posterior_Wild-capture-ghg.RData"

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