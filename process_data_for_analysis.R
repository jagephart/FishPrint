# Clean compiled LCA data needed for all models

rm(list = ls())
library(tidyverse)
library(rstan)
library(data.table)
library(countrycode) # part of clean.lca
library(bayesplot) # for mcmc_areas_ridges
library(shinystan)
library(brms)
library(tidybayes)

# Mac
datadir <- "/Volumes/jgephart/BFA Environment 2/Data"
outdir <- "/Volumes/jgephart/BFA Environment 2/Outputs"
# Windows
# datadir <- "K:/BFA Environment 2/Data"
# outdir <- "K:BFA Environment 2/Outputs"

lca_dat <- read.csv(file.path(datadir, "LCA_compiled_20201109.csv"), fileEncoding="UTF-8-BOM") #fileEncoding needed when reading in file from windows computer (suppresses BOM hidden characters)
source("Functions.R")

# Clean LCA data
lca_dat_clean <- clean.lca(LCA_data = lca_dat)

# Rebuild FAO fish production from zip file
fishstat_dat <- rebuild_fish("/Volumes/jgephart/FishStatR/Data/Production-Global/ZippedFiles/GlobalProduction_2019.1.0.zip")

# Classify species into taxa groupings
# Use ISSCAAP grouping (in FAO data) to help with classification
lca_dat_clean_groups <- add_taxa_group(lca_dat_clean, fishstat_dat)

# Output clean data with groups (later will read this back in to join with model predictions)
write.csv(lca_dat_clean_groups, file.path(datadir, "lca_clean_with_groups.csv"), row.names = FALSE)

# Output taxa groupings and sample sizes:
#data.frame(table(lca_dat_clean_groups$taxa_group_name))
#data.frame(table(lca_dat_clean_groups$taxa)) # abbreviated version of taxa_group_name for writing models
#lca_dat_clean_groups %>% select(taxa_group_name, clean_sci_name) %>% unique() %>% arrange(taxa_group_name)
#write.csv(data.frame(table(lca_dat_clean_groups$taxa_group_name)), file.path(outdir, "taxa_group_sample_size.csv"))
#write.csv(lca_dat_clean_groups %>% select(taxa_group_name, clean_sci_name) %>% unique() %>% arrange(taxa_group_name), file.path(outdir, "taxa_group_composition.csv"))
