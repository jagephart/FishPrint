# Bayesian estimation of proportions of each feed component (soy, other crops, FMFOs, and other animal)

rm(list=ls())
library(tidyverse)
library(rstan)
library(taxize)
library(data.table)
library(countrycode) # part of clean.lca
library(bayesplot) # for mcmc_areas_ridges

# Mac
datadir <- "/Volumes/jgephart/BFA Environment 2/Data"
outdir <- "/Volumes/jgephart/BFA Environment 2/Outputs"
# Windows
# datadir <- "K:/BFA Environment 2/Data"
# outdir <- "K:BFA Environment 2/Outputs"

lca_dat <- read.csv(file.path(datadir, "LCA_compiled_20201006.csv"), fileEncoding="UTF-8-BOM") #fileEncoding needed when reading in file from windows computer (suppresses BOM hidden characters)
source("Functions.R")

lca_dat_clean <- clean.lca(LCA_data = lca_dat)

######################################################################################################
# Model 1: Remove all NAs


lca_dat_no_na <- lca_dat_clean %>%
  select(clean_sci_name, Feed_soy_percent, Feed_othercrops_percent, Feed_FMFO_percent, Feed_animal_percent) %>%
  # NOTE multinomial-dirchlet model requires all elements > 0 (change some to 0.001 for now?)
  # mutate(feed_soy_new = if_else(Feed_soy_percent == 0, true = 0.01, false = Feed_soy_percent),
  #        feed_crops_new = if_else(Feed_othercrops_percent == 0, true = 0.01, false = Feed_othercrops_percent),
  #        feed_fmfo_new = if_else(Feed_FMFO_percent == 0, true = 0.01, false = Feed_FMFO_percent),
  #        feed_animal_new = if_else(Feed_animal_percent == 0, true = 0.01, false = Feed_animal_percent)) %>%
  # # Renomoralize values so they sum to 1
  # mutate(new_sum = rowSums(select(., contains("new")))) %>%
  # mutate(feed_soy_new = feed_soy_new / new_sum,
  #        feed_crops_new = feed_crops_new / new_sum,
  #        feed_fmfo_new = feed_fmfo_new / new_sum,
  #        feed_animal_new = feed_animal_new / new_sum) %>%
  filter(is.na(Feed_soy_percent)==FALSE)

stan_pooled <- 'data {
  int<lower=0> n;  // number of observations
  int<lower = 1> k; // number of trials
  real<lower = 0> alpha;
}
generated quantities {
  vector[k] theta = dirichlet_rng(rep_vector(alpha, k));
}
model {
  theta[n] ~ dirichlet(alpha)
}'

no_missing_mod <- stan_model(model_code = stan_pooled, verbose = TRUE)
# Note: For Windows, apparently OK to ignore this warning message:
# Warning message:
#   In system(paste(CXX, ARGS), ignore.stdout = TRUE, ignore.stderr = TRUE) :
#   'C:/rtools40/usr/mingw_/bin/g++' not found

# Fit model:
fit_grouped <- sampling(object = grouped_mod, data = list())

