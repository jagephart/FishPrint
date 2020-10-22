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
  mutate(feed_soy_new = if_else(Feed_soy_percent == 0, true = 0.0001, false = Feed_soy_percent),
         feed_crops_new = if_else(Feed_othercrops_percent == 0, true = 0.0001, false = Feed_othercrops_percent),
         feed_fmfo_new = if_else(Feed_FMFO_percent == 0, true = 0.0001, false = Feed_FMFO_percent),
         feed_animal_new = if_else(Feed_animal_percent == 0, true = 0.0001, false = Feed_animal_percent)) %>%
  # Renomoralize values so they sum to 1
  mutate(sum = rowSums(select(., contains("new")))) %>%
  mutate(feed_soy_new = feed_soy_new / sum,
         feed_crops_new = feed_crops_new / sum,
         feed_fmfo_new = feed_fmfo_new / sum,
         feed_animal_new = feed_animal_new / sum) %>%
  filter(is.na(Feed_soy_percent)==FALSE)



k = 4
n = 3
feed_weights <- lca_dat_no_na %>%
  filter(clean_sci_name == "Oncorhynchus mykiss") %>%
  select(contains("new")) %>%
  as.matrix()

# Try to get dirichlet to work with just one set of studies: Oncorhynchus mykiss

# note: dirichlet_rng is just a random number generator:
# rep_vector(x, m) creates a column consisting of m copies of x
# generated quantities {
#   vector[k] theta = dirichlet_rng(rep_vector(alpha, k));
# }

# OLD CODE: only estimating alpha (the shape parameters for the dirichlet)
# notice tighter distributions for alpha 1 and alpha 4 - due to replicate data values for soy feed and animal feed?
# stan_pooled <- 'data {
#   int<lower=0> n;  // number of observations
#   int<lower=1> k; // number of feed types
#   simplex[k] feed_weights[n]; // array of feed weights simplexes
# }
# parameters {
#   vector<lower=0>[k] alpha;
# }
# model {
#   for (i in 1:n) {
#     feed_weights[i] ~ dirichlet(alpha);
#   }
# }'

# Estimate feed component proportions for a single species
stan_pooled <- 'data {
  int<lower=0> n;  // number of observations
  int<lower=1> k; // number of feed types
  simplex[k] feed_weights[n]; // array of feed weights simplexes
}
parameters {
  vector<lower=0>[k] alpha;
  simplex[k] theta;
}
model {
  for (i in 1:n) {
    feed_weights[i] ~ dirichlet(alpha); // estimate vector of alphas based on the data of feed weights
  }
  theta ~ dirichlet(alpha); // now, estimate feed weights based on the vector of alphas
}'

no_missing_mod <- stan_model(model_code = stan_pooled, verbose = TRUE)
# Note: For Windows, apparently OK to ignore this warning message:
# Warning message:
#   In system(paste(CXX, ARGS), ignore.stdout = TRUE, ignore.stderr = TRUE) :
#   'C:/rtools40/usr/mingw_/bin/g++' not found

# Fit model:
fit_pooled <- sampling(object = no_missing_mod, data = list(n = n,
                                                            k = k,
                                                            feed_weights = feed_weights))
print(fit_pooled)

distribution_pooled <- as.matrix(fit_pooled)
p_alpha <- mcmc_areas_ridges(distribution_pooled,
                       pars = vars(contains("alpha")),
                       prob = 0.8)
p_alpha
p_theta <- mcmc_areas_ridges(distribution_pooled,
                             pars = vars(contains("theta")),
                             prob = 0.8)
p_theta

