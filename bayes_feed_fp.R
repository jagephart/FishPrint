# Combine bayes_prop_feed and bayes_fcr to calculate footprint


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

lca_dat_no_zeroes <- clean.lca(LCA_data = lca_dat) %>%
  #select(clean_sci_name, Feed_soy_percent, Feed_othercrops_percent, Feed_FMFO_percent, Feed_animal_percent, taxa_group_name) %>%
  # NOTE multinomial-dirchlet model requires all elements > 0 (change some to 0.001 for now?)
  mutate(feed_soy_new = if_else(Feed_soy_percent == 0, true = 0.01, false = Feed_soy_percent),
         feed_crops_new = if_else(Feed_othercrops_percent == 0, true = 0.01, false = Feed_othercrops_percent),
         feed_fmfo_new = if_else(Feed_FMFO_percent == 0, true = 0.01, false = Feed_FMFO_percent),
         feed_animal_new = if_else(Feed_animal_percent == 0, true = 0.01, false = Feed_animal_percent)) %>%
  # Renomoralize values so they sum to 1
  mutate(sum = rowSums(select(., contains("new")))) %>%
  mutate(feed_soy_new = feed_soy_new / sum,
         feed_crops_new = feed_crops_new / sum,
         feed_fmfo_new = feed_fmfo_new / sum,
         feed_animal_new = feed_animal_new / sum)

fp_dat <- read.csv(file.path(datadir, "Feed_FP_raw.csv"))
fp_clean <- clean.feedFP(fp_dat)

######################################################################################################
# Model 1: Remove all NAs - estimate feed footprint for Oncorhynchus mykiss

# Remove NAs
lca_dat_no_na <- lca_dat_no_zeroes %>%
  filter(clean_sci_name == "Oncorhynchus mykiss") %>%
  filter(is.na(Feed_soy_percent)==FALSE) 

# Set data for model
# for FCR model:
x <- lca_dat_no_na$FCR

# for Feed proportion model:
k = 4
n = 3
feed_weights <- lca_dat_no_na %>%
  select(contains("new")) %>%
  as.matrix()

# for final foot print calculation
fp_dat <- fp_clean %>%
  filter(Category != "Energy") %>%
  select(FP, Category, FP_val) 

# ORDER: Animal, crop, FMFO, soy
fp_carbon_dat <- fp_dat %>%
  filter(FP == "Carbon") %>%
  select(-FP) %>%
  pivot_wider(names_from = Category, values_from = FP_val) %>%
  as.matrix() %>%
  c()

# note: dirichlet_rng is just a random number generator:
# rep_vector(x, m) creates a column consisting of m copies of x
# generated quantities {
#   vector[k] theta = dirichlet_rng(rep_vector(alpha, k));
# }

# Estimate feed component proportions for a single species
stan_pooled <- 'data {
  int<lower=0> n;  // number of observations
  vector[n] x; // data
  int<lower=1> k; // number of feed types
  simplex[k] feed_weights[n]; // array of feed weights simplexes
  vector[k] fp_carbon_dat;
}
parameters {
  // FCR model:
  real<lower=0> mu;
  real<lower=0> sigma;
  // Feed proportion model:
  vector<lower=0>[k] alpha;
  simplex[k] theta;
}
model {
  
  x ~ normal(mu, sigma); // note: stan interprets second param as standard deviation

  for (i in 1:n) {
    feed_weights[i] ~ dirichlet(alpha); // estimate vector of alphas based on the data of feed weights
  }
  theta ~ dirichlet(alpha); // now, estimate feed weights based on the vector of alphas
}
generated quantities {
  real species_carbon_footprint;
  vector[k] weighted_feed_footprint;
  real total_feed_footprint;
  
  weighted_feed_footprint = fp_carbon_dat .* theta;
  total_feed_footprint = sum(weighted_feed_footprint);
  species_carbon_footprint = mu * total_feed_footprint;
}'

no_missing_mod <- stan_model(model_code = stan_pooled, verbose = TRUE)
# Note: For Windows, apparently OK to ignore this warning message:
# Warning message:
#   In system(paste(CXX, ARGS), ignore.stdout = TRUE, ignore.stderr = TRUE) :
#   'C:/rtools40/usr/mingw_/bin/g++' not found

# Fit model:
fit_pooled <- sampling(object = no_missing_mod, data = list(n = n,
                                                            x = x,
                                                            k = k,
                                                            feed_weights = feed_weights,
                                                            fp_carbon_dat = fp_carbon_dat),
                       iter = 10000, cores = 4,
                       control = list(adapt_delta = 0.99))
print(fit_pooled)

feeds <- c("soy", "crops", "fmfo", "animal")
feed_key <- data.frame(weighted_feed_footprint = paste("weighted_feed_footprint[", feeds, "]", sep = ""))

fit_pooled_clean <- fit_pooled
names(fit_pooled_clean)[grep(names(fit_pooled_clean), pattern = "weighted_feed_footprint")] <- feed_key$weighted_feed_footprint

distribution_pooled <- as.matrix(fit_pooled_clean)

# FIX IT - replace parameter names and add plots for other parameters
plot_theme <- theme(axis.text=element_text(size=14, color = "black"))

p_footprint <- mcmc_areas_ridges(distribution_pooled,
                             pars = vars(contains("footprint")),
                             prob = 0.8) + 
  ggtitle("Oncorhynchus mykiss full feed footprint model", "with 80% credible intervals") +
  plot_theme +
  xlim(0, 10)

p_footprint
ggsave(filename = file.path(outdir, "bayes-example_trout_feed-proportion_alphas.png"), width = 11, height = 8.5)
