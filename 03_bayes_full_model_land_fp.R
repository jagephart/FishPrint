# Get Bayesian mean of sci and taxa-level land impacts 

###################### REMINDER FOR LAND: 
# ON-FARM impacts only apply to SYSTEM == PONDS or RECIRCULATING and TANKS
# OFF-FARM impacts EXCLUDE wherever FCR = 0 and bivalves and plants

# Step 0: Set directories, load packages and data
rm(list = ls())

# Libraries for processing and analyses
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

# Load Data
lca_dat <- read.csv(file.path(datadir, "lca-dat-imputed-vars_rep-n-farms.csv"), fileEncoding="UTF-8-BOM") #fileEncoding needed when reading in file from windows computer (suppresses BOM hidden characters)

# Load feed impact constants:
# Choose allocation method
# Choose Impact.category
# IMPORTANT - multiply all values by 1000 to convert to kg CO2 per tonne (currently in kg CO2 per kg)
impact <- "Global warming potential" # i.e., Carbon impact
set_allocation <- "Mass"

fp_dat <- read.csv(file.path(datadir, "20201217_weighted_feed_fp.csv")) %>%
  filter(Allocation == set_allocation) %>%
  mutate(ave_stressor_per_tonne = ave_stressor * 1000)

# ORDER OF feed_weights data vector: soy, crops, fmfo, animal
head(feed_weights)
# ORDER of footprint data vector should match this:
set_fp_order <- c("Soy", "Crop", "Fishery", "Livestock")

fp_constant <- fp_dat %>%
  filter(Impact.category == impact) %>% 
  arrange(match(Input.type, set_fp_order)) %>% # Match index and arrange by custom order
  select(ave_stressor_per_tonne) %>%
  as.matrix() %>%
  c()


# Set data
N = nrow(land_footprint_dat)
N_SCI <- length(unique(land_footprint_dat$sci))
n_to_sci <- land_footprint_dat$sci
n_to_tx <- land_footprint_dat$tx
N_TX <- length(unique(land_footprint_dat$tx))
sci_to_tx <- land_footprint_dat %>%
  select(sci, tx) %>%
  unique() %>%
  pull(tx)
yield <- land_footprint_dat$yield

stan_data <- list(N = N,
                  N_SCI = N_SCI, 
                  n_to_sci = n_to_sci,
                  N_TX = N_TX,
                  #n_to_tx = n_to_tx,
                  sci_to_tx = sci_to_tx,
                  yield = yield)

# Estimate foot print for all scientific names and taxa groups (removed the "all-seafood" level for simplicity)
# GAMMA distribution hierarchical model
# stan_no_na <- 'data {
#   // data for gamma model for FCR
#   int<lower=0> N;  // number of observations
#   vector<lower=0>[N] yield; // data
#   int N_TX; // number of taxa groups
#   int N_SCI; // number of scientific names
#   int n_to_sci[N]; // sciname index
#   int sci_to_tx[N_SCI]; // taxa-group indices
# }
# parameters {
#   vector<lower=0>[N_TX] tx_mu;
#   vector<lower=0>[N_SCI] sci_mu;
#   real<lower=0> tx_sigma; // only need to define sigmas if using option 1
#   real<lower=0> sci_sigma; 
#   //vector<lower=0>[N_SCI] sci_shape; // define shape here if not transforming (using option 2 below)
# }
# transformed parameters {
#   // define transofrmed params for gamma model for FCRs
#   vector<lower=0>[N_SCI] sci_shape;
#   vector<lower=0>[N_SCI] sci_rate;
#   vector<lower=0>[N_TX] tx_shape;
#   vector<lower=0>[N_TX] tx_rate;
# 
#   // reparamaterize gamma to get mu and sigma; defining these here instead of the model section allows us to see these parameters in the output
#   // option 1: mean and variance
#   for (n_tx in 1:N_TX){
#     tx_shape[n_tx] = square(tx_mu[n_tx]) / square(tx_sigma);
#     tx_rate[n_tx] = tx_mu[n_tx] / square(tx_sigma);
#   }
#   for (n_sci in 1:N_SCI){
#     sci_shape[n_sci] = square(sci_mu[n_sci]) / square(sci_sigma);
#     sci_rate[n_sci] = sci_mu[n_sci] / square(sci_sigma);
#   }
#   
#   // option 2: rate = shape / mean
#   //for (n_tx in 1:N_TX){
#   //  tx_rate[n_tx] = tx_shape[n_tx] / tx_mu[n_tx];
#   //}
#   //for (n_sci in 1:N_SCI){
#   //  sci_rate[n_sci] = sci_shape[n_sci] / sci_mu[n_sci];
#   //}
#   
# }
# model {
#   // define priors for gamma model
#   // Put priors on mu and sigma (instead of shape and rate) since this is more intuitive:
#   //tx_mu ~ uniform(0, 100); // note: uniform(0,100) for all of these doesnt help much with convergence
#   //sci_mu ~ uniform(0, 100);
#   //tx_sigma ~ uniform(0, 100); // only need sigmas if calculating shape and rate with mu and sigma
#   //sci_sigma ~ uniform(0, 100);
#   //tx_sigma ~ cauchy(.05, .01);
#   // try adding weakly informative priors on shape and/or rate parameters
#   // shape priors require target += notation: values come outputs of shape parameters in non-hierarchical model - did not help
#   //target += cauchy_lpdf(sci_shape[1] | 16, 8);
#   //target += cauchy_lpdf(sci_shape[2] | 200, 100);
#   //target += cauchy_lpdf(sci_shape[3] | 8, 4);
#   //target += cauchy_lpdf(sci_shape[4] | 0.3, 0.15);
#   //target += cauchy_lpdf(sci_shape[5] | 7, 3.5);
#   //target += cauchy_lpdf(sci_shape[6] | 0.4, 0.2);
#   //target += cauchy_lpdf(sci_shape[7] | 0.4, 0.2);
#   //target += cauchy_lpdf(sci_shape[8] | 14, 7);
#   //target += cauchy_lpdf(sci_shape[9] | 1, 0.5);
#   //target += cauchy_lpdf(sci_shape[10] | 1, 0.5);
#   
#   // FIX IT - create generated quantities section for land = 1 / yield
#   // likelihood
#   // gamma model sci-name and taxa-level
#   for (n in 1:N){
#     yield[n] ~ gamma(sci_shape[n_to_sci[n]], sci_rate[n_to_sci[n]]);
#   }
#   for (n_sci in 1:N_SCI){
#     sci_mu[n_sci] ~ gamma(tx_shape[sci_to_tx[n_sci]], tx_rate[sci_to_tx[n_sci]]);
#   }
# }'

# NORMAL distribution hierarchical model
stan_no_na <- 'data {
  // data for gamma model for FCR
  int<lower=0> N;  // number of observations
  vector<lower=0>[N] yield; // data
  int N_TX; // number of taxa groups
  int N_SCI; // number of scientific names
  int n_to_sci[N]; // sciname index
  int sci_to_tx[N_SCI]; // taxa-group indices
}
parameters {
  vector<lower=0>[N_TX] tx_mu;
  vector<lower=0>[N_SCI] sci_mu;
  //vector<lower=0>[N_TX] tx_sigma;
  //vector<lower=0>[N_SCI] sci_sigma; 
  real<lower=0> tx_sigma;
  real<lower=0> sci_sigma; 
}
model {
  // define priors for gamma model
  // Put priors on mu and sigma (instead of shape and rate) since this is more intuitive:
  //tx_mu ~ uniform(0, 100); // note: uniform(0,100) for all of these doesnt help much with convergence
  //sci_mu ~ uniform(0, 100);
  //tx_sigma ~ uniform(0, 100); 
  //sci_sigma ~ uniform(0, 100);
  //tx_sigma ~ cauchy(.05, .01);

  
  // likelihood
  // normal model sci-name and taxa-level
  for (n in 1:N){
    yield[n] ~ normal(sci_mu[n_to_sci[n]], sci_sigma);
  }

  for (n_sci in 1:N_SCI){
    sci_mu[n_sci] ~ normal(tx_mu[sci_to_tx[n_sci]], tx_sigma);
  }
}
generated quantities {
  // Land = 1 / Yield - model converges better if yield ~ normal() instead of 1/yield ~ normal()
  vector<lower=0>[N_TX] tx_land;
  vector<lower=0>[N_SCI] sci_land;
  
  // Calculations
  for (n_tx in 1:N_TX) {
    tx_land[n_tx] = 1 / tx_mu[n_tx];
  }
  for (n_sci in 1:N_SCI) {
    sci_land[n_sci] = 1 / sci_mu[n_sci];
  }
  
}'

no_na_mod <- stan_model(model_code = stan_no_na)
# Note: For Windows, apparently OK to ignore this warning message:
# Warning message:
#   In system(paste(CXX, ARGS), ignore.stdout = TRUE, ignore.stderr = TRUE) :
#   'C:/rtools40/usr/mingw_/bin/g++' not found

# Fit model:
# Set seed while testing
fit_no_na <- sampling(object = no_na_mod, data = stan_data, cores = 4, seed = "11729", iter = 10000, control = list(adapt_delta = 0.99))
#fit_no_na <- sampling(object = no_na_mod, data = stan_data, cores = 4, iter = 10000, control = list(adapt_delta = 0.99))
summary(fit_no_na)$summary

#launch_shinystan(fit_no_na)

######################################################################################################
# PLOT RESULTS
plot_theme <- theme(title = element_text(size = 18),
                    axis.title.x = element_text(size = 16),
                    axis.text=element_text(size=14, color = "black"))

# Key for naming sci and taxa levels
# Get full taxa group names back
index_key <- land_footprint_dat %>%
  select(clean_sci_name, sci, taxa, tx, n_in_sci, n_in_taxa) %>%
  unique() %>%
  mutate(taxa = as.character(taxa),
         full_taxa_name = case_when(taxa == "fresh_crust" ~ "freshwater crustaceans",
                                    taxa == "misc_diad" ~ "misc diadromous fishes",
                                    taxa == "misc_fresh" ~ "misc freshwater fishes",
                                    taxa == "misc_marine" ~ "misc marine fishes",
                                    TRUE ~ taxa),
         taxa = as.factor(taxa),
         full_taxa_name = as.factor(full_taxa_name))

# Use tidybayes + ggdist for finer control of aes mapping (instead of bayesplots) 
get_variables(fit_no_na)

# Sci-level land footprints as point intervals:
fit_no_na %>%
  spread_draws(sci_land[sci]) %>%
  median_qi(.width = 0.8) %>%
  left_join(index_key, by = "sci") %>% # Join with index key to get sci and taxa names
  mutate(clean_sci_name = fct_reorder(clean_sci_name, as.character(full_taxa_name))) %>%
  ggplot(aes(y = clean_sci_name, x = sci_land, xmin = .lower, xmax = .upper, color = full_taxa_name)) +
  geom_pointinterval() +
  theme_classic() + 
  plot_theme + 
  labs(x = "hectares per tonne", y = "", title = "Land", color = "taxa group")
ggsave(filename = file.path(outdir, "plot_land-footprint_sci-level.png"), width = 11, height = 8.5)

# Taxa-level land footprints as densities: 
# FIX IT - DENSITIES ARE SUPER STEEP (see just trout for example) - if this is true for final dataset, just plot as point-interval
fit_no_na %>%
  spread_draws(tx_land[tx]) %>%
  #median_qi(.width = 0.8) %>% # need at least 2 points to select a bandwidth automatically
  left_join(index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  #filter(tx_land > 100)
  filter(full_taxa_name == "freshwater crustaceans") %>%
  ggplot(aes(y = full_taxa_name, x = tx_land)) +
  stat_halfeye(aes(slab_fill = full_taxa_name)) +
  theme_classic() + 
  plot_theme + 
  theme(legend.position = "none") +
  labs(x = "hectares per tonne", y = "", title = "Land")
ggsave(filename = file.path(outdir, "plot_land-footprint_taxa-level.png"), width = 11, height = 8.5)


# OPTION 2: Taxa-level plots with color themes:
# If we want to mimic bayesplot color schemes, can get hexadecimal colors and input manually to stat_halfeye aesthetics
color_scheme_get("blue")
color_scheme_get("green")

