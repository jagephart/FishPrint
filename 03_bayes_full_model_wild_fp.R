# Calculate greenhouse gas footprint for wild caught fisheries

rm(list = ls())
library(tidyverse)
library(rstan)
library(bayesplot) # for mcmc_areas_ridges
library(shinystan)
library(brms)
library(tidybayes)

# Mac
datadir <- "/Volumes/jgephart/BFA Environment 2/Data"
outdir <- "/Volumes/jgephart/BFA Environment 2/Outputs"

# FIX IT - figure out weightings - easiest would be to pass a single weighting to each species' ghg
wild_dat <- read.csv(file.path(datadir, "fisheries_fuel_use.csv")) %>% as_tibble()

wild_dat_new_weights <- wild_dat %>%
  filter(species_group != "Finfish") %>%
  # Remove mixed gear and nei observations
  filter(!str_detect(pattern = " nei", species)) %>%
  filter(gear != "Other, Mixed, or Unknown") %>%
  # Remove observations with 0 gear, species, or consumption weighting
  filter(gear_weighting > 0 & species_weighting > 0 & consumption_weighting > 0) %>%
  # Re-weight gear within species
  group_by(species) %>%
  mutate(gear_weights_new = gear_weighting/sum(gear_weighting)) %>%
  ungroup() %>%
  # Re-weight species and consumption within species_groups
  group_by(species_group) %>%
  mutate(species_weights_new = species_weighting/sum(species_weighting),
         consumption_weights_new = consumption_weighting/sum(consumption_weighting)) %>%
  ungroup() %>%
  # Calculate overall weights and re-weight
  mutate(prod_of_weights = gear_weights_new * species_weights_new * consumption_weights_new) %>%
  group_by(species_group) %>%
  mutate(overall_weights = prod_of_weights/sum(prod_of_weights))

# LEFT OFF HERE - Find Bayesian mean of ghg column with species nested within species_group, then apply overall_weights

jg_results <- read.csv(file.path(datadir, "fisheries_fuel_use.csv")) %>%
# Remove mixed gear and nei observations
  filter(!str_detect(pattern = " nei", species)) %>%
  filter(gear != "Other, Mixed, or Unknown") %>%
  # Remove observations with 0 gear, species, or consumption weighting
  filter(gear_weighting > 0 & species_weighting > 0 & consumption_weighting > 0) %>%
  # Re-weight gear within each species
  group_by(species_group, species) %>%
  mutate(gear_weighting_new = gear_weighting/sum(gear_weighting)) %>%
  # Create species gear-weighted means
  summarise(species_ghg_kg_t = sum(ghg*gear_weighting_new), 
            species_weighting = mean(species_weighting), 
            consumption_weighting = mean(consumption_weighting)) %>%
  # Re-weight species and consumption within taxa group
  ungroup() %>%
  group_by(species_group) %>%
  mutate(species_consumption_weighting = (species_weighting*consumption_weighting)/sum(species_weighting*consumption_weighting)) %>%
  summarise(ghg_kg_t = sum(species_ghg_kg_t*species_consumption_weighting))

# FOR NOW, no weightings and pooling observations up to taxa (i.e., species_group) level and skip scientific name level
# Set data
N = nrow(wild_dat)
N_SCI <- length(unique(wild_dat$species))
n_to_sci <- wild_dat$sci
n_to_tx <- wild_dat$tx
N_TX <- length(unique(wild_dat$tx))
sci_to_tx <- wild_dat %>%
  select(sci, tx) %>%
  unique() %>%
  pull(tx)
ghg <- wild_dat$ghg

stan_data <- list(N = N,
                  N_TX = N_TX,
                  n_to_tx = n_to_tx,
                  ghg = ghg)


# GAMMA distribution model no levels
stan_no_na <- 'data {
  // data for gamma model
  int<lower=0> N;  // number of observations
  vector<lower=0>[N] ghg; // data
  int n_to_tx[N]; // tx index
  int N_TX; // number of taxa groups
}
parameters {
  vector<lower=0>[N_TX] tx_mu;
  //real<lower=0> tx_sigma;
  vector<lower=0>[N_TX] tx_shape;
}
transformed parameters {
  // define transfomed params for gamma model
  // vector<lower=0>[N_TX] tx_shape;
  vector<lower=0>[N_TX] tx_rate;

  // reparamaterize gamma to get mu and sigma; defining these here instead of the model section allows us to see these parameters in the output
  //for (n_tx in 1:N_TX){
  //  tx_shape[n_tx] = square(tx_mu[n_tx]) ./ square(tx_sigma);
  //  tx_rate[n_tx] = tx_mu[n_tx] ./ square(tx_sigma);
  //}
  // reparamaterize option 2:
  for (n_tx in 1:N_TX){
    tx_rate[n_tx] = tx_shape[n_tx] / tx_mu[n_tx];
  }
}
model {
  // define priors for gamma model
  // Put priors on mu and sigma:
  //tx_mu ~ uniform(0, 100); // 
  //sci_mu ~ uniform(0, 100);
  //tx_sigma ~ uniform(0, 100);
  tx_shape ~ cauchy(0, 5);

  // likelihood
  // gamma model sci-name and taxa-level
  for (n in 1:N){
    ghg[n] ~ gamma(tx_shape[n_to_tx[n]], tx_rate[n_to_tx[n]]);
  }
}'

# NORMAL distribution model no levels
# stan_no_na <- 'data {
#   // data for gamma model
#   int<lower=0> N;  // number of observations
#   vector<lower=0>[N] ghg; // data
#   int n_to_tx[N]; // tx index
#   int N_TX; // number of taxa groups
# }
# parameters {
#   vector<lower=0>[N_TX] tx_mu;
#   real<lower=0> tx_sigma;
# }
# model {
#   // define priors for gamma model
#   // Put priors on mu and sigma:
#   //tx_mu ~ uniform(0, 100); //
#   //sci_mu ~ uniform(0, 100);
#   //tx_sigma ~ uniform(0, 100);
# 
#   // likelihood
#   // normal model sci-name and taxa-level
#   for (n in 1:N){
#     ghg[n] ~ normal(tx_mu[n_to_tx[n]], tx_sigma);
#   }
# }'


no_na_mod <- stan_model(model_code = stan_no_na)
# Note: For Windows, apparently OK to ignore this warning message:
# Warning message:
#   In system(paste(CXX, ARGS), ignore.stdout = TRUE, ignore.stderr = TRUE) :
#   'C:/rtools40/usr/mingw_/bin/g++' not found

# Fit model:
# Set seed while testing
fit_no_na <- sampling(object = no_na_mod, data = stan_data, cores = 4, seed = "11729", iter = 5000, control = list(adapt_delta = 0.99))
summary(fit_no_na)$summary

#launch_shinystan(fit_no_na)

######################################################################################################
# PLOTS
plot_theme <- theme(title = element_text(size = 18),
                    axis.title.x = element_text(size = 16),
                    axis.text=element_text(size=14, color = "black"))

var_title = "kg CO2-eq per tonne"
groups_title = ""
plot_title = "Wild catch greenhouse gas emissions"

# FIRST PLOT DATA
ggplot(data = wild_dat, aes(x = ghg, y = species_group)) +
  #geom_violin() +
  geom_boxplot(aes(x = ghg, y = species_group)) +
  #stat_halfeye() +
  theme_classic() +
  plot_theme +
  labs(y = groups_title, x = var_title, title = plot_title)
ggsave(filename = file.path(outdir, "plot_wild-catch-ghg-data_tx-level.png"), width = 11, height = 8.5)
  
# PLOT RESULTS
# Key for naming sci and taxa levels
index_key <- wild_dat %>%
  select(species_group, tx) %>%
  unique()

# Use tidybayes + ggdist for finer control of aes mapping (instead of bayesplots) 
get_variables(fit_no_na)

# Tx-level footprints as point intervals:
fit_no_na %>%
  spread_draws(tx_mu[tx]) %>%
  median_qi(.width = 0.95) %>%
  left_join(index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  ggplot(aes(y = species_group, x = tx_mu, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() +
  theme_classic() + 
  plot_theme + 
  labs(x = var_title, y = groups_title, title = plot_title)
ggsave(filename = file.path(outdir, "plot_wild-catch-ghg-footprint_tx-level.png"), width = 11, height = 8.5)

# Option: Same plot with halfeye densities
# fit_no_na %>%
#   spread_draws(tx_mu[tx]) %>%
#   #median_qi(.width = 0.8) %>%
#   left_join(index_key, by = "tx") %>%
#   ggplot(aes(y = species_group, x = tx_mu)) +
#   stat_halfeye() +
#   theme_classic() + 
#   plot_theme + 
#   labs(x = var_title, y = groups_title, title = plot_title)
# ggsave(filename = file.path(outdir, "plot_wild-catch-ghg-footprint_tx-level_densities.png"), width = 11, height = 8.5)