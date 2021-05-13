# Author: Kelvin Gorospe
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

# GHG data
wild_dat <- read.csv(file.path(datadir, "fisheries_fuel_use.csv")) %>% as_tibble()

# Apply gear, species, and consumption weighting
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
  mutate(overall_weights = prod_of_weights/sum(prod_of_weights)) %>%
  ungroup() %>%
  # ADD sci and taxa-level indices
  mutate(clean_sci_name = as.factor(species),
         sci = as.numeric(clean_sci_name),
         taxa = as.factor(species_group),
         tx = as.numeric(taxa)) %>%
  # LAST FORMATING STEP - always arrange by clean_sci_name
  arrange(species)

# OPTION: Apply EDIBLE PORTIONS adjustment
wild_edible <- read.csv(file.path(datadir, "capture_edible_CFs.csv"))
wild_dat_new_weights <- wild_dat_new_weights %>% # Join and apply edible portions weightings
  left_join(wild_edible, by = c("taxa" = "full_taxa_name")) %>% 
  mutate(ghg = ghg * 1/(edible_mean/100))


# Get priors on taxa-level GHG
source("Functions.R")
priors_csv <- clean_wild_priors("Priors - Capture.csv") %>%
  select(taxa = Group.name, Mean.GHG) %>%
  filter(taxa %in% unique(wild_dat_new_weights$taxa)) %>%
  right_join(wild_dat_new_weights %>% select(taxa) %>% unique(), by = "taxa") %>%
  arrange(taxa)
  
# Format priors for STAN
# can't pass NAs into STAN - drop NAs but keep track of vector positions
prior_vec_index <- which(is.na(priors_csv$Mean.GHG)==FALSE)
priors <- priors_csv$Mean.GHG[prior_vec_index]

#################################################################
# Set data, indices, and weights for STAN

# DATA
ghg <- wild_dat_new_weights$ghg

# INDICES
N = nrow(wild_dat_new_weights)
N_SCI <- length(unique(wild_dat_new_weights$species))
n_to_sci <- wild_dat_new_weights$sci
n_to_tx <- wild_dat_new_weights$tx
N_TX <- length(unique(wild_dat_new_weights$tx))
sci_to_tx <- wild_dat_new_weights %>%
  select(sci, tx) %>%
  unique() %>%
  pull(tx)

# WEIGHTS:
# Get sci-level weightings for generating taxa-level quantities:
# Overall_weights column was calculated per row, before normalizing within a taxa group
# To get sci-level weightings, sum within species
# IMPORTANT arrange by clean_sci_name so that the order matches data
sci_w <- wild_dat_new_weights %>%
  group_by(clean_sci_name, taxa) %>%
  summarise(sci_weights = sum(overall_weights)) %>%
  ungroup() %>%
  arrange(clean_sci_name) %>%
  pull(sci_weights)
  
# Need the following for generating ragged array in STAN - indexing vector of sci_mu's based on their taxa identity (each slice will have a different length)
where_tx <- order(sci_to_tx) # i.e., give the positions in sci_to_tx in order of their taxa level (where in sci_to_tx is taxa level 1, followed by taxa level 2, etc)
# How many sci are in each taxa - need this to declare length of each sci_mu vector
n_sci_in_tx <- wild_dat_new_weights %>%
  select(sci, tx) %>%
  unique() %>%
  group_by(tx) %>%
  mutate(n_sci_in_tx = n()) %>%
  ungroup() %>%
  arrange(tx) %>%
  select(tx, n_sci_in_tx) %>%
  unique() %>%
  pull(n_sci_in_tx)
slice_where_tx <- cumsum(n_sci_in_tx) # These are the breaks in where_tx corresponding to each taxa level - need this to split up where_tx
slice_where_tx <- c(0, slice_where_tx)

# Set data for stan:
# REMINDER RE: PRIORS vs NO PRIORS - make sure STAN code below for defining/applying priors is allowed to run or commented out as needed
# NO PRIORS
# stan_data <- list(N = N,
#                   N_SCI = N_SCI,
#                   n_to_sci = n_to_sci,
#                   N_TX = N_TX,
#                   sci_to_tx = sci_to_tx,
#                   ghg = ghg,
#                   sci_w = sci_w,
#                   where_tx = where_tx,
#                   n_sci_in_tx = n_sci_in_tx,
#                   slice_where_tx = slice_where_tx)

# WITH PRIORS
# REMINDER RE: PRIORS vs NO PRIORS - make sure STAN code below for defining/applying priors is allowed to run or commented out as needed
stan_data <- list(N = N,
                  N_SCI = N_SCI,
                  n_to_sci = n_to_sci,
                  N_TX = N_TX,
                  sci_to_tx = sci_to_tx,
                  ghg = ghg,
                  sci_w = sci_w,
                  where_tx = where_tx,
                  n_sci_in_tx = n_sci_in_tx,
                  slice_where_tx = slice_where_tx,
                  priors = priors,
                  prior_vec_index = prior_vec_index)

# NORMAL DISTRIBUTION model - fed and non-fed
stan_no_na <- 'data {
  // indices
  int<lower=0> N;  // number of observations
  int N_TX; // number of taxa groups
  int N_SCI; // number of scientific names
  int n_to_sci[N]; // sciname index
  int sci_to_tx[N_SCI]; // taxa-group indices
  
  // data 
  vector<lower=0>[N] ghg; // data
  
  // PRIORS
  vector[8] priors;
  int prior_vec_index[8];
  
  // indices for slicing vectors for calculating weighted means
  vector<lower=0>[N_SCI] sci_w; // sci-level production weights
  int where_tx[N_SCI]; // order sci_to_tx by taxa-level
  int n_sci_in_tx[N_TX]; // number of sci in each taxa-level, ordered by taxa-level
  int slice_where_tx[N_TX + 1]; // breaks in where_tx by taxa-level
}
parameters {
  // On farm model
  vector<lower=0>[N_TX] tx_mu_ghg;
  vector<lower=0>[N_SCI] sci_mu_ghg;
  real<lower=0> tx_sigma_ghg;
  real<lower=0> sci_sigma_ghg;
}
model {
  // PRIORS
  tx_mu_ghg[prior_vec_index] ~ normal(priors, 1000);

  // example priors for dirichlet model for feed proportions
  // sci_phi defined as sci_phi[n_to_sci][K]
  // sci_phi[2][1] ~ normal(0.13, 5); // mean for Oncorhynhchus mykiss soy feed

  // weak priors on sigma
  tx_sigma_ghg ~ cauchy(0, 10000);
  sci_sigma_ghg ~ cauchy(0, 10000);

  // normal model for sci and taxa-level on-farm footrpint
  ghg ~ normal(sci_mu_ghg[n_to_sci], sci_sigma_ghg);
  sci_mu_ghg ~ normal(tx_mu_ghg[sci_to_tx], tx_sigma_ghg);
}
generated quantities {
  // Declare vectors for weightings
  vector[N_TX] tx_ghg_w;
  vector[N_SCI] sci_ghg_w;
 
  // Apply weightings
  sci_ghg_w = sci_mu_ghg .* sci_w; // WEIGHTED sci-level
  
  for (n_tx in 1:N_TX) {
    vector[n_sci_in_tx[n_tx]] sci_w_vec; // declare vector of sci_ghg in taxa-level n_tx
    sci_w_vec = sci_ghg_w[where_tx[slice_where_tx[n_tx]+1:slice_where_tx[n_tx+1]]]; // get all the sci_mu in taxa-level n_tx
    tx_ghg_w[n_tx] = sum(sci_w_vec); // sum sci_w_vec to get WEIGHTED tx-level outputs
  }
  
}'

no_na_mod <- stan_model(model_code = stan_no_na)
# Note: For Windows, apparently OK to ignore this warning message:
# Warning message:
#   In system(paste(CXX, ARGS), ignore.stdout = TRUE, ignore.stderr = TRUE) :
#   'C:/rtools40/usr/mingw_/bin/g++' not found

# Fit model:
# Set seed while testing
fit_no_na <- sampling(object = no_na_mod, 
                      data = stan_data, 
                      cores = 4, 
                      seed = "11729", 
                      iter = 2500, 
                      control = list(adapt_delta = 0.99, max_treedepth = 15))
#fit_no_na <- sampling(object = no_na_mod, data = stan_data, cores = 4, iter = 5000, control = list(adapt_delta = 0.99))
summary(fit_no_na)$summary

######################################################################################################
# STEP 3: OUTPUT RESULTS

###########################################################
# RESTARTING POINT
rm(list=ls()[!(ls() %in% c("datadir", "outdir", "impact", "set_allocation",
                           "wild_dat_new_weights", "fit_no_na"))])
save.image(file = file.path(outdir, paste(Sys.Date(), "_full-model-posterior_Wild-capture-ghg.RData", sep = "")))

###########################################################

# SET THEME
x <- seq(0, 1, length.out = 16)
base_color <- "#364F6B"
show_col(seq_gradient_pal(base_color, "white")(x)) # Get hexadecimals for other colors
interval_palette <- c("#9EA8B7", "#6A7A90", "#364F6B") # Order: light to dark
sci_plot_theme <- theme(title = element_text(size = 18),
                        axis.title.x = element_text(size = 16),
                        axis.text=element_text(size=10, color = "black"))
tx_plot_theme <- list(theme(title = element_text(size = 20),
                       axis.title.x = element_text(size = 20),
                       axis.text=element_text(size=20, color = "black"),
                       legend.position = "none"),
                      scale_color_manual(values = interval_palette))

units_for_plot = "kg CO2-eq per tonne"

# Key for naming taxa levels
# Get full taxa group names back
tx_index_key <- wild_dat_new_weights %>%
  group_by(clean_sci_name) %>%
  mutate(n_obs = n()) %>%
  ungroup() %>%
  select(taxa, tx) %>%
  unique() %>%
  arrange(taxa) 

# Key for naming sci levels
# Get full taxa group names back
sci_index_key <- wild_dat_new_weights %>%
  group_by(clean_sci_name) %>%
  mutate(n_obs = n()) %>%
  ungroup() %>%
  select(clean_sci_name, sci, taxa, tx, n_obs) %>%
  unique() %>%
  arrange(taxa)

# Use tidybayes + ggdist for finer control of aes mapping (instead of bayesplots) 
get_variables(fit_no_na)

######################################################################################################
# PLOT final outputs

# WEIGHTED taxa-level GHG
fit_no_na %>%
  spread_draws(tx_ghg_w[tx]) %>%
  median_qi(tx_ghg_w, .width = c(0.95, 0.8, 0.5)) %>%
  left_join(tx_index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  ggplot(aes(y = taxa, x = tx_ghg_w)) +
  geom_interval(aes(xmin = .lower, xmax = .upper)) +
  theme_classic() + 
  #coord_cartesian(xlim = c(0, 12500)) +
  tx_plot_theme + 
  labs(x = units_for_plot, y = "", title = "")
ggsave(filename = file.path(outdir, "plot_WILD-GHG-TAXA-LEVEL-WEIGHTED.png"), width = 11, height = 8.5)


# Same but as CSV output
fit_no_na %>%
  spread_draws(tx_ghg_w[tx]) %>%
  median_qi(tx_ghg_w, .width = c(0.95, 0.8, 0.5)) %>%
  left_join(tx_index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  rename(total_stressor = tx_ghg_w) %>%
  write.csv(file = file.path(outdir, "summary_WILD-GHG-TAXA-LEVEL-WEIGHTED.csv"), row.names = FALSE)

# UNWEIGHTED taxa-level GHG
fit_no_na %>%
  spread_draws(tx_mu_ghg[tx]) %>%
  median_qi(tx_mu_ghg, .width = c(0.95, 0.8, 0.5)) %>%
  left_join(tx_index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  # SET plants and bivalves to 0
  ggplot(aes(y = taxa, x = tx_mu_ghg)) +
  geom_interval(aes(xmin = .lower, xmax = .upper)) +
  theme_classic() + 
  tx_plot_theme + 
  labs(x = units_for_plot, y = "", title = "")
ggsave(filename = file.path(outdir, "plot_WILD-GHG-TAXA-LEVEL-UNWEIGHTED.png"), width = 11, height = 8.5)

# WEIGHTED sci-level GHG
fit_no_na %>%
  spread_draws(sci_ghg_w[sci]) %>%
  median_qi(sci_ghg_w, .width = c(0.95, 0.8, 0.5)) %>%
  left_join(sci_index_key, by = "sci") %>% # Join with index key to get sci and taxa names
  ggplot(aes(y = paste(taxa, clean_sci_name, sep = ""), x = sci_ghg_w)) +
  geom_interval(aes(xmin = .lower, xmax = .upper)) +
  scale_color_brewer() +
  theme_classic() + 
  sci_plot_theme + 
  labs(x = units_for_plot, y = "", title = "")
ggsave(filename = file.path(outdir, "plot_WILD-GHG-SCI-LEVEL-WEIGHTED.png"), width = 11, height = 8.5)

# UNWEIGHTED taxa-level GHG
fit_no_na %>%
  spread_draws(sci_mu_ghg[sci]) %>%
  median_qi(sci_mu_ghg, .width = c(0.95, 0.8, 0.5)) %>%
  left_join(sci_index_key, by = "sci") %>% # Join with index key to get sci and taxa names
  ggplot(aes(y = paste(taxa, clean_sci_name, sep = ""), x = sci_mu_ghg)) +
  geom_interval(aes(xmin = .lower, xmax = .upper)) +
  scale_color_brewer() +
  theme_classic() + 
  sci_plot_theme + 
  labs(x = units_for_plot, y = "", title = "")
ggsave(filename = file.path(outdir, "plot_WILD-GHG-SCI-LEVEL-UNWEIGHTED.png"), width = 11, height = 8.5)
