# Get Bayesian mean of sci and taxa-level land and water impacts 

###################### REMINDER FOR LAND: 
# ON-FARM impacts only apply to SYSTEM == PONDS or RECIRCULATING and TANKS
# OFF-FARM impacts EXCLUDE wherever FCR = 0 and bivalves and plants

###################### REMINDER FOR WATER: 
# CALCULATION for ON-FARM impacts: MULTIPLY mean_evap_mm / 1000 (to get to m2) * LAND * grow out period in days (means per taxa group)/365
# i.e., just multiplying LAND with constants within the "generated quantities" section

# ON-FARM IMPACTS only apply to:
# Freshwater taxa only - i.e., oth_carp, catfish, hypoph_carp, tilapia, trout, fresh crust
# System = Pond or Recirculating Tanks (Land already deals with filtering this out)
# The rest get a ZERO for ON-FARM IMPACTS


######################
# Step 0: Set directories, load packages
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

######################
# Step 1: Set Data, indices, weights, constants


# Load Data
lca_full_dat <- read.csv(file.path(datadir, "lca-dat-imputed-vars_rep-n-farms.csv"), fileEncoding="UTF-8-BOM") %>%

# First to model data
lca_model_dat <- lca_full_dat %>%
  filter(is.na(fcr)==FALSE) %>% # NA here means that data point was not imputed in step 2 - remove these from model
  select(study_id, clean_sci_name, taxa, intensity, system, 
         feed_soy, feed_crops, feed_fmfo, feed_animal, 
         fcr, 
         yield = Yield_m2_per_t) %>%
  # TEST ON "UNIFIED" Dataset first - system == ponds or recirculating tanks & FCR != 0 - later need to add conditionals to STAN
  filter(fcr != 0)  %>%
  filter(system %in% c("Ponds", "Recirculating and tanks"))
  
# Load on-farm evaporative water loss (NOAA data - country-level mean of monthly climatological means 1981-2010)
evap_clim <- read.csv(file.path(datadir, "20201222_clim_summarise_by_country.csv")) %>%
  mutate(iso3c = countrycode(admin, origin = "country.name", destination = "iso3c")) %>%
  select(-X) %>%
  drop_na()

# FIX IT - make sure this works with; just copied over from water model
# Merge evap dat with Countries in lca data, only need mean_evap_mm and ID columns: study_id, Country, iso3c
# Then merge with feed data to get full model dataset
evap_dat <- lca_dat_clean_groups %>%
  left_join(evap_clim, by = "iso3c") %>%
  # Deal with some countries manually: Calculate mean for Germany, Denmark; Assign global mean to N/A
  mutate(mean_evap_mm = case_when(Country == "Germany, Denmark" ~ evap_clim %>% filter(iso3c %in% c("DEU", "DNK")) %>% pull(mean_evap_mm) %>% mean(),
                                  Country == "Scotland" ~ evap_clim %>% filter(iso3c == "GBR") %>% pull(mean_evap_mm),
                                  Country == "N/A" ~ mean(evap_clim$mean_evap_mm),
                                  TRUE ~ mean_evap_mm)) %>%
  select(study_id, Country, iso3c, clean_sci_name, taxa, intensity = Intensity, system = Production_system_group, mean_evap_mm) %>%
  # Merge with feed data for full model dataset
  full_join(feed_dat_merge, by = intersect(names(feed_dat_merge), names(.))) %>%
  full_join(fcr_dat_merge, by = intersect(names(feed_dat_merge), names(fcr_dat_merge)))

# VARIABLE-SPECIFIC DATA:

# For on_farm ghg
farm_c <- lca_full_dat$total_ghg
# For FCR model:
fcr <- carbon_footprint_dat$fcr
# For Feed proportion model:
K = 4
feed_weights <- carbon_footprint_dat %>%
  select(feed_soy, feed_crops, feed_fmfo, feed_animal) %>%
  as.matrix()
# Also needed for feed proportion dirichlet: counts per sci name and counts per taxa group (also included as data in the model):
sci_kappa <- carbon_footprint_dat %>% 
  group_by(sci) %>% 
  summarise(n_obs = n()) %>%
  ungroup() %>%
  arrange(sci) %>%
  pull(n_obs)
tx_kappa <- carbon_footprint_dat %>% 
  group_by(tx) %>% 
  summarise(n_obs = n()) %>%
  ungroup() %>%
  arrange(tx) %>%
  pull(n_obs)

# INDICES:
N = nrow(carbon_footprint_dat)
N_SCI <- length(unique(carbon_footprint_dat$sci))
n_to_sci <- carbon_footprint_dat$sci
n_to_tx <- carbon_footprint_dat$tx
N_TX <- length(unique(carbon_footprint_dat$tx))
sci_to_tx <- carbon_footprint_dat %>%
  select(sci, tx) %>%
  unique() %>%
  pull(tx)

# FEED IMPACT CONSTANTS:
# Choose allocation method
# Choose CARBON for Impact.category
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

# WEIGHTS:
# Get sci-level weightings for generating taxa-level quantities:
# IMPORTANT arrange by clean_sci_name so that the order matches data
sci_prod_weights <- read.csv(file.path(datadir, "aqua_prod_weightings.csv")) %>%
  arrange(clean_sci_name)
# Drop sci-names not found in the data
sci_prod_weights <- sci_prod_weights %>%
  filter(clean_sci_name %in% carbon_footprint_dat$clean_sci_name)
# Check that the order of sci names in both weights and data are the same (sum = 0)
sum(sci_prod_weights$clean_sci_name != unique(carbon_footprint_dat$clean_sci_name))
sci_w <- sci_prod_weights$weighting

# Need the following for generating ragged array in STAN - indexing vector of sci_mu's based on their taxa identity (each slice will have a different length)
where_tx <- order(sci_to_tx) # i.e., give the positions in sci_to_tx in order of their taxa level (where in sci_to_tx is taxa level 1, followed by taxa level 2, etc)
# How many sci are in each taxa - need this to declare length of each sci_mu vector
n_sci_in_tx <- carbon_footprint_dat %>%
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
stan_data <- list(N = N,
                  N_SCI = N_SCI, 
                  n_to_sci = n_to_sci,
                  N_TX = N_TX,
                  sci_to_tx = sci_to_tx,
                  fcr = fcr,
                  K = K,
                  feed_weights = feed_weights,
                  farm_c = farm_c,
                  sci_kappa = sci_kappa,
                  tx_kappa = tx_kappa,
                  fp_constant = fp_constant,
                  sci_w = sci_w,
                  where_tx = where_tx,
                  n_sci_in_tx = n_sci_in_tx,
                  slice_where_tx = slice_where_tx)

# NORMAL distribution hierarchical model
stan_no_na <- 'data {
  // indices
  int<lower=0> N;  // number of observations
  int N_TX; // number of taxa groups
  int N_SCI; // number of scientific names
  int n_to_sci[N]; // sciname index
  int sci_to_tx[N_SCI]; // taxa-group indices

  // data for feed footrpint
  vector<lower=0>[N] fcr; // fcr data
  int K; // number of feed types
  simplex[K] feed_weights[N]; // array of observed feed weights simplexes
  int sci_kappa[N_SCI]; // number of observations per sci-name
  int tx_kappa[N_TX]; // number of observations per taxa group
  
  // data for on-farm footrpint
  vector<lower=0>[N] farm; // data

  // constants to apply to feed footrpint
  //vector[K] fp_constant;
  
  // data for slicing vectors for calculating weighted means
  vector<lower=0>[N_SCI] sci_w; // sci-level production weights
  int where_tx[N_SCI]; // order sci_to_tx by taxa-level
  int n_sci_in_tx[N_TX]; // number of sci in each taxa-level, ordered by taxa-level
  int slice_where_tx[N_TX + 1]; // breaks in where_tx by taxa-level
}
parameters {
  // FCR model
  vector<lower=0>[N_TX] tx_mu_fcr;
  vector<lower=0>[N_SCI] sci_mu_fcr;
  real<lower=0> tx_sigma_fcr;
  real<lower=0> sci_sigma_fcr;
  //vector<lower=0>[N_TX] tx_sigma_fcr;
  //vector<lower=0>[N_SCI] sci_sigma_fcr;
  
  // Feed proportion model:
  simplex[K] sci_theta[N_SCI]; // vectors of estimated sci-level feed weight simplexes
  simplex[K] tx_theta[N_TX];
  
  // On farm model
  vector<lower=0>[N_TX] tx_mu_farm;
  vector<lower=0>[N_SCI] sci_mu_farm;
  real<lower=0> tx_sigma_farm;
  real<lower=0> sci_sigma_farm;
}
transformed parameters {
  // define params for dirichlet model for feed proportions
  vector<lower=0>[K] sci_alpha[N_SCI];
  vector<lower=0>[K] tx_alpha[N_TX];

  // dirichlet model reparameterization
  // reparameterize alphas as a vector of means (phi) and counts (kappas)
  // theta is expected value of mean feed weights
  // kappa is strength of the prior measured in number of prior observations (minus K)
  for (n_tx in 1:N_TX) {
    tx_alpha[n_tx] = tx_kappa[n_tx] * tx_theta[n_tx];
  }
  for (n_sci in 1:N_SCI) {
    sci_alpha[n_sci] = sci_kappa[n_sci] * sci_theta[n_sci];
  }
}
model {
  // example priors for gamma model for FCRs
  // Put priors on mu and sigma (instead of shape and rate) since this is more intuitive:
  //tx_mu ~ uniform(0, 100);
  //tx_sigma ~ uniform(0, 100);

  // example priors for dirichlet model for feed proportions
  // sci_phi defined as sci_phi[n_to_sci][K]
  // sci_phi[2][1] ~ normal(0.13, 5); // mean for Oncorhynhchus mykiss soy feed

  // weak priors on sigma
  tx_sigma_fcr ~ cauchy(0, 1);
  sci_sigma_fcr ~ cauchy(0, 1);
  tx_sigma_farm ~ cauchy(0, 10);
  sci_sigma_farm ~ cauchy(0, 10);

  // likelihood
  // normal model sci-name and taxa-level for FCR
  fcr ~ normal(sci_mu_fcr[n_to_sci], sci_sigma_fcr);
  sci_mu_fcr ~ normal(tx_mu_fcr[sci_to_tx], tx_sigma_fcr);

  // dirichlet model for feed proportions
  // use data to estimate sci-level dirichlet shape param (alphas)
  // LOOPING feed proportions:
  for (n in 1:N){
    feed_weights[n] ~ dirichlet(to_vector(sci_alpha[n_to_sci[n]]));
  }
  // use sci-level thetas to estimate taxa-level dirichlet shape param (alphas)
  for (n_sci in 1:N_SCI){
    sci_theta[n_sci] ~ dirichlet(to_vector(tx_alpha[sci_to_tx[n_sci]]));
  }
  
  // normal model for sci and taxa-level on-farm footrpint
  farm ~ normal(sci_mu_farm[n_to_sci], sci_sigma_farm);
  sci_mu_farm ~ normal(tx_mu_farm[sci_to_tx], tx_sigma_farm);
}
generated quantities {
  // Declare vectors for weightings
  vector[N_TX] tx_feed_fp;
  vector[N_SCI] sci_feed_fp;
  vector[N_SCI] sci_total_fp; // unweighted sci-level
  vector[N_TX] tx_total_fp; //
  vector[N_SCI] sci_total_fp_w; // weighted sci-level
  vector[N_TX] tx_total_fp_w; // weighted tx-level

  // Feed footrpint calculations
  for (n_tx in 1:N_TX) {
    tx_feed_fp[n_tx] = tx_mu_feed[n_tx] * sum(fp_constant .* tx_theta[n_tx]);
  }
  for (n_sci in 1:N_SCI) {
    sci_feed_fp[n_sci] = sci_mu_feed[n_sci] * sum(fp_constant .* sci_theta[n_sci]);
  }

  // Sum off farm (feed) and on farm footprints (UNWEIGHTED sci and tx-level outputs)
  for (n_sci in 1:N_SCI) {
    sci_total_fp[n_sci] = sci_feed_fp[n_sci] + sci_mu_farm[n_sci];
  }
  for (n_tx in 1:N_TX) {
    tx_total_fp[n_tx] = tx_feed_fp[n_tx] + tx_mu_farm[n_tx];
  }

  // Apply weightings
  //sci_total_fp_w = sci_total_fp .* sci_w; // WEIGHTED sci-level outputs

  //for (n_tx in 1:N_TX){
  //  vector[n_sci_in_tx[n_tx]] sci_mu_w_vec; // declare vector of sci_mu in taxa-level n_tx
  //  sci_mu_w_vec = sci_total_fp_w[where_tx[slice_where_tx[n_tx]+1:slice_where_tx[n_tx+1]]]; // get all the sci_mu in taxa-level n_tx
  //  tx_total_fp_w[n_tx] = sum(sci_mu_w_vec); // sum sci_mu_w_vec to get WEIGHTED tx-level outputs
  //}
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

