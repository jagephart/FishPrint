# Get Bayesian means of on-farm (non-feed) and off-farm (feed) CARBON impacts 

###################### REMINDER FOR CARBON: 
# OFF-FARM impacts should be given a ZERO when FCR = 0 and for all bivalves and plants
# ON-FARM carbon footprint calculated as:
# (Electricity use * country-specific GHG of electricity) + (Diesel * GHG of diesel) + (Petrol * GHG of petrol) + (Natural gas * GHG of natural gas)

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
# Step 1: Load and format data

# Load Data
lca_full_dat <- read.csv(file.path(datadir, "lca-dat-imputed-vars_rep-n-farms.csv"), fileEncoding="UTF-8-BOM")

# Get farm-associated carbon footprints data
# Add iso3c
electric_fp_dat <- read.csv(file.path(datadir, "electricity_GWP.csv")) %>%
  mutate(iso3c = countrycode(Country, origin = "country.name", destination = "iso3c"))
other_energy_fp_dat <- read.csv(file.path(datadir, "energy_carriers_impact_factors.csv"))

# units for diesel and petrol constants are in kg CO2-eq / L
diesel_fp <- other_energy_fp_dat %>% filter(Impact.category == "Global warming potential" & Input == "Diesel") %>% pull(Value)
petrol_fp <- other_energy_fp_dat %>% filter(Impact.category == "Global warming potential" & Input == "Petrol") %>% pull(Value)

# NOTE: natural gas is in kg CO2-eq / m3, need to convert to kg CO2-eq / L - multiply by 0.001 m3 / L
natgas_fp <- other_energy_fp_dat %>% filter(Impact.category == "Global warming potential" & Input == "Natural gas") %>% pull(Value) * 0.001

# LEFT OFF HERE: model converges when only running for fed species (i.e., filter(fcr != 0))
# FIX IT - testing model first on all fed species; then need to modify STAN to include these
# Format data for model:
lca_model_dat <- lca_full_dat %>%
  select(study_id, iso3c, clean_sci_name, taxa, intensity, system, 
         feed_soy, feed_crops, feed_fmfo, feed_animal, 
         fcr, 
         electric = Electricity_kwh,
         diesel = Diesel_L,
         petrol = Petrol_L,
         natgas = NaturalGas_L) %>%
  drop_na() %>% # NAs in these variables means they were not able to be imputed in step 2 
  #filter(fcr != 0)  %>% # FIRST TEST MODEL ON ALL FED SPECIES 
  # Multiply energy input by their Carbon (i.e., GHG) footprint, then sum across energy inputs
  # (Electricity use * country-specific GHG of electricity) + (Diesel * GHG of diesel) + (Petrol * GHG of petrol) + (Natural gas * GHG of natural gas)
  # Calculate electriciy GHG footprint
  left_join(electric_fp_dat %>% select(-Country), by = "iso3c") %>% 
  mutate(GWP_perkWh_kgCO2eq = if_else(is.na(GWP_perkWh_kgCO2eq), true = mean(GWP_perkWh_kgCO2eq, na.rm = TRUE), false = GWP_perkWh_kgCO2eq)) %>% # for studies with no country data, just use the average across countries
  mutate(electric_ghg = electric * GWP_perkWh_kgCO2eq) %>%
  # Calculate diesel GHG footprint
  mutate(diesel_ghg = diesel * diesel_fp) %>%
  # Calculate petrol GHG footprint
  mutate(petrol_ghg = petrol * petrol_fp) %>%
  # Calculate natural gas GHG footprint
  mutate(natgas_ghg = natgas * natgas_fp) %>%
  # Calculate sum total of GHG footprint
  mutate(total_ghg = electric_ghg + diesel_ghg + petrol_ghg + natgas_ghg) %>%
  arrange(clean_sci_name, taxa) %>%
  mutate(clean_sci_name = as.factor(clean_sci_name),
         sci = as.numeric(clean_sci_name),
         taxa = as.factor(taxa),
         tx = as.numeric(taxa),
         fed = if_else(fcr == 0, true = 0, false = 1)) 

# Set data, indices, constants, weights for STAN

# VARIABLE-SPECIFIC DATA:

# For on_farm ghg
farm <- lca_model_dat$total_ghg
# For FCR model:
fcr <- lca_model_dat$fcr
# For Feed proportion model:
K = 4
feed_weights <- lca_model_dat %>%
  select(feed_soy, feed_crops, feed_fmfo, feed_animal) %>%
  as.matrix()
# Also needed for feed proportion dirichlet: counts per sci name and counts per taxa group (also included as data in the model):
sci_kappa <- lca_model_dat %>% 
  group_by(sci) %>% 
  summarise(n_obs = n()) %>%
  ungroup() %>%
  arrange(sci) %>%
  pull(n_obs)
tx_kappa <- lca_model_dat %>% 
  group_by(tx) %>% 
  summarise(n_obs = n()) %>%
  ungroup() %>%
  arrange(tx) %>%
  pull(n_obs)

# INDICES:
N = nrow(lca_model_dat)
N_SCI <- length(unique(lca_model_dat$sci))
n_to_sci <- lca_model_dat$sci
n_to_tx <- lca_model_dat$tx
N_TX <- length(unique(lca_model_dat$tx))
sci_to_tx <- lca_model_dat %>%
  select(sci, tx) %>%
  unique() %>%
  pull(tx)
fed <- lca_model_dat$fed

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
  filter(clean_sci_name %in% lca_model_dat$clean_sci_name)
# Check that the order of sci names in both weights and data are the same (sum = 0)
sum(sci_prod_weights$clean_sci_name != unique(lca_model_dat$clean_sci_name))
sci_w <- sci_prod_weights$weighting

# Need the following for generating ragged array in STAN - indexing vector of sci_mu's based on their taxa identity (each slice will have a different length)
where_tx <- order(sci_to_tx) # i.e., give the positions in sci_to_tx in order of their taxa level (where in sci_to_tx is taxa level 1, followed by taxa level 2, etc)
# How many sci are in each taxa - need this to declare length of each sci_mu vector
n_sci_in_tx <- lca_model_dat %>%
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
                  farm = farm,
                  sci_kappa = sci_kappa,
                  tx_kappa = tx_kappa,
                  fp_constant = fp_constant,
                  sci_w = sci_w,
                  where_tx = where_tx,
                  n_sci_in_tx = n_sci_in_tx,
                  slice_where_tx = slice_where_tx)

# LEFT OFF HERE
# NORMAL DISTRIBUTION model - fed and non-fed
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
  vector[K] fp_constant;
  
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
    tx_feed_fp[n_tx] = tx_mu_fcr[n_tx] * sum(fp_constant .* tx_theta[n_tx]);
  }
  for (n_sci in 1:N_SCI) {
    sci_feed_fp[n_sci] = sci_mu_fcr[n_sci] * sum(fp_constant .* sci_theta[n_sci]);
  }

  // Sum off farm (feed) and on farm footprints (UNWEIGHTED sci and tx-level outputs)
  for (n_sci in 1:N_SCI) {
    sci_total_fp[n_sci] = sci_feed_fp[n_sci] + sci_mu_farm[n_sci];
  }
  for (n_tx in 1:N_TX) {
    tx_total_fp[n_tx] = tx_feed_fp[n_tx] + tx_mu_farm[n_tx];
  }

  // Apply weightings
  sci_total_fp_w = sci_total_fp .* sci_w; // WEIGHTED sci-level outputs

  for (n_tx in 1:N_TX){
    vector[n_sci_in_tx[n_tx]] sci_mu_w_vec; // declare vector of sci_mu in taxa-level n_tx
    sci_mu_w_vec = sci_total_fp_w[where_tx[slice_where_tx[n_tx]+1:slice_where_tx[n_tx+1]]]; // get all the sci_mu in taxa-level n_tx
    tx_total_fp_w[n_tx] = sum(sci_mu_w_vec); // sum sci_mu_w_vec to get WEIGHTED tx-level outputs
  }
}'


# NORMAL distribution hierarchical model - only fed species
# stan_no_na <- 'data {
#   // indices
#   int<lower=0> N;  // number of observations
#   int N_TX; // number of taxa groups
#   int N_SCI; // number of scientific names
#   int n_to_sci[N]; // sciname index
#   int sci_to_tx[N_SCI]; // taxa-group indices
# 
#   // data for feed footrpint
#   vector<lower=0>[N] fcr; // fcr data
#   int K; // number of feed types
#   simplex[K] feed_weights[N]; // array of observed feed weights simplexes
#   int sci_kappa[N_SCI]; // number of observations per sci-name
#   int tx_kappa[N_TX]; // number of observations per taxa group
#   
#   // data for on-farm footrpint
#   vector<lower=0>[N] farm; // data
# 
#   // constants to apply to feed footrpint
#   vector[K] fp_constant;
#   
#   // data for slicing vectors for calculating weighted means
#   vector<lower=0>[N_SCI] sci_w; // sci-level production weights
#   int where_tx[N_SCI]; // order sci_to_tx by taxa-level
#   int n_sci_in_tx[N_TX]; // number of sci in each taxa-level, ordered by taxa-level
#   int slice_where_tx[N_TX + 1]; // breaks in where_tx by taxa-level
# }
# parameters {
#   // FCR model
#   vector<lower=0>[N_TX] tx_mu_fcr;
#   vector<lower=0>[N_SCI] sci_mu_fcr;
#   real<lower=0> tx_sigma_fcr;
#   real<lower=0> sci_sigma_fcr;
#   //vector<lower=0>[N_TX] tx_sigma_fcr;
#   //vector<lower=0>[N_SCI] sci_sigma_fcr;
#   
#   // Feed proportion model:
#   simplex[K] sci_theta[N_SCI]; // vectors of estimated sci-level feed weight simplexes
#   simplex[K] tx_theta[N_TX];
#   
#   // On farm model
#   vector<lower=0>[N_TX] tx_mu_farm;
#   vector<lower=0>[N_SCI] sci_mu_farm;
#   real<lower=0> tx_sigma_farm;
#   real<lower=0> sci_sigma_farm;
# }
# transformed parameters {
#   // define params for dirichlet model for feed proportions
#   vector<lower=0>[K] sci_alpha[N_SCI];
#   vector<lower=0>[K] tx_alpha[N_TX];
# 
#   // dirichlet model reparameterization
#   // reparameterize alphas as a vector of means (phi) and counts (kappas)
#   // theta is expected value of mean feed weights
#   // kappa is strength of the prior measured in number of prior observations (minus K)
#   for (n_tx in 1:N_TX) {
#     tx_alpha[n_tx] = tx_kappa[n_tx] * tx_theta[n_tx];
#   }
#   for (n_sci in 1:N_SCI) {
#     sci_alpha[n_sci] = sci_kappa[n_sci] * sci_theta[n_sci];
#   }
# }
# model {
#   // example priors for gamma model for FCRs
#   // Put priors on mu and sigma (instead of shape and rate) since this is more intuitive:
#   //tx_mu ~ uniform(0, 100);
#   //tx_sigma ~ uniform(0, 100);
# 
#   // example priors for dirichlet model for feed proportions
#   // sci_phi defined as sci_phi[n_to_sci][K]
#   // sci_phi[2][1] ~ normal(0.13, 5); // mean for Oncorhynhchus mykiss soy feed
# 
#   // weak priors on sigma
#   tx_sigma_fcr ~ cauchy(0, 1);
#   sci_sigma_fcr ~ cauchy(0, 1);
#   tx_sigma_farm ~ cauchy(0, 10);
#   sci_sigma_farm ~ cauchy(0, 10);
# 
#   // likelihood
#   // normal model sci-name and taxa-level for FCR
#   fcr ~ normal(sci_mu_fcr[n_to_sci], sci_sigma_fcr);
#   sci_mu_fcr ~ normal(tx_mu_fcr[sci_to_tx], tx_sigma_fcr);
# 
#   // dirichlet model for feed proportions
#   // use data to estimate sci-level dirichlet shape param (alphas)
#   // LOOPING feed proportions:
#   for (n in 1:N){
#     feed_weights[n] ~ dirichlet(to_vector(sci_alpha[n_to_sci[n]]));
#   }
#   // use sci-level thetas to estimate taxa-level dirichlet shape param (alphas)
#   for (n_sci in 1:N_SCI){
#     sci_theta[n_sci] ~ dirichlet(to_vector(tx_alpha[sci_to_tx[n_sci]]));
#   }
#   
#   // normal model for sci and taxa-level on-farm footrpint
#   farm ~ normal(sci_mu_farm[n_to_sci], sci_sigma_farm);
#   sci_mu_farm ~ normal(tx_mu_farm[sci_to_tx], tx_sigma_farm);
# }
# generated quantities {
#   // Declare vectors for weightings
#   vector[N_TX] tx_feed_fp;
#   vector[N_SCI] sci_feed_fp;
#   vector[N_SCI] sci_total_fp; // unweighted sci-level
#   vector[N_TX] tx_total_fp; //
#   vector[N_SCI] sci_total_fp_w; // weighted sci-level
#   vector[N_TX] tx_total_fp_w; // weighted tx-level
# 
#   // Feed footrpint calculations
#   for (n_tx in 1:N_TX) {
#     tx_feed_fp[n_tx] = tx_mu_fcr[n_tx] * sum(fp_constant .* tx_theta[n_tx]);
#   }
#   for (n_sci in 1:N_SCI) {
#     sci_feed_fp[n_sci] = sci_mu_fcr[n_sci] * sum(fp_constant .* sci_theta[n_sci]);
#   }
# 
#   // Sum off farm (feed) and on farm footprints (UNWEIGHTED sci and tx-level outputs)
#   for (n_sci in 1:N_SCI) {
#     sci_total_fp[n_sci] = sci_feed_fp[n_sci] + sci_mu_farm[n_sci];
#   }
#   for (n_tx in 1:N_TX) {
#     tx_total_fp[n_tx] = tx_feed_fp[n_tx] + tx_mu_farm[n_tx];
#   }
# 
#   // Apply weightings
#   sci_total_fp_w = sci_total_fp .* sci_w; // WEIGHTED sci-level outputs
# 
#   for (n_tx in 1:N_TX){
#     vector[n_sci_in_tx[n_tx]] sci_mu_w_vec; // declare vector of sci_mu in taxa-level n_tx
#     sci_mu_w_vec = sci_total_fp_w[where_tx[slice_where_tx[n_tx]+1:slice_where_tx[n_tx+1]]]; // get all the sci_mu in taxa-level n_tx
#     tx_total_fp_w[n_tx] = sum(sci_mu_w_vec); // sum sci_mu_w_vec to get WEIGHTED tx-level outputs
#   }
# }'

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

#launch_shinystan(fit_no_na)
######################################################################################################
# RESTARTING POINT
# FIX IT - which objects to clear before saving?
# rm(list=ls()[!(ls() %in% c("datadir", "outdir", 
#                            "lca_dat_clean_groups", "feed_model_dat_categories",
#                            "full_feed_dat", "full_fcr_dat", "feed_footprint_dat", "fit_no_na"))])
#save.image(file.path(outdir, paste(Sys.Date(), "_full-model_", impact, "_", set_allocation, "-allocation_all-data-prior-to-plotting.RData", sep = "")))

# datadir <- "/Volumes/jgephart/BFA Environment 2/Data"
# outdir <- "/Volumes/jgephart/BFA Environment 2/Outputs"
# load(file.path(outdir, "<file-name>.RData"))
# set_allocation <- "Mass"
######################################################################################################
# PLOT RESULTS

# SET THEME
sci_plot_theme <- theme(title = element_text(size = 18),
                    axis.title.x = element_text(size = 16),
                    axis.text=element_text(size=14, color = "black"))
tx_plot_theme <- theme(title = element_text(size = 20),
                        axis.title.x = element_text(size = 20),
                        axis.text=element_text(size=20, color = "black"))

# Key for naming sci and taxa levels
# Get full taxa group names back
index_key <- carbon_footprint_dat %>%
  group_by(clean_sci_name) %>%
  mutate(n_obs = n()) %>%
  ungroup() %>%
  select(clean_sci_name, sci, taxa, tx, n_obs) %>%
  unique() %>%
  arrange(taxa) %>%
  mutate(taxa = as.character(taxa),
         full_taxa_name = case_when(taxa == "hypoph_carp" ~ "bighead/silverhead carp",
                                    taxa == "misc_marine" ~ "misc marine fishes",
                                    taxa == "misc_fresh" ~ "misc freshwater fishes",
                                    taxa == "misc_diad" ~ "misc diadromous fishes",
                                    taxa == "fresh_crust" ~ "freshwater crustaceans",
                                    TRUE ~ taxa),
         taxa = as.factor(taxa),
         full_taxa_name = as.factor(full_taxa_name))

# Key for naming [sci, feed] and [taxa, feed] levels
# Order of feeds: soy, crops, fmfo, animal
tx_feed_key <- carbon_footprint_dat %>%
  select(contains(c("taxa", "tx", "soy", "crops", "fmfo", "animal"))) %>%
  pivot_longer(cols = contains(c("soy", "crops", "fmfo", "animal")), names_to = "feed") %>%
  select(-value) %>%
  unique() %>%
  mutate(feed_index = case_when(str_detect(feed, "soy") ~ 1,
                                str_detect(feed, "crops") ~ 2,
                                str_detect(feed, "fmfo") ~ 3,
                                str_detect(feed, "animal") ~ 4)) %>%
  # Clean feed names
  mutate(feed = gsub(feed, pattern = "feed_", replacement = "")) %>%
  mutate(index = paste("[", tx, ",", feed_index, "]", sep = "")) %>%
  mutate(tx_theta_name = paste("theta[", taxa, ", ", feed, "]", sep = "")) %>%
  mutate(tx_alpha_name = paste("alpha[", taxa, ", ", feed, "]", sep = "")) %>%
  # IMPORTANT before replaceing param names: ARRANGE BY FEED, THEN TAXA NAME TO MATCH HOW NAMES ARE ARRANGED IN STANFIT OBJECT
  arrange(feed_index, tx)

sci_feed_key <- carbon_footprint_dat %>%
  select(contains(c("clean_sci_name", "taxa", "sci", "soy", "crops", "fmfo", "animal"))) %>%
  pivot_longer(cols = contains(c("soy", "crops", "fmfo", "animal")), names_to = "feed") %>%
  select(-value) %>%
  unique() %>%
  mutate(feed_index = case_when(str_detect(feed, "soy") ~ 1,
                                str_detect(feed, "crops") ~ 2,
                                str_detect(feed, "fmfo") ~ 3,
                                str_detect(feed, "animal") ~ 4)) %>%
  # Clean feed names
  mutate(feed = gsub(feed, pattern = "feed_", replacement = "")) %>%
  mutate(index = paste("[", sci, ",", feed_index, "]", sep = "")) %>%
  mutate(sci_theta_name = paste("theta[", clean_sci_name, ", ", feed, "]", sep = "")) %>%
  mutate(sci_alpha_name = paste("alpha[", clean_sci_name, ", ", feed, "]", sep = "")) %>%
  # IMPORTANT before replaceing param names: ARRANGE BY FEED, THEN TAXA NAME TO MATCH HOW NAMES ARE ARRANGED IN STANFIT OBJECT
  arrange(feed_index, sci)


# Use tidybayes + ggdist for finer control of aes mapping (instead of bayesplots) 
get_variables(fit_no_na)

######################################################################################################
# PLOT final outputs (total on + off farm impacts)

######################################################################################################
# PLOT other intermediate-level calculations
######################################################################################################
# WEIGHTED OUTPUTS:
# Sci-level feed footprints as point intervals:
# Carbon
fit_no_na %>%
  spread_draws(sci_total_fp_w[sci]) %>%
  median_qi(.width = 0.95) %>%
  left_join(index_key, by = "sci") %>% # Join with index key to get sci and taxa names
  mutate(clean_sci_name = fct_reorder(clean_sci_name, as.character(full_taxa_name))) %>%
  ggplot(aes(y = clean_sci_name, x = sci_total_fp_w, xmin = .lower, xmax = .upper, color = full_taxa_name)) +
  geom_pointinterval() +
  theme_classic() + 
  sci_plot_theme + 
  labs(x = "kg CO2-eq per tonne", y = "", title = "Carbon", color = "taxa group")
ggsave(filename = file.path(outdir, paste("plot_feed-impact_carbon_", set_allocation, "-allocation_sci-level-WEIGHTED.png", sep = "")), width = 11, height = 8.5)

# Taxa-level feed footprints as point intervals:
# Carbon
fit_no_na %>%
  spread_draws(tx_total_fp_w[tx]) %>%
  median_qi(.width = 0.95) %>%
  left_join(index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  ggplot(aes(y = full_taxa_name, x = tx_total_fp_w, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() +
  #stat_halfeye(aes(slab_fill = full_taxa_name)) +
  theme_classic() + 
  tx_plot_theme + 
  theme(legend.position = "none") +
  labs(x = "kg CO2-eq per tonne", y = "", title = "Carbon")
ggsave(filename = file.path(outdir, paste("plot_feed-impact_carbon_", set_allocation, "-allocation_taxa-level-WEIGHTED.png", sep = "")), width = 11, height = 8.5)


# UNWEIGHTED OUTPUTS:
# Sci-level feed footprints as point intervals:
# Carbon
fit_no_na %>%
  spread_draws(sci_feed_fp[sci]) %>%
  median_qi(.width = 0.95) %>%
  left_join(index_key, by = "sci") %>% # Join with index key to get sci and taxa names
  mutate(clean_sci_name = fct_reorder(clean_sci_name, as.character(full_taxa_name))) %>%
  ggplot(aes(y = clean_sci_name, x = sci_feed_fp, xmin = .lower, xmax = .upper, color = full_taxa_name)) +
  geom_pointinterval() +
  theme_classic() + 
  sci_plot_theme + 
  labs(x = "kg CO2-eq per tonne", y = "", title = "Carbon", color = "taxa group")
ggsave(filename = file.path(outdir, paste("plot_feed-impact_carbon_", set_allocation, "-allocation_sci-level.png", sep = "")), width = 11, height = 8.5)

# Taxa-level feed footprints as point intervals:
# Carbon
fit_no_na %>%
  spread_draws(tx_c_fp[tx]) %>%
  median_qi(.width = 0.95) %>%
  left_join(index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  ggplot(aes(y = full_taxa_name, x = tx_c_fp, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() +
  #stat_halfeye(aes(slab_fill = full_taxa_name)) +
  theme_classic() + 
  tx_plot_theme + 
  theme(legend.position = "none") +
  labs(x = "kg CO2-eq per tonne", y = "", title = "Carbon")
ggsave(filename = file.path(outdir, paste("plot_feed-impact_carbon_", set_allocation, "-allocation_taxa-level.png", sep = "")), width = 11, height = 8.5)

# FIX IT - fct_reorder not working for feed proportions?
# Sci-level theta (feed proportions)
# Make separate plot for each feed component
feed_component <- c("soy", "crops", "fmfo", "animal")
for (i in 1:length(feed_component)){
  plot_dat <- fit_no_na %>%
    spread_draws(sci_theta[sci, feed_index]) %>%
    median_qi(.width = 0.95) %>%
    left_join(sci_feed_key, by = c("sci" = "sci", "feed_index" = "feed_index")) %>% # Join with index key to get sci and taxa names
    filter(feed == feed_component[i]) %>%
    mutate(clean_sci_name = fct_reorder(clean_sci_name, as.character(taxa)))
    #mutate(clean_sci_name = fct_reorder(clean_sci_name, sci_theta))
    p <- ggplot(plot_dat, aes(y = clean_sci_name, x = sci_theta, xmin = .lower, xmax = .upper, color = taxa)) +
      geom_pointinterval() +
      #stat_halfeye(aes(slab_fill = full_taxa_name)) +
      coord_cartesian(xlim = c(0, 1)) +
      theme_classic() + 
      sci_plot_theme + 
      labs(x = "Feed Proportion", y = "", title = feed_component[i], color = "taxa group")
    print(p)
    ggsave(filename = file.path(outdir, paste("plot_feed-proportion_", feed_component[i], "_sci-level.png", sep = "")), width = 11, height = 8.5)
}

# Taxa-level theta (feed proportions)
# Make separate plot for each feed component
feed_component <- c("soy", "crops", "fmfo", "animal")
for (i in 1:length(feed_component)){
  plot_dat <- fit_no_na %>%
    spread_draws(tx_theta[tx, feed_index]) %>%
    median_qi(.width = 0.95) %>%
    left_join(tx_feed_key, by = c("tx" = "tx", "feed_index" = "feed_index")) %>% # Join with index key to get sci and taxa names
    filter(feed == feed_component[i])
  p <- ggplot(plot_dat, aes(y = taxa, x = tx_theta, xmin = .lower, xmax = .upper)) +
    geom_pointinterval() +
    #stat_halfeye(aes(slab_fill = full_taxa_name)) +
    coord_cartesian(xlim = c(0, 1)) +
    theme_classic() + 
    tx_plot_theme + 
    theme(legend.position = "none") +
    labs(x = "Feed Proportion", y = "", title = feed_component[i])
  print(p)
  ggsave(filename = file.path(outdir, paste("plot_feed-proportion_", feed_component[i], "_taxa-level.png", sep = "")), width = 11, height = 8.5)
}

######################################################################################################

# Taxa-level FCR
fit_no_na %>%
  spread_draws(tx_mu_fcr[tx]) %>%
  median_qi(.width = 0.95) %>%
  left_join(index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  ggplot(aes(y = full_taxa_name, x = tx_mu_fcr, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() +
  #stat_halfeye(aes(slab_fill = full_taxa_name)) +
  theme_classic() + 
  tx_plot_theme + 
  theme(legend.position = "none") +
  labs(x = "", y = "", title = "FCR")
ggsave(filename = file.path(outdir, paste("plot_fcr_taxa-level.png", sep = "")), width = 11, height = 8.5)




# OPTION 2: Taxa-level plots with color themes:
# If we want to mimic bayesplot color schemes, can get hexadecimal colors and input manually to stat_halfeye aesthetics
color_scheme_get("blue")
color_scheme_get("green")





