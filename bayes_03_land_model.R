# Author: Kelvin Gorospe
# Get Bayesian sci and taxa-level means of on-farm (non-feed) and off-farm (feed) LAND impacts 

###################### REMINDER FOR LAND: 
# ON-FARM impacts only apply to SYSTEM == PONDS or RECIRCULATING and TANKS (everything else is zero)
# OFF-FARM impacts EXCLUDE wherever FCR = 0 and bivalves and plants
# CALCULATION for ON-FARM land = yield

######################
# STEP 0: SET DIRECTORIES, LOAD PACKAGES
rm(list = ls())

# Libraries for processing and analyses
library(tidyverse)
library(rstan)
library(data.table)
library(countrycode)
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
# STEP 1: LOAD AND FORMAT DATA

# Load Data
#lca_full_dat <- read.csv(file.path(outdir, "lca-dat-imputed-vars_rep-sqrt-n-farms_live-weight.csv"), fileEncoding="UTF-8-BOM")
lca_full_dat <- read.csv(file.path(outdir, "lca-dat-imputed-vars_rep-sqrt-n-farms_edible-weight.csv"), fileEncoding="UTF-8-BOM")

# Format data for model:
lca_model_dat <- lca_full_dat %>%
  select(study_id, iso3c, clean_sci_name, taxa, intensity, system, 
         feed_soy, feed_crops, feed_fmfo, feed_animal, 
         fcr, 
         yield = Yield_m2_per_t) %>%
  # Option 1: TEST ON "UNIFIED" Dataset first: only system == ponds or recirculating tanks & FCR != 0 - later need to add conditionals to STAN
  # filter(fcr != 0)  %>%
  # filter(system %in% c("Ponds", "Recirculating and tanks")) %>%
  # OPTION 2: INCLUDE FED AND NON-FED SPECIES BUT mutate feed proportions to be the average within it's clean_sci_name; otherwise, give it an arbitrary simplex (0.25 per component) to avoid STAN error for simplexes that don't sum to 1
  group_by(clean_sci_name) %>%
  mutate(feed_soy = if_else(fcr==0, true = mean(feed_soy), false = feed_soy),
         feed_crops = if_else(fcr==0, true = mean(feed_crops), false = feed_crops),
         feed_fmfo = if_else(fcr==0, true = mean(feed_fmfo), false = feed_fmfo),
         feed_animal = if_else(fcr==0, true = mean(feed_animal), false = feed_animal)) %>%
  ungroup() %>%
  # Need to re-normalize everything (using the mean value does not guarantee that proportion sums to 1)
  mutate(feed_sum = feed_soy + feed_crops + feed_fmfo + feed_animal) %>%
  mutate(feed_soy = if_else(feed_sum != 0, true = feed_soy / feed_sum, false = feed_soy), # need conditional feed_sum != 0 - otherwise species with FCR == 0 and no sci-level average will renormalize using denominator of 0
         feed_crops = if_else(feed_sum != 0, true = feed_crops / feed_sum, false = feed_crops),
         feed_fmfo = if_else(feed_sum != 0, true = feed_fmfo / feed_sum, false = feed_fmfo),
         feed_animal = if_else(feed_sum != 0, true = feed_animal / feed_sum, false = feed_animal)) %>%
  # Give the remaining things 0.25
  mutate(feed_soy = if_else(fcr == 0 & feed_soy == 0 & feed_crops == 0 & feed_fmfo == 0 & feed_animal == 0, true = 0.25, false = feed_soy),
         feed_crops = if_else(fcr == 0 & feed_crops == 0 & feed_fmfo == 0 & feed_animal == 0, true = 0.25, false = feed_crops),
         feed_fmfo = if_else(fcr == 0 & feed_fmfo == 0 & feed_animal == 0, true = 0.25, false = feed_fmfo),
         feed_animal = if_else(fcr == 0 & feed_animal == 0, true = 0.25, false = feed_animal)) %>%
  # Set yield = 0 (including where yield is NA) for studies that are not Ponds or Recirculating/tanks
  mutate(yield = if_else(system %in% c("Ponds", "Recirculating and tanks")==FALSE, true = 0, false = yield)) %>%
  drop_na() %>% # After fixing NA's in yield, drop all other NAs
  # LAST FORMATING STEP - always arrange by clean_sci_name
  arrange(clean_sci_name) %>%
  mutate(clean_sci_name = as.factor(clean_sci_name),
         sci = as.numeric(clean_sci_name),
         taxa = as.factor(taxa),
         tx = as.numeric(taxa)) 

# Get priors on taxa-level FCR
# Can ignore warning: NAs introduced by coercion (inserts NAs for blank cells)
source("Functions.R")
priors_csv <- clean_priors("Priors - Aquaculture.csv") %>%
  select(contains(c("taxa", "FCR"))) %>%
  arrange(taxa) # Arrange by taxa so that index matches tx in lca_model_dat

# Format priors for STAN
# can't pass NAs into STAN - drop NAs but keep track of vector positions
prior_vec_index <- which(is.na(priors_csv$Ave.FCR)==FALSE)
priors <- priors_csv$Ave.FCR[prior_vec_index]

#####################
# Set data, indices, constants, weights for STAN

# VARIABLE-SPECIFIC DATA:

# For on farm land
land <- lca_model_dat$yield
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

# WEIGHTS:
# Get sci-level weightings for generating taxa-level quantities:
# IMPORTANT arrange by clean_sci_name so that the order matches data
sci_prod_weights <- read.csv(file.path(outdir, "aqua_prod_weightings.csv")) %>%
  arrange(clean_sci_name)
# Drop sci-names not found in the data
sci_prod_weights <- sci_prod_weights %>%
  filter(clean_sci_name %in% lca_model_dat$clean_sci_name)
# Check that the order of sci names in both weights and data are the same (sum = 0)
sum(sci_prod_weights$clean_sci_name != unique(lca_model_dat$clean_sci_name))
sci_w <- sci_prod_weights$prod_weighting

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

# FEED IMPACT CONSTANTS: need to include both land and water constants in the model

# Choose allocation method
#set_allocation <- "Mass"
#set_allocation <- "Gross energy content"
set_allocation <- "Economic"

fp_dat <- read.csv(file.path(outdir, "weighted_feed_fp.csv")) %>%
  filter(Allocation == set_allocation) %>%
  mutate(ave_stressor_per_tonne = ave_stressor * 1000) # Multiply by 1000 to convert to kg CO2 per tonne

# ORDER OF feed_weights data vector: soy, crops, fmfo, animal
head(feed_weights)
# ORDER of footprint data vector should match this:
set_fp_order <- c("Soy", "Crop", "Fishery", "Livestock")

land_feed_fp <- fp_dat %>%
  filter(Impact.category == "Land use") %>% 
  arrange(match(Input.type, set_fp_order)) %>% # Match index and arrange by custom order
  select(ave_stressor_per_tonne) %>%
  as.matrix() %>%
  c()

######################################################################################################
# STEP 2: RUN STAN MODEL

# Set data for stan:
# REMINDER RE: PRIORS vs NO PRIORS - make sure STAN code below for defining/applying priors is allowed to run or commented out as needed
# SPECIFIC FOR LAND MODEL: When NOT using priors, add weak priors on sigma_land (see STAN code)
# NO PRIORS
# stan_data <- list(N = N,
#                   N_SCI = N_SCI,
#                   n_to_sci = n_to_sci,
#                   N_TX = N_TX,
#                   sci_to_tx = sci_to_tx,
#                   fcr = fcr,
#                   K = K,
#                   feed_weights = feed_weights,
#                   land = land,
#                   sci_kappa = sci_kappa,
#                   tx_kappa = tx_kappa,
#                   land_feed_fp = land_feed_fp,
#                   sci_w = sci_w,
#                   where_tx = where_tx,
#                   n_sci_in_tx = n_sci_in_tx,
#                   slice_where_tx = slice_where_tx)

# WITH FCR PRIORS
# REMINDER RE: PRIORS vs NO PRIORS - make sure STAN code below for defining/applying priors is allowed to run or commented out as needed
# REMINDER: Attempted using priors on both FCR and LAND - model does not converge, stick with priors on just FCR
stan_data <- list(N = N,
                  N_SCI = N_SCI,
                  n_to_sci = n_to_sci,
                  N_TX = N_TX,
                  sci_to_tx = sci_to_tx,
                  fcr = fcr,
                  K = K,
                  feed_weights = feed_weights,
                  land = land,
                  sci_kappa = sci_kappa,
                  tx_kappa = tx_kappa,
                  land_feed_fp = land_feed_fp,
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

  // data for feed footrpint
  vector<lower=0>[N] fcr; // fcr data
  int K; // number of feed types
  simplex[K] feed_weights[N]; // array of observed feed weights simplexes
  int sci_kappa[N_SCI]; // number of observations per sci-name
  int tx_kappa[N_TX]; // number of observations per taxa group
  
  // PRIORS
  vector[11] priors; // priors on FCR
  int prior_vec_index[11];
  
  // data for on-farm footrpint
  vector<lower=0>[N] land; // data

  // constants to apply to feed footrpint
  vector[K] land_feed_fp;
  
  // data for slicing vectors for calculating weighted means
  vector<lower=0>[N_SCI] sci_w; // sci-level production weights
  int where_tx[N_SCI]; // order sci_to_tx by taxa-level
  int n_sci_in_tx[N_TX]; // number of sci in each taxa-level, ordered by taxa-level
  int slice_where_tx[N_TX + 1]; // breaks in where_tx by taxa-level
}
parameters {
  // FCR model
  vector[N_TX] tx_mu_fcr; // putting lower=0 bounds will cause mu_fcr to skew positive when zero (eg, plants, bivalves)
  vector[N_SCI] sci_mu_fcr;
  real<lower=0> tx_sigma_fcr;
  real<lower=0> sci_sigma_fcr;
  
  // Feed proportion model:
  simplex[K] sci_theta[N_SCI]; // vectors of estimated sci-level feed weight simplexes
  simplex[K] tx_theta[N_TX];
  
  // On farm model
  vector[N_TX] tx_land_farm; // putting lower=0 bounds will cause tx_land_farm to skew positive when zero
  vector[N_SCI] sci_land_farm;
  real<lower=0> tx_sigma_land;
  real<lower=0> sci_sigma_land;
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
  // PRIORS
  tx_mu_fcr[prior_vec_index] ~ normal(priors, 1); // priors on FCR

  // weak priors on sigma
  tx_sigma_fcr ~ cauchy(0, 1);
  sci_sigma_fcr ~ cauchy(0, 1);
  //tx_sigma_land ~ cauchy(0, 10000); // only need priors on sigma_land when NOT using priors on mu_fcr
  //sci_sigma_land ~ cauchy(0, 10000);

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
  
  // normal model for sci and taxa-level on-farm land footrpint
  land ~ normal(sci_land_farm[n_to_sci], sci_sigma_land);
  sci_land_farm ~ normal(tx_land_farm[sci_to_tx], tx_sigma_land);
}
generated quantities {
  // Declare vectors for weightings
  vector[N_TX] tx_land_feed;
  vector[N_SCI] sci_land_feed;
  vector[N_SCI] sci_land_total; // unweighted
  vector[N_TX] tx_land_total; //
  vector[N_TX] tx_land_total_w; // weighted
  vector[N_TX] tx_land_feed_w;
  vector[N_TX] tx_land_farm_w;
  vector[N_SCI] sci_land_feed_w;
  vector[N_SCI] sci_land_farm_w;

  // Feed LAND footprint calculations
  for (n_tx in 1:N_TX) {
    tx_land_feed[n_tx] = tx_mu_fcr[n_tx] * sum(land_feed_fp .* tx_theta[n_tx]);
  }
  for (n_sci in 1:N_SCI) {
    sci_land_feed[n_sci] = sci_mu_fcr[n_sci] * sum(land_feed_fp .* sci_theta[n_sci]);
  }

  // Sum off farm (feed) and on farm footprints (UNWEIGHTED sci and tx-level outputs)
  for (n_sci in 1:N_SCI) {
    sci_land_total[n_sci] = sci_land_feed[n_sci] + sci_land_farm[n_sci];
  }
  for (n_tx in 1:N_TX) {
    tx_land_total[n_tx] = tx_land_feed[n_tx] + tx_land_farm[n_tx];
  }

  // Apply weightings to LAND
  sci_land_feed_w = sci_land_feed .* sci_w; // WEIGHTED sci-level off-farm impacts
  sci_land_farm_w = sci_land_farm .* sci_w; // WEIGHTED sci-level on-farm impacts

  for (n_tx in 1:N_TX) {
    vector[n_sci_in_tx[n_tx]] land_feed_vec; // declare vector of sci_feed in taxa-level n_tx
    vector[n_sci_in_tx[n_tx]] land_farm_vec; // declare vector of sci_farm in taxa-level n_tx

    land_feed_vec = sci_land_feed_w[where_tx[slice_where_tx[n_tx]+1:slice_where_tx[n_tx+1]]]; // get all the sci_mu in taxa-level n_tx
    tx_land_feed_w[n_tx] = sum(land_feed_vec); // sum land_feed_vec to get WEIGHTED tx-level outputs

    land_farm_vec = sci_land_farm_w[where_tx[slice_where_tx[n_tx]+1:slice_where_tx[n_tx+1]]]; // get all the sci_mu in taxa-level n_tx
    tx_land_farm_w[n_tx] = sum(land_farm_vec); // sum land_farm_vec to get WEIGHTED tx-level outputs
  }

  // Sum of weighted on and off-farm impacts
  for (n_tx in 1:N_TX){
    tx_land_total_w[n_tx] = tx_land_feed_w[n_tx] + tx_land_farm_w[n_tx];
  }
}'

no_na_mod <- stan_model(model_code = stan_no_na)
# Note: For Windows, apparently OK to ignore this warning message:
# Warning message:
#   In system(paste(CXX, ARGS), ignore.stdout = TRUE, ignore.stderr = TRUE) :
#   'C:/rtools40/usr/mingw_/bin/g++' not found

# Fit model:
# Set seed while testing
start_sampling <- Sys.time()
fit_no_na <- sampling(object = no_na_mod, 
                      data = stan_data, 
                      cores = 4, 
                      iter = 2000, 
                      seed = "11729",
                      control = list(adapt_delta = 0.99, max_treedepth = 15))
#fit_no_na <- sampling(object = no_na_mod, data = stan_data, cores = 4, iter = 5000, control = list(adapt_delta = 0.99))
end_sampling <- Sys.time()
end_sampling - start_sampling
summary(fit_no_na)$summary
# Index specific parameters by name:
#(summary(fit_no_na)$summary)[rownames(summary(fit_no_na)$summary) %in% c("tx_sigma_land", "sci_sigma_land", "tx_sigma_fcr", "sci_sigma_fcr"),]

#launch_shinystan(fit_no_na)

######################################################################################################
# STEP 3: OUTPUT RESULTS

impact <- "Land Use"

###########################################################
# RESTARTING POINT
rm(list=ls()[!(ls() %in% c("datadir", "outdir", "impact", "set_allocation",
                           "lca_model_dat", "fit_no_na"))])
save.image(file = file.path(outdir, paste(Sys.Date(), "_full-model-posterior_", impact, "_", set_allocation, "-allocation.RData", sep = "")))

###########################################################
# PLOT final LAND outputs (off-farm, on-farm, and total impacts)

# SET THEME
x <- seq(0, 1, length.out = 16)
base_color <- "#57D182"
library(scales)
show_col(seq_gradient_pal(base_color, "white")(x)) # Get hexadecimals for other colors
# COMPLEMENTARY GREEN (Azote guidelines)
interval_palette <- c("#B7EBC4", "#8BDEA3", "#57D182")
# PRIMARY GREEN:
#interval_palette <- c("#D9EAB2", "#C2DD86", "#A9D158") # Order: light to dark
full_taxa_name_order <- c("plants", "bivalves", "shrimp", "misc marine fishes", "milkfish", "salmon", "misc diadromous fishes", "trout", "tilapia", "catfish", "misc carps", "bighead/silverhead carp")

sci_plot_theme <- theme(title = element_text(size = 18),
                        axis.title.x = element_text(size = 16),
                        axis.text=element_text(size=14, color = "black"))
tx_plot_theme <- list(theme(title = element_text(size = 20),
                       axis.title.x = element_text(size = 20),
                       axis.text=element_text(size=20, color = "black"),
                       legend.position = "none"),
                      scale_color_manual(values = interval_palette))

# Set units:
units_for_plot = bquote('m'^2*'a per tonne')

# Key for naming taxa levels
# Get full taxa group names back
tx_index_key <- lca_model_dat %>%
  group_by(clean_sci_name) %>%
  mutate(n_obs = n()) %>%
  ungroup() %>%
  select(taxa, tx) %>%
  unique() %>%
  arrange(taxa) %>%
  mutate(taxa = as.character(taxa),
         full_taxa_name = case_when(taxa == "hypoph_carp" ~ "bighead/silverhead carp",
                                    taxa == "misc_marine" ~ "misc marine fishes",
                                    taxa == "misc_fresh" ~ "misc freshwater fishes",
                                    taxa == "misc_diad" ~ "misc diadromous fishes",
                                    taxa == "oth_carp" ~ "misc carps",
                                    taxa == "fresh_crust" ~ "freshwater crustaceans",
                                    TRUE ~ taxa),
         taxa = as.factor(taxa),
         full_taxa_name = as.factor(full_taxa_name))

# Key for naming sci levels
# Get full taxa group names back
sci_index_key <- lca_model_dat %>%
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
                                    taxa == "oth_carp" ~ "misc carps",
                                    taxa == "fresh_crust" ~ "freshwater crustaceans",
                                    TRUE ~ taxa),
         taxa = as.factor(taxa),
         full_taxa_name = as.factor(full_taxa_name))

# Key for naming [sci, feed] and [taxa, feed] levels
# Order of feeds: soy, crops, fmfo, animal
tx_feed_key <- lca_model_dat %>%
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
  arrange(feed_index, tx) %>%
  mutate(taxa = as.character(taxa),
         full_taxa_name = case_when(taxa == "hypoph_carp" ~ "bighead/silverhead carp",
                                    taxa == "misc_marine" ~ "misc marine fishes",
                                    taxa == "misc_fresh" ~ "misc freshwater fishes",
                                    taxa == "misc_diad" ~ "misc diadromous fishes",
                                    taxa == "oth_carp" ~ "misc carps",
                                    taxa == "fresh_crust" ~ "freshwater crustaceans",
                                    TRUE ~ taxa),
         taxa = as.factor(taxa),
         full_taxa_name = as.factor(full_taxa_name))

sci_feed_key <- lca_model_dat %>%
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
  arrange(feed_index, sci) %>%
  mutate(taxa = as.character(taxa),
         full_taxa_name = case_when(taxa == "hypoph_carp" ~ "bighead/silverhead carp",
                                    taxa == "misc_marine" ~ "misc marine fishes",
                                    taxa == "misc_fresh" ~ "misc freshwater fishes",
                                    taxa == "misc_diad" ~ "misc diadromous fishes",
                                    taxa == "oth_carp" ~ "misc carps",
                                    taxa == "fresh_crust" ~ "freshwater crustaceans",
                                    TRUE ~ taxa),
         taxa = as.factor(taxa),
         full_taxa_name = as.factor(full_taxa_name))


# Use tidybayes + ggdist for finer control of aes mapping (instead of bayesplots) 
get_variables(fit_no_na)

###########################################################
## WEIGHTED (only taxa level is relevant - i.e., sum of weighted sci-levels)
###########################################################
# Mean off-farm (feed) impact taxa-level
fit_no_na %>%
  spread_draws(tx_land_feed_w[tx]) %>%
  median_qi(.width = c(0.95, 0.8, 0.5)) %>%
  left_join(tx_index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  # SET plants and bivalves to 0
  mutate(tx_land_feed_w = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = tx_land_feed_w),
         .lower = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = .lower),
         .upper = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = .upper)) %>%
  # REORDER taxa axis
  mutate(full_taxa_name = fct_relevel(full_taxa_name, full_taxa_name_order)) %>%
  #mutate(full_taxa_name = fct_reorder(full_taxa_name, tx_feed_fp_w)) %>%
  ggplot(aes(y = full_taxa_name, x = tx_land_feed_w)) +
  geom_interval(aes(xmin = .lower, xmax = .upper)) +
  geom_point(x = 0, y = "bivalves") +
  geom_point(x = 0, y = "plants") +
  theme_classic() + 
  tx_plot_theme + 
  labs(x = units_for_plot, y = "", title = "Mean off-farm (feed) impact")
ggsave(filename = file.path(outdir, paste("plot_", impact, "_", set_allocation, "-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.png", sep = "")), width = 11, height = 8.5)
#ggsave(filename = file.path(outdir, "plot_geom_pointinterval.png"))

# Same but as CSV output
fit_no_na %>%
  spread_draws(tx_land_feed_w[tx]) %>%
  median_qi(.width = 0.95) %>%
  left_join(tx_index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  mutate(tx_land_feed_w = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = tx_land_feed_w),
         .lower = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = .lower),
         .upper = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = .upper)) %>%
  rename(off_farm = tx_land_feed_w) %>%
  write.csv(file = file.path(outdir, paste("summary_", impact, "_", set_allocation, "-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.csv", sep = "")), row.names = FALSE)

# Mean on-farm impact taxa-level:
# First - set all taxa that have no PONDS or RECIRCULATING and TANKS to be ZERO, all the rest are true distributions
tx_null_farm <- lca_model_dat %>%
  group_by(taxa) %>%
  mutate(ponds = if_else(system %in% c("Ponds", "Recirculating and tanks"), true = 1, false = 0),
         total_ponds = sum(ponds),
         n_sci = n()) %>%
  ungroup() %>%
  filter(total_ponds == 0) %>%
  pull(taxa) %>%
  unique()

fit_no_na %>%
  spread_draws(tx_land_farm_w[tx]) %>%
  median_qi(.width = c(0.95, 0.8, 0.5)) %>%
  left_join(tx_index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  # SET tx_null_farm to 0
  mutate(tx_land_farm_w = if_else(taxa %in% tx_null_farm, true = 0, false = tx_land_farm_w),
         .lower = if_else(taxa %in% tx_null_farm, true = 0, false = .lower),
         .upper = if_else(taxa %in% tx_null_farm, true = 0, false = .upper)) %>%
  # REORDER taxa axis
  mutate(full_taxa_name = fct_relevel(full_taxa_name, full_taxa_name_order)) %>%
  #mutate(full_taxa_name = fct_reorder(full_taxa_name, tx_feed_fp_w)) %>%
  ggplot(aes(y = full_taxa_name, x = tx_land_farm_w)) +
  geom_interval(aes(xmin = .lower, xmax = .upper)) +
  geom_point(x = 0, y = "bivalves") +
  geom_point(x = 0, y = "plants") +
  theme_classic() + 
  tx_plot_theme + 
  labs(x = units_for_plot, y = "", title = "Mean on-farm impact")
ggsave(filename = file.path(outdir, paste("plot_", impact, "_", set_allocation, "-allocation_ON-FARM-TAXA-LEVEL-WEIGHTED.png", sep = "")), width = 11, height = 8.5)

# Same but as CSV output
fit_no_na %>%
  spread_draws(tx_land_farm_w[tx]) %>%
  median_qi(.width = 0.95) %>%
  left_join(tx_index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  # SET tx_null_farm to 0
  mutate(tx_land_farm_w = if_else(taxa %in% tx_null_farm, true = 0, false = tx_land_farm_w),
         .lower = if_else(taxa %in% tx_null_farm, true = 0, false = .lower),
         .upper = if_else(taxa %in% tx_null_farm, true = 0, false = .upper)) %>%
  rename(on_farm = tx_land_farm_w) %>%
  write.csv(file = file.path(outdir, paste("summary_", impact, "_", set_allocation, "-allocation_ON-FARM-TAXA-LEVEL-WEIGHTED.csv", sep = "")), row.names = FALSE)

# Mean total impact taxa-level
fit_no_na %>%
  spread_draws(tx_land_total_w[tx]) %>%
  median_qi(.width = c(0.95, 0.8, 0.5)) %>%
  left_join(tx_index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  # Set plant and bivalves distributions to 0 (both are 0 for on and off farm)
  mutate(tx_land_total_w = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = tx_land_total_w),
         .lower = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = .lower),
         .upper = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = .upper)) %>%
  # REORDER taxa axis
  mutate(full_taxa_name = fct_relevel(full_taxa_name, full_taxa_name_order)) %>%
  #mutate(full_taxa_name = fct_reorder(full_taxa_name, tx_feed_fp_w)) %>%
  ggplot(., aes(y = full_taxa_name, x = tx_land_total_w)) +
  geom_interval(aes(xmin = .lower, xmax = .upper)) +
  geom_point(x = 0, y = "bivalves") +
  geom_point(x = 0, y = "plants") +
  theme_classic() + 
  tx_plot_theme + 
  labs(x = units_for_plot, y = "", title = "Total (on and off-farm) impact", color = "taxa group")
ggsave(filename = file.path(outdir, paste("plot_", impact, "_", set_allocation, "-allocation_TOTAL-IMPACT-TAXA-LEVEL-WEIGHTED.png", sep = "")), width = 11, height = 8.5)

# Same but as CSV output
fit_no_na %>%
  spread_draws(tx_land_total_w[tx]) %>%
  median_qi(.width = 0.95) %>%
  left_join(tx_index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  # Set plant and bivalves distributions to 0 (both are 0 for on and off farm)
  mutate(tx_land_total_w = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = tx_land_total_w),
         .lower = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = .lower),
         .upper = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = .upper)) %>%
  rename(total_stressor = tx_land_total_w) %>%
  write.csv(file = file.path(outdir, paste("summary_", impact, "_", set_allocation, "-allocation_TOTAL-IMPACT-TAXA-LEVEL-WEIGHTED.csv", sep = "")), row.names = FALSE)


###########################################################
## PLOT WEIGHTED SCI-LEVEL ESTIMATES ANYWAY TO COMPARE WITH UNWEIGHTED ESTIMATES
# Mean 
fit_no_na %>%
  spread_draws(sci_land_feed_w[sci]) %>%
  median_qi(.width = 0.95) %>%
  left_join(sci_index_key, by = "sci") %>% # Join with index key to get sci and taxa names
  mutate(clean_sci_name = fct_reorder(clean_sci_name, as.character(full_taxa_name))) %>%
  ggplot(aes(y = clean_sci_name, x = sci_land_feed_w, xmin = .lower, xmax = .upper, color = full_taxa_name)) +
  geom_pointinterval() +
  theme_classic() + 
  sci_plot_theme + 
  labs(x = units_for_plot, y = "", title = "Weighted sci-level off-farm impacts", color = "taxa group")
ggsave(filename = file.path(outdir, paste("plot_", impact, "_", set_allocation, "-allocation_OFF-FARM-SCI-LEVEL-WEIGHTED.png", sep = "")), width = 11, height = 8.5)


fit_no_na %>%
  spread_draws(sci_land_farm_w[sci]) %>%
  median_qi(.width = 0.95) %>%
  left_join(sci_index_key, by = "sci") %>% # Join with index key to get sci and taxa names
  mutate(clean_sci_name = fct_reorder(clean_sci_name, as.character(full_taxa_name))) %>%
  ggplot(aes(y = clean_sci_name, x = sci_land_farm_w, xmin = .lower, xmax = .upper, color = full_taxa_name)) +
  geom_pointinterval() +
  theme_classic() + 
  sci_plot_theme + 
  labs(x = units_for_plot, y = "", title = "Weighted sci-level on-farm impacts", color = "taxa group")
ggsave(filename = file.path(outdir, paste("plot_", impact, "_", set_allocation, "-allocation_ON-FARM-SCI-LEVEL-WEIGHTED.png", sep = "")), width = 11, height = 8.5)


###########################################################
## UNWEIGHTED
###########################################################
# Mean off-farm (feed) impact sci-level
fit_no_na %>%
  spread_draws(sci_land_feed[sci]) %>%
  median_qi(.width = 0.95) %>%
  left_join(sci_index_key, by = "sci") %>% # Join with index key to get sci and taxa names
  # SET plants and bivalves to 0
  mutate(sci_land_feed = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = sci_land_feed),
         .lower = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = .lower),
         .upper = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = .upper)) %>%
  mutate(clean_sci_name = fct_reorder(clean_sci_name, as.character(full_taxa_name))) %>%
  ggplot(aes(y = clean_sci_name, x = sci_land_feed, xmin = .lower, xmax = .upper, color = full_taxa_name)) +
  geom_pointinterval() +
  theme_classic() + 
  sci_plot_theme + 
  labs(x = units_for_plot, y = "", title = "Mean off-farm (feed) impact", color = "taxa group")
ggsave(filename = file.path(outdir, paste("plot_", impact, "_", set_allocation, "-allocation_OFF-FARM-SCI-LEVEL-UNWEIGHTED.png", sep = "")), width = 11, height = 8.5)

# Taxa-level
fit_no_na %>%
  spread_draws(tx_land_feed[tx]) %>%
  median_qi(.width = c(0.95, 0.8, 0.5)) %>%
  left_join(tx_index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  # SET plants and bivalves to 0
  mutate(tx_land_feed = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = tx_land_feed),
         .lower = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = .lower),
         .upper = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = .upper)) %>%
  # REORDER taxa axis
  mutate(full_taxa_name = fct_relevel(full_taxa_name, full_taxa_name_order)) %>%
  #mutate(full_taxa_name = fct_reorder(full_taxa_name, tx_feed_fp_w)) %>%
  ggplot(aes(y = full_taxa_name, x = tx_land_feed)) +
  geom_interval(aes(xmin = .lower, xmax = .upper)) +
  geom_point(x = 0, y = "bivalves") +
  geom_point(x = 0, y = "plants") +
  theme_classic() + 
  tx_plot_theme + 
  labs(x = units_for_plot, y = "", title = "Mean off-farm (feed) impact")
ggsave(filename = file.path(outdir, paste("plot_", impact, "_", set_allocation, "-allocation_OFF-FARM-TAXA-LEVEL-UNWEIGHTED.png", sep = "")), width = 11, height = 8.5)

###########################################################
# Mean on-farm impact sci-level
# NOTE: Mean on-farm estimates are bounded by zero as declared in STAN model (FCR was left unbounded to allow for mean = 0, so off-farm is also unabounded)

# First, set all sci names that have no Ponds or Recirculating/tanks to 0
sci_null_farm <- lca_model_dat %>%
  group_by(clean_sci_name) %>%
  mutate(ponds = if_else(system %in% c("Ponds", "Recirculating and tanks"), true = 1, false = 0),
         total_ponds = sum(ponds),
         n_sci = n()) %>%
  ungroup() %>%
  filter(total_ponds == 0) %>%
  pull(clean_sci_name) %>%
  unique()

fit_no_na %>%
  spread_draws(sci_land_farm[sci]) %>%
  median_qi(.width = 0.95) %>%
  left_join(sci_index_key, by = "sci") %>% # Join with index key to get sci and taxa names
  # SET sci_null_farm to 0
  mutate(sci_land_farm = if_else(clean_sci_name %in% sci_null_farm, true = 0, false = sci_land_farm),
         .lower = if_else(clean_sci_name %in% sci_null_farm, true = 0, false = .lower),
         .upper = if_else(clean_sci_name %in% sci_null_farm, true = 0, false = .upper)) %>%
  mutate(clean_sci_name = fct_reorder(clean_sci_name, as.character(full_taxa_name))) %>%
  ggplot(aes(y = clean_sci_name, x = sci_land_farm, xmin = .lower, xmax = .upper, color = full_taxa_name)) +
  geom_pointinterval() +
  theme_classic() + 
  sci_plot_theme + 
  labs(x = units_for_plot, y = "", title = "Mean on-farm impact", color = "taxa group")
ggsave(filename = file.path(outdir, paste("plot_", impact, "_", set_allocation, "-allocation_ON-FARM-SCI-LEVEL-UNWEIGHTED.png", sep = "")), width = 11, height = 8.5)

# Taxa-level
fit_no_na %>%
  spread_draws(tx_land_farm[tx]) %>%
  median_qi(.width = c(0.95, 0.8, 0.5)) %>%
  left_join(tx_index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  # SET tx_null_farm to 0
  mutate(tx_land_farm = if_else(taxa %in% tx_null_farm, true = 0, false = tx_land_farm),
         .lower = if_else(taxa %in% tx_null_farm, true = 0, false = .lower),
         .upper = if_else(taxa %in% tx_null_farm, true = 0, false = .upper)) %>%
  # REORDER taxa axis
  mutate(full_taxa_name = fct_relevel(full_taxa_name, full_taxa_name_order)) %>%
  #mutate(full_taxa_name = fct_reorder(full_taxa_name, tx_feed_fp_w)) %>%
  ggplot(aes(y = full_taxa_name, x = tx_land_farm)) +
  geom_interval(aes(xmin = .lower, xmax = .upper)) +
  geom_point(x = 0, y = "bivalves") +
  geom_point(x = 0, y = "plants") +
  theme_classic() + 
  tx_plot_theme + 
  labs(x = units_for_plot, y = "", title = "Mean on-farm impact", color = "taxa group")
ggsave(filename = file.path(outdir, paste("plot_", impact, "_", set_allocation, "-allocation_ON-FARM-TAXA-LEVEL-UNWEIGHTED.png", sep = "")), width = 11, height = 8.5)

###########################################################
# Mean unweighted total (on + off-farm) impact sci-level
fit_no_na %>%
  spread_draws(sci_land_total[sci]) %>%
  median_qi(.width = 0.95) %>%
  left_join(sci_index_key, by = "sci") %>% # Join with index key to get sci and taxa names
  # Set plant and bivalves distributions to 0 (both are 0 for on and off farm)
  mutate(sci_land_total = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = sci_land_total),
         .lower = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = .lower),
         .upper = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = .upper)) %>%
  mutate(clean_sci_name = fct_reorder(clean_sci_name, as.character(full_taxa_name))) %>%
  ggplot(aes(y = clean_sci_name, x = sci_land_total, xmin = .lower, xmax = .upper, color = full_taxa_name)) +
  geom_pointinterval() +
  theme_classic() + 
  sci_plot_theme + 
  labs(x = units_for_plot, y = "", title = "Total (on and off-farm) impact", color = "taxa group") 
ggsave(filename = file.path(outdir, paste("plot_", impact, "_", set_allocation, "-allocation_TOTAL-IMPACT-SCI-LEVEL-UNWEIGHTED.png", sep = "")), width = 11, height = 8.5)

# Taxa-level
fit_no_na %>%
  spread_draws(tx_land_total[tx]) %>%
  median_qi(.width = c(0.95, 0.8, 0.5)) %>%
  left_join(tx_index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  # Set plant and bivalves distributions to 0 (both are 0 for on and off farm)
  mutate(tx_land_total = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = tx_land_total),
         .lower = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = .lower),
         .upper = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = .upper)) %>%
  # REORDER taxa axis
  mutate(full_taxa_name = fct_relevel(full_taxa_name, full_taxa_name_order)) %>%
  #mutate(full_taxa_name = fct_reorder(full_taxa_name, tx_feed_fp_w)) %>%
  ggplot(aes(y = full_taxa_name, x = tx_land_total)) +
  geom_interval(aes(xmin = .lower, xmax = .upper)) +
  geom_point(x = 0, y = "bivalves") +
  geom_point(x = 0, y = "plants") +
  theme_classic() + 
  tx_plot_theme + 
  labs(x = units_for_plot, y = "", title = "Total (on and off-farm) impact", color = "taxa group")
ggsave(filename = file.path(outdir, paste("plot_", impact, "_", set_allocation, "-allocation_TOTAL-IMPACT-TAXA-LEVEL-UNWEIGHTED.png", sep = "")), width = 11, height = 8.5)

######################################################################################################