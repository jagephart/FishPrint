# Get Bayesian sci and taxa-level means of on-farm (non-feed) and off-farm (feed) CARBON impacts 

###################### REMINDER FOR CARBON: 
# OFF-FARM impacts should be ZERO when FCR = 0 (i.e., for full taxa groups bivalves and plants, and other sci-names within other taxa groups)
# ON-FARM carbon footprint calculated as:
# (Electricity use * country-specific GHG of electricity) + (Diesel * GHG of diesel) + (Petrol * GHG of petrol) + (Natural gas * GHG of natural gas)

######################
# STEP 0: SET DIRECTORIES, LOAD PACKAGES
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
# STEP 1: LOAD AND FORMAT DATA

# Load Data
lca_full_dat <- read.csv(file.path(datadir, "2021-01-06_lca-dat-imputed-vars_rep-sqrt-n-farms.csv"), fileEncoding="UTF-8-BOM")

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

# Format data for model:
lca_model_dat <- lca_full_dat %>%
  select(study_id, iso3c, clean_sci_name, taxa, intensity, system, 
         feed_soy, feed_crops, feed_fmfo, feed_animal, 
         fcr, 
         electric = Electricity_kwh,
         diesel = Diesel_L,
         petrol = Petrol_L,
         natgas = NaturalGas_L) %>%
  drop_na() %>% # just in case
  # OPTION 1: TEST MODEL ON SPECIES THAT ARE FED 
  #filter(fcr != 0)  %>% 
  # OPTION 2: INCLUDE FED AND NON-FED SPECIES BUT mutate feed proportions to be the average within it's clean_sci_name; otherwise, give it an arbitrary simplex (0.25 per component) to avoid STAN error for simplexes that don't sum to 1
  # FIX IT - in terms of on-farm footprint, this is OK because feed proportions are multiplied by FCR == 0 so on-farm footprint for these studies will be 0
  # BUT this will affect the pooled sci and taxa level feed proportions since they enter as 0.25
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
  # LAST FORMATING STEP - always arrange by clean_sci_name
  arrange(clean_sci_name) %>%
  mutate(clean_sci_name = as.factor(clean_sci_name),
         sci = as.numeric(clean_sci_name),
         taxa = as.factor(taxa),
         tx = as.numeric(taxa)) 



# FEED IMPACT CONSTANTS:
# Loop through entire code by allocation_method (get new feed impact constants and recreate everything after lca_model_dat)
#for (i in c("Mass", "Gross energy content", "Economic")){
# FIX IT: Delete next line, use complete list for loop
#for (i in c("Gross energy content", "Economic")){
#set_allocation <- i

# Or Choose allocation method manually
set_allocation <- "Mass"
#set_allocation <- "Gross energy content"
#set_allocation <- "Economic"

# Choose CARBON for Impact.category
impact <- "Global warming potential" # i.e., Carbon impact

fp_dat <- read.csv(file.path(datadir, "weighted_feed_fp.csv")) %>%
  filter(Allocation == set_allocation) %>%
  mutate(ave_stressor_per_tonne = ave_stressor * 1000) # Multiply by 1000 to convert to kg CO2 per tonne

# ORDER OF feed_weights data vector: soy, crops, fmfo, animal
#head(feed_weights) # created later
# ORDER of footprint data vector should match this:
set_fp_order <- c("Soy", "Crop", "Fishery", "Livestock")

fp_constant <- fp_dat %>%
  filter(Impact.category == impact) %>% 
  arrange(match(Input.type, set_fp_order)) %>% # Match index and arrange by custom order
  select(ave_stressor_per_tonne) %>%
  as.matrix() %>%
  c()

# Get priors on taxa-level FCR
# Can ignore warning: NAs introduced by coercion (inserts NAs for blank cells)
source("Functions.R")
priors_csv <- clean_priors("Priors - Nonfeed.csv") %>%
  select(contains(c("taxa", "FCR"))) %>%
  arrange(taxa) # Arrange by taxa so that index matches tx in lca_model_dat

# Format priors for STAN
# can't pass NAs into STAN - drop NAs but keep track of vector positions
prior_vec_index <- which(is.na(priors_csv$Ave.FCR)==FALSE)
priors <- priors_csv$Ave.FCR[prior_vec_index]
# priors_1 <- priors_csv$Ave.FCR[1]
# priors_4 <- priors_csv$Ave.FCR[4]
# priors_6_12 <- priors_csv$Ave.FCR[6:12]

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


######################################################################################################
# STEP 2: RUN STAN MODEL

# Set data for stan:
# NO PRIORS
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

# WITH PRIORS
# stan_data <- list(N = N,
#                   N_SCI = N_SCI,
#                   n_to_sci = n_to_sci,
#                   N_TX = N_TX,
#                   sci_to_tx = sci_to_tx,
#                   fcr = fcr,
#                   K = K,
#                   feed_weights = feed_weights,
#                   farm = farm,
#                   sci_kappa = sci_kappa,
#                   tx_kappa = tx_kappa,
#                   fp_constant = fp_constant,
#                   sci_w = sci_w,
#                   where_tx = where_tx,
#                   n_sci_in_tx = n_sci_in_tx,
#                   slice_where_tx = slice_where_tx,
#                   priors = priors,
#                   prior_vec_index = prior_vec_index)


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
  //vector[11] priors;
  //int prior_vec_index[11];
  
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
  vector[N_TX] tx_mu_fcr; // putting lower=0 bounds will cause mu_fcr to skew positive when zero (eg, plants, bivalves)
  vector[N_SCI] sci_mu_fcr;
  real<lower=0> tx_sigma_fcr;
  real<lower=0> sci_sigma_fcr;

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
  // PRIORS
  //tx_mu_fcr[prior_vec_index] ~ normal(priors, 1);
  
  // example priors for dirichlet model for feed proportions
  // sci_phi defined as sci_phi[n_to_sci][K]
  // sci_phi[2][1] ~ normal(0.13, 5); // mean for Oncorhynhchus mykiss soy feed

  // weak priors on sigma
  tx_sigma_fcr ~ cauchy(0, 1);
  sci_sigma_fcr ~ cauchy(0, 1);
  tx_sigma_farm ~ cauchy(0, 1000);
  sci_sigma_farm ~ cauchy(0, 10000);

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
  vector[N_SCI] sci_total_fp; // unweighted
  vector[N_TX] tx_total_fp; //
  //vector[N_SCI] sci_total_fp_w; // only need this if applying weights directly to the total impact
  vector[N_TX] tx_total_fp_w; // weighted
  vector[N_TX] tx_feed_fp_w;
  vector[N_TX] tx_farm_fp_w;
  vector[N_SCI] sci_feed_fp_w;
  vector[N_SCI] sci_mu_farm_w;
  
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
  sci_feed_fp_w = sci_feed_fp .* sci_w; // WEIGHTED sci-level off-farm impacts
  sci_mu_farm_w = sci_mu_farm .* sci_w; // WEIGHTED sci-level on-farm impacts
  
  for (n_tx in 1:N_TX) {
    vector[n_sci_in_tx[n_tx]] sci_feed_w_vec; // declare vector of sci_feed in taxa-level n_tx
    vector[n_sci_in_tx[n_tx]] sci_farm_w_vec; // declare vector of sci_farm in taxa-level n_tx
    
    sci_feed_w_vec = sci_feed_fp_w[where_tx[slice_where_tx[n_tx]+1:slice_where_tx[n_tx+1]]]; // get all the sci_mu in taxa-level n_tx
    tx_feed_fp_w[n_tx] = sum(sci_feed_w_vec); // sum sci_feed_w_vec to get WEIGHTED tx-level outputs
    
    sci_farm_w_vec = sci_mu_farm_w[where_tx[slice_where_tx[n_tx]+1:slice_where_tx[n_tx+1]]]; // get all the sci_mu in taxa-level n_tx
    tx_farm_fp_w[n_tx] = sum(sci_farm_w_vec); // sum sci_farm_w_vec to get WEIGHTED tx-level outputs
  }
  
  // Sum of weighted on and off-farm impacts
  for (n_tx in 1:N_TX){
    tx_total_fp_w[n_tx] = tx_feed_fp_w[n_tx] + tx_farm_fp_w[n_tx];
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
                      seed = "11729", 
                      iter = 2000, 
                      control = list(adapt_delta = 0.99, max_treedepth = 15))
#fit_no_na <- sampling(object = no_na_mod, data = stan_data, cores = 4, iter = 5000, control = list(adapt_delta = 0.99))
end_sampling <- Sys.time()
end_sampling - start_sampling
summary(fit_no_na)$summary

#launch_shinystan(fit_no_na)
######################################################################################################
# RESTARTING POINT
rm(list=ls()[!(ls() %in% c("datadir", "outdir", "impact", "set_allocation",
                            "lca_model_dat", "fit_no_na"))])
save.image(file = file.path(outdir, paste(Sys.Date(), "_full-model-posterior_", impact, "_", set_allocation, "-allocation.RData", sep = "")))

######################################################################################################
# STEP 3: OUTPUT RESULTS

# SET THEME
x <- seq(0, 1, length.out = 16)
base_color <- "#364F6B"
library(scales)
show_col(seq_gradient_pal(base_color, "white")(x)) # Get hexadecimals for other colors
interval_palette <- c("#9EA8B7", "#6A7A90", "#364F6B") # Order: light to dark
full_taxa_name_order <- c("plants", "bivalves", "shrimp", "misc marine fishes", "milkfish", "salmon", "misc diadromous fishes", "trout", "tilapia", "catfish", "misc carps", "bighead/silverhead carp")
sci_plot_theme <- theme(title = element_text(size = 18),
                    axis.title.x = element_text(size = 16),
                    axis.text=element_text(size=14, color = "black"))
tx_plot_theme <- list(theme(title = element_text(size = 20),
                        axis.title.x = element_text(size = 20),
                        axis.text=element_text(size=20, color = "black"),
                       legend.position = "none"),
                      scale_color_manual(values = interval_palette))

#scale_color_manual(values = c("red", "green", "blue")) + # This works


# Set units:
if (impact == "Global warming potential") {
  units_for_plot = "kg CO2-eq per tonne"
} else if (impact == "Freshwater eutrophication") {
  units_for_plot = "kg P-eq per tonne"
} else if (impact == "Marine eutrophication") {
  units_for_plot = "kg N-eq per tonne"
} else if (impact == "Land use") {
  units_for_plot = bquote('m'^2~'a per tonne')
} else if (impact == "Water consumption") {
  units_for_plot = bquote('m'^3~'per tonne')
}

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

######################################################################################################
# Output plots and summaries (off-farm, on-farm, and total impacts)

###########################################################
## WEIGHTED (only taxa level is relevant - i.e., sum of weighted sci-levels)
###########################################################
# Mean off-farm (feed) impact taxa-level
fit_no_na %>%
  spread_draws(tx_feed_fp_w[tx]) %>%
  median_qi(.width = c(0.95, 0.8, 0.5)) %>%
  left_join(tx_index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  mutate(tx_feed_fp_w = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = tx_feed_fp_w),
         .lower = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = .lower),
         .upper = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = .upper)) %>%
  # REORDER taxa axis
  mutate(full_taxa_name = fct_relevel(full_taxa_name, full_taxa_name_order)) %>%
  #mutate(full_taxa_name = fct_reorder(full_taxa_name, tx_feed_fp_w)) %>%
  ggplot(aes(y = full_taxa_name, x = tx_feed_fp_w)) +
  geom_interval(aes(xmin = .lower, xmax = .upper)) +
  geom_point(x = c(0), y = c("bivalves")) +
  geom_point(x = 0, y = "plants") +
  theme_classic() + 
  tx_plot_theme + 
  labs(x = units_for_plot, y = "", title = "Mean off-farm (feed) impact")
ggsave(filename = file.path(outdir, paste("plot_", impact, "_", set_allocation, "-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.png", sep = "")), width = 11, height = 8.5)

# Same but as CSV output
fit_no_na %>%
  spread_draws(tx_feed_fp_w[tx]) %>%
  median_qi(.width = 0.95) %>%
  left_join(tx_index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  mutate(tx_feed_fp_w = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = tx_feed_fp_w),
         .lower = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = .lower),
         .upper = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = .upper)) %>%
  rename(off_farm = tx_feed_fp_w) %>%
  write.csv(file = file.path(outdir, paste("summary_", impact, "_", set_allocation, "-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.csv", sep = "")), row.names = FALSE)

# Mean on-farm impact taxa-level
fit_no_na %>%
  spread_draws(tx_farm_fp_w[tx]) %>%
  median_qi(.width = c(0.95, 0.8, 0.5)) %>%
  left_join(tx_index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  # REORDER taxa axis
  mutate(full_taxa_name = fct_relevel(full_taxa_name, full_taxa_name_order)) %>%
  #mutate(full_taxa_name = fct_reorder(full_taxa_name, tx_farm_fp_w)) %>%
  ggplot(aes(y = full_taxa_name, x = tx_farm_fp_w)) +
  geom_interval(aes(xmin = .lower, xmax = .upper)) +
  theme_classic() + 
  tx_plot_theme + 
  labs(x = units_for_plot, y = "", title = "Mean on-farm impact")
ggsave(filename = file.path(outdir, paste("plot_", impact, "_", set_allocation, "-allocation_ON-FARM-TAXA-LEVEL-WEIGHTED.png", sep = "")), width = 11, height = 8.5)

# Same but as CSV output
fit_no_na %>%
  spread_draws(tx_farm_fp_w[tx]) %>%
  median_qi(.width = 0.95) %>%
  left_join(tx_index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  rename(on_farm = tx_farm_fp_w) %>%
  write.csv(file = file.path(outdir, paste("summary_", impact, "_", set_allocation, "-allocation_ON-FARM-TAXA-LEVEL-WEIGHTED.csv", sep = "")), row.names = FALSE)

# Mean total impact taxa-level
fit_no_na %>%
  spread_draws(tx_total_fp_w[tx]) %>%
  median_qi(.width = c(0.95, 0.8, 0.5)) %>%
  left_join(tx_index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  # Filter out tx_total_fp_w for plants and bivalves
  filter(taxa %in% c("bivalves", "plants")==FALSE) %>%
  # Substitute on_farm data as total impact for plants and bivalves
  bind_rows(fit_no_na %>% 
              spread_draws(tx_farm_fp_w[tx]) %>% 
              median_qi(.width = c(0.95, 0.8, 0.5)) %>% 
              left_join(tx_index_key, by = "tx") %>% 
              filter(taxa %in% c("bivalves", "plants")) %>% 
              rename(tx_total_fp_w = tx_farm_fp_w)) %>%
  # REORDER taxa axis
  mutate(full_taxa_name = fct_relevel(full_taxa_name, full_taxa_name_order)) %>%
  #mutate(full_taxa_name = fct_reorder(full_taxa_name, tx_feed_fp_w)) %>%
  ggplot(aes(y = full_taxa_name, x = tx_total_fp_w, xmin = .lower, xmax = .upper)) +
  geom_interval(aes(xmin = .lower, xmax = .upper)) +
  theme_classic() + 
  tx_plot_theme + 
  labs(x = units_for_plot, y = "", title = "Total (on and off-farm) impact", color = "taxa group")
ggsave(filename = file.path(outdir, paste("plot_", impact, "_", set_allocation, "-allocation_TOTAL-IMPACT-TAXA-LEVEL-WEIGHTED.png", sep = "")), width = 11, height = 8.5)

# Same but as CSV output
fit_no_na %>%
  spread_draws(tx_total_fp_w[tx]) %>%
  median_qi(.width = 0.95) %>%
  left_join(tx_index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  # Set .lower limit of plants and bivalves to be 0
  mutate(.lower = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = .lower)) %>% 
  rename(total_stressor = tx_total_fp_w) %>%
  write.csv(file = file.path(outdir, paste("summary_", impact, "_", set_allocation, "-allocation_TOTAL-IMPACT-TAXA-LEVEL-WEIGHTED.csv", sep = "")), row.names = FALSE)

###########################################################
## PLOT WEIGHTED SCI-LEVEL ESTIMATES ANYWAY TO COMPARE WITH UNWEIGHTED ESTIMATES
# Mean 
fit_no_na %>%
  spread_draws(sci_feed_fp_w[sci]) %>%
  median_qi(.width = 0.95) %>%
  left_join(sci_index_key, by = "sci") %>% # Join with index key to get sci and taxa names
  mutate(clean_sci_name = fct_reorder(clean_sci_name, as.character(full_taxa_name))) %>%
  ggplot(aes(y = clean_sci_name, x = sci_feed_fp_w, xmin = .lower, xmax = .upper, color = full_taxa_name)) +
  geom_pointinterval() +
  theme_classic() + 
  sci_plot_theme + 
  labs(x = units_for_plot, y = "", title = "Weighted sci-level off-farm impacts", color = "taxa group")
ggsave(filename = file.path(outdir, paste("plot_", impact, "_", set_allocation, "-allocation_OFF-FARM-SCI-LEVEL-WEIGHTED.png", sep = "")), width = 11, height = 8.5)


fit_no_na %>%
  spread_draws(sci_mu_farm_w[sci]) %>%
  median_qi(.width = 0.95) %>%
  left_join(sci_index_key, by = "sci") %>% # Join with index key to get sci and taxa names
  mutate(clean_sci_name = fct_reorder(clean_sci_name, as.character(full_taxa_name))) %>%
  ggplot(aes(y = clean_sci_name, x = sci_mu_farm_w, xmin = .lower, xmax = .upper, color = full_taxa_name)) +
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
  spread_draws(sci_feed_fp[sci]) %>%
  median_qi(.width = 0.95) %>%
  left_join(sci_index_key, by = "sci") %>% # Join with index key to get sci and taxa names
  # SET plants and bivalves to 0
  mutate(sci_feed_fp = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = sci_feed_fp),
         .lower = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = .lower),
         .upper = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = .upper)) %>%
  mutate(clean_sci_name = fct_reorder(clean_sci_name, as.character(full_taxa_name))) %>%
  ggplot(aes(y = clean_sci_name, x = sci_feed_fp, xmin = .lower, xmax = .upper, color = full_taxa_name)) +
  geom_pointinterval() +
  theme_classic() + 
  sci_plot_theme + 
  labs(x = units_for_plot, y = "", title = "Mean off-farm (feed) impact", color = "taxa group")
ggsave(filename = file.path(outdir, paste("plot_", impact, "_", set_allocation, "-allocation_OFF-FARM-SCI-LEVEL-UNWEIGHTED.png", sep = "")), width = 11, height = 8.5)

# Taxa-level
fit_no_na %>%
  spread_draws(tx_feed_fp[tx]) %>%
  median_qi(.width = c(0.95, 0.8, 0.5)) %>%
  left_join(tx_index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  # SET plants and bivalves to 0
  mutate(tx_feed_fp = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = tx_feed_fp),
         .lower = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = .lower),
         .upper = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = .upper)) %>%
  # REORDER taxa axis
  mutate(full_taxa_name = fct_relevel(full_taxa_name, full_taxa_name_order)) %>%
  #mutate(full_taxa_name = fct_reorder(full_taxa_name, tx_feed_fp_w)) %>%
  ggplot(aes(y = full_taxa_name, x = tx_feed_fp)) +
  geom_interval(aes(xmin = .lower, xmax = .upper)) +
  geom_point(x = c(0), y = c("bivalves")) +
  geom_point(x = 0, y = "plants") +
  theme_classic() + 
  tx_plot_theme + 
  labs(x = units_for_plot, y = "", title = "Mean off-farm (feed) impact")
ggsave(filename = file.path(outdir, paste("plot_", impact, "_", set_allocation, "-allocation_OFF-FARM-TAXA-LEVEL-UNWEIGHTED.png", sep = "")), width = 11, height = 8.5)

###########################################################
# Mean on-farm impact sci-level
# NOTE: Mean on-farm estimates are bounded by zero as declared in STAN model (FCR was left unbounded to allow for mean = 0, so off-farm is also unabounded)
fit_no_na %>%
  spread_draws(sci_mu_farm[sci]) %>%
  median_qi(.width = 0.95) %>%
  left_join(sci_index_key, by = "sci") %>% # Join with index key to get sci and taxa names
  mutate(clean_sci_name = fct_reorder(clean_sci_name, as.character(full_taxa_name))) %>%
  ggplot(aes(y = clean_sci_name, x = sci_mu_farm, xmin = .lower, xmax = .upper, color = full_taxa_name)) +
  geom_pointinterval() +
  theme_classic() + 
  sci_plot_theme + 
  labs(x = units_for_plot, y = "", title = "Mean on-farm impact", color = "taxa group")
ggsave(filename = file.path(outdir, paste("plot_", impact, "_", set_allocation, "-allocation_ON-FARM-SCI-LEVEL-UNWEIGHTED.png", sep = "")), width = 11, height = 8.5)

# Taxa-level
fit_no_na %>%
  spread_draws(tx_mu_farm[tx]) %>%
  median_qi(.width = c(0.95, 0.8, 0.5)) %>%
  left_join(tx_index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  # REORDER taxa axis
  mutate(full_taxa_name = fct_relevel(full_taxa_name, full_taxa_name_order)) %>%
  #mutate(full_taxa_name = fct_reorder(full_taxa_name, tx_feed_fp_w)) %>%
  ggplot(aes(y = full_taxa_name, x = tx_mu_farm)) +
  geom_interval(aes(xmin = .lower, xmax = .upper)) +
  theme_classic() + 
  tx_plot_theme + 
  labs(x = units_for_plot, y = "", title = "Mean on-farm impact", color = "taxa group")
ggsave(filename = file.path(outdir, paste("plot_", impact, "_", set_allocation, "-allocation_ON-FARM-TAXA-LEVEL-UNWEIGHTED.png", sep = "")), width = 11, height = 8.5)

###########################################################
# Mean unweighted total (on + off-farm) impact sci-level
fit_no_na %>%
  spread_draws(sci_total_fp[sci]) %>%
  median_qi(.width = 0.95) %>%
  left_join(sci_index_key, by = "sci") %>% # Join with index key to get sci and taxa names
  mutate(clean_sci_name = fct_reorder(clean_sci_name, as.character(full_taxa_name))) %>%
  ggplot(aes(y = clean_sci_name, x = sci_total_fp, xmin = .lower, xmax = .upper, color = full_taxa_name)) +
  geom_pointinterval() +
  theme_classic() + 
  sci_plot_theme + 
  labs(x = units_for_plot, y = "", title = "Total (on and off-farm) impact", color = "taxa group") 
ggsave(filename = file.path(outdir, paste("plot_", impact, "_", set_allocation, "-allocation_TOTAL-IMPACT-SCI-LEVEL-UNWEIGHTED.png", sep = "")), width = 11, height = 8.5)

# Taxa-level
fit_no_na %>%
  spread_draws(tx_total_fp[tx]) %>%
  median_qi(.width = c(0.95, 0.8, 0.5)) %>%
  left_join(tx_index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  # REORDER taxa axis
  mutate(full_taxa_name = fct_relevel(full_taxa_name, full_taxa_name_order)) %>%
  #mutate(full_taxa_name = fct_reorder(full_taxa_name, tx_feed_fp_w)) %>%
  ggplot(aes(y = full_taxa_name, x = tx_total_fp)) +
  geom_interval(aes(xmin = .lower, xmax = .upper)) +
  theme_classic() + 
  tx_plot_theme + 
  labs(x = units_for_plot, y = "", title = "Total (on and off-farm) impact", color = "taxa group")
ggsave(filename = file.path(outdir, paste("plot_", impact, "_", set_allocation, "-allocation_TOTAL-IMPACT-TAXA-LEVEL-UNWEIGHTED.png", sep = "")), width = 11, height = 8.5)

rm(fit_no_na) # Clear fit-no-na before restarting loop

#} # End loop by allocation method: for (i in c("Mass", "Gross energy content", "Economic")){

# BEFORE CLEARING WORKSPACE, run 04_plot_fcr_and_feed at least once (outputs will be the same for all models)


