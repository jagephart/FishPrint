# Get Bayesian sci and taxa-level means of on-farm (non-feed) and off-farm (feed) NITROGEN and PHOSPHORUS impacts 

###################### REMINDER FOR NITROGEN and PHOSPHORUS: 
# Run once with impact <- "Freshwater eutrophication", feed_element <- "P", and fish_element <- "P_t_liveweight_t" - i.e., Phosphorus impact and 
# And once with impact <- "Marine eutrophication", feed_element <- "N", and fish_element <- "N_t_liveweight_t" -  i.e., Nitrogen impact
# OFF-FARM impacts should be ZERO when FCR = 0 (i.e., for full taxa groups bivalves and plants, and other sci-names within other taxa groups)
# ON-FARM carbon footprint calculated as: FCR * Nitrogen (or Phosphorus) content of feed - Nitrogen (or Phosphorus) content of fish (matched by clean_sci_name)

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

# SELECT BETWEEN PHOSPHORUS VS NITROGEN MODEL
# option 1: Phosphorus model
# fish_element <- "P_t_liveweight_t"
# impact <- "Freshwater eutrophication"
# feed_element <- "P"
# interval_palette <- c("#FFB4C4", "#FF86A4", "#FC5185") # Order: light to dark

# option 2: Nitrogen model:
fish_element <- "N_t_liveweight_t"
impact <- "Marine eutrophication"
feed_element <- "N"
interval_palette <- c("#FFD5A9", "#FFBD79", "#FFA647") # COMPLEMENTARY COLOR - use orange instead of yellow

# Load Data
#lca_full_dat <- read.csv(file.path(datadir, "2021-05-05_lca-dat-imputed-vars_rep-sqrt-n-farms_live-weight.csv"), fileEncoding="UTF-8-BOM")
lca_full_dat <- read.csv(file.path(datadir, "2021-05-05_lca-dat-imputed-vars_rep-sqrt-n-farms_edible-weight.csv"), fileEncoding="UTF-8-BOM")

# fish_content_dat <- read.csv(file.path(datadir, "fish_NP_clean.csv")) %>%
#   select(clean_sci_name, !!sym(fish_element))

lca_model_dat <- lca_full_dat %>%
  #left_join(fish_content_dat, by = "clean_sci_name") %>%
  select(study_id, iso3c, clean_sci_name, taxa, intensity, system, 
         feed_soy, feed_crops, feed_fmfo, feed_animal, 
         fcr,
         !!sym(fish_element)) %>%
  drop_na() %>%  # just in case
  # OPTION 1: TEST MODEL ON SPECIES THAT ARE FED 
  #filter(fcr != 0)  %>% 
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
  # LAST FORMATING STEP - always arrange by clean_sci_name
  arrange(clean_sci_name) %>%
  mutate(clean_sci_name = as.factor(clean_sci_name),
         sci = as.numeric(clean_sci_name),
         taxa = as.factor(taxa),
         tx = as.numeric(taxa)) 

# CHOOSE FEED IMPACT CONSTANT:
set_allocation <- "Mass"
#set_allocation <- "Gross energy content"
#set_allocation <- "Economic"

# IMPORTANT - multiply all values by 1000 to convert to kg CO2 per tonne (currently in kg CO2 per kg)
  fp_dat <- read.csv(file.path(datadir, "weighted_feed_fp.csv")) %>%
    filter(Allocation == set_allocation) %>%
    mutate(ave_stressor_per_tonne = ave_stressor * 1000)
  
  # ORDER OF feed_weights data vector: soy, crops, fmfo, animal
  #head(feed_weights)
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


#####################################################################
# Set data, indices, constants, weights for STAN

# VARIABLE-SPECIFIC DATA:

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

# FISH CONTENT CONSTANTS
fish_content <- lca_model_dat %>%
  select(clean_sci_name, !!sym(fish_element)) %>%
  arrange(clean_sci_name) %>%
  unique() %>%
  pull(!!sym(fish_element))

# FEED CONTENT CONSTANTS (for on-farm model):
source("Functions.R") 
feed_NP <-clean_feedNutrition(feedNutrition_data = read.csv(file.path(datadir, "United-States-Canadian-Tables-of-Feed-1982-pages-68-921-with_CrudeProtein.csv"),
                                                            stringsAsFactors = FALSE)) %>%
  mutate(feed_type = case_when(
    (ingredient == "Soy") ~ "soy",
    (ingredient == "Crop") ~ "crops",
    (ingredient == "Fishery") ~ "fmfo",
    (ingredient == "Animal by-products") ~ "animal"
  )) %>%
  select(-c("ingredient", "sd")) %>%
  pivot_wider(names_from = "element", values_from = "value")

# Choose feed element (N or P)
set_feed_content_order <- c("soy", "crops", "fmfo", "animal")

# Divide by 100 to get the correct units (N/P feed content data are in percentages)
feed_content <- feed_NP %>%
  arrange(match(feed_type, set_feed_content_order)) %>%
  pull(!!sym(feed_element)) / 100

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
# REMINDER RE: PRIORS vs NO PRIORS - make sure STAN code below for defining/applying priors is allowed to run or commented out as needed
# NO PRIORS
# stan_data <- list(N = N,
#                   N_SCI = N_SCI,
#                   n_to_sci = n_to_sci,
#                   N_TX = N_TX,
#                   sci_to_tx = sci_to_tx,
#                   fcr = fcr,
#                   K = K,
#                   feed_weights = feed_weights,
#                   sci_kappa = sci_kappa,
#                   tx_kappa = tx_kappa,
#                   fp_constant = fp_constant,
#                   fish_content = fish_content,
#                   feed_content = feed_content,
#                   sci_w = sci_w,
#                   where_tx = where_tx,
#                   n_sci_in_tx = n_sci_in_tx,
#                   slice_where_tx = slice_where_tx)

#WITH PRIORS
# REMINDER RE: PRIORS vs NO PRIORS - make sure STAN code below for defining/applying priors is allowed to run or commented out as needed
stan_data <- list(N = N,
                  N_SCI = N_SCI,
                  n_to_sci = n_to_sci,
                  N_TX = N_TX,
                  sci_to_tx = sci_to_tx,
                  fcr = fcr,
                  K = K,
                  feed_weights = feed_weights,
                  sci_kappa = sci_kappa,
                  tx_kappa = tx_kappa,
                  fp_constant = fp_constant,
                  fish_content = fish_content,
                  feed_content = feed_content,
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
  vector[11] priors;
  int prior_vec_index[11];

  // constants to apply to feed footrpint
  vector[K] fp_constant;
  
  // constants for on-farm footrpint
  vector<lower=0>[N_SCI] fish_content;
  vector[K] feed_content;
  
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
  tx_mu_fcr[prior_vec_index] ~ normal(priors, 1);

  // weak priors on sigma
  tx_sigma_fcr ~ cauchy(0, 1);
  sci_sigma_fcr ~ cauchy(0, 1);
  
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
  
}
generated quantities {
  // Declare variables for weightings
  vector[N_TX] tx_feed_fp;
  vector[N_SCI] sci_feed_fp;
  vector[N_SCI] sci_total_fp; // unweighted
  //vector[N_TX] tx_total_fp; // not applicable for this model
  vector[N_TX] tx_total_fp_w; // weighted
  vector[N_TX] tx_feed_fp_w;
  vector[N_TX] tx_farm_fp_w;
  vector[N_SCI] sci_feed_fp_w;
  vector[N_SCI] sci_mu_farm_w;
    
  // Declare variables for on farm (unweighted) model
  vector[N_SCI] sci_mu_farm;

  // Feed footrpint calculations
  for (n_tx in 1:N_TX) {
    tx_feed_fp[n_tx] = tx_mu_fcr[n_tx] * sum(fp_constant .* tx_theta[n_tx]);
  }
  for (n_sci in 1:N_SCI) {
    sci_feed_fp[n_sci] = sci_mu_fcr[n_sci] * sum(fp_constant .* sci_theta[n_sci]);
  }
  
  // Multiply by 1000 to get to correct units (kg / tonne)
  // On farm footprint calculations - code only calculates unweighted on-farm for sci-level (no fish_content data at taxa level)
  for (n_sci in 1:N_SCI) {
    sci_mu_farm[n_sci] = 1000 * (sci_mu_fcr[n_sci] * sum(feed_content .* sci_theta[n_sci]) - fish_content[n_sci]);
  }

  // Sum off farm (feed) and on farm footprints (UNWEIGHTED sci-level only, no unweighted estimate for taxa level)
  for (n_sci in 1:N_SCI) {
    sci_total_fp[n_sci] = sci_feed_fp[n_sci] + sci_mu_farm[n_sci];
  }

  // Apply weightings
  
  // Apply individually to sci_feed_fp and sci_mu_farm and sum to get sci_total_fp
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
                      iter = 2500, 
                      control = list(adapt_delta = 0.99, max_treedepth = 15))
#fit_no_na <- sampling(object = no_na_mod, data = stan_data, cores = 4, iter = 5000, control = list(adapt_delta = 0.99))
end_sampling <- Sys.time()
end_sampling - start_sampling
summary(fit_no_na)$summary

#launch_shinystan(fit_no_na)

###########################################################
# STEP 3: OUTPUT RESULTS
###########################################################
# RESTARTING POINT
rm(list=ls()[!(ls() %in% c("datadir", "outdir", "impact", "set_allocation",
                           "lca_model_dat", "fit_no_na", "interval_palette"))])
save.image(file = file.path(outdir, paste(Sys.Date(), "_full-model-posterior_", impact, "_", set_allocation, "-allocation.RData", sep = "")))

###########################################################

# SET THEME
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
# PLOT final outputs (off-farm, on-farm, and total impacts)

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
  geom_point(x = 0, y = "bivalves") +
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
  #mutate(full_taxa_name = fct_reorder(full_taxa_name, tx_feed_fp_w)) %>%
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
  # Set .lower limit of plants and bivalves to be 0
  #mutate(.lower = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = .lower)) %>% 
  # REORDER taxa axis
  mutate(full_taxa_name = fct_relevel(full_taxa_name, full_taxa_name_order)) %>%
  #mutate(full_taxa_name = fct_reorder(full_taxa_name, tx_feed_fp_w)) %>%
  ggplot(aes(y = full_taxa_name, x = tx_total_fp_w)) +
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
  geom_point(x = c(0), y = c("bivalves")) +
  geom_point(x = 0, y = "plants") +
  geom_interval(aes(xmin = .lower, xmax = .upper)) +
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

# NA: TAXA-LEVEL UNWEIGHTED ON-FARM IMPACTS NOT CALCULATED BY THIS MODEL (no taxa-level fish content info to calculate this in the "generated quantities" section)

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

# NA: TAXA-LEVEL UNWEIGHTED ON-FARM IMPACTS NOT CALCULATED BY THIS MODEL (no taxa-level fish content info to calculate this in the "generated quantities" section)

######################################################################################################
# NEXT: Before clearing workspace, use 04_plot_common_outputs.R - plot other intermediate-level calculations (these are universally shared among all models)



