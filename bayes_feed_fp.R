# Combine bayes_prop_feed and bayes_fcr to calculate footprint

rm(list=ls())
library(tidyverse)
library(rstan)
library(taxize)
library(data.table)
library(countrycode) # part of clean.lca
library(bayesplot) # for mcmc_areas_ridges
library(shinystan)

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
# Model 2: Remove all NAs - estimate feed footprint for all sci names

# Remove NAs
lca_dat_no_na <- lca_dat_no_zeroes %>%
  filter(is.na(Feed_soy_percent)==FALSE) %>%
  filter(is.na(FCR) == FALSE) %>%
  filter(FCR != 0) %>%
  mutate(clean_sci_name = as.factor(clean_sci_name),
         sci = as.numeric(clean_sci_name),
         taxa_group_name = as.factor(taxa_group_name),
         tx = as.numeric(taxa_group_name)) %>%
  select(clean_sci_name, sci, taxa_group_name, tx, FCR, Feed_soy_percent, Feed_othercrops_percent, Feed_FMFO_percent, Feed_animal_percent) %>%
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
         feed_animal_new = feed_animal_new / sum) %>%
  select(clean_sci_name, sci, taxa_group_name, tx, FCR, contains("new"))

# Try for just two scientific names:
lca_dat_no_na <- lca_dat_no_na %>%
  filter(clean_sci_name %in% c("Oncorhynchus mykiss", "Salmo salar"))

# BOX PLOTS OF DATA:
# Theme for ALL PLOTS (including mcmc plots)
plot_theme <- theme(axis.text=element_text(size=14, color = "black"))

# FCR:
plot_fcr <- lca_dat_no_na %>%
  select(clean_sci_name, FCR) %>%
  #mutate(clean_sci_name = paste(clean_sci_name, row_number(), sep = "")) %>%
  pivot_longer(cols = FCR)
ggplot(data = plot_fcr, aes(x = clean_sci_name, y = value)) +
  geom_boxplot() +
  theme_classic() +
  plot_theme +
  labs(title = "Boxplots of FCRs",
       x = "",
       y = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#ggsave(file.path(outdir, "boxplot_fcr_no_na.png"), height = 8, width = 11.5)


# feed proportion:
plot_feed_prop <- lca_dat_no_na %>%
  select(clean_sci_name, soy = feed_soy_new, crops = feed_crops_new, fmfo = feed_fmfo_new, animal = feed_animal_new) %>%
  pivot_longer(cols = soy:animal)



feed_vars <- c("soy", "crops", "fmfo", "animal")
for (i in 1:length(feed_vars)) {
  p <- ggplot(data = plot_feed_prop %>% filter(name == feed_vars[i]), aes(x = clean_sci_name, y = value)) +
    geom_boxplot() +
    theme_classic() +
    plot_theme +
    labs(title = paste("Boxplots of ", feed_vars[i], " feed proportions", sep = ""),
         x = "",
         y = "")  +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p)
  #ggsave(file.path(outdir, "boxplot_feed-prop_no_na.png"), height = 8, width = 11.5)
}

# Set data for model

# overall model:
N = nrow(lca_dat_no_na)
N_SCI <- length(unique(lca_dat_no_na$sci))
sci <- lca_dat_no_na$sci
N_TX <- length(unique(lca_dat_no_na$tx))
tx <- lca_dat_no_na$tx

# for FCR model:
x <- lca_dat_no_na$FCR

# for Feed proportion model:
K = 4
feed_weights <- lca_dat_no_na %>%
  select(contains("new")) %>%
  as.matrix()
  
# Get counts per sci name and counts per taxa group (also included as data in the model):
sci_kappa <- lca_dat_no_na %>% 
  select(contains(c("new", "sci", "obs"))) %>%
  group_by(sci) %>% 
  summarise(n_obs = n()) %>%
  ungroup() %>%
  arrange(sci) %>%
  pull(n_obs)

tx_kappa <- lca_dat_no_na %>% 
  select(contains(c("new", "tx", "obs"))) %>%
  group_by(tx) %>% 
  summarise(n_obs = n()) %>%
  ungroup() %>%
  arrange(tx) %>%
  pull(n_obs)

# data (constants) for final foot print calculation 
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

fp_nitrogen_dat <- fp_dat %>%
  filter(FP == "Nitrogen") %>%
  select(-FP) %>%
  pivot_wider(names_from = Category, values_from = FP_val) %>%
  as.matrix() %>%
  c()

fp_phosphorus_dat <- fp_dat %>%
  filter(FP == "Phosphorus") %>%
  select(-FP) %>%
  pivot_wider(names_from = Category, values_from = FP_val) %>%
  as.matrix() %>%
  c()

fp_land_dat <- fp_dat %>%
  filter(FP == "Land") %>%
  select(-FP) %>%
  pivot_wider(names_from = Category, values_from = FP_val) %>%
  as.matrix() %>%
  c()

fp_water_dat <- fp_dat %>%
  filter(FP == "Water") %>%
  select(-FP) %>%
  pivot_wider(names_from = Category, values_from = FP_val) %>%
  as.matrix() %>%
  c()

# Prior information
# Mean feed proportions per sci-name
sci_phi_mean <- lca_dat_no_na %>% 
  select(contains(c("new", "sci", "clean_sci_name"))) %>%
  group_by(clean_sci_name, sci) %>% 
  summarise(across(contains("new"), mean)) %>%
  ungroup() %>%
  arrange(sci)

# Mean feed proportions per taxa group
tx_phi_mean <- lca_dat_no_na %>% 
  select(contains(c("new", "tx", "taxa_group_name"))) %>%
  group_by(taxa_group_name, tx) %>% 
  summarise(across(contains("new"), mean)) %>%
  ungroup() %>%
  arrange(tx)


stan_data <- list(N = N,
                  N_SCI = N_SCI, 
                  sci = sci,
                  N_TX = N_TX,
                  tx = tx,
                  x = x,
                  K = K,
                  feed_weights = feed_weights,
                  sci_kappa = sci_kappa,
                  tx_kappa = tx_kappa,
                  fp_carbon_dat = fp_carbon_dat,
                  fp_nitrogen_dat = fp_nitrogen_dat,
                  fp_phosphorus_dat = fp_phosphorus_dat,
                  fp_land_dat = fp_land_dat,
                  fp_water_dat = fp_water_dat)

# Estimate foot print for all scientific names without NAs
stan_no_na <- 'data {
  // data for gamma model for FCR
  int<lower=0> N;  // number of observations
  vector<lower=0>[N] x; // data
  int N_TX; // number of taxa groups
  int tx[N]; // taxa group index (ordered by unique sci index)
  int N_SCI; // number of scientific names
  int sci[N]; // sciname index
  
  // data for dirichlet model for feed
  int K; // number of feed types
  simplex[K] feed_weights[N]; // array of observed feed weights simplexes
  int sci_kappa[N_SCI]; // number of observations per sci-name
  int tx_kappa[N_TX]; // number of observations per taxa group
  
  // constants for generated quantities
  vector[K] fp_carbon_dat;
  vector[K] fp_nitrogen_dat;
  vector[K] fp_phosphorus_dat;
  vector[K] fp_land_dat;
  vector[K] fp_water_dat;
}
parameters {
  // FCR model:
  real<lower=0> mu;
  real<lower=0> sigma;
  vector<lower=0>[N_TX] tx_mu;
  real<lower=0> tx_sigma;
  vector<lower=0>[N_SCI] sci_mu;
  real<lower=0> sci_sigma;
  
  // Feed proportion model:
  simplex[K] sci_phi[N_SCI];
  simplex[K] sci_theta[N_SCI]; // vectors of estimated sci-level feed weight simplexes
  simplex[K] tx_phi[N_TX];
  simplex[K] tx_theta[N_TX];
  simplex[K] phi;
  simplex[K] theta;
  
  // Params for the dirichlet priors:
  // real<lower=0> sigma_1;
  // real<lower=0> sigma_2;
}
transformed parameters {
  // define transofrmed params for gamma model for FCRs
  real shape;
  real rate; 
  vector[N_SCI] sci_shape;
  vector[N_SCI] sci_rate;
  vector[N_TX] tx_shape;
  vector[N_TX] tx_rate;
  
  // define params for dirichlet model for feed proportions
  vector<lower=0>[K] sci_alpha[N_SCI];
  vector<lower=0>[K] tx_alpha[N_TX];
  vector<lower=0>[K] alpha;
  
  // reparamaterize gamma to get mu and sigma; defining these here instead of the model section allows us to see these parameters in the output
  // global-level
  shape = square(mu) / square(sigma);
  rate = mu / square(sigma);
  // sci name and taxa group levels
  for (n in 1:N){
    tx_shape[tx[n]] = square(tx_mu[tx[n]]) ./ square(tx_sigma);
    tx_rate[tx[n]] = tx_mu[tx[n]] ./ square(tx_sigma);
    sci_shape[sci[n]] = square(sci_mu[sci[n]]) ./ square(sci_sigma);
    sci_rate[sci[n]] = sci_mu[sci[n]] ./ square(sci_sigma);
  }
  
  // reparameterize alphas as a vector of means (phi) and counts (kappas)
  // phi is expected value of theta (mean feed weights)
  // kappa is strength of the prior measured in number of prior observations (minus K)
  alpha = N * phi;
  for (n_tx in 1:N_TX) {
    tx_alpha[n_tx] = tx_kappa[n_tx] * tx_phi[n_tx];
  }    
  
  for (n_sci in 1:N_SCI) {
    sci_alpha[n_sci] = sci_kappa[n_sci] * sci_phi[n_sci];
  }
}
model {
  // define priors for gamma model for FCRs
  // Put priors on mu and sigma (instead of shape and rate) since this is more intuitive:
  mu ~ uniform(0, 10); // note: uniform(0,100) for all of these doesnt help much with convergence
  sigma ~ uniform(0, 10);
  tx_mu ~ uniform(0, 10);
  tx_sigma ~ uniform(0, 10);
  sci_mu ~ uniform(0, 10);
  sci_sigma ~ uniform(0, 10);
  
  // define priors for dirichlet model for feed proportions
  // sci_phi defined as sci_phi[sci][K]
  
  // option 1: define feed proportion priors as lower upper bounds
  // sci_phi[2][1] ~ uniform(0.1, 0.2); // hypothetical lower and upper bounds for Oncorhynchus mykiss soy 
  // sci_phi[6][1] ~ uniform(0.05, 0.2); // lower upper bounds fo Salmo salar soy
  
  // option 2: define feed proportions as means (need to define sigmas in parameters block: real<lower=0> sigma_1, sigma_2 etc;)
  // sci_phi[2][1] ~ normal(0.13, sigma_1); // mean for Oncorhynhchus mykiss soy feed 
  // sci_phi[6][1] ~ normal(0.7, sigma_2); // mean for Salmo salar soy feed
  // tx_phi[6][1] ~ normal(0.13, sigma_1); // mean for taxa-level trout soy feed 
  // tx_phi[4][1] ~ normal(0.7, sigma_2); // mean for taxa-level salmon/char soy feed
  
  // likelihood
  
  // global-level for dirichlet
  tx_mu ~ gamma(shape, rate);
  
  for (n in 1:N){
    // gamma model sci-name and taxa-level (global level is part of transformed params)
    sci_mu[sci[n]] ~ gamma(tx_shape[tx[n]], tx_rate[tx[n]]);
    x[n] ~ gamma(sci_shape[sci[n]], sci_rate[sci[n]]);
    
    // dirichlet model sci-name and taxa-level
    tx_phi[tx[n]] ~ dirichlet(alpha);
    sci_phi[sci[n]] ~ dirichlet(to_vector(tx_alpha[tx[n]]));
    feed_weights[n] ~ dirichlet(to_vector(sci_alpha[sci[n]])); 
  }

  // dirichlet model - estimate feed weights based on estimated alphas
  // global level estimates
  theta ~ dirichlet(to_vector(alpha));
  // taxa level estimates
  for (n_tx in 1:N_TX) {
    tx_theta[n_tx] ~ dirichlet(to_vector(tx_alpha[n_tx]));
  }
  // sci-name level estimates
  for (n_sci in 1:N_SCI) {
    sci_theta[n_sci] ~ dirichlet(to_vector(sci_alpha[n_sci]));
  }
}'

# Add generated quantities later
# generated quantities {
#   // Carbon
#   real<lower=0> species_carbon_footprint;
#   vector[K] feed_carbon_footprint;
#   real total_feed_carbon_footprint;
#   // Nitrogen
#   real<lower=0> species_nitrogen_footprint;
#   vector[K] feed_nitrogen_footprint;
#   real total_feed_nitrogen_footprint;
#   // Phosphorus
#   real<lower=0> species_phosphorus_footprint;
#   vector[K] feed_phosphorus_footprint;
#   real total_feed_phosphorus_footprint;
#   // Land
#   real<lower=0> species_land_footprint;
#   vector[K] feed_land_footprint;
#   real total_feed_land_footprint;
#   // Water
#   real<lower=0> species_water_footprint;
#   vector[K] feed_water_footprint;
#   real total_feed_water_footprint;
#   
#   // Calculations
#   feed_carbon_footprint = fp_carbon_dat .* theta;
#   total_feed_carbon_footprint = sum(feed_carbon_footprint);
#   species_carbon_footprint = mu * total_feed_carbon_footprint;
# 
#   feed_nitrogen_footprint = fp_nitrogen_dat .* theta;
#   total_feed_nitrogen_footprint = sum(feed_nitrogen_footprint);
#   species_nitrogen_footprint = mu * total_feed_nitrogen_footprint;
# 
#   feed_phosphorus_footprint = fp_phosphorus_dat .* theta;
#   total_feed_phosphorus_footprint = sum(feed_phosphorus_footprint);
#   species_phosphorus_footprint = mu * total_feed_phosphorus_footprint;
#   
#   feed_land_footprint = fp_land_dat .* theta;
#   total_feed_land_footprint = sum(feed_land_footprint);
#   species_land_footprint = mu * total_feed_land_footprint;
#   
#   feed_water_footprint = fp_water_dat .* theta;
#   total_feed_water_footprint = sum(feed_water_footprint);
#   species_water_footprint = mu * total_feed_water_footprint;
# }'

no_na_mod <- stan_model(model_code = stan_no_na, verbose = TRUE)
# Note: For Windows, apparently OK to ignore this warning message:
# Warning message:
#   In system(paste(CXX, ARGS), ignore.stdout = TRUE, ignore.stderr = TRUE) :
#   'C:/rtools40/usr/mingw_/bin/g++' not found

# Fit model:
# Set seed while testing
fit_no_na <- sampling(object = no_na_mod, data = stan_data, cores = 4, seed = "11729")
                       #,iter = 10000, cores = 4,
                       #control = list(adapt_delta = 0.99))
print(fit_no_na)

launch_shinystan(fit_no_na)

######################################################################################################
# Model 1: Remove all NAs - estimate feed footprint for Oncorhynchus mykiss

# Remove NAs
lca_dat_no_na <- lca_dat_no_zeroes %>%
  filter(clean_sci_name == "Oncorhynchus mykiss") %>%
  filter(is.na(Feed_soy_percent)==FALSE) 


# BOX PLOTS OF DATA:

# Theme for ALL PLOTS (including mcmc plots)
plot_theme <- theme(axis.text=element_text(size=14, color = "black"))

# FCR:
plot_fcr <- lca_dat_no_na %>%
  select(clean_sci_name, FCR) %>%
  mutate(clean_sci_name = paste(clean_sci_name, row_number(), sep = "")) %>%
  pivot_longer(cols = FCR)
ggplot(data = plot_fcr, aes(x = name, y = value)) +
  geom_boxplot() +
  theme_classic() +
  plot_theme +
  labs(title = "Boxplots of FCRs for Oncorhynchus mykiss",
       x = "",
       y = "")
ggsave(file.path(outdir, "boxplot_fcr_trout.png"), height = 8, width = 11.5)


# feed proportion:
plot_feed_prop <- lca_dat_no_na %>%
  select(clean_sci_name, soy = feed_soy_new, crops = feed_crops_new, fmfo = feed_fmfo_new, animal = feed_animal_new) %>%
  mutate(clean_sci_name = paste(clean_sci_name, row_number(), sep = "")) %>%
  pivot_longer(cols = soy:animal)

ggplot(data = plot_feed_prop, aes(x = name, y = value)) +
  geom_boxplot() +
  theme_classic() +
  plot_theme +
  labs(title = "Boxplots of feed proportions for Oncorhynchus mykiss",
       x = "",
       y = "")
ggsave(file.path(outdir, "boxplot_feed-prop_trout.png"), height = 8, width = 11.5)

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

fp_nitrogen_dat <- fp_dat %>%
  filter(FP == "Nitrogen") %>%
  select(-FP) %>%
  pivot_wider(names_from = Category, values_from = FP_val) %>%
  as.matrix() %>%
  c()

fp_phosphorus_dat <- fp_dat %>%
  filter(FP == "Phosphorus") %>%
  select(-FP) %>%
  pivot_wider(names_from = Category, values_from = FP_val) %>%
  as.matrix() %>%
  c()

fp_land_dat <- fp_dat %>%
  filter(FP == "Land") %>%
  select(-FP) %>%
  pivot_wider(names_from = Category, values_from = FP_val) %>%
  as.matrix() %>%
  c()

fp_water_dat <- fp_dat %>%
  filter(FP == "Water") %>%
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
  vector[k] fp_nitrogen_dat;
  vector[k] fp_phosphorus_dat;
  vector[k] fp_land_dat;
  vector[k] fp_water_dat;
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
  // Carbon
  real<lower=0> species_carbon_footprint;
  vector[k] feed_carbon_footprint;
  real total_feed_carbon_footprint;
  // Nitrogen
  real<lower=0> species_nitrogen_footprint;
  vector[k] feed_nitrogen_footprint;
  real total_feed_nitrogen_footprint;
  // Phosphorus
  real<lower=0> species_phosphorus_footprint;
  vector[k] feed_phosphorus_footprint;
  real total_feed_phosphorus_footprint;
  // Land
  real<lower=0> species_land_footprint;
  vector[k] feed_land_footprint;
  real total_feed_land_footprint;
  // Water
  real<lower=0> species_water_footprint;
  vector[k] feed_water_footprint;
  real total_feed_water_footprint;
  
  // Calculations
  feed_carbon_footprint = fp_carbon_dat .* theta;
  total_feed_carbon_footprint = sum(feed_carbon_footprint);
  species_carbon_footprint = mu * total_feed_carbon_footprint;

  feed_nitrogen_footprint = fp_nitrogen_dat .* theta;
  total_feed_nitrogen_footprint = sum(feed_nitrogen_footprint);
  species_nitrogen_footprint = mu * total_feed_nitrogen_footprint;

  feed_phosphorus_footprint = fp_phosphorus_dat .* theta;
  total_feed_phosphorus_footprint = sum(feed_phosphorus_footprint);
  species_phosphorus_footprint = mu * total_feed_phosphorus_footprint;
  
  feed_land_footprint = fp_land_dat .* theta;
  total_feed_land_footprint = sum(feed_land_footprint);
  species_land_footprint = mu * total_feed_land_footprint;
  
  feed_water_footprint = fp_water_dat .* theta;
  total_feed_water_footprint = sum(feed_water_footprint);
  species_water_footprint = mu * total_feed_water_footprint;
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
                                                            fp_carbon_dat = fp_carbon_dat,
                                                            fp_nitrogen_dat = fp_nitrogen_dat,
                                                            fp_phosphorus_dat = fp_phosphorus_dat,
                                                            fp_land_dat = fp_land_dat,
                                                            fp_water_dat = fp_water_dat))#,
                       # iter = 10000, cores = 4,
                       # control = list(adapt_delta = 0.99))
print(fit_pooled)

launch_shinystan(fit_pooled)

feeds <- c("soy", "crops", "fmfo", "animal")
feed_key <- data.frame(carbon_footprint = paste("carbon_footprint[", feeds, "]", sep = ""),
                       nitrogen_footprint = paste("nitrogen_footprint[", feeds, "]", sep = ""),
                       phosphorus_footprint = paste("phosphorus_footprint[", feeds, "]", sep = ""),
                       land_footprint = paste("land_footprint[", feeds, "]", sep = ""),
                       water_footprint = paste("water_footprint[", feeds, "]", sep = ""))

fit_pooled_clean <- fit_pooled
names(fit_pooled_clean)[grep(names(fit_pooled_clean), pattern = "feed_carbon_footprint\\[")] <- feed_key$carbon_footprint
names(fit_pooled_clean)[grep(names(fit_pooled_clean), pattern = "feed_nitrogen_footprint\\[")] <- feed_key$nitrogen_footprint
names(fit_pooled_clean)[grep(names(fit_pooled_clean), pattern = "feed_phosphorus_footprint\\[")] <- feed_key$phosphorus_footprint
names(fit_pooled_clean)[grep(names(fit_pooled_clean), pattern = "feed_land_footprint\\[")] <- feed_key$land_footprint
names(fit_pooled_clean)[grep(names(fit_pooled_clean), pattern = "feed_water_footprint\\[")] <- feed_key$water_footprint

distribution_pooled <- as.matrix(fit_pooled_clean)

# FIX IT - replace parameter names and add plots for other parameters
p_footprint <- mcmc_areas(distribution_pooled,
                          pars = vars(contains("carbon_footprint")),
                          prob = 0.8,
                          prob_outer = 0.9,
                          area_method = "scaled height",
                          point_est = "median") + 
  ggtitle("Oncorhynchus mykiss full feed footprint model", "with 80% credible intervals") +
  xlim(0, 10) +
  plot_theme 

p_footprint
ggsave(filename = file.path(outdir, "bayes-example_trout_carbon-feed-footprint.png"), width = 11, height = 8.5)


p_footprint <- mcmc_areas(distribution_pooled,
                          pars = vars(contains("nitrogen_footprint")),
                          prob = 0.8,
                          prob_outer = 0.9,
                          area_method = "scaled height",
                          point_est = "median") + 
  ggtitle("Oncorhynchus mykiss full feed footprint model", "with 80% credible intervals") +
  xlim(0, 0.01) +
  plot_theme 

p_footprint
ggsave(filename = file.path(outdir, "bayes-example_trout_nitrogen-feed-footprint.png"), width = 11, height = 8.5)

p_footprint <- mcmc_areas(distribution_pooled,
                          pars = vars(contains("phosphorus_footprint")),
                          prob = 0.8,
                          prob_outer = 0.9,
                          area_method = "scaled height",
                          point_est = "median") + 
  ggtitle("Oncorhynchus mykiss full feed footprint model", "with 80% credible intervals") +
  xlim(0, 0.001) +
  plot_theme 

p_footprint
ggsave(filename = file.path(outdir, "bayes-example_trout_phosphorus-feed-footprint.png"), width = 11, height = 8.5)

p_footprint <- mcmc_areas(distribution_pooled,
                          pars = vars(contains("land_footprint")),
                          prob = 0.8,
                          prob_outer = 0.9,
                          area_method = "scaled height",
                          point_est = "median") + 
  ggtitle("Oncorhynchus mykiss full feed footprint model", "with 80% credible intervals") +
  xlim(0, 5) +
  plot_theme 

p_footprint
ggsave(filename = file.path(outdir, "bayes-example_trout_land-feed-footprint.png"), width = 11, height = 8.5)

p_footprint <- mcmc_areas(distribution_pooled,
                          pars = vars(contains("water_footprint")),
                          prob = 0.8,
                          prob_outer = 0.9,
                          area_method = "scaled height",
                          point_est = "median") + 
  ggtitle("Oncorhynchus mykiss full feed footprint model", "with 80% credible intervals") +
  xlim(0, 0.2) +
  plot_theme 

p_footprint
ggsave(filename = file.path(outdir, "bayes-example_trout_water-feed-footprint.png"), width = 11, height = 8.5)
