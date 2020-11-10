# Bayesian estimation of proportions of each feed component (soy, other crops, FMFOs, and other animal)

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

lca_dat <- read.csv(file.path(datadir, "LCA_compiled_20201109.csv"), fileEncoding="UTF-8-BOM") #fileEncoding needed when reading in file from windows computer (suppresses BOM hidden characters)
source("Functions.R")

lca_dat_no_zeroes <- clean.lca(LCA_data = lca_dat) %>%
  select(clean_sci_name, Feed_soy_percent, Feed_othercrops_percent, Feed_FMFO_percent, Feed_animal_percent, taxa_group_name) %>%
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
  select(clean_sci_name, taxa_group_name, contains("new"))

# Set the FINAL value to be no less than 0.01
lca_dat_no_zeroes <- clean.lca(LCA_data = lca_dat) %>%
  select(clean_sci_name, Feed_soy_percent, Feed_othercrops_percent, Feed_FMFO_percent, Feed_animal_percent, taxa_group_name) %>%
  # NOTE multinomial-dirchlet model requires all elements > 0 (change some to 0.001 for now?)
  mutate(feed_soy_new = if_else(Feed_soy_percent == 0, true = 0.0105, false = Feed_soy_percent),
         feed_crops_new = if_else(Feed_othercrops_percent == 0, true = 0.0105, false = Feed_othercrops_percent),
         feed_fmfo_new = if_else(Feed_FMFO_percent == 0, true = 0.0105, false = Feed_FMFO_percent),
         feed_animal_new = if_else(Feed_animal_percent == 0, true = 0.0105, false = Feed_animal_percent)) %>%
  # Renomoralize values so they sum to 1
  mutate(sum = rowSums(select(., contains("new")))) %>%
  mutate(feed_soy_new = feed_soy_new / sum,
         feed_crops_new = feed_crops_new / sum,
         feed_fmfo_new = feed_fmfo_new / sum,
         feed_animal_new = feed_animal_new / sum) %>%
  select(clean_sci_name, taxa_group_name, contains("new"))

######################################################################################################
# Model 1: Remove all NAs - estimate proportion feed for a set of studies of one species

# Remove NAs
lca_dat_no_na <- lca_dat_no_zeroes %>%
  filter(is.na(Feed_soy_percent)==FALSE) 

# Try to get dirichlet to work with just one set of studies: Oncorhynchus mykiss
# Set data for model:
k = 4
n = 3
feed_weights <- lca_dat_no_na %>%
  filter(clean_sci_name == "Oncorhynchus mykiss") %>%
  select(contains("new")) %>%
  as.matrix()

# note: dirichlet_rng is just a random number generator:
# rep_vector(x, m) creates a column consisting of m copies of x
# generated quantities {
#   vector[k] theta = dirichlet_rng(rep_vector(alpha, k));
# }

# Estimate feed component proportions for a single species
stan_pooled <- 'data {
  int<lower=0> n;  // number of observations
  int<lower=1> k; // number of feed types
  simplex[k] feed_weights[n]; // array of feed weights simplexes
}
parameters {
  vector<lower=0>[k] alpha;
  simplex[k] theta;
}
model {
  for (i in 1:n) {
    feed_weights[i] ~ dirichlet(alpha); // estimate vector of alphas based on the data of feed weights
  }
  theta ~ dirichlet(alpha); // now, estimate feed weights based on the vector of alphas
}'


no_missing_mod <- stan_model(model_code = stan_pooled, verbose = TRUE)
# Note: For Windows, apparently OK to ignore this warning message:
# Warning message:
#   In system(paste(CXX, ARGS), ignore.stdout = TRUE, ignore.stderr = TRUE) :
#   'C:/rtools40/usr/mingw_/bin/g++' not found

# Fit model:
fit_pooled <- sampling(object = no_missing_mod, data = list(n = n,
                                                            k = k,
                                                            feed_weights = feed_weights))
print(fit_pooled)

feeds <- c("soy", "crops", "fmfo", "animal")
feed_key <- data.frame(alpha_param = paste("alpha[", feeds, "]", sep = ""),
                       theta_param = paste("theta[", feeds, "]", sep = ""))

fit_pooled_clean <- fit_pooled
names(fit_pooled_clean)[grep(names(fit_pooled_clean), pattern = "alpha")] <- feed_key$alpha_param
names(fit_pooled_clean)[grep(names(fit_pooled_clean), pattern = "theta")] <- feed_key$theta_param


distribution_pooled <- as.matrix(fit_pooled_clean)

plot_theme <- theme(axis.text=element_text(size=14, color = "black"))

p_alpha <- mcmc_areas_ridges(distribution_pooled,
                       pars = vars(contains("alpha")),
                       prob = 0.8) + 
  ggtitle("Oncorhynchus mykiss feed proportion model", "with 80% credible intervals") +
  plot_theme

p_alpha
ggsave(filename = file.path(outdir, "bayes-example_trout_feed-proportion_alphas.png"), width = 11, height = 8.5)

p_theta <- mcmc_areas_ridges(distribution_pooled,
                             pars = vars(contains("theta")),
                             prob = 0.8) + 
  ggtitle("Oncorhynchus mykiss feed proportion model", "with 80% credible intervals") +
  plot_theme

p_theta
ggsave(filename = file.path(outdir, "bayes-example_trout_feed-proportion_thetas.png"), width = 11, height = 8.5)


######################################################################################################
# Model 2: Remove all NAs - estimate proportion feed for groups of scientific names in the dataset (but no hierarchies)

lca_dat_no_na <- lca_dat_no_zeroes %>%
  filter(is.na(feed_soy_new)==FALSE) 

#lca_groups <- lca_dat_no_na %>%
  #filter(clean_sci_name %in% c("Oncorhynchus mykiss", "Salmo salar")) %>% # converges
  #filter(clean_sci_name %in% c("Oncorhynchus mykiss", "Salmo salar", "Macrobrachium amazonicum")) %>% # creates divergent transitions
  #filter(clean_sci_name %in% c("Oncorhynchus mykiss", "Salmo salar", "Oreochromis niloticus")) %>% # converges
  #filter(clean_sci_name %in% c("Oncorhynchus mykiss", "Salmo salar", "Pangasius")) %>% # creates divergent transitions
  #filter(clean_sci_name %in% c("Oncorhynchus mykiss", "Salmo salar", "Penaeus monodon")) %>% # creates divergent transitions
  #filter(clean_sci_name %in% c("Penaeus monodon", "Salmo salar")) %>% # creates divergent transitions
  # # Add indices for each sci-name
  # mutate(clean_sci_name = as.factor(clean_sci_name),
  #        sci = as.numeric(clean_sci_name))

# lca_groups <- lca_dat_no_na %>%
#   filter(clean_sci_name %in% c("Macrobrachium amazonicum", "Penaeus monodon"))  %>%
#   # Add indices for each sci-name
#   mutate(clean_sci_name = as.factor(clean_sci_name),
#          sci = as.numeric(clean_sci_name))

# Now that alpha and theta are vectorized, can include all groups
# lca_groups <- lca_dat_no_na %>%
#   # Add indices for each sci-name
#   mutate(clean_sci_name = as.factor(clean_sci_name),
#          sci = as.numeric(clean_sci_name))

# Try including groups with only n>1; also remove Thunnus orientalis since both data points are identical (effectively n = 1)
lca_groups <- lca_dat_no_na %>%
  group_by(clean_sci_name) %>%
  mutate(n_sci = n()) %>%
  ungroup() %>%
  filter(n_sci > 1) %>%
  filter(clean_sci_name != "Thunnus orientalis") %>%
  # Add indices for each sci-name
  mutate(clean_sci_name = as.factor(clean_sci_name),
         sci = as.numeric(clean_sci_name))

feed_vars <- c("feed_soy_new", "feed_crops_new", "feed_fmfo_new", "feed_animal_new")
for (i in 1:length(feed_vars)) {
  p <- ggplot(lca_groups, aes(x = clean_sci_name, y = !!sym(feed_vars[i]))) +
    geom_boxplot() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 16)) + 
    labs(title = "Boxplots of feed proportion by scientific name")
  print(p)
}

# Set data for model:
feed_weights <- lca_groups %>%
  select(contains("new")) %>%
  as.matrix()
k = 4
n = nrow(feed_weights)
n_sci = length(unique(lca_groups$sci))
sci = lca_groups$sci

# SIMULATE FAKE DATA TO TEST MODEL
# library(MCMCpack)
# samp_1 <- rdirichlet(n = 10, alpha = c(1,1,1,1))
# samp_2 <- rdirichlet(n = 10, alpha = c(10, 1, 1, 1))
# feed_weights <- rbind(samp_1, samp_2)
# k = 4
# n = nrow(feed_weights)
# n_sci = 2
# sci = c(rep(1, n/2), rep (2, n/2))

# Vectorize over alpha and theta
stan_pooled <- 'data {
  int n;  // number of observations
  int k; // number of feed types
  int n_sci;
  simplex[k] feed_weights[n]; // array of observed feed weights simplexes
  int sci[n]; // sci-name indices
}
parameters {
  vector<lower=0>[k] alpha[n_sci]; // vector of dirichlet priors, one for each sci name
  simplex[k] theta[n_sci]; // vector of estimated sci-level feed weight simplexes;
}
model {
  // priors on alpha
  //for (m in 1:k){
  //  alpha[n_sci][m] ~ uniform(0.1, 10);
  //}

  
  for (i in 1:n) {
    feed_weights[i] ~ dirichlet(to_vector(alpha[sci[i]]));
  }
  // now, estimate feed weights based on the vector of alphas
  for (j in 1:n_sci) {
    theta[j] ~ dirichlet(to_vector(alpha[j]));
  }
}'

# # Translated and scaled simplex:
# From: https://mc-stan.org/docs/2_21/stan-users-guide/parameterizing-centered-vectors.html
# stan_pooled <- 'data {
#   int n;  // number of observations
#   int k; // number of feed types
#   int n_sci;
#   simplex[k] feed_weights[n]; // array of observed feed weights simplexes
#   int sci[n]; // sci-name indices
# }
# parameters {
#   vector<lower=0>[k] alpha[n_sci]; // vector of dirichlet priors, one for each sci name
#   simplex[k] theta_raw[n_sci]; // vector of estimated sci-level feed weight simplexes;
#   real theta_scale[n_sci];
# }
# transformed parameters {
#   vector[k] theta;
#   for (j in 1:n_sci) {
#     theta = theta_scale[j] * (theta_raw[j] - inv(k));
#   }
# 
# }
# model {
#   for (i in 1:n) {
#     feed_weights[i] ~ dirichlet(to_vector(alpha[sci[i]]));
#   }
#   // now, estimate feed weights based on the vector of alphas
#   for (j in 1:n_sci) {
#     theta_raw[j] ~ dirichlet(to_vector(alpha[j]));
#   }
# }'


no_missing_mod <- stan_model(model_code = stan_pooled, verbose = TRUE)
# Note: For Windows, apparently OK to ignore this warning message:
# Warning message:
#   In system(paste(CXX, ARGS), ignore.stdout = TRUE, ignore.stderr = TRUE) :
#   'C:/rtools40/usr/mingw_/bin/g++' not found

# Fit model:
fit_grouped <- sampling(object = no_missing_mod, data = list(n = n,
                                                             k = k,
                                                             feed_weights = feed_weights,
                                                             n_sci = n_sci,
                                                             sci = sci),
                        cores = 4, seed = "11729")
                        #cores = 4, iter = 10000) # iter = 10000
                        # control = list(adapt_delta = 0.99)) # address divergent transitions by increasing delta, i.e., take smaller steps
print(fit_grouped)

# Format of parameters is: theta[sci_name, feed]
sci_feed_key <- lca_groups %>%
  select(contains(c("clean_sci_name", "new", "sci"))) %>%
  pivot_longer(cols = contains("new"), names_to = "feed") %>%
  select(-value) %>%
  unique() %>%
  mutate(feed_index = case_when(str_detect(feed, "soy") ~ 1,
                                str_detect(feed, "crops") ~ 2,
                                str_detect(feed, "fmfo") ~ 3,
                                str_detect(feed, "animal") ~ 4)) %>%
  # Clean feed names
  mutate(feed = gsub(feed, pattern = "feed_", replacement = "")) %>%
  mutate(feed = gsub(feed, pattern = "_new", replacement = "")) %>%
  mutate(param_name = paste("[", sci, ",", feed_index, "]", sep = "")) %>%
  mutate(alpha_param_name = paste("alpha", clean_sci_name, feed, sep = "-")) %>%
  mutate(theta_param_name = paste("theta", clean_sci_name, feed, sep = "-")) %>%
  # IMPORTANT before replaceing param names: ARRANGE BY FEED, THEN SCIENTIFIC NAME TO MATCH HOW NAMES ARE ARRANGED IN STANFIT OBJECT
  arrange(feed_index, sci)

# Replace param names; first copy to fit_grouped_clean to avoid having to re-run sampling as a result of doing something wrong to fit_grouped
fit_grouped_clean <- fit_grouped
names(fit_grouped_clean)[grep(names(fit_grouped_clean), pattern = "alpha\\[")] <- sci_feed_key$alpha_param_name
names(fit_grouped_clean)[grep(names(fit_grouped_clean), pattern = "theta\\[")] <- sci_feed_key$theta_param_name

distribution_grouped <- as.matrix(fit_grouped_clean)
p_alpha <- mcmc_areas(distribution_grouped,
                      pars = vars(contains("alpha")),
                      prob = 0.8,
                      area_method = "scaled height") + ggtitle("")
p_alpha
p_theta <- mcmc_areas(distribution_grouped,
                      pars = vars(contains("theta")),
                      prob = 0.8,
                      area_method = "scaled height") + ggtitle("")
p_theta


######################################################################################################
# Model 2.1: Same as model 2 but with informative priors:
# Remove all NAs - estimate proportion feed for just two scientific names in the dataset

lca_dat_no_na <- lca_dat_no_zeroes %>%
  filter(is.na(feed_soy_new)==FALSE) 

lca_groups <- lca_dat_no_na %>%
  filter(clean_sci_name %in% c("Oncorhynchus mykiss", "Salmo salar")) %>%
  #filter(clean_sci_name %in% c("Oncorhynchus mykiss", "Salmo salar", "Macrobrachium amazonicum")) %>% # creates divergent transitions
  #filter(clean_sci_name %in% c("Oncorhynchus mykiss", "Salmo salar", "Oreochromis niloticus")) %>% # converges
  #filter(clean_sci_name %in% c("Oncorhynchus mykiss", "Salmo salar", "Pangasius")) %>% # creates divergent transitions
  #filter(clean_sci_name %in% c("Oncorhynchus mykiss", "Salmo salar", "Penaeus monodon")) %>% # creates divergent transitions
  #filter(clean_sci_name %in% c("Penaeus monodon", "Salmo salar")) %>% # creates divergent transitions
  # Add indices for each sci-name
  mutate(clean_sci_name = as.factor(clean_sci_name),
         sci = as.numeric(clean_sci_name))

# lca_groups <- lca_dat_no_na %>%
#   filter(clean_sci_name %in% c("Macrobrachium amazonicum", "Penaeus monodon"))  %>%
#   # Add indices for each sci-name
#   mutate(clean_sci_name = as.factor(clean_sci_name),
#          sci = as.numeric(clean_sci_name))

# Now that alpha and theta are vectorized, can include all groups
# lca_groups <- lca_dat_no_na %>%
#   # Add indices for each sci-name
#   mutate(clean_sci_name = as.factor(clean_sci_name),
#          sci = as.numeric(clean_sci_name))

# Set data for model:
feed_weights <- lca_groups %>%
  select(contains("new")) %>%
  as.matrix()
k = 4
n = nrow(feed_weights)
n_sci = length(unique(lca_groups$sci))
sci = lca_groups$sci
# Get the mean observations across all sci-names
phi_mean <- lca_groups %>% 
  group_by(sci) %>% 
  summarise(across(where(is.numeric), mean),
            n_obs = n()) %>%
  ungroup() %>%
  select(contains(c("new", "sci", "obs"))) %>%
  arrange(sci)
phi <- phi_mean %>%
  select(contains("new")) %>%
  as.matrix() 
kappa <- phi_mean %>% pull(n_obs) + k



# OLD CODE for just two groups:
# Prior specification following: https://mc-stan.org/docs/2_18/stan-users-guide/reparameterizations.html
# stan_pooled <- 'data {
#   int<lower=0> n;  // number of observations
#   int<lower=1> k; // number of feed types
#   simplex[k] feed_weights[n]; // array of observed feed weights simplexes
#   int<lower=1, upper=2> sci[n]; // sci-name indices
#   simplex[k] phi_1;
#   simplex[k] phi_2;
#   real<lower=0> kappa_1;
#   real<lower=0> kappa_2;
# }
# parameters {
#   simplex[k] theta_1; // vector of estimated sci-level feed weight simplexes; FIX IT - alpha must be a vector but should be able to set up theta as a list of vectors (just like feed_weights)
#   simplex[k] theta_2;
# }
# transformed parameters {
#   // reparameterize alpha distributions as a vector of means and counts
#   // phi is expected value of theta (mean feed weights)
#   // strength of the prior measured in number of prior observations
#   vector[k] alpha_1 = kappa_1 * phi_1;
#   vector[k] alpha_2 = kappa_2 * phi_2;
# }
# model {
#   for (i in 1:n) {
#     if (sci[i]==1){
#       feed_weights[i] ~ dirichlet(alpha_1); // estimate vector of alphas based on the data of feed weights
#     }
#     if (sci[i]==2){
#       feed_weights[i] ~ dirichlet(alpha_2); // estimate vector of alphas based on the data of feed weights
#     }
#   }
#   // now, estimate feed weights based on the vector of alphas
#   theta_1 ~ dirichlet(alpha_1);
#   theta_2 ~ dirichlet(alpha_2);
# }'


# This code vectorizes over alpha and theta, allowing all groups to be estiamted
# this stan_data list passes phi in as data
# stan_data = list(n = n,
#                  k = k,
#                  feed_weights = feed_weights,
#                  n_sci = n_sci,
#                  sci = sci,
#                  phi = phi,
#                  kappa = kappa)

# Code that passes priors in as data
# stan_pooled <- 'data {
#   int n;  // number of observations
#   int k; // number of feed types
#   int n_sci; // number of sci names
#   simplex[k] feed_weights[n]; // array of observed feed weights simplexes
#   int sci[n]; // sci-name indices
#   simplex[k] phi[n_sci];
#   int kappa[n_sci];
# }
# parameters {
#   // alpha parameter now moved into transformed parameter section
#   simplex[k] theta[n_sci]; // vectors of estimated sci-level feed weight simplexes;
# }
# transformed parameters {
#   // reparameterize alpha distributions as a vector of means and counts
#   // phi is expected value of theta (mean feed weights)
#   // kappa is strength of the prior measured in number of prior observations (minus K)
#   vector<lower=0>[k] alpha[n_sci];
#   for (m in 1:n) {
#     alpha[sci[m]] = kappa[sci[m]] * phi[sci[m]];
#   }
# }
# model {
# 
#   for (i in 1:n) {
#     feed_weights[i] ~ dirichlet(to_vector(alpha[sci[i]]));
#     // theta[sci[i]] ~ dirichlet(to_vector(alpha[sci[i]])); // this has problems converging here
#   }
#   // now, estimate feed weights based on the vector of alphas
#   for (j in 1:n_sci) {
#     theta[j] ~ dirichlet(to_vector(alpha[j]));
#   }
# }'


# NEW CODE: instead of passing phi in as data, pass it as a parameter with a distribution
# Appears like model is only valid when only one element in the phi simplex (per scientific name) is given a prior

# this stan_data list only defines kappa (not phi) as data
stan_data = list(n = n,
                 k = k,
                 feed_weights = feed_weights,
                 n_sci = n_sci,
                 sci = sci,
                 kappa = kappa)


stan_pooled <- 'data {
  int n;  // number of observations
  int k; // number of feed types
  int n_sci; // number of sci names
  simplex[k] feed_weights[n]; // array of observed feed weights simplexes
  int sci[n]; // sci-name indices
  int kappa[n_sci];
}
parameters {
  // alpha parameter now moved into transformed parameter section
  simplex[k] phi[n_sci];
  simplex[k] theta[n_sci]; // vectors of estimated sci-level feed weight simplexes
  
  // sigma parameters for mean priors
  real<lower=0> sigma_1;
  // real<lower=0> sigma_2;
}
transformed parameters {
  // reparameterize alpha distributions as a vector of means and counts
  // phi is expected value of theta (mean feed weights)
  // kappa is strength of the prior measured in number of prior observations (minus K)
  vector<lower=0>[k] alpha[n_sci];
  for (m in 1:n) {
    alpha[sci[m]] = kappa[sci[m]] * phi[sci[m]];
  }
}
model {
  // priors on specific phi
  // phi defined as phi[sci][k]
  
  // option 1: define feed proportion priors as lower upper bounds (but can only give a prior for one element per simplex - i.e., priors on phi[6][1] and phi[6][2] causes error probably because elements within a simplex are constrained?)
  // phi[sci][k] ~ uniform(0.1, 0.2); // example prior on lower and upper bounds
  
  // option 2: define feed proportions as means (need to define sigmas in parameters block: real<lower=0> sigma_1, sigma_2 etc; etc;)
  // phi[sci][k] ~ normal(0.13, sigma_1); // example prior on mean
  sigma_1 ~ uniform(0, 10);
  phi[1][1] ~ normal(0.13, sigma_1);

  for (i in 1:n) {
    feed_weights[i] ~ dirichlet(to_vector(alpha[sci[i]]));
    // theta[sci[i]] ~ dirichlet(to_vector(alpha[sci[i]])); // this has problems converging here
  }
  // now, estimate feed weights based on the vector of alphas
  for (j in 1:n_sci) {
    theta[j] ~ dirichlet(to_vector(alpha[j]));
  }
}'


no_missing_mod <- stan_model(model_code = stan_pooled, verbose = TRUE)
# Note: For Windows, apparently OK to ignore this warning message:
# Warning message:
#   In system(paste(CXX, ARGS), ignore.stdout = TRUE, ignore.stderr = TRUE) :
#   'C:/rtools40/usr/mingw_/bin/g++' not found

# RUNS but gives warning about divergent transitions
# Fit model:
fit_grouped <- sampling(object = no_missing_mod, data = stan_data, cores = 4, seed = "11729")
#,cores = 4, iter = 10000,
#control = list(adapt_delta = 0.99)) # address divergent transitions by increasing delta, i.e., take smaller steps
print(fit_grouped)

launch_shinystan(fit_grouped)

# Format of parameters is: theta[sci_name, feed]
sci_feed_key <- lca_groups %>%
  select(contains(c("clean_sci_name", "new", "sci"))) %>%
  pivot_longer(cols = contains("new"), names_to = "feed") %>%
  select(-value) %>%
  unique() %>%
  mutate(feed_index = case_when(str_detect(feed, "soy") ~ 1,
                                str_detect(feed, "crops") ~ 2,
                                str_detect(feed, "fmfo") ~ 3,
                                str_detect(feed, "animal") ~ 4)) %>%
  # Clean feed names
  mutate(feed = gsub(feed, pattern = "feed_", replacement = "")) %>%
  mutate(feed = gsub(feed, pattern = "_new", replacement = "")) %>%
  mutate(param_name = paste("[", sci, ",", feed_index, "]", sep = "")) %>%
  mutate(alpha_param_name = paste("alpha", clean_sci_name, feed, sep = "-")) %>%
  mutate(theta_param_name = paste("theta", clean_sci_name, feed, sep = "-")) %>%
  # IMPORTANT before replaceing param names: ARRANGE BY FEED, THEN SCIENTIFIC NAME TO MATCH HOW NAMES ARE ARRANGED IN STANFIT OBJECT
  arrange(feed_index, sci)

# Replace param names; first copy to fit_grouped_clean to avoid having to re-run sampling as a result of doing something wrong to fit_grouped
fit_grouped_clean <- fit_grouped
names(fit_grouped_clean)[grep(names(fit_grouped_clean), pattern = "alpha\\[")] <- sci_feed_key$alpha_param_name
names(fit_grouped_clean)[grep(names(fit_grouped_clean), pattern = "theta\\[")] <- sci_feed_key$theta_param_name

distribution_grouped <- as.matrix(fit_grouped_clean)
p_alpha <- mcmc_areas(distribution_grouped,
                      pars = vars(contains("alpha")),
                      prob = 0.8,
                      area_method = "scaled height")
p_alpha
p_theta <- mcmc_areas(distribution_grouped,
                      pars = vars(contains("theta")),
                      prob = 0.8,
                      area_method = "scaled height")
p_theta

######################################################################################################
# Model 3: Still no NAs, make into three-level model
# Remove all NAs - estimate proportion feed for just two scientific names in the dataset

lca_dat_no_na <- lca_dat_no_zeroes %>%
  filter(is.na(feed_soy_new)==FALSE)

# Keep all data, but Add indices
# lca_groups <- lca_dat_no_na %>%
#   # Add indices for each sci-name
#   mutate(clean_sci_name = as.factor(clean_sci_name),
#          sci = as.numeric(clean_sci_name),
#          taxa_group_name = as.factor(taxa_group_name),
#          tx = as.numeric(taxa_group_name)) %>%
#   arrange(sci)

# Test a smaller dataset (just salmon/char and marine shrimp - i.e., two taxa levels + overall level)
lca_groups <- lca_dat_no_na %>%
  filter(taxa_group_name %in% c("salmon/char", "marine shrimp")) %>%
  mutate(clean_sci_name = as.factor(clean_sci_name),
         sci = as.numeric(clean_sci_name),
         taxa_group_name = as.factor(taxa_group_name),
         tx = as.numeric(taxa_group_name)) %>%
  arrange(sci)


# # Try analyzing only groups with  n>1; also remove Thunnus orientalis since both data points are identical (effectively n = 1)
# lca_groups <- lca_dat_no_na %>%
#   group_by(clean_sci_name) %>%
#   mutate(n_sci = n()) %>%
#   ungroup() %>%
#   filter(n_sci > 1) %>%
#   filter(clean_sci_name != "Thunnus orientalis") %>%
#   # Add indices for each sci-name
#   mutate(clean_sci_name = as.factor(clean_sci_name),
#          sci = as.numeric(clean_sci_name),
#          taxa_group_name = as.factor(taxa_group_name),
#          tx = as.numeric(taxa_group_name)) %>%
#   arrange(sci)

# Set data for model:
feed_weights <- lca_groups %>%
  select(contains("new")) %>%
  as.matrix()
K = 4
N = nrow(feed_weights)
N_SCI = length(unique(lca_groups$sci))
N_TX = length(unique(lca_groups$tx))
sci = lca_groups$sci
tx = lca_groups$tx

# Get counts per sci name and counts per taxa group (also included as data in the model):
sci_kappa <- lca_groups %>% 
  select(contains(c("new", "sci", "obs"))) %>%
  group_by(sci) %>% 
  summarise(n_obs = n()) %>%
  ungroup() %>%
  arrange(sci) %>%
  pull(n_obs)

tx_kappa <- lca_groups %>% 
  select(contains(c("new", "tx", "obs"))) %>%
  group_by(tx) %>% 
  summarise(n_obs = n()) %>%
  ungroup() %>%
  arrange(tx) %>%
  pull(n_obs)

# For priors, get the mean of observations per sci-name
sci_phi_mean <- lca_groups %>% 
  select(contains(c("new", "sci"))) %>%
  group_by(sci) %>% 
  summarise(across(contains("new"), mean)) %>%
  ungroup() %>%
  arrange(sci)

# Get mean observations per taxa group
tx_phi_mean <- lca_groups %>% 
  select(contains(c("new", "tx"))) %>%
  group_by(tx) %>% 
  summarise(across(contains("new"), mean)) %>%
  ungroup() %>%
  arrange(tx)

# Three-level model with no priors:
stan_data = list(N = N,
                 K = K,
                 feed_weights = feed_weights,
                 N_SCI = N_SCI,
                 N_TX = N_TX,
                 sci = sci,
                 tx = tx)

# stan_pooled <- 'data {
#   int N;  // number of total observations
#   int K; // number of feed types
#   int N_SCI; // number of sci names
#   int N_TX; // number of taxa groups
#   simplex[K] feed_weights[N]; // array of observed feed weights simplexes
#   int sci[N]; // sci-name indices
#   int tx[N]; // taxa-group indices
# }
# parameters {
#   vector<lower=0>[K] sci_alpha[N_SCI]; // vector of dirichlet priors, one for each sci name
#   simplex[K] sci_theta[N_SCI]; // vectors of estimated sci-level feed weight simplexes
#   vector<lower=0>[K] tx_alpha[N_TX];
#   simplex[K] tx_theta[N_TX];
#   vector[K] alpha;
#   simplex[K] theta;
# }
# model {
# 
#   // likelihood
#   for (n in 1:N) {
#     tx_theta[tx[n]] ~ dirichlet(alpha);
#     sci_theta[sci[n]] ~ dirichlet(to_vector(tx_alpha[tx[n]]));
#     feed_weights[n] ~ dirichlet(to_vector(sci_alpha[sci[n]]));
#   }
#   // now, estimate feed weights based on the vector of alphas
#   theta ~ dirichlet(to_vector(alpha));
#   for (n_tx in 1:N_TX) {
#     tx_theta[n_tx] ~ dirichlet(to_vector(tx_alpha[n_tx]));
#   }
#   for (n_sci in 1:N_SCI) {
#     sci_theta[n_sci] ~ dirichlet(to_vector(sci_alpha[n_sci]));
#   }
# }'

# Three level model (no priors), but to help with convergence, try offsetting and scaling simplex:
# From: https://mc-stan.org/docs/2_21/stan-users-guide/parameterizing-centered-vectors.html
stan_data = list(N = N,
                 K = K,
                 feed_weights = feed_weights,
                 N_SCI = N_SCI,
                 N_TX = N_TX,
                 sci = sci,
                 tx = tx)

stan_pooled <- 'data {
  int N;  // number of total observations
  int K; // number of feed types
  int N_SCI; // number of sci names
  int N_TX; // number of taxa groups
  simplex[K] feed_weights[N]; // array of observed feed weights simplexes
  int sci[N]; // sci-name indices
  int tx[N]; // taxa-group indices
}
parameters {
  vector<lower=0>[K] sci_alpha[N_SCI]; // vector of dirichlet priors, one for each sci name (alpha is not a simplex)
  simplex[K] sci_theta_raw[N_SCI]; // vectors of estimated sci-level feed weight simplexes
  vector<lower=0>[K] tx_alpha[N_TX];
  simplex[K] tx_theta_raw[N_TX];
  vector<lower=0>[K] alpha;
  simplex[K] theta_raw;

  // scaling parameters
  real sci_theta_scale[N_SCI]; // vectors of estimated sci-level feed weight simplexes
  real tx_theta_scale[N_TX];
  real theta_scale;
}
transformed parameters {
  vector[K] sci_theta[N_SCI];
  vector[K] tx_theta[N_TX];
  vector[K] theta;

  for (n_sci in 1:N_SCI){
    sci_theta[N_SCI] = sci_theta_scale[N_SCI] * (sci_theta_raw[N_SCI] - inv(K));
  }

  for (n_tx in 1:N_TX) {
    tx_theta[N_TX] = tx_theta_scale[N_TX] * (tx_theta_raw[N_TX] - inv(K));
  }
  theta = theta_scale * (theta_raw - inv(K));

}
model {

  // likelihood
  for (n in 1:N) {
    tx_theta_raw[tx[n]] ~ dirichlet(alpha);
    sci_theta_raw[sci[n]] ~ dirichlet(to_vector(tx_alpha[tx[n]]));
    feed_weights[n] ~ dirichlet(to_vector(sci_alpha[sci[n]]));
  }
  // now, estimate feed weights based on the vector of alphas
  theta_raw ~ dirichlet(to_vector(alpha));
  for (n_tx in 1:N_TX) {
    tx_theta_raw[n_tx] ~ dirichlet(to_vector(tx_alpha[n_tx]));
  }
  for (n_sci in 1:N_SCI) {
    sci_theta_raw[n_sci] ~ dirichlet(to_vector(sci_alpha[n_sci]));
  }
}'

# Three-level model with priors
# Appears like model is only valid when only one element in the phi simplex (per scientific name) is given a prior
# this stan_data list only defines kappa (not phi) as data
# stan_data = list(N = N,
#                  K = K,
#                  feed_weights = feed_weights,
#                  N_SCI = N_SCI,
#                  N_TX = N_TX,
#                  sci = sci,
#                  tx = tx,
#                  sci_kappa = sci_kappa,
#                  tx_kappa = tx_kappa)
# 
# stan_pooled <- 'data {
#   int N;  // number of total observations
#   int K; // number of feed types
#   int N_SCI; // number of sci names
#   int N_TX; // number of taxa groups
#   simplex[K] feed_weights[N]; // array of observed feed weights simplexes
#   int sci[N]; // sci-name indices
#   int tx[N]; // taxa-group indices
#   int sci_kappa[N_SCI]; // number of observations per sci-name
#   int tx_kappa[N_TX]; // number of observations per taxa group
# }
# parameters {
#   simplex[K] sci_phi[N_SCI]; // vectors of sci-level feed weight priors
#   simplex[K] sci_theta[N_SCI]; // vectors of estimated sci-level feed weight simplexes
#   simplex[K] tx_phi[N_TX];
#   simplex[K] tx_theta[N_TX];
#   simplex[K] phi;
#   simplex[K] theta;
# }
# transformed parameters {
#   // define params
#   vector<lower=0>[K] sci_alpha[N_SCI];
#   vector<lower=0>[K] tx_alpha[N_TX];
#   vector<lower=0>[K] alpha;
#   
#   // reparameterize alphas as a vector of means (phi) and counts (kappas)
#   // phi is expected value of theta (mean feed weights)
#   // kappa is strength of the prior measured in number of prior observations (minus K)
#   alpha = N * phi;
#   for (n_tx in 1:N_TX) {
#     tx_alpha[n_tx] = tx_kappa[n_tx] * tx_phi[n_tx];
#   }    
#   
#   for (n_sci in 1:N_SCI) {
#     sci_alpha[n_sci] = sci_kappa[n_sci] * sci_phi[n_sci];
#   }
# }
# model {
#   // priors on specific phi
#   // sci_phi defined as sci_phi[sci][K]
#   
#   // option 1: define feed proportion priors as lower upper bounds
#   sci_phi[2][1] ~ uniform(0.1, 0.2); // hypothetical lower and upper bounds for Oncorhynchus mykiss soy 
#   
#   // option 2: define feed proportions as means (need to define sigmas in parameters block: real<lower=0> sigma_1, sigma_2 etc; etc;)
#   // sci_phi[2][2] ~ normal(0.13, sigma_1); // mean for Oncorhynhchus mykiss soy feed 
# 
#   // likelihood
#   for (n in 1:N) {
#     tx_phi[tx[n]] ~ dirichlet(alpha);
#     sci_phi[sci[n]] ~ dirichlet(to_vector(tx_alpha[tx[n]]));
#     feed_weights[n] ~ dirichlet(to_vector(sci_alpha[sci[n]])); 
#   }
#   // now, estimate feed weights based on the vector of alphas
#   theta ~ dirichlet(to_vector(alpha));
#   for (n_tx in 1:N_TX) {
#     tx_theta[n_tx] ~ dirichlet(to_vector(tx_alpha[n_tx]));
#   }
#   for (n_sci in 1:N_SCI) {
#     sci_theta[n_sci] ~ dirichlet(to_vector(sci_alpha[n_sci]));
#   }
# }'

no_missing_mod <- stan_model(model_code = stan_pooled, verbose = TRUE)
# Note: For Windows, apparently OK to ignore this warning message:
# Warning message:
#   In system(paste(CXX, ARGS), ignore.stdout = TRUE, ignore.stderr = TRUE) :
#   'C:/rtools40/usr/mingw_/bin/g++' not found

# RUNS but gives warning about divergent transitions
# Fit model:
fit_grouped <- sampling(object = no_missing_mod, data = stan_data, cores = 4, seed = "11729")
# Increasing adapt_delta decreases the divergences but doesn't get rid of them
# fit_grouped <- sampling(object = no_missing_mod, data = stan_data, cores = 4, seed = "11720", control = list(adapt_delta = 0.9))
# fit_grouped <- sampling(object = no_missing_mod, data = stan_data, cores = 4, seed = "11720", control = list(adapt_delta = 0.99))


launch_shinystan(fit_grouped)

#,cores = 4, iter = 10000,
#control = list(adapt_delta = 0.99)) # address divergent transitions by increasing delta, i.e., take smaller steps
print(fit_grouped)

# Create formatted names for sci-name level
# Format of parameters is: theta[sci_name, feed]
sci_feed_key <- lca_groups %>%
  select(contains(c("clean_sci_name", "new", "sci"))) %>%
  pivot_longer(cols = contains("new"), names_to = "feed") %>%
  select(-value) %>%
  unique() %>%
  mutate(feed_index = case_when(str_detect(feed, "soy") ~ 1,
                                str_detect(feed, "crops") ~ 2,
                                str_detect(feed, "fmfo") ~ 3,
                                str_detect(feed, "animal") ~ 4)) %>%
  # Clean feed names
  mutate(feed = gsub(feed, pattern = "feed_", replacement = "")) %>%
  mutate(feed = gsub(feed, pattern = "_new", replacement = "")) %>%
  mutate(index = paste("[", sci, ",", feed_index, "]", sep = "")) %>%
  mutate(sci_phi_param_name = paste("phi[", clean_sci_name, ", ", feed, "]", sep = "")) %>%
  mutate(sci_theta_param_name = paste("theta[", clean_sci_name, ", ", feed, "]", sep = "")) %>%
  mutate(sci_alpha_param_name = paste("alpha[", clean_sci_name, ", ", feed, "]", sep = "")) %>%
  # IMPORTANT before replaceing param names: ARRANGE BY FEED, THEN SCIENTIFIC NAME TO MATCH HOW NAMES ARE ARRANGED IN STANFIT OBJECT
  arrange(feed_index, sci)

# Create formatted names for taxa-group level
# Format of parameters is: theta[sci_name, feed]
tx_feed_key <- lca_groups %>%
  select(contains(c("taxa_group_name", "new", "tx"))) %>%
  pivot_longer(cols = contains("new"), names_to = "feed") %>%
  select(-value) %>%
  unique() %>%
  mutate(feed_index = case_when(str_detect(feed, "soy") ~ 1,
                                str_detect(feed, "crops") ~ 2,
                                str_detect(feed, "fmfo") ~ 3,
                                str_detect(feed, "animal") ~ 4)) %>%
  # Clean feed names
  mutate(feed = gsub(feed, pattern = "feed_", replacement = "")) %>%
  mutate(feed = gsub(feed, pattern = "_new", replacement = "")) %>%
  mutate(index = paste("[", tx, ",", feed_index, "]", sep = "")) %>%
  mutate(tx_phi_param_name = paste("phi[", taxa_group_name, ", ", feed, "]", sep = "")) %>%
  mutate(tx_theta_param_name = paste("theta[", taxa_group_name, ", ", feed, "]", sep = "")) %>%
  mutate(tx_alpha_param_name = paste("alpha[", taxa_group_name, ", ", feed, "]", sep = "")) %>%
  # IMPORTANT before replaceing param names: ARRANGE BY FEED, THEN TAXA NAME TO MATCH HOW NAMES ARE ARRANGED IN STANFIT OBJECT
  arrange(feed_index, tx)

overall_feed_key <- lca_groups %>%
  select(contains("new")) %>%
  pivot_longer(cols = contains("new"), names_to = "feed") %>%
  select(-value) %>%
  unique() %>%
  mutate(feed_index = case_when(str_detect(feed, "soy") ~ 1,
                                str_detect(feed, "crops") ~ 2,
                                str_detect(feed, "fmfo") ~ 3,
                                str_detect(feed, "animal") ~ 4)) %>%
  # Clean feed names
  mutate(feed = gsub(feed, pattern = "feed_", replacement = "")) %>%
  mutate(feed = gsub(feed, pattern = "_new", replacement = "")) %>%
  mutate(overall_phi_param_name = paste("phi[overall, ", feed, "]", sep = "")) %>%
  mutate(overall_theta_param_name = paste("theta[overall, ", feed, "]", sep = "")) %>%
  mutate(overall_alpha_param_name = paste("alpha[overall, ", feed, "]", sep = "")) %>%
  # IMPORTANT before replaceing param names: ARRANGE BY FEED TO MATCH HOW NAMES ARE ARRANGED IN STANFIT OBJECT
  arrange(feed_index)

# Replace param names; first copy to fit_grouped_clean to avoid having to re-run sampling as a result of doing something wrong to fit_grouped
fit_grouped_clean <- fit_grouped
# Sci-Level
names(fit_grouped_clean)[grep(names(fit_grouped_clean), pattern = "sci_phi")] <- sci_feed_key$sci_phi_param_name
names(fit_grouped_clean)[grep(names(fit_grouped_clean), pattern = "sci_alpha")] <- sci_feed_key$sci_alpha_param_name
names(fit_grouped_clean)[grep(names(fit_grouped_clean), pattern = "sci_theta")] <- sci_feed_key$sci_theta_param_name
# Taxa-level
names(fit_grouped_clean)[grep(names(fit_grouped_clean), pattern = "tx_phi")] <- tx_feed_key$tx_phi_param_name
names(fit_grouped_clean)[grep(names(fit_grouped_clean), pattern = "tx_alpha")] <- tx_feed_key$tx_alpha_param_name
names(fit_grouped_clean)[grep(names(fit_grouped_clean), pattern = "tx_theta")] <- tx_feed_key$tx_theta_param_name
# Global-level
names(fit_grouped_clean)[grep(names(fit_grouped_clean), pattern = "phi\\[[1-4]")] <- overall_feed_key$overall_phi_param_name
names(fit_grouped_clean)[grep(names(fit_grouped_clean), pattern = "alpha\\[[1-4]")] <- overall_feed_key$overall_alpha_param_name
names(fit_grouped_clean)[grep(names(fit_grouped_clean), pattern = "theta\\[[1-4]")] <- overall_feed_key$overall_theta_param_name

distribution_grouped <- as.matrix(fit_grouped_clean)
p_alpha <- mcmc_areas(distribution_grouped,
                      pars = vars(contains("alpha")),
                      prob = 0.8,
                      area_method = "scaled height")
p_alpha
p_theta <- mcmc_areas(distribution_grouped,
                      pars = vars(contains("theta")),
                      prob = 0.8,
                      area_method = "scaled height")
p_theta


