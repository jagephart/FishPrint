# Bayesian estimation of proportions of each feed component (soy, other crops, FMFOs, and other animal)

rm(list=ls())
library(tidyverse)
library(rstan)
library(taxize)
library(data.table)
library(countrycode) # part of clean.lca
library(bayesplot) # for mcmc_areas_ridges

# Mac
datadir <- "/Volumes/jgephart/BFA Environment 2/Data"
outdir <- "/Volumes/jgephart/BFA Environment 2/Outputs"
# Windows
# datadir <- "K:/BFA Environment 2/Data"
# outdir <- "K:BFA Environment 2/Outputs"

lca_dat <- read.csv(file.path(datadir, "LCA_compiled_20201006.csv"), fileEncoding="UTF-8-BOM") #fileEncoding needed when reading in file from windows computer (suppresses BOM hidden characters)
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
         feed_animal_new = feed_animal_new / sum)

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
  filter(is.na(Feed_soy_percent)==FALSE) 

lca_groups <- lca_dat_no_na %>%
  filter(clean_sci_name %in% c("Oncorhynchus mykiss", "Salmo salar")) %>%
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

# OLD CODE: create separate alpha's for each group; alpha must be a vector but should be able to set up theta as a list of vectors (just like feed_weights)
# Similar to: https://www.alexpghayes.com/blog/some-things-ive-learned-about-stan/
# Estimate feed component proportions for two species
# stan_pooled <- 'data {
#   int<lower=0> n;  // number of observations
#   int<lower=1> k; // number of feed types
#   simplex[k] feed_weights[n]; // array of observed feed weights simplexes
#   int<lower=1, upper=2> sci[n]; // sci-name indices
# }
# parameters {
#   vector<lower=0>[k] alpha_1; // alpha MUST be a vector, otherwise warning: vector ~ dirichlet(vector);
#   simplex[k] theta_1; // vector of estimated sci-level feed weight simplexes; FIX IT - alpha must be a vector but should be able to set up theta as a list of vectors (just like feed_weights)
#   vector<lower=0>[k] alpha_2; //
#   simplex[k] theta_2;
# }
# model {
# 
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
  for (i in 1:n) {
    feed_weights[i] ~ dirichlet(to_vector(alpha[sci[i]]));
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

# Fit model:
fit_grouped <- sampling(object = no_missing_mod, data = list(n = n,
                                                             k = k,
                                                             feed_weights = feed_weights,
                                                             n_sci = n_sci,
                                                             sci = sci),
                        cores = 4)
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
  filter(is.na(Feed_soy_percent)==FALSE) 

# lca_groups <- lca_dat_no_na %>%
#   filter(clean_sci_name %in% c("Oncorhynchus mykiss", "Salmo salar")) %>%
#   # Add indices for each sci-name
#   mutate(clean_sci_name = as.factor(clean_sci_name),
#          sci = as.numeric(clean_sci_name))

# lca_groups <- lca_dat_no_na %>%
#   filter(clean_sci_name %in% c("Macrobrachium amazonicum", "Penaeus monodon"))  %>%
#   # Add indices for each sci-name
#   mutate(clean_sci_name = as.factor(clean_sci_name),
#          sci = as.numeric(clean_sci_name))

# Now that alpha and theta are vectorized, can include all groups
lca_groups <- lca_dat_no_na %>%
  # Add indices for each sci-name
  mutate(clean_sci_name = as.factor(clean_sci_name),
         sci = as.numeric(clean_sci_name))

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

stan_data = list(n = n,
                 k = k,
                 feed_weights = feed_weights,
                 n_sci = n_sci,
                 sci = sci,
                 phi = phi,
                 kappa = kappa)

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

# LEFT OFF HERE - posteriors for the alphas all point-estimates because these are all deterministic, how to give it a proper prior DISTRIBUTION?
# NEW CODE - try to vectorize over alpha and theta:
stan_pooled <- 'data {
  int n;  // number of observations
  int k; // number of feed types
  int n_sci; // number of sci names
  simplex[k] feed_weights[n]; // array of observed feed weights simplexes
  int sci[n]; // sci-name indices
  simplex[k] phi[n_sci];
  int kappa[n_sci];
}
parameters {
  // alpha parameter now moved into transformed parameter section
  simplex[k] theta[n_sci]; // vectors of estimated sci-level feed weight simplexes;
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
fit_grouped <- sampling(object = no_missing_mod, data = stan_data, cores = 4)
#,cores = 4, iter = 10000,
#control = list(adapt_delta = 0.99)) # address divergent transitions by increasing delta, i.e., take smaller steps
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
                      area_method = "scaled height")
p_alpha
p_theta <- mcmc_areas(distribution_grouped,
                      pars = vars(contains("theta")),
                      prob = 0.8,
                      area_method = "scaled height")
p_theta

