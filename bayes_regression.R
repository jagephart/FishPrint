# Bayesian regression to estimate variables based on taxa-group, intensity, and production system

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

lca_dat_clean <- clean.lca(LCA_data = lca_dat)

######################################################################################################
# Model 1: Remove ALL NA's (predictors and response), and estimate feed conversion ratio for all scientific names, non-hierarchical

# If desired, replicate data based on clean_sample_size column:
# lca_dat_clean <- rep_data(lca_dat_clean)
# IMPORTANT for convergence: filter out FCR == 0

# Remove NAs (which also removes bivalves) and 0's
lca_dat_groups <- lca_dat_clean %>%
  filter(is.na(FCR) == FALSE) %>% 
  filter(FCR != 0) %>%
  mutate(taxa = as.factor(taxa_group_name),
         tx = as.numeric(taxa),
         intensity = as.factor(Intensity),
         its = as.numeric(intensity),
         system = as.factor(Production_system_group),
         ps = as.numeric(system)) %>%
  select(clean_sci_name, FCR, taxa, tx, intensity, its, system, ps) %>%
  drop_na() # This last drop NA filters any tx, its, or ps NAs

# Set data:
N = nrow(lca_dat_groups)
y = lca_dat_groups$FCR
X <- model.matrix(object = ~taxa + intensity + system, 
                  data = lca_dat_groups %>% select(taxa, intensity, system))
K = ncol(X) # number of predictors

stan_data <- list(N = N, y = y, K = K, X = X)

# FIX IT - is canonical inverse link more appropriate for categorical variables?

# FIX IT - reparamterize on shape and mean?
# I found it is better to use
# alpha = shape;
# beta = shape/mu;
# You can keep your definition of mu the same, but instead of using phi as a precision term, we specify a shape parameter. I believe this parameterization is what rstanarm uses.

stan_one_level <- 'data {
  int N; // number of observations
  real y[N]; // response
  int K; // number of predictors, i.e., number of columns in the design matrix X
  matrix [N, K] X; // design matrix X
}
parameters {
  vector[K] beta; // regression coefficients
  real sigma; // standard deviation
}
transformed parameters {
  vector[N] mu; // expected values
  vector[N] shape; // shape parameter for gamma
  vector[N] rate; // rate parameter for gamma 
  
  // reparameterize shape and rate
  mu = exp(X * beta); // log-link since shape and rate params are positive
  shape = mu .* mu / square(sigma);
  rate = mu / square(sigma);
}
model {
  // priors
  // beta[1] ~ cauchy(0,10); // prior for intercept following Gelman 2008 
  
  //for (k in 2:K){
  //  beta[K] ~ cauchy(0, 2.5); // prior for the slopes
  //}
  
  // likelihood
  y ~ gamma(shape, rate);
}
generated quantities {
  // simulate data by drawing from the posterior
  vector[N] y_rep;
  for (n in 1:N){
  y_rep[n] = gamma_rng(shape[n], rate[n]);
  }
}'

# Compile
simple_mod <- stan_model(model_code = stan_one_level, verbose = TRUE)
# Note: For Windows, apparently OK to ignore this warning message:
# Warning message:
#   In system(paste(CXX, ARGS), ignore.stdout = TRUE, ignore.stderr = TRUE) :
#   'C:/rtools40/usr/mingw_/bin/g++' not found

# Fit model:
fit_simple <- sampling(object = simple_mod, data = stan_data)

print(fit_simple)

# Response variable scinames
sci_index_key <- lca_dat_groups %>%
  select(clean_sci_name, intensity, system) %>%
  arrange(clean_sci_name) %>%
  mutate(sci_param_name = paste("sp_mu[", clean_sci_name, "]", sep =""),
         shape_param_name = paste("sp_shape[", clean_sci_name, "]", sep =""),
         rate_param_name = paste("sp_rate[", clean_sci_name, "]", sep = ""))
#write.csv(sp_index_key, file.path(outdir, "sp_info.csv"), row.names = FALSE)

# LEFT OFF HERE
# Clean up beta parameter names
beta_index_key <- colnames(X) %>% str_remove(pattern = "taxa|intensity|system|\\(|\\)")
  


# Replace param names
names(fit_grouped)[grep(names(fit_grouped), pattern = "sp_mu")] <- sp_index_key$mu_param_name
names(fit_grouped)[grep(names(fit_grouped), pattern = "sp_shape")] <- sp_index_key$shape_param_name
names(fit_grouped)[grep(names(fit_grouped), pattern = "sp_rate")] <- sp_index_key$rate_param_name