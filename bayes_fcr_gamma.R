# Bayesian estimation of FCR (with gamma distribution):

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

# Remove bivalves from FCR analysis
lca_dat_clean <- lca_dat_clean %>%
  filter(taxa_group_name != "bivalves")

######################################################################################################
# Model 1: Remove NA's, and estimate group-level feed conversion ratio for Nile tilapia, Oreochromis niloticus (species with the most FCR data)

# Gamma distribution model
stan_pooled <- 'data {
  int<lower=0> n;  // number of observations
  vector[n] x; // data
}
parameters {
  real<lower=0> mu;
  real<lower=0> sigma;
}
model {
  // reparamaterize gamma based on mu and sigma; sigma is standard deviation
  real shape = square(mu) / square(sigma);
  real rate = mu / square(sigma);
  
  x ~ gamma(shape, rate);
  // target += gamma_lpdf(x | shape, rate); // alternative notation

}'


# Fit model:
fit_pooled <- stan(model_code = stan_pooled, data = list(x = x, n = n))
#iter = 10000, warmup = 1000, chain = 3, cores = 3)
# default: chains = 4, iter = 2000, warmup = floor(iter/2)

print(fit_pooled)
# Note: lp__ is the sum of the vector of log probabilities (but after removing any constant scale factors, making it not useful for model comparison)
# https://www.jax.org/news-and-insights/jax-blog/2015/october/lp-in-stan-output#:~:text=Therefore%2C%20%E2%80%9Clp%E2%80%9D%20is%20actually,useful%20for%20model%20comparison%20purposes.

# Diagnostics
stan_trace(fit_pooled)
#stan_trace(fit_pooled, pars = c('mu'))
#stan_trace(fit_pooled, pars = c('sigma'))

distribution_grouped <- as.matrix(fit_pooled)

plot_theme <- theme(axis.text=element_text(size=14, color = "black"))

p <- mcmc_areas_ridges(distribution_grouped,
                       pars = vars(contains(c("mu", "sigma"))),
                       prob = 0.8) +
  ggtitle("Oncorhynchus mykiss FCR model", "with 80% credible intervals") +
  plot_theme

p 
ggsave(file.path(outdir, "bayes-example_trout_fcr-gamma-target.png"), height = 8.5, width = 11)

######################################################################################################
# Model 2: Remove ALL NA's, and estimate group-level feed conversion ratio for all species and a global feed conversion ratio (mu)

# NOTE: Stan does not support NA's in data - must be modeled explicitly
# Estimate group-level means
lca_dat_groups <- lca_dat_clean %>%
  filter(is.na(FCR) == FALSE) %>%
  mutate(clean_sci_name = as.factor(clean_sci_name),
         sp= as.numeric(clean_sci_name)) %>%
  select(clean_sci_name, FCR, sp)

# Boxplot of data
ggplot(data = lca_dat_groups, aes(x = clean_sci_name, y = FCR)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 16)) + 
  labs(title = "Boxplots of species-level FCR")
#ggsave(file.path(outdir, "plot_boxplot_FCR-by-species.png"), height = 8, width = 11.5)

x <- lca_dat_groups$FCR
n <- nrow(lca_dat_groups)
j <- length(unique(lca_dat_groups$sp))
sp <- lca_dat_groups$sp

stan_grouped <- 'data {
  int<lower=0> n;  // number of observations
  vector<lower=0>[n] x; // data
  int j; // number of species
  int sp[n]; // species indicators
}
parameters {
  real<lower=0> mu;
  real<lower=0> sigma;
  vector<lower=0>[j] sp_mu;
  real<lower=0> sp_sigma;
}
transformed parameters {
  // reparamaterize sci-name level gamma to get sci-name mu and sigma
  vector[j] sp_shape = square(sp_mu) ./ square(sp_sigma);
  vector[j] sp_rate = sp_mu ./ square(sp_sigma);
  // reparamaterize global gamma to get global mu and sigma
  real shape = square(mu) / square(sigma);
  real rate = mu / square(sigma);
}
model {
  // Put priors on shape and rate params, otherwise will sample negative
  // shape ~ uniform(0,100);
  // rate ~ uniform(0,100);
  // sp_shape[sp] ~ uniform(0, 100); // shape and rate are also always positive 
  // sp_rate[sp] ~ uniform(0, 100);
  
  // Put priors on mu and sigma since this is more intuitive:
  mu ~ uniform(0, 10);
  sigma ~ uniform(0, 100);
  sp_mu[sp] ~ uniform(0, 10);
  sp_sigma ~ uniform(0, 100);
  
  // likelihood
  // target += gamma_lpdf(sp_mu | shape, rate); 
  sp_mu ~ gamma(shape, rate);
  //target += gamma_lpdf(x | sp_shape[sp], sp_rate[sp]); 
  x ~ gamma(sp_shape[sp], sp_rate[sp]); // equivalent notation to target += but faster
  
}'

# Compile
no_missing_mod <- stan_model(model_code = stan_grouped, verbose = TRUE)
# Note: For Windows, apparently OK to ignore this warning message:
# Warning message:
#   In system(paste(CXX, ARGS), ignore.stdout = TRUE, ignore.stderr = TRUE) :
#   'C:/rtools40/usr/mingw_/bin/g++' not found

# Fit model:
fit_grouped <- sampling(object = no_missing_mod, data = list(x = x,
                                                             n = n,
                                                             j = j,
                                                             sp = sp),
                        iter = 10000, warmup = 1000, chain = 4, cores = 4)

print(fit_grouped)

# How to interpret sp index numbers:
# What are the sample sizes per group:
sp_index_key <- lca_dat_groups %>%
  group_by(clean_sci_name) %>%
  mutate(n_sci_name = n()) %>%
  ungroup() %>%
  select(sp, clean_sci_name, n_sci_name) %>%
  arrange(sp) %>%
  unique() %>%
  mutate(mu_param_name = paste("sp_mu[", clean_sci_name, "]", sep =""),
         shape_param_name = paste("sp_shape[", clean_sci_name, "]", sep =""),
         rate_param_name = paste("sp_rate[", clean_sci_name, "]", sep = ""))
#write.csv(sp_index_key, file.path(outdir, "sp_info.csv"), row.names = FALSE)

# Replace param names
names(fit_grouped)[grep(names(fit_grouped), pattern = "sp_mu")] <- sp_index_key$mu_param_name
names(fit_grouped)[grep(names(fit_grouped), pattern = "sp_shape")] <- sp_index_key$shape_param_name
names(fit_grouped)[grep(names(fit_grouped), pattern = "sp_rate")] <- sp_index_key$rate_param_name

# Diagnostics
#stan_trace(fit_grouped) # prints first 10 parameters
#stan_trace(fit_pooled, pars = c('mu'))
#stan_trace(fit_pooled, pars = c('sigma'))

# Make example MCMC plots of mu and sigma
# posterior_grouped <- rstan::extract(fit_grouped, inc_warmup = TRUE, permuted = FALSE)
# color_scheme_set("mix-blue-pink")
# p <- mcmc_trace(posterior_grouped,  pars = c("mu", "sigma"), n_warmup = 1000,
#                 facet_args = list(nrow = 2, labeller = label_parsed))
# p + facet_text(size = 15)

#ggsave(file.path(outdir, "plot_mcmc-plot_FCR-by-species.png"), height = 8, width = 11.5)

# Make plots of Posterior distributions with 80% credible intervals
distribution_grouped <- as.matrix(fit_grouped)
p2 <- mcmc_areas(distribution_grouped,
                 pars = vars(starts_with("mu")|contains("sp_mu")),
                 prob = 0.8,
                 prob_outer = 0.9,
                 area_method = "scaled height")

p2 + ggtitle("Posterior distributions", "with 80% credible intervals")
#ggsave(file.path(outdir, "plot_post-distribution-plot_FCR-by-species.png"), height = 8, width = 11.5)

p2 <- mcmc_areas_ridges(distribution_grouped,
                        pars = vars(contains("sp_shape")),
                        prob = 0.8,
                        prob_outer = 0.9)
p2 + ggtitle("Posterior distributions", "with 80% credible intervals")

p2 <- mcmc_areas_ridges(distribution_grouped,
                        pars = vars(contains("sp_rate")),
                        prob = 0.8,
                        prob_outer = 0.9)

p2 + ggtitle("Posterior distributions", "with 80% credible intervals")


######################################################################################################
# LEFT OFF HERE:
# Model 3: Gamma distribution with 3 levels (sci-name, taxa-group, all-seafood): start with just sci-name and taxa-group for now