# Bayesian estimation of Fish Foot Prints:

rm(list=ls())
library(tidyverse)
library(rstan)
library(taxsize)
library(data.table)

datadir <- "/Volumes/jgephart/BFA Environment 2/Data"
outdir <- "/Volumes/jgephart/BFA Environment 2/Outputs"
lca_dat <- read.csv(file.path(datadir, "LCA_compiled.csv"))

# Clean lca_dat for species matching:
lca_dat_clean <- lca_dat %>% 
  mutate(clean_sci_name = case_when(str_detect(Species.scientific.name, "spp") ~ str_replace(Species.scientific.name, pattern = " spp\\.| spp", replacement = ""),
                                    TRUE ~ Species.scientific.name)) %>%
  select(Species.scientific.name, clean_sci_name)

# SESYNC Bayesian course: https://cchecastaldo.github.io/BayesianShortCourse/Syllabus.html
# START OVER WITH THIS TUTORIAL: https://cran.r-project.org/web/packages/bridgesampling/vignettes/bridgesampling_example_stan.html

# Model 1: Remove NA's, and estimate group-level feed conversion ratio for Nile tilapia, Oreochromis niloticus (species with the most FCR data)
lca_dat_simple <- lca_dat_clean %>%
  filter(is.na(FCR) == FALSE) %>%
  filter(Species.scientific.name=="Oreochromis niloticus")

# Each observation y is normally distributed with corresponding mean theta, and known variance, sigma2
# Each theta is drawn from a normal group-level distribution with mean mu and variance tau2

x <- lca_dat_simple$FCR
n <- nrow(lca_dat_simple)

# Use a half-cauchy distribution (weakly informative) prior on sigma:
stan_pooled <- 'data {
  int<lower=0> n;  // number of observations
  vector[n] x; // data
}
parameters {
  real<lower=0> mu;
  real<lower=0> sigma;
}
model {
  // priors
  sigma ~ cauchy(0, 5);
  // notice: no prior on mu; any param with no prior is given a uniform

  // likelihood
  x ~ normal(mu, sigma);

}'

# Fit model:
fit_pooled <- stan(model_code = stan_pooled, data = list(x = x, n = n),
                   iter = 10000, warmup = 1000, chain = 3, cores = 3)
# default: chains = 4, iter = 2000, warmup = floor(iter/2)

print(fit_pooled)
# Note: lp__ is the sum of the vector of log probabilities (but after removing any constant scale factors, making it not useful for model comparison)
# https://www.jax.org/news-and-insights/jax-blog/2015/october/lp-in-stan-output#:~:text=Therefore%2C%20%E2%80%9Clp%E2%80%9D%20is%20actually,useful%20for%20model%20comparison%20purposes.

# Diagnostics
#stan_trace(fit_pooled)
stan_trace(fit_pooled, pars = c('mu'))
stan_trace(fit_pooled, pars = c('sigma'))

################################################################################################################
# Model sigma with an inv_gamma prior:
stan_pooled_2 <- 'data {
  int<lower=0> n;  // number of observations
  vector[n] x; // data
}
parameters {
  real<lower=0> mu;
  real<lower=0> sigma2;
  real<lower=0> alpha;
  real<lower=0> beta;
}
transformed parameters {
real<lower=0> sigma;
sigma = sqrt(sigma2);
}
model {
  // priors
  sigma ~ inv_gamma(alpha, beta);

  // likelihood
  x ~ normal(mu, sigma);
}'
# Use the inverse-gamma (instead of just the gamma distribution) to model the variance because of conjugacy: 
# https://stats.stackexchange.com/questions/350924/why-do-we-use-inverse-gamma-as-prior-on-variance-when-empirical-variance-is-gam

# Fit model:
fit_pooled_2 <- stan(model_code = stan_pooled_2, data = list(x = x,
                                                         n = n),
                     iter = 50000, warmup = 1000, chain = 3, cores = 3)
# Can't get this to converge

################################################################################################################
# Calculate group-level means for all groups using half-cauchy model:

# Model 2: Remove NA's, and estimate group-level feed conversion ratio for all species
lca_dat_groups <- lca_dat_clean %>%
  filter(is.na(FCR) == FALSE) %>%
  filter(Species.scientific.name != "") %>%
  mutate(Species.scientific.name = as.factor(Species.scientific.name),
         grp = as.numeric(Species.scientific.name)) %>%
  select(Species.scientific.name, FCR, grp)

# FIX IT - blanks for scientific names

ggplot(data = lca_dat_groups, aes(x = Species.scientific.name, y = FCR)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 16)) + 
  labs(title = "Boxplots of species-level FCR")
ggsave(file.path(outdir, "plot_boxplot_FCR-by-species.png"), height = 8, width = 11.5)


x <- lca_dat_groups$FCR
n <- nrow(lca_dat_groups)
j <- length(unique(lca_dat_groups$grp))
grp <- lca_dat_groups$grp


stan_grouped <- 'data {
  int<lower=0> n;  // number of observations
  vector[n] x; // data
  int j; // number of groups
  int grp[n]; // group indicators
}
parameters {
  real<lower=0> mu;
  real<lower=0> sigma;
  vector[j] grp_mu;
  real<lower=0> grp_sigma;
}

model {
  // priors
  grp_sigma ~ cauchy(0, 5);
  // note: because priors are optional, not sure whether I should bother giving grp_sigma a cauchy distribution

  // likelihood
  grp_mu ~ normal(mu, grp_sigma);
  x ~ normal(grp_mu[grp], sigma);

}'

# Fit model:
fit_grouped <- stan(model_code = stan_grouped, data = list(x = x,
                                                           n = n,
                                                           j = j,
                                                           grp = grp))

print(fit_grouped)

# Diagnostics
stan_trace(fit_grouped) # prints first 10 parameters
#stan_trace(fit_pooled, pars = c('mu'))
#stan_trace(fit_pooled, pars = c('sigma'))

################################################################################################################
# Include NA's while estimating group-level means using half-cauchy model:

# FIX IT - missing scientific names
# Model 2: Remove NA's, and estimate group-level feed conversion ratio for all species

# Use package taxize to get higher classification levels for each species

# Use WORMS database, outputs are more simpler, fewer types of ranks
classify_ncbi <- classification("Penaeus monodon", db = "ncbi")
classify_ncbi[[1]]
classify_worms <- classification("Penaeus monodon", db = "worms")
classify_worms[[1]]

## LEFT OFF HERE - make sure code runs below using new clean_sci_name column (removed "spp")
# Get full species list
species_list <- unique(lca_dat_clean$clean_sci_name)
# TEST: species_list <- c("Penaeus monodon", "Oreochromis niloticus")

# Use classification function to get higher ranks for all species
classify_worms <- lapply(species_list, classification, db = "worms")

classify_test <- data.frame(no_results = unlist(lapply(classify_worms, is.na))) %>%
  filter(no_results == TRUE)

# Function to reformat ranks into columns
format_worms <- function(classify_worms) {
  higher_ranks <- classify_worms[[1]] %>%
    filter(rank %in% c("Class", "Subclass", "Superorder", "Order", "Suborder", "Superfamily", "Family", "Subfamily", "Genus", "Species")) %>%
    select(name, rank) %>%
    pivot_wider(names_from = rank, values_from = name)
}

# Reformat all classification function outputs into columns and bind into single data table
worms_cols <- lapply(classify_worms, format_worms)
worms_cols_dt <- data.table::rbindlist(worms_cols, use.names = TRUE, fill = TRUE)

# Join with lca_dat_clean
lca_with_ranks <- lca_dat_clean %>%
  left_join(worms_cols_dt, by = c("clean_sci_name" = "Species"))

# Find another grouping variable besides species.scientific.name
lca_dat_clean %>%
  filter(Species.scientific.name != "") %>%
  mutate(Species.scientific.name = as.factor(Species.scientific.name),
         grp = as.numeric(Species.scientific.name)) %>%
  select(Species.common.name, Source, Species.scientific.name, grp, Strain, Product) %>%
  group_by(Species.scientific.name) %>%
  mutate(n = n()) %>%
  arrange(Species.common.name)

