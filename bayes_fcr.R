# Bayesian estimation of Fish Foot Prints:

rm(list=ls())
library(tidyverse)
library(rstan)
library(taxize)
library(data.table)

datadir <- "/Volumes/jgephart/BFA Environment 2/Data"
outdir <- "/Volumes/jgephart/BFA Environment 2/Outputs"
lca_dat <- read.csv(file.path(datadir, "LCA_compiled.csv"))

# FIX IT - blanks for scientific names; add a taxanomic name that can be recognized by taxize
# FIX IT - remove bivalves from FCR analysis
# Clean lca_dat:
# Clean species names
# Remove unnecessary columns
lca_dat_clean <- lca_dat %>% 
  mutate(clean_sci_name = case_when(str_detect(Species.scientific.name, "spp") ~ str_replace(Species.scientific.name, pattern = " spp\\.| spp", replacement = ""),
                                    Species.scientific.name == "Morone chrysops x M. saxatilis" ~ "Morone",
                                    TRUE ~ Species.scientific.name)) %>%
  filter(Species.scientific.name != "") %>%
  mutate(FCR = case_when(str_detect(Species.scientific.name, "Thunnus") ~ FCR/5,
                         TRUE ~ FCR)) %>%
  select(-c(Date_entered, Description, Note_on_system, Notes, Person_entering, Product, Production_system, Sample_size, SeaWEED.ID, Source, Specific_location, Strain))

# SESYNC Bayesian course: https://cchecastaldo.github.io/BayesianShortCourse/Syllabus.html
# START OVER WITH THIS TUTORIAL: https://cran.r-project.org/web/packages/bridgesampling/vignettes/bridgesampling_example_stan.html

# Model 1: Remove NA's, and estimate group-level feed conversion ratio for Nile tilapia, Oreochromis niloticus (species with the most FCR data)
lca_dat_simple <- lca_dat_clean %>%
  filter(is.na(FCR) == FALSE) %>%
  filter(clean_sci_name == "Oreochromis niloticus")

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
  // sigma ~ cauchy(0, 5); // if we want to put a prior on sigma
  // notice: no prior on mu; any param with no prior is given a uniform

  // likelihood
  x ~ normal(mu, sigma); // note: stan interprets second param as standard deviation

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
# Model 1a: Calculate variance as a "transformed parameter"
stan_pooled_1a <- 'data {
  int<lower=0> n;  // number of observations
  vector[n] x; // data
}
parameters {
  real<lower=0> mu;
  real<lower=0> sigma2;
}
transformed parameters {
real<lower=0> sigma;
sigma = sqrt(sigma2);
}
model {
  // priors
  // sigma ~ cauchy(0, 5);

  // likelihood
  x ~ normal(mu, sigma);
}'
# Use the inverse-gamma (instead of just the gamma distribution) to model the variance because of conjugacy: 
# https://stats.stackexchange.com/questions/350924/why-do-we-use-inverse-gamma-as-prior-on-variance-when-empirical-variance-is-gam

# Fit model:
fit_pooled_1a <- stan(model_code = stan_pooled_1a, data = list(x = x,
                                                             n = n),
                     iter = 50000, warmup = 1000, chain = 3, cores = 3)
print(fit_pooled_1a)

################################################################################################################
# Model 2: Remove ALL NA's, and estimate group-level feed conversion ratio for all species
# NOTE: Stan does not support NA's in data - must be modeled explicitly
# Estimate group-level means
lca_dat_groups <- lca_dat_clean %>%
  filter(is.na(FCR) == FALSE) %>%
  mutate(clean_sci_name = as.factor(clean_sci_name),
         sp= as.numeric(clean_sci_name)) %>%
  select(clean_sci_name, FCR, sp)

# Remove "farms" that are NA AND have no other duplicate species
# lca_dat_groups <- lca_dat_clean %>%
#   group_by(clean_sci_name) %>%
#   mutate(n_farms = n()) %>%
#   ungroup() %>%
#   filter((is.na(FCR) & n_farms == 1)==FALSE) #%>%
#   mutate(clean_sci_name = as.factor(clean_sci_name),
#          sp = as.numeric(clean_sci_name)) %>%
#   select(clean_sci_name, FCR, sp)


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
  vector[n] x; // data
  int j; // number of species
  int sp[n]; // species indicators
}
parameters {
  real<lower=0> mu;
  real<lower=0> sp_sigma;
  vector[j] sp_mu;
  real<lower=0> sigma;
}

model {
  // priors
  // sigma ~ cauchy(0, 5);
  // note: because priors are optional, not sure whether I should bother giving sigma a cauchy distribution

  // likelihood
  sp_mu ~ normal(mu, sigma);
  x ~ normal(sp_mu[sp], sp_sigma);

}'

# Fit model:
fit_grouped <- stan(model_code = stan_grouped, data = list(x = x,
                                                           n = n,
                                                           j = j,
                                                           sp = sp))

print(fit_grouped)

# Diagnostics
stan_trace(fit_grouped) # prints first 10 parameters
#stan_trace(fit_pooled, pars = c('mu'))
#stan_trace(fit_pooled, pars = c('sigma'))


################################################################################################################
# Model 2a: Include NA's while estimating group-level means

# First, examine entries that are NA's and have no other con-specifics: how to group these?
lca_dat_clean %>%
  group_by(clean_sci_name) %>%
  mutate(n_farms = n()) %>%
  ungroup() %>%
  filter(is.na(FCR) & n_farms == 1) %>%
  select(clean_sci_name)

# STOP HERE: There are no single-entry missing FCR data at a species-level that we need to estimate using grouped hierarchies
# At some point, may need to build in hierarchies for this model so that we are able to output final foot print estimates for different levels (e.g., family, genus, etc), but wait until we have a discussion about this



