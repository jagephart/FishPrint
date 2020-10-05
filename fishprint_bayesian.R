# Bayesian estimation of Fish Foot Prints:

rm(list=ls())
library(tidyverse)
library(rstan)

datadir <- "/Volumes/jgephart/BFA Environment 2/Data"
outdir <- "/Volumes/jgephart/BFA Environment 2/Outputs"
lca_dat <- read.csv(file.path(datadir, "LCA_compiled.csv"))

# SESYNC Bayesian course: https://cchecastaldo.github.io/BayesianShortCourse/Syllabus.html
# START OVER WITH THIS TUTORIAL: https://cran.r-project.org/web/packages/bridgesampling/vignettes/bridgesampling_example_stan.html

# Model 1: Remove NA's, and estimate group-level feed conversion ratio for Nile tilapia, Oreochromis niloticus (species with the most FCR data)
lca_dat_simple <- lca_dat %>%
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
  // no need to specify prior on mu if flat

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
  sigma ~ cauchy(0, 5);

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
# Try group-levels for half-cauchy model:

# Model 2: Remove NA's, and estimate group-level feed conversion ratio for all species
lca_dat_groups <- lca_dat %>%
  filter(is.na(FCR) == FALSE) %>%
  filter(Species.scientific.name!="") %>%
  mutate(Species.scientific.name = as.factor(Species.scientific.name),
         grp = as.numeric(Species.scientific.name))

# FIX IT - blanks for scientific names

ggplot(data = lca_dat_groups, aes(x = Species.scientific.name, y = FCR)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 16)) + 
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
  vector[j] grp_mu
}

model {
  // priors
  
  sigma ~ cauchy(0, 5)

  // likelihood
  grp_mu ~ normal(mu, grp_sigma)
  x ~ normal(grp_mu[grp], sigma);

}'

# Fit model:
fit_grouped <- stan(model_code = stan_grouped, data = list(x = x,
                                                           n = n,
                                                           j = j,
                                                           grp = grp))


# When ready to calculate footprint, create section for "transformed parameters":
# https://stats.stackexchange.com/questions/375696/what-is-the-purpose-of-transformed-variables-in-stan

# Watch Gelman intro in Stan: https://www.youtube.com/watch?v=T1gYvX5c2sM&t=27s
# https://stats.stackexchange.com/questions/375696/what-is-the-purpose-of-transformed-variables-in-stan



# NEXT: STAN code for mean by groups: https://discourse.mc-stan.org/t/mean-by-groups/1300
# Also: https://discourse.mc-stan.org/t/simple-nested-normal-dist-model-help/5429


