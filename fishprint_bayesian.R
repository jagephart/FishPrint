# Bayesian estimation of Fish Foot Prints:

rm(list=ls())
library(tidyverse)
library(rstan)

datadir <- "/Volumes/jgephart/BFA Environment 2/Data"
lca_dat <- read.csv(file.path(datadir, "LCA_compiled.csv"))

# START OVER WITH THIS TUTORIAL: https://cran.r-project.org/web/packages/bridgesampling/vignettes/bridgesampling_example_stan.html

# Model 1: estimate group-level feed conversion ratio for Nile tilapia, Oreochromis niloticus (species with the most FCR data)
lca_dat_clean <- lca_dat %>%
  filter(is.na(FCR) == FALSE) %>%
  filter(Species.scientific.name=="Oreochromis niloticus")

# Each observation y is normally distributed with corresponding mean theta, and known variance, sigma2
# Each theta is drawn from a normal group-level distribution with mean mu and variance tau2

y <- lca_dat_clean$FCR
n <- nrow(lca_dat_clean)

# The theta's are then drawn from a group-level normal distribution with mean mu, and standard deviation, eta



# Each LAMBDA is then drawn from the group-level normal distribution with mean, mu and standard deviation, sigma
mu <- mean(lca_dat_clean$FCR)
sigma_x <- sd(lca_dat_clean$FCR)

stan_group_mean <- 'data {
int<lower=1> n; // number of observations
vector[n] x; // observations
}
parameters {

}
model {
target += normal_lpdf(x | mu, sigma_x)
}
'



# And the group level variance, tau is drawn from an inverse-gamma distribution (with shape params alpha, beta)
# Use the inverse-gamma (instead of just the gamma distribution) to model the variance because of conjugacy: 
# https://stats.stackexchange.com/questions/350924/why-do-we-use-inverse-gamma-as-prior-on-variance-when-empirical-variance-is-gam

# NEXT: STAN code for mean by groups: https://discourse.mc-stan.org/t/mean-by-groups/1300


