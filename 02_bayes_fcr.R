# Bayesian estimation of FCR:

rm(list=ls())
library(tidyverse)
library(rstan)
library(taxize)
library(data.table)
library(bayesplot) # for mcmc_areas_ridges

# Mac
datadir <- "/Volumes/jgephart/BFA Environment 2/Data"
outdir <- "/Volumes/jgephart/BFA Environment 2/Outputs"
# Windows
# datadir <- "K:/BFA Environment 2/Data"
# outdir <- "K:BFA Environment 2/Outputs"
lca_dat_clean <- read.csv(file.path(datadir, "lca_clean_with_ranks.csv"))

# FIX IT - remove bivalves from FCR analysis

######################################################################################################
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
stan_trace(fit_pooled)
#stan_trace(fit_pooled, pars = c('mu'))
#stan_trace(fit_pooled, pars = c('sigma'))

######################################################################################################
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

######################################################################################################
# Model 2: Remove ALL NA's, and estimate group-level feed conversion ratio for all species and a global feed conversion ratio (mu)
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
ggsave(file.path(outdir, "plot_boxplot_FCR-by-species.png"), height = 8, width = 11.5)


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
                                                           sp = sp),
                    iter = 50000, warmup = 1000, chain = 3, cores = 3)

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
  mutate(param_name = paste("sp_mu[", clean_sci_name, "]", sep =""))
write.csv(sp_index_key, file.path(outdir, "sp_info.csv"), row.names = FALSE)

# Replace param names
names(fit_grouped)[grep(names(fit_grouped), pattern = "sp_mu")] <- sp_index_key$param_name

# Diagnostics
stan_trace(fit_grouped) # prints first 10 parameters
#stan_trace(fit_pooled, pars = c('mu'))
#stan_trace(fit_pooled, pars = c('sigma'))

# Make example MCMC plots of mu and sigma
posterior_grouped <- rstan::extract(fit_grouped, inc_warmup = TRUE, permuted = FALSE)
color_scheme_set("mix-blue-pink")
p <- mcmc_trace(posterior_grouped,  pars = c("mu", "sigma"), n_warmup = 1000,
                facet_args = list(nrow = 2, labeller = label_parsed))
p + facet_text(size = 15)

ggsave(file.path(outdir, "plot_mcmc-plot_FCR-by-species.png"), height = 8, width = 11.5)

# Make plots of Posterior distributions with 80% credible intervals
distribution_grouped <- as.matrix(fit_grouped)
p2 <- mcmc_areas_ridges(distribution_grouped,
           pars = vars(contains("mu")),
           prob = 0.8)

p2 + ggtitle("Posterior distributions", "with 80% credible intervals")
ggsave(file.path(outdir, "plot_post-distribution-plot_FCR-by-species.png"), height = 8, width = 11.5)


######################################################################################################
# Model 2a: Include NA's (but only for studies that have OTHER studies of the same taxa) and estimate group-level and global-level mu

# First, remove studies that have missing FCR data and no other studies of the same taxa
lca_dat_with_missing <- lca_dat_clean %>%
  group_by(clean_sci_name) %>%
  mutate(n_study_with_data = sum(is.na(FCR)==FALSE)) %>% # how many studies per species that have data
  ungroup() %>%
  filter(is.na(FCR)==FALSE | is.na(FCR) & n_study_with_data > 0) %>% # only keep studies if they have FCR data OR if FCR == NA and there is at least one study with data
  mutate(clean_sci_name = as.factor(clean_sci_name),
         sp = as.numeric(clean_sci_name)) %>%
  select(clean_sci_name, FCR, sp, n_study_with_data) %>%
  arrange(clean_sci_name)


# Warning message because now there are NA's in the data frame
ggplot(data = lca_dat_with_missing, aes(x = clean_sci_name, y = FCR)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 16)) + 
  labs(title = "Boxplots of species-level FCR")
#ggsave(file.path(outdir, "plot_boxplot_FCR-by-species.png"), height = 8, width = 11.5)

# Separate out the missing data from the observed data:
lca_dat_observed <- lca_dat_with_missing %>%
  filter(is.na(FCR)==FALSE)

lca_dat_na <- lca_dat_with_missing %>%
  filter(is.na(FCR))

x_obs <- lca_dat_observed$FCR
x_mis <- lca_dat_na$FCR
n_obs <- length(x_obs)
n_mis <- length(x_mis)
j <- length(unique(lca_dat_observed$sp))
sp_obs <- lca_dat_observed$sp
sp_mis <- lca_dat_na$sp

stan_grouped <- 'data {
  int<lower=0> n_obs;  // number of observations
  int<lower=0> n_mis;  // number of missing observations
  vector[n_obs] x_obs; // data
  int j; // number of species
  int sp_obs[n_obs]; // species indicators for observed data
  int sp_mis[n_mis]; // species indicators for NAs
}
parameters {
  real<lower=0> mu;
  real<lower=0> sp_sigma;
  vector[j] sp_mu;
  real<lower=0> sigma;
  real x_mis[n_mis]; // missing data are treated as parameters
}

model {
  // priors
  // sigma ~ cauchy(0, 5);
  // note: because priors are optional, not sure whether I should bother giving sigma a cauchy distribution

  // likelihood
  sp_mu ~ normal(mu, sigma);
  x_obs ~ normal(sp_mu[sp_obs], sp_sigma);
  x_mis ~ normal(sp_mu[sp_mis], sp_sigma);

}'

grouped_mod <- stan_model(model_code = stan_grouped, verbose = TRUE)
# Note: For Windows, apparently OK to ignore this warning message:
# Warning message:
#   In system(paste(CXX, ARGS), ignore.stdout = TRUE, ignore.stderr = TRUE) :
#   'C:/rtools40/usr/mingw_/bin/g++' not found

# Fit model:
  fit_grouped <- sampling(object = grouped_mod, data = list(x_obs = x_obs,
                                                           x_mis = x_mis,
                                                           n_obs = n_obs,
                                                           n_mis = n_mis,
                                                           j = j,
                                                           sp_obs = sp_obs,
                                                           sp_mis = sp_mis))

print(fit_grouped)


# How to interpret sp index numbers:
# What are the sample sizes per group:
sp_index_key <- lca_dat_with_missing %>%
  group_by(clean_sci_name) %>%
  mutate(n_sci_name = n()) %>%
  ungroup() %>%
  select(sp, clean_sci_name, n_sci_name) %>%
  arrange(sp) %>%
  unique() %>%
  mutate(param_name = paste("sp_mu[", clean_sci_name, "]", sep =""))
#write.csv(sp_index_key, file.path(outdir, "sp_info.csv"), row.names = FALSE)

# How to interpret x_mis index numbers:
x_mis_key <- lca_dat_with_missing %>%
  filter(is.na(FCR)) %>%
  group_by(clean_sci_name) %>%
  mutate(sci_name_index = row_number()) %>%
  mutate(param_name = paste("x_mis[", clean_sci_name, " ", sci_name_index, "]", sep = ""))

# Replace param names; first copy to fit_grouped_clean so as not to overwrite original sampling output
fit_grouped_clean <- fit_grouped
names(fit_grouped_clean)[grep(names(fit_grouped_clean), pattern = "sp_mu")] <- sp_index_key$param_name
names(fit_grouped_clean)[grep(names(fit_grouped_clean), pattern = "x_mis")] <- x_mis_key$param_name


# Make plots of Posterior distributions with 80% credible intervals
distribution_grouped <- as.matrix(fit_grouped_clean)
p_mu <- mcmc_areas_ridges(distribution_grouped,
                        pars = vars(contains("mu")),
                        prob = 0.8)
p_mu + ggtitle("Posterior distributions", "with 80% credible intervals")

p_mis <- mcmc_areas_ridges(distribution_grouped,
                       pars = vars(contains("x_mis")),
                       prob = 0.8)
p_mis + ggtitle("Posterior distributions", "with 80% credible intervals")

p_sigma <- mcmc_areas_ridges(distribution_grouped,
                             pars = vars(contains("sigma")),
                             prob = 0.8)
p_sigma + ggtitle("Posterior distributions", "with 80% credible intervals")


######################################################################################################
# Model 3: Include NA's and estimate FCR for studies from lower levels - i.e., model more than just species-level and global-level

# Which group-levels should be modeled?
# First, get the sci_name of studies that have missing FCR data and no other studies of the same taxa
missing_dat <- lca_dat_clean %>%
  group_by(clean_sci_name) %>%
  mutate(n_study_with_data = sum(is.na(FCR)==FALSE)) %>%
  ungroup() %>%
  filter(is.na(FCR) & n_study_with_data == 0) %>%
  select(clean_sci_name, sci_name_rank, n_study_with_data) %>%
  unique()

# missing_dat
# A tibble: 3 x 1
# clean_sci_name
# <chr>         
# 1 Actinopterygii 
# 2 Macrobrachium 
# 3 Brachyura     
# 4 Mytilus edulis species                       
# 5 Mytilus edulis species                       
# 6 Mytilus edulis species                       

# Test missing_dat and create a vector (drop_taxa) of sci_names that have no other studies from lower classification levels
# e.g., the missing FCR data for the study on genus Macrobrachium can come from the species-level studies of Macrobrachium spp.
# but the missing data for infraorder = Brachyura has no other studies of lower-level taxa to draw this info from
drop_taxa <- NULL
for (i in 1:nrow(missing_dat)){
  taxa_na <- missing_dat[i,]$clean_sci_name
  rank_search <- sym(missing_dat[i,]$sci_name_rank)
  lca_dat_rank <- lca_dat_clean %>%
    filter(!!rank_search == taxa_na) %>%
    filter(is.na(FCR)==FALSE)
  if (nrow(lca_dat_rank) == 0){
    print(paste(taxa_na, " has no overlapping taxa with observed data", sep = ""))
    drop_taxa <- append(drop_taxa, values = taxa_na)
  }
}

# FIX IT - how to deal with Brachyura and Mytilus edulis
# options:
# remove them, but downstream analysis will have to use a higher-level (e.g., Suborder for infraorder Brachyura) when calculating posterior samples specific to Brachyura
# keep them in, see if the model still runs?

# 
# Create species level and other higher group levels
lca_dat_groups_full <- lca_dat_clean %>%
  filter(clean_sci_name %in% drop_taxa == FALSE) %>%
  mutate(across(where(is.character), as.factor),
         sp = as.numeric(species),
         gen = as.numeric(genus),
         fam = as.numeric(family),
         ord = as.numeric(order),
         sbord = as.numeric(suborder),
         spcl = as.numeric(superclass)) %>%
  droplevels()

# To keep code simple for now, just model species, genus (to get Macrobrachium), superclass (for Actinopterygii), and global (across all taxa)

# Data that enter at the "x" level (individual studies that are species)
x_obs <- lca_dat_groups_full %>% filter(is.na(sp) == FALSE & is.na(FCR) == FALSE & sci_name_rank == "species") %>% select(clean_sci_name, FCR, sp)
#x_mis <- lca_dat_groups_full %>% filter(is.na(sp) == FALSE & is.na(FCR) == TRUE & sci_name_rank == "species") %>% select(clean_sci_name, FCR) %>% unique() # FIX IT - leave these out? don't need to estimate these?

x_gen_obs <- lca_dat_groups_full %>% filter(is.na(gen) == FALSE & is.na(FCR) == FALSE & sci_name_rank == "genus") %>% select(clean_sci_name, FCR) %>% unique()
x_gen_mis <- lca_dat_groups_full %>% filter(is.na(gen) == FALSE & is.na(FCR) == TRUE & sci_name_rank == "genus") %>% select(clean_sci_name, FCR) %>% unique()

# FIX IT -leave these out since there's no observed data and we're estimating from lower levels:
#x_spcl_obs <- lca_dat_groups_full %>% filter(is.na(gen) == FALSE & is.na(FCR) == FALSE & sci_name_rank == "superclass") %>% select(clean_sci_name, FCR) %>% unique() 
#x_spcl_mis <- lca_dat_groups_full %>% filter(is.na(gen) == FALSE & is.na(FCR) == TRUE & sci_name_rank == "superclass") %>% select(clean_sci_name, FCR) %>% unique() 

n_x_obs <- length(x_obs)
n_gen_obs <- length(x_gen_obs)
n_gen_mis <- length(x_gen_mis)
j <- length(unique(lca_dat_observed$sp))
sp_obs <- lca_dat_observed$sp
sp_mis <- lca_dat_na$sp

stan_grouped <- 'data {
  int<lower=0> n_obs;  // number of observations
  int<lower=0> n_mis;  // number of missing observations
  vector[n_obs] x_obs; // data
  int j; // number of species
  int sp_obs[n_obs]; // species indicators for observed data
  int sp_mis[n_mis]; // species indicators for NAs
}
parameters {
  real<lower=0> mu;
  real<lower=0> sp_sigma;
  vector[j] sp_mu;
  real<lower=0> sigma;
  real x_mis[n_mis]; // missing data are treated as parameters
}

model {
  // priors
  // sigma ~ cauchy(0, 5);
  // note: because priors are optional, not sure whether I should bother giving sigma a cauchy distribution

  // likelihood
  spcl_mu[spcl_index] ~ normal(mu, sigma);
  gen_mu[gen_obs] ~ normal(spcl_mu[spcl_index], spcl_sigma);
  gen_mu[gen_mis] ~ normal(spcl_mu[spcl_index], spcl_sigma);
  sp_mu[sp_obs] ~ normal(gen_mu[gen_obs], gen_sigma);
  sp_mu[sp_obs] ~ normal(gen_mu[gen_mis], gen_sigma);
  x_obs ~ normal(sp_mu[sp_obs], sp_sigma);

}'