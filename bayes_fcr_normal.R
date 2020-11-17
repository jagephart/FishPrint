# Bayesian estimation of FCR (with normal distribution):

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

lca_dat <- read.csv(file.path(datadir, "LCA_compiled_20201109.csv"), fileEncoding="UTF-8-BOM") #fileEncoding needed when reading in file from windows computer (suppresses BOM hidden characters)
source("Functions.R")

lca_dat_clean <- clean.lca(LCA_data = lca_dat)

# Remove bivalves from FCR analysis
lca_dat_clean <- lca_dat_clean %>%
  filter(taxa_group_name != "bivalves")
  
######################################################################################################
# Model 1: Remove NA's, and estimate group-level feed conversion ratio for Nile tilapia, Oreochromis niloticus (species with the most FCR data)
lca_dat_simple <- lca_dat_clean %>%
  filter(is.na(FCR) == FALSE) %>%
  filter(clean_sci_name == "Oncorhynchus mykiss")
  #filter(clean_sci_name == "Oreochromis niloticus")


# Each observation y is normally distributed with corresponding mean theta, and known variance, sigma2
# Each theta is drawn from a normal group-level distribution with mean mu and variance tau2

x <- lca_dat_simple$FCR
n <- nrow(lca_dat_simple)

# Normal distribution model:
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
#ggsave(file.path(outdir, "bayes-example_trout_fcr-gamma-target.png"), height = 8.5, width = 11)

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
#ggsave(file.path(outdir, "plot_boxplot_FCR-by-species.png"), height = 8, width = 11.5)


x <- lca_dat_groups$FCR
n <- nrow(lca_dat_groups)
j <- length(unique(lca_dat_groups$sp))
sp <- lca_dat_groups$sp

# Normal distribution model:
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
# Model 2a: Include NA's (but only for studies that have OTHER studies of the same taxa) and estimate sciname-level and allseafood-level mu

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
# Model 3: Include all NA's and include three hierarchical levels: sciname, taxagroup, and all-seafood

# Are there any scinames that don't have any other data to draw from at the taxagroup level
missing_dat <- lca_dat_clean %>%
  group_by(taxa_group_name) %>%
  mutate(n_study_with_data = sum(is.na(FCR)==FALSE)) %>%
  ungroup() %>%
  filter(is.na(FCR) & n_study_with_data == 0) %>%
  select(clean_sci_name, taxa_group_name) %>%
  unique()

# OLD CODE FOR TESTING taxonomic overlap of missing_data with rest of dataset:
# missing_dat
# A tibble: 4 x 2
# clean_sci_name sci_name_rank
# <chr>          <chr>        
# 1 Actinopterygii superclass   
# 2 Macrobrachium  genus        
# 3 Brachyura      infraorder   
# 4 Mytilus edulis species                      

# Test missing_dat and create a vector (drop_taxa) of sci_names that have no other studies from lower classification levels
# e.g., the missing FCR data for the study on genus Macrobrachium can come from the species-level studies of Macrobrachium spp.
# but the missing data for infraorder = Brachyura has no other studies of lower-level taxa to draw this info from
# drop_taxa <- NULL
# for (i in 1:nrow(missing_dat)){
#   taxa_na <- missing_dat[i,]$clean_sci_name
#   rank_search <- sym(missing_dat[i,]$sci_name_rank)
#   lca_dat_rank <- lca_dat_clean %>%
#     filter(!!rank_search == taxa_na) %>%
#     filter(is.na(FCR)==FALSE)
#   if (nrow(lca_dat_rank) == 0){
#     print(paste(taxa_na, " ", rank_search, " has no overlapping taxa with observed data", sep = ""))
#     drop_taxa <- append(drop_taxa, values = taxa_na)
#   }
# }

# Create sciname and taxagroup levels
lca_dat_groups_full <- lca_dat_clean %>%
  mutate(across(where(is.character), as.factor),
         sci_level = as.numeric(clean_sci_name),
         grp_level = as.numeric(taxa_group_name)) %>%
  droplevels() %>%
  arrange(clean_sci_name)

name_assignments_key <- lca_dat_groups_full %>%
  select(Scientific.Name, Common.Name, clean_sci_name, taxa_group_name) %>%
  unique()
#write.csv(name_assignments_key, file.path(outdir, "name_assignments_key.csv"), row.names = FALSE)

groupings_and_sample_size_key <- lca_dat_groups_full %>%
  group_by(clean_sci_name) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  select(clean_sci_name, taxa_group_name, n) %>%
  unique()
#write.csv(groupings_and_sample_size_key, file.path(outdir, "sample_size_key.csv"), row.names = FALSE)

# Data "x" that enter at the clean_sci_name level (individual studies, regardless of classification ranking)
x_obs_dat <- lca_dat_groups_full %>% filter(is.na(FCR) == FALSE) %>% select(clean_sci_name, FCR, sci_level, taxa_group_name, grp_level)
x_mis_dat <- lca_dat_groups_full %>% filter(is.na(FCR) == TRUE) %>% select(clean_sci_name, FCR, sci_level, taxa_group_name, grp_level) 
x_obs <- x_obs_dat$FCR
x_mis <- x_mis_dat$FCR
sci_obs <- x_obs_dat$sci_level
sci_mis <- x_mis_dat$sci_level
n_obs <- length(x_obs)
n_mis <- length(x_mis)
n_sci <- length(unique(lca_dat_groups_full$sci_level))

## Create a vector of group-level indicators where j-th element gives group ID for sci-name ID j
grp_obs_key <- unique(x_obs_dat[c("sci_level", "grp_level")])[,"grp_level"]
grp_mis_key <- unique(x_mis_dat[c("sci_level", "grp_level")])[,"grp_level"]
n_grp_obs <- length(grp_obs_key)
n_grp_mis <- length(grp_mis_key)
n_grp <- length(unique(lca_dat_groups_full$grp_level))

#ORIGINAL grouped model (RUNS, but still trying to get to converge)
stan_grouped <- 'data {
  int<lower=0> n_obs;  // number of observations
  int<lower=0> n_mis;  // number of missing observations
  vector[n_obs] x_obs; // sciname-level data
  int sci_obs[n_obs]; // sciname indicators for observed data
  int sci_mis[n_mis]; // sciname indicators for missing data
  int n_sci; // number of total scinames
  int<lower=0> n_grp_obs; // number of groups in observed data
  int<lower=0> n_grp_mis; // number of groups in missing data
  int grp_obs_key[n_grp_obs]; // group indicators for observed data
  int grp_mis_key[n_grp_mis]; // group indicators for missing data
  int n_grp; // number of total group names
}
parameters {
  real<lower=0> mu;
  real<lower=0> sigma;
  vector[n_grp] grp_mu;
  real<lower=0> grp_sigma;
  vector[n_sci] sci_mu;
  real<lower=0> sci_sigma;
  real x_mis[n_mis]; // missing data are treated as parameters

}

model {
  // priors
  // sigma ~ cauchy(0, 5);
  // note: because priors are optional, not sure whether I should bother giving sigma a cauchy distribution

  // likelihood
  grp_mu ~ normal(mu, sigma);
  for (i in 1:n_grp_obs){
    sci_mu ~ normal(grp_mu[grp_obs_key[i]], grp_sigma);
  }
  for (j in 1:n_grp_mis){
    sci_mu ~ normal(grp_mu[grp_mis_key[j]], grp_sigma);
  }
  x_obs ~ normal(sci_mu[sci_obs], sci_sigma);
  x_mis ~ normal(sci_mu[sci_mis], sci_sigma);

}'

# Compile model
grouped_mod <- stan_model(model_code = stan_grouped, verbose = TRUE)
# Note: For Windows, apparently OK to ignore this warning message:
# Warning message:
#   In system(paste(CXX, ARGS), ignore.stdout = TRUE, ignore.stderr = TRUE) :
#   'C:/rtools40/usr/mingw_/bin/g++' not found
# See forum discussion: https://discourse.mc-stan.org/t/rstan-on-windows/16673/60

# Fit model:
fit_grouped <- sampling(object = grouped_mod, data = list(x_obs = x_obs,
                                                          x_mis = x_mis,
                                                          n_obs = n_obs,
                                                          n_mis = n_mis,
                                                          n_sci = n_sci,
                                                          sci_obs = sci_obs,
                                                          sci_mis = sci_mis,
                                                          grp_obs_key = grp_obs_key,
                                                          grp_mis_key = grp_mis_key,
                                                          n_grp_obs = n_grp_obs,
                                                          n_grp_mis = n_grp_mis,
                                                          n_grp = n_grp),
                        iter = 50000, warmup = 1000, chain = 3, cores = 1)
# IF STILL CRASHING, TRY RUNNING ON ONE CORE
print(fit_grouped)

# How to interpret sci index numbers:
sci_index_key <- lca_dat_groups_full %>%
  select(clean_sci_name, sci_level) %>%
  arrange(sci_level) %>%
  unique() %>%
  mutate(sci_param_name = paste("sci_mu[", clean_sci_name, "]", sep = ""))
#write.csv(sp_index_key, file.path(outdir, "sp_info.csv"), row.names = FALSE)

# How to interpret grp index numbers:
grp_index_key <- lca_dat_groups_full %>%
  select(taxa_group_name, grp_level) %>%
  arrange(grp_level) %>%
  unique() %>%
  mutate(grp_param_name = paste("grp_mu[", taxa_group_name, "]", sep = ""))

# How to interpret x_mis index numbers:
x_mis_key <- lca_dat_groups_full %>%
  filter(is.na(FCR)) %>%
  group_by(clean_sci_name) %>%
  mutate(sci_name_index = row_number()) %>%
  mutate(x_mis_param_name = paste("x_mis[", clean_sci_name, " ", sci_name_index, " ", taxa_group_name, "]", sep = "")) %>%
  select(clean_sci_name, sci_name_index, x_mis_param_name)

# Replace param names; first copy to fit_grouped_clean to avoid having to re-run sampling as a result of doing something wrong to fit_grouped
fit_grouped_clean <- fit_grouped
names(fit_grouped_clean)[grep(names(fit_grouped_clean), pattern = "sci_mu")] <- sci_index_key$sci_param_name
names(fit_grouped_clean)[grep(names(fit_grouped_clean), pattern = "grp_mu")] <- grp_index_key$grp_param_name
names(fit_grouped_clean)[grep(names(fit_grouped_clean), pattern = "x_mis")] <- x_mis_key$x_mis_param_name


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
