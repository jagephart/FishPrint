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

lca_dat <- read.csv(file.path(datadir, "LCA_compiled_20201109.csv"), fileEncoding="UTF-8-BOM") #fileEncoding needed when reading in file from windows computer (suppresses BOM hidden characters)
source("Functions.R")

lca_dat_clean <- clean.lca(LCA_data = lca_dat)

# Rebuild FAO fish production from zip file
fishstat_dat <- rebuild_fish("/Volumes/jgephart/FishStatR/Data/Production-Global/ZippedFiles/GlobalProduction_2019.1.0.zip")

# Classify species into taxa groupings
# Use ISSCAAP grouping (in FAO data) to help with classification
lca_dat_clean_groups <- add_taxa_group(lca_dat_clean, fishstat_dat)

# Output taxa groupings and sample sizes:
#data.frame(table(lca_dat_clean_groups$taxa_group_name))
#data.frame(table(lca_dat_clean_groups$taxa)) # abbreviated version of taxa_group_name for writing models
#lca_dat_clean_groups %>% select(taxa_group_name, clean_sci_name) %>% unique() %>% arrange(taxa_group_name)
#write.csv(data.frame(table(lca_dat_clean_groups$taxa_group_name)), file.path(outdir, "taxa_group_sample_size.csv"))
#write.csv(lca_dat_clean_groups %>% select(taxa_group_name, clean_sci_name) %>% unique() %>% arrange(taxa_group_name), file.path(outdir, "taxa_group_composition.csv"))

# CREATE STUDY ID COLUMN - use this for rejoining outputs from multiple regression models back together
# Select relevant data columns and arrange by categorical info
# Need FCR to identify which species aren't fed (FCR == 0)
lca_categories <- lca_dat_clean_groups %>%
  select(FCR, contains("new"), clean_sci_name, taxa, intensity = Intensity, system = Production_system_group) %>%
  arrange(clean_sci_name, taxa, intensity, system) %>%
  filter(taxa %in% c("mussel")==FALSE) %>% # Remove taxa that don't belong in FCR/feed analysis - mussels
  mutate(study_id = row_number())


######################################################################################################
# Model 1: Remove NA's, and estimate group-level feed conversion ratio for Nile tilapia, Oreochromis niloticus (species with the most FCR data)
lca_dat_simple <- lca_dat_clean %>%
  filter(is.na(FCR) == FALSE) %>%
  filter(clean_sci_name == "Oncorhynchus mykiss")
  #filter(clean_sci_name == "Dicentrarchus labrax") # Note: doesn't run for n = 1

# Set data
x <- lca_dat_simple$FCR
n <- nrow(lca_dat_simple)

# Gamma distribution model
stan_simple <- 'data {
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

# Compile
simple_mod <- stan_model(model_code = stan_simple, verbose = TRUE)
# Note: For Windows, apparently OK to ignore this warning message:
# Warning message:
#   In system(paste(CXX, ARGS), ignore.stdout = TRUE, ignore.stderr = TRUE) :
#   'C:/rtools40/usr/mingw_/bin/g++' not found

# Fit model:
fit_simple <- sampling(object = simple_mod, data = list(x = x,
                                                        n = n))

print(fit_simple)

launch_shinystan(fit_simple)

# Note: lp__ is the sum of the vector of log probabilities (but after removing any constant scale factors, making it not useful for model comparison)
# https://www.jax.org/news-and-insights/jax-blog/2015/october/lp-in-stan-output#:~:text=Therefore%2C%20%E2%80%9Clp%E2%80%9D%20is%20actually,useful%20for%20model%20comparison%20purposes.

# Diagnostics
stan_trace(fit_simple)


distribution_simple <- as.matrix(fit_simple)

plot_theme <- theme(axis.text=element_text(size=14, color = "black"))

p <- mcmc_areas_ridges(distribution_simple,
                       pars = vars(contains(c("mu", "sigma"))),
                       prob = 0.8) +
  ggtitle("Oncorhynchus mykiss FCR model", "with 80% credible intervals") +
  plot_theme

p 
#ggsave(file.path(outdir, "bayes-example_trout_fcr-gamma-target.png"), height = 8.5, width = 11)

######################################################################################################
# Model 2: Remove ALL NA's, and estimate group-level feed conversion ratio for all species and a global feed conversion ratio (mu)

# If desired, replicate data based on clean_sample_size column:
# lca_dat_clean <- rep_data(lca_dat_clean)
# IMPORTANT for convergence: filter out FCR == 0

lca_dat_groups <- lca_dat_clean %>%
  filter(is.na(FCR) == FALSE) %>% 
  filter(FCR != 0) %>%
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

stan_data <- list(x = x, n = n, j = j, sp = sp)

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
  // reparamaterize global gamma to get global mu and sigma
  real shape; 
  real rate;
  vector[j] sp_shape;
  vector[j] sp_rate;
  
  // global level
  shape = square(mu) / square(sigma);
  rate = mu / square(sigma);
  
  // sci-name level
  for (k in 1:n){
    sp_shape[sp[k]] = square(sp_mu[sp[k]]) ./ square(sp_sigma);
    sp_rate[sp[k]] = sp_mu[sp[k]] ./ square(sp_sigma);
  }
}
model {
  // Put priors on mu and sigma since this is more intuitive:
  //mu ~ uniform(0.05, 100);
  //sigma ~ uniform(0, 100);
  //sp_sigma ~ uniform(0, 100);
  
  // Specific prior on Salmonidae
  sp_mu[19] ~ uniform(6, 10);

  // likelihood
  // target += gamma_lpdf(sp_mu | shape, rate); // alternative notation
  sp_mu ~ gamma(shape, rate);
  for (i in 1:n){
    x[i] ~ gamma(sp_shape[sp[i]], sp_rate[sp[i]]); 
  }

}'

# Compile
no_missing_mod <- stan_model(model_code = stan_grouped, verbose = TRUE)
# Note: For Windows, apparently OK to ignore this warning message:
# Warning message:
#   In system(paste(CXX, ARGS), ignore.stdout = TRUE, ignore.stderr = TRUE) :
#   'C:/rtools40/usr/mingw_/bin/g++' not found

# Fit model:
fit_grouped <- sampling(object = no_missing_mod, data = stan_data,
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
                 area_method = "scaled height")
p2 + ggtitle("Posterior distributions", "with 80% credible intervals")
#ggsave(file.path(outdir, "plot_post-distribution-plot_FCR-by-species.png"), height = 8, width = 11.5)

p2 <- mcmc_areas(distribution_grouped,
                 pars = vars(contains("sp_shape")),
                 prob = 0.8,
                 area_method = "scaled height")
p2 + ggtitle("Posterior distributions", "with 80% credible intervals")

p2 <- mcmc_areas(distribution_grouped,
                 pars = vars(contains("sp_rate")),
                 prob = 0.8,
                 area_method = "scaled height")
p2 + ggtitle("Posterior distributions", "with 80% credible intervals")


######################################################################################################
# Model 3: Still no NA's, but Gamma distribution with 2-3 levels (sci-name, taxa-group, all-seafood)

######################################################################################################
# Model 3.1: Two-levels


# Test on smaller subset
lca_dat_groups <- lca_dat_clean %>%
  #filter(taxa_group_name %in% c("salmon/char")) %>%
  filter(taxa_group_name %in% c("salmon/char", "marine shrimp")) %>%
  filter(is.na(FCR) == FALSE) %>% 
  filter(FCR != 0) %>%
  mutate(taxa_group_name = as.factor(taxa_group_name),
         tx = as.numeric(taxa_group_name)) %>%
  select(clean_sci_name, taxa_group_name, tx, FCR) %>%
  arrange(tx)



# Run analysis on all non-NA data
# lca_dat_groups <- lca_dat_clean %>%
#   filter(is.na(FCR) == FALSE) %>% 
#   filter(FCR != 0) %>%
#   mutate(taxa_group_name = as.factor(taxa_group_name),
#          tx = as.numeric(taxa_group_name)) %>%
#   select(clean_sci_name, taxa_group_name, tx, FCR) %>%
# arrange(sci)

# Boxplot of data
ggplot(data = lca_dat_groups, aes(x = clean_sci_name, y = FCR)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 16)) + 
  labs(title = "Boxplots of FCR by scientific name")
#ggsave(file.path(outdir, "plot_boxplot_FCR-by-species.png"), height = 8, width = 11.5)

ggplot(data = lca_dat_groups, aes(x = taxa_group_name, y = FCR)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 16)) + 
  labs(title = "Boxplots of FCR by taxa group")

groupings_and_sample_size_key <- lca_dat_groups %>%
  group_by(taxa_group_name) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  select(taxa_group_name, n) %>%
  unique()

x <- lca_dat_groups$FCR
N <- nrow(lca_dat_groups)
N_TX <- length(unique(lca_dat_groups$tx))
tx <- lca_dat_groups$tx

stan_data <- list(x = x, N = N, N_TX = N_TX, tx = tx)

stan_grouped <- 'data {
  int<lower=0> N;  // number of observations
  vector<lower=0>[N] x; // data
  int N_TX; // number of taxa groups
  int tx[N]; // taxa group index
}
parameters {
  real<lower=0> mu;
  real<lower=0> sigma;
  vector<lower=0>[N_TX] tx_mu;
  real<lower=0> tx_sigma;
}
transformed parameters {
  // reparamaterize gamma to get mu and sigma; defining these here instead of the model section allows us to see these parameters in the output
  real shape;
  real rate;
  vector[N_TX] tx_shape;
  vector[N_TX] tx_rate;

  // global-level
  shape = square(mu) / square(sigma);
  rate = mu / square(sigma);
  
  // taxa group level
  for (n_tx in 1:N_TX){
    tx_shape[n_tx] = square(tx_mu[n_tx]) ./ square(tx_sigma);
    tx_rate[n_tx] = tx_mu[n_tx] ./ square(tx_sigma);
  }
}
model {
  // Put priors on mu and sigma (instead of shape and rate) since this is more intuitive:
  mu ~ uniform(0, 10);
  sigma ~ uniform(0, 10);
  tx_mu ~ uniform(0, 10);
  tx_sigma ~ uniform(0, 10);

  // likelihood
  for (n in 1:N){
    x[n] ~ gamma(tx_shape[tx[n]], tx_rate[tx[n]]); // equivalent notation to target += but faster
  }
  
  for (n_tx in 1:N_TX){
    tx_mu[n_tx] ~ gamma(shape, rate);
  }
}'


# OLD CODE: combined loop
# stan_grouped <- 'data {
#   int<lower=0> N;  // number of observations
#   vector<lower=0>[N] x; // data
#   int N_TX; // number of taxa groups
#   int tx[N]; // taxa group index (ordered by unique sci index)
#   int N_SCI; // number of scientific names
#   int sci[N]; // sciname index
# }
# parameters {
#   real<lower=0> mu;
#   real<lower=0> sigma;
#   vector<lower=0>[N_TX] tx_mu;
#   real<lower=0> tx_sigma;
#   vector<lower=0>[N_SCI] sci_mu;
#   real<lower=0> sci_sigma;
# }
# transformed parameters {
#   // reparamaterize gamma to get mu and sigma; defining these here instead of the model section allows us to see these parameters in the output
#   real shape;
#   real rate; 
#   vector[N_SCI] sci_shape;
#   vector[N_SCI] sci_rate;
#   vector[N_TX] tx_shape;
#   vector[N_TX] tx_rate;
#   
#   // global-level
#   shape = square(mu) / square(sigma);
#   rate = mu / square(sigma);
#   // sci name and taxa group levels
#   for (j in 1:N){
#     tx_shape[tx[j]] = square(tx_mu[tx[j]]) ./ square(tx_sigma);
#     tx_rate[tx[j]] = tx_mu[tx[j]] ./ square(tx_sigma);
#     sci_shape[sci[j]] = square(sci_mu[sci[j]]) ./ square(sci_sigma);
#     sci_rate[sci[j]] = sci_mu[sci[j]] ./ square(sci_sigma);
#   }
# }
# model {
#   // Put priors on mu and sigma (instead of shape and rate) since this is more intuitive:
#   //mu ~ uniform(0.05, 100);
#   //sigma ~ uniform(0, 100);
#   //sci_sigma ~ uniform(0, 100);
#   
#   // Specific prior for Salmonidae
#   // sci_mu[19] ~ uniform(6, 10);
#   
#   // likelihood
#   tx_mu ~ gamma(shape, rate);
#   for (i in 1:N){
#     sci_mu[sci[i]] ~ gamma(tx_shape[tx[i]], tx_rate[tx[i]]);
#     x[i] ~ gamma(sci_shape[sci[i]], sci_rate[sci[i]]); // equivalent notation to target += but faster
#   }
# }'

# Compile
no_missing_mod <- stan_model(model_code = stan_grouped, verbose = TRUE)
# Note: For Windows, apparently OK to ignore this warning message:
# Warning message:
#   In system(paste(CXX, ARGS), ignore.stdout = TRUE, ignore.stderr = TRUE) :
#   'C:/rtools40/usr/mingw_/bin/g++' not found

# Fit model:
fit_grouped <- sampling(object = no_missing_mod, data = stan_data, cores = 4)
#iter = 10000)

print(fit_grouped)

######################################################################################################
# Model 3.2: Three-levels
# If desired, replicate data based on clean_sample_size column:
# lca_dat_clean <- rep_data(lca_dat_clean)
# IMPORTANT for convergence: filter out FCR == 0



# Test on smaller subset
lca_dat_groups <- lca_categories %>%
  #filter(taxa_group_name %in% c("salmon/char")) %>%
  filter(taxa %in% c("salmon", "shrimp")) %>%
  filter(is.na(FCR) == FALSE) %>% 
  filter(FCR != 0) %>%
  mutate(clean_sci_name = as.factor(clean_sci_name),
         sci = as.numeric(clean_sci_name),
         taxa = as.factor(taxa),
         tx = as.numeric(taxa)) %>%
  select(clean_sci_name, sci, taxa, tx, FCR) %>%
  arrange(sci)

# Boxplot of data
ggplot(data = lca_dat_groups, aes(x = clean_sci_name, y = FCR)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 16)) + 
  labs(title = "Boxplots of FCR by scientific name")
#ggsave(file.path(outdir, "plot_boxplot_FCR-by-species.png"), height = 8, width = 11.5)

ggplot(data = lca_dat_groups, aes(x = taxa, y = FCR)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 16)) + 
  labs(title = "Boxplots of FCR by taxa group")

groupings_and_sample_size_key <- lca_dat_groups %>%
  group_by(clean_sci_name) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  select(clean_sci_name, taxa, n) %>%
  unique()

x <- lca_dat_groups$FCR
N <- nrow(lca_dat_groups)
N_SCI <- length(unique(lca_dat_groups$sci))
n_to_sci <- lca_dat_groups$sci
N_TX <- length(unique(lca_dat_groups$tx))
sci_to_tx <- lca_dat_groups %>%
  select(sci, tx) %>%
  unique() %>%
  pull(tx)
  
#tx <- unique(lca_dat_groups[c("sci", "tx")])[,"tx"]

stan_data <- list(x = x, N = N, N_SCI = N_SCI, n_to_sci = n_to_sci, N_TX = N_TX, sci_to_tx = sci_to_tx)

stan_grouped <- 'data {
  int<lower=0> N;  // number of observations
  vector<lower=0>[N] x; // data
  int N_TX; // number of taxa groups
  int N_SCI; // number of scientific names
  int sci_to_tx[N_SCI]; // taxa group index
  int n_to_sci[N]; // sciname index
}
parameters {
  real<lower=0> mu;
  real<lower=0> sigma;
  vector<lower=0>[N_TX] tx_mu;
  real<lower=0> tx_sigma;
  vector<lower=0>[N_SCI] sci_mu;
  real<lower=0> sci_sigma;
}
transformed parameters {
  // reparamaterize gamma to get mu and sigma; defining these here instead of the model section allows us to see these parameters in the output
  real shape;
  real rate;
  vector[N_SCI] sci_shape;
  vector[N_SCI] sci_rate;
  vector[N_TX] tx_shape;
  vector[N_TX] tx_rate;

  // global-level
  shape = square(mu) / square(sigma);
  rate = mu / square(sigma);
  
  // taxa group level
  for (n_tx in 1:N_TX){
    tx_shape[n_tx] = square(tx_mu[n_tx]) ./ square(tx_sigma);
    tx_rate[n_tx] = tx_mu[n_tx] ./ square(tx_sigma);
  }
  
  // tx_mu = sum(sci_mu * weight_production) // possibly how to do weighting by production
  
  // sci name level
  for (n_sci in 1:N_SCI){
    sci_shape[n_sci] = square(sci_mu[n_sci]) ./ square(sci_sigma);
    sci_rate[n_sci] = sci_mu[n_sci] ./ square(sci_sigma);
  }
}
model {
  // Put priors on mu and sigma (instead of shape and rate) since this is more intuitive:
  mu ~ uniform(0, 100);
  sigma ~ uniform(0, 100);
  tx_mu ~ uniform(0, 100);
  tx_sigma ~ uniform(0, 100);
  sci_mu ~ uniform(0, 100);
  
  // Model only converges when prior given to Salmonidae 
  //for (i in 1:N_SCI){
  //  if (i == 6) 
  //  sci_mu[i] ~ uniform(6, 10);
  //  else
  //  sci_mu[i] ~ uniform(0, 10);
  //}
  sci_sigma ~ uniform(0, 100);

  // likelihood
  for (n in 1:N){
    x[n] ~ gamma(sci_shape[n_to_sci[n]], sci_rate[n_to_sci[n]]); // equivalent notation to target += but faster
  }
  
  for (n_sci in 1:N_SCI){
    sci_mu[n_sci] ~ gamma(tx_shape[sci_to_tx[n_sci]], tx_rate[sci_to_tx[n_sci]]);
  }
  
  for (n_tx in 1:N_TX){
    tx_mu[n_tx] ~ gamma(shape, rate);
  }
}'


# Compile
no_missing_mod <- stan_model(model_code = stan_grouped)
# Note: For Windows, apparently OK to ignore this warning message:
# Warning message:
#   In system(paste(CXX, ARGS), ignore.stdout = TRUE, ignore.stderr = TRUE) :
#   'C:/rtools40/usr/mingw_/bin/g++' not found

# Fit model:
fit_grouped <- sampling(object = no_missing_mod, data = stan_data, cores = 4,
                        control = list(adapt_delta = 0.99))
                        #iter = 10000)

print(fit_grouped)

# How to interpret sci index numbers:
# What are the sample sizes per group:
sci_index_key <- lca_dat_groups %>%
  group_by(clean_sci_name) %>%
  mutate(n_sci_name = n()) %>%
  ungroup() %>%
  select(sci, clean_sci_name, n_sci_name) %>%
  arrange(sci) %>%
  unique() %>%
  mutate(mu_param_name = paste("sci_mu[", clean_sci_name, "]", sep =""),
         shape_param_name = paste("sci_shape[", clean_sci_name, "]", sep =""),
         rate_param_name = paste("sci_rate[", clean_sci_name, "]", sep = ""))
#write.csv(sp_index_key, file.path(outdir, "sp_info.csv"), row.names = FALSE)

taxa_index_key <- lca_dat_groups %>%
  group_by(taxa_group_name) %>%
  mutate(n_taxa_name = n()) %>%
  ungroup() %>%
  select(tx, taxa_group_name, n_taxa_name) %>%
  arrange(tx) %>%
  unique() %>%
  mutate(mu_param_name = paste("tx_mu[", taxa_group_name, "]", sep =""),
         shape_param_name = paste("tx_shape[", taxa_group_name, "]", sep =""),
         rate_param_name = paste("tx_rate[", taxa_group_name, "]", sep = ""))

# Replace param names
names(fit_grouped)[grep(names(fit_grouped), pattern = "sci_mu")] <- sci_index_key$mu_param_name
names(fit_grouped)[grep(names(fit_grouped), pattern = "sci_shape")] <- sci_index_key$shape_param_name
names(fit_grouped)[grep(names(fit_grouped), pattern = "sci_rate")] <- sci_index_key$rate_param_name

names(fit_grouped)[grep(names(fit_grouped), pattern = "tx_mu")] <- taxa_index_key$mu_param_name
names(fit_grouped)[grep(names(fit_grouped), pattern = "tx_shape")] <- taxa_index_key$shape_param_name
names(fit_grouped)[grep(names(fit_grouped), pattern = "tx_rate")] <- taxa_index_key$rate_param_name

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

# Sci-name level:
p2 <- mcmc_areas(distribution_grouped,
                 pars = vars(starts_with("mu")|contains("sci_mu")),
                 prob = 0.8,
                 area_method = "scaled height")
p2 + ggtitle("Posterior distributions", "with 80% credible intervals")
#ggsave(file.path(outdir, "plot_post-distribution-plot_FCR-by-species.png"), height = 8, width = 11.5)

p2 <- mcmc_areas(distribution_grouped,
                 pars = vars(contains("sci_shape")),
                 prob = 0.8,
                 area_method = "scaled height")
p2 + ggtitle("Posterior distributions", "with 80% credible intervals")

p2 <- mcmc_areas(distribution_grouped,
                        pars = vars(contains("sci_rate")),
                        prob = 0.8,
                        area_method = "scaled height")
p2 + ggtitle("Posterior distributions", "with 80% credible intervals")

# Taxa-group level
p2 <- mcmc_areas(distribution_grouped,
                 pars = vars(starts_with("mu")|contains("tx_mu")),
                 prob = 0.8,
                 area_method = "scaled height")
p2 + ggtitle("Posterior distributions", "with 80% credible intervals")
#ggsave(file.path(outdir, "plot_post-distribution-plot_FCR-by-species.png"), height = 8, width = 11.5)

p2 <- mcmc_areas(distribution_grouped,
                 pars = vars(contains("tx_shape")),
                 prob = 0.8,
                 area_method = "scaled height")

p2 <- mcmc_areas(distribution_grouped,
                 pars = vars(contains("tx_rate")),
                 prob = 0.8,
                 area_method = "scaled height")

p2 + ggtitle("Posterior distributions", "with 80% credible intervals")
#ggsave(file.path(outdir, "plot_post-distribution-plot_FCR-by-species.png"), height = 8, width = 11.5)