# Combine bayes_prop_feed and bayes_fcr to calculate footprint

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
  #select(clean_sci_name, Feed_soy_percent, Feed_othercrops_percent, Feed_FMFO_percent, Feed_animal_percent, taxa_group_name) %>%
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

fp_dat <- read.csv(file.path(datadir, "Feed_FP_raw.csv"))
fp_clean <- clean.feedFP(fp_dat)

######################################################################################################
# Model 2: Remove all NAs - estimate feed footprint for all sci names

# Remove NAs
lca_dat_no_na <- lca_dat_no_zeroes %>%
  filter(is.na(Feed_soy_percent)==FALSE) 

# BOX PLOTS OF DATA:

# Theme for ALL PLOTS (including mcmc plots)
plot_theme <- theme(axis.text=element_text(size=14, color = "black"))

# FCR:
plot_fcr <- lca_dat_no_na %>%
  select(clean_sci_name, FCR) %>%
  #mutate(clean_sci_name = paste(clean_sci_name, row_number(), sep = "")) %>%
  pivot_longer(cols = FCR)
ggplot(data = plot_fcr, aes(x = clean_sci_name, y = value)) +
  geom_boxplot() +
  theme_classic() +
  plot_theme +
  labs(title = "Boxplots of FCRs",
       x = "",
       y = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#ggsave(file.path(outdir, "boxplot_fcr_no_na.png"), height = 8, width = 11.5)


# feed proportion:
plot_feed_prop <- lca_dat_no_na %>%
  select(clean_sci_name, soy = feed_soy_new, crops = feed_crops_new, fmfo = feed_fmfo_new, animal = feed_animal_new) %>%
  pivot_longer(cols = soy:animal)



feed_vars <- c("soy", "crops", "fmfo", "animal")
for (i in 1:length(feed_vars)) {
  p <- ggplot(data = plot_feed_prop %>% filter(name == feed_vars[i]), aes(x = clean_sci_name, y = value)) +
    geom_boxplot() +
    theme_classic() +
    plot_theme +
    labs(title = paste("Boxplots of ", feed_vars[i], " feed proportions", sep = ""),
         x = "",
         y = "")  +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p)
  #ggsave(file.path(outdir, "boxplot_feed-prop_no_na.png"), height = 8, width = 11.5)
}

# STILL NEED TO REDO THIS SECTION
# Set data for model
# for FCR model:
x <- lca_dat_no_na$FCR

# for Feed proportion model:
k = 4
n = 3
feed_weights <- lca_dat_no_na %>%
  select(contains("new")) %>%
  as.matrix()

# for final foot print calculation
fp_dat <- fp_clean %>%
  filter(Category != "Energy") %>%
  select(FP, Category, FP_val) 

# ORDER: Animal, crop, FMFO, soy
fp_carbon_dat <- fp_dat %>%
  filter(FP == "Carbon") %>%
  select(-FP) %>%
  pivot_wider(names_from = Category, values_from = FP_val) %>%
  as.matrix() %>%
  c()

fp_nitrogen_dat <- fp_dat %>%
  filter(FP == "Nitrogen") %>%
  select(-FP) %>%
  pivot_wider(names_from = Category, values_from = FP_val) %>%
  as.matrix() %>%
  c()

fp_phosphorus_dat <- fp_dat %>%
  filter(FP == "Phosphorus") %>%
  select(-FP) %>%
  pivot_wider(names_from = Category, values_from = FP_val) %>%
  as.matrix() %>%
  c()

fp_land_dat <- fp_dat %>%
  filter(FP == "Land") %>%
  select(-FP) %>%
  pivot_wider(names_from = Category, values_from = FP_val) %>%
  as.matrix() %>%
  c()

fp_water_dat <- fp_dat %>%
  filter(FP == "Water") %>%
  select(-FP) %>%
  pivot_wider(names_from = Category, values_from = FP_val) %>%
  as.matrix() %>%
  c()

# Estimate foot print for all scientific names without NAs
stan_pooled <- 'data {
  int<lower=0> n;  // number of observations
  vector[n] x; // data
  int<lower=1> k; // number of feed types
  simplex[k] feed_weights[n]; // array of feed weights simplexes
  vector[k] fp_carbon_dat;
  vector[k] fp_nitrogen_dat;
  vector[k] fp_phosphorus_dat;
  vector[k] fp_land_dat;
  vector[k] fp_water_dat;
}
parameters {
  // FCR model:
  real<lower=0> mu;
  real<lower=0> sigma;
  // Feed proportion model:
  vector<lower=0>[k] alpha;
  simplex[k] theta;
}
model {
  
  x ~ normal(mu, sigma); // note: stan interprets second param as standard deviation

  for (i in 1:n) {
    feed_weights[i] ~ dirichlet(alpha); // estimate vector of alphas based on the data of feed weights
  }
  theta ~ dirichlet(alpha); // now, estimate feed weights based on the vector of alphas
}
generated quantities {
  // Carbon
  real<lower=0> species_carbon_footprint;
  vector[k] feed_carbon_footprint;
  real total_feed_carbon_footprint;
  // Nitrogen
  real<lower=0> species_nitrogen_footprint;
  vector[k] feed_nitrogen_footprint;
  real total_feed_nitrogen_footprint;
  // Phosphorus
  real<lower=0> species_phosphorus_footprint;
  vector[k] feed_phosphorus_footprint;
  real total_feed_phosphorus_footprint;
  // Land
  real<lower=0> species_land_footprint;
  vector[k] feed_land_footprint;
  real total_feed_land_footprint;
  // Water
  real<lower=0> species_water_footprint;
  vector[k] feed_water_footprint;
  real total_feed_water_footprint;
  
  // Calculations
  feed_carbon_footprint = fp_carbon_dat .* theta;
  total_feed_carbon_footprint = sum(feed_carbon_footprint);
  species_carbon_footprint = mu * total_feed_carbon_footprint;

  feed_nitrogen_footprint = fp_nitrogen_dat .* theta;
  total_feed_nitrogen_footprint = sum(feed_nitrogen_footprint);
  species_nitrogen_footprint = mu * total_feed_nitrogen_footprint;

  feed_phosphorus_footprint = fp_phosphorus_dat .* theta;
  total_feed_phosphorus_footprint = sum(feed_phosphorus_footprint);
  species_phosphorus_footprint = mu * total_feed_phosphorus_footprint;
  
  feed_land_footprint = fp_land_dat .* theta;
  total_feed_land_footprint = sum(feed_land_footprint);
  species_land_footprint = mu * total_feed_land_footprint;
  
  feed_water_footprint = fp_water_dat .* theta;
  total_feed_water_footprint = sum(feed_water_footprint);
  species_water_footprint = mu * total_feed_water_footprint;
}'


######################################################################################################
# Model 1: Remove all NAs - estimate feed footprint for Oncorhynchus mykiss

# Remove NAs
lca_dat_no_na <- lca_dat_no_zeroes %>%
  filter(clean_sci_name == "Oncorhynchus mykiss") %>%
  filter(is.na(Feed_soy_percent)==FALSE) 


# BOX PLOTS OF DATA:

# Theme for ALL PLOTS (including mcmc plots)
plot_theme <- theme(axis.text=element_text(size=14, color = "black"))

# FCR:
plot_fcr <- lca_dat_no_na %>%
  select(clean_sci_name, FCR) %>%
  mutate(clean_sci_name = paste(clean_sci_name, row_number(), sep = "")) %>%
  pivot_longer(cols = FCR)
ggplot(data = plot_fcr, aes(x = name, y = value)) +
  geom_boxplot() +
  theme_classic() +
  plot_theme +
  labs(title = "Boxplots of FCRs for Oncorhynchus mykiss",
       x = "",
       y = "")
ggsave(file.path(outdir, "boxplot_fcr_trout.png"), height = 8, width = 11.5)


# feed proportion:
plot_feed_prop <- lca_dat_no_na %>%
  select(clean_sci_name, soy = feed_soy_new, crops = feed_crops_new, fmfo = feed_fmfo_new, animal = feed_animal_new) %>%
  mutate(clean_sci_name = paste(clean_sci_name, row_number(), sep = "")) %>%
  pivot_longer(cols = soy:animal)

ggplot(data = plot_feed_prop, aes(x = name, y = value)) +
  geom_boxplot() +
  theme_classic() +
  plot_theme +
  labs(title = "Boxplots of feed proportions for Oncorhynchus mykiss",
       x = "",
       y = "")
ggsave(file.path(outdir, "boxplot_feed-prop_trout.png"), height = 8, width = 11.5)

# Set data for model
# for FCR model:
x <- lca_dat_no_na$FCR

# for Feed proportion model:
k = 4
n = 3
feed_weights <- lca_dat_no_na %>%
  select(contains("new")) %>%
  as.matrix()

# for final foot print calculation
fp_dat <- fp_clean %>%
  filter(Category != "Energy") %>%
  select(FP, Category, FP_val) 

# ORDER: Animal, crop, FMFO, soy
fp_carbon_dat <- fp_dat %>%
  filter(FP == "Carbon") %>%
  select(-FP) %>%
  pivot_wider(names_from = Category, values_from = FP_val) %>%
  as.matrix() %>%
  c()

fp_nitrogen_dat <- fp_dat %>%
  filter(FP == "Nitrogen") %>%
  select(-FP) %>%
  pivot_wider(names_from = Category, values_from = FP_val) %>%
  as.matrix() %>%
  c()

fp_phosphorus_dat <- fp_dat %>%
  filter(FP == "Phosphorus") %>%
  select(-FP) %>%
  pivot_wider(names_from = Category, values_from = FP_val) %>%
  as.matrix() %>%
  c()

fp_land_dat <- fp_dat %>%
  filter(FP == "Land") %>%
  select(-FP) %>%
  pivot_wider(names_from = Category, values_from = FP_val) %>%
  as.matrix() %>%
  c()

fp_water_dat <- fp_dat %>%
  filter(FP == "Water") %>%
  select(-FP) %>%
  pivot_wider(names_from = Category, values_from = FP_val) %>%
  as.matrix() %>%
  c()

# note: dirichlet_rng is just a random number generator:
# rep_vector(x, m) creates a column consisting of m copies of x
# generated quantities {
#   vector[k] theta = dirichlet_rng(rep_vector(alpha, k));
# }

# Estimate feed component proportions for a single species
stan_pooled <- 'data {
  int<lower=0> n;  // number of observations
  vector[n] x; // data
  int<lower=1> k; // number of feed types
  simplex[k] feed_weights[n]; // array of feed weights simplexes
  vector[k] fp_carbon_dat;
  vector[k] fp_nitrogen_dat;
  vector[k] fp_phosphorus_dat;
  vector[k] fp_land_dat;
  vector[k] fp_water_dat;
}
parameters {
  // FCR model:
  real<lower=0> mu;
  real<lower=0> sigma;
  // Feed proportion model:
  vector<lower=0>[k] alpha;
  simplex[k] theta;
}
model {
  
  x ~ normal(mu, sigma); // note: stan interprets second param as standard deviation

  for (i in 1:n) {
    feed_weights[i] ~ dirichlet(alpha); // estimate vector of alphas based on the data of feed weights
  }
  theta ~ dirichlet(alpha); // now, estimate feed weights based on the vector of alphas
}
generated quantities {
  // Carbon
  real<lower=0> species_carbon_footprint;
  vector[k] feed_carbon_footprint;
  real total_feed_carbon_footprint;
  // Nitrogen
  real<lower=0> species_nitrogen_footprint;
  vector[k] feed_nitrogen_footprint;
  real total_feed_nitrogen_footprint;
  // Phosphorus
  real<lower=0> species_phosphorus_footprint;
  vector[k] feed_phosphorus_footprint;
  real total_feed_phosphorus_footprint;
  // Land
  real<lower=0> species_land_footprint;
  vector[k] feed_land_footprint;
  real total_feed_land_footprint;
  // Water
  real<lower=0> species_water_footprint;
  vector[k] feed_water_footprint;
  real total_feed_water_footprint;
  
  // Calculations
  feed_carbon_footprint = fp_carbon_dat .* theta;
  total_feed_carbon_footprint = sum(feed_carbon_footprint);
  species_carbon_footprint = mu * total_feed_carbon_footprint;

  feed_nitrogen_footprint = fp_nitrogen_dat .* theta;
  total_feed_nitrogen_footprint = sum(feed_nitrogen_footprint);
  species_nitrogen_footprint = mu * total_feed_nitrogen_footprint;

  feed_phosphorus_footprint = fp_phosphorus_dat .* theta;
  total_feed_phosphorus_footprint = sum(feed_phosphorus_footprint);
  species_phosphorus_footprint = mu * total_feed_phosphorus_footprint;
  
  feed_land_footprint = fp_land_dat .* theta;
  total_feed_land_footprint = sum(feed_land_footprint);
  species_land_footprint = mu * total_feed_land_footprint;
  
  feed_water_footprint = fp_water_dat .* theta;
  total_feed_water_footprint = sum(feed_water_footprint);
  species_water_footprint = mu * total_feed_water_footprint;
}'

no_missing_mod <- stan_model(model_code = stan_pooled, verbose = TRUE)
# Note: For Windows, apparently OK to ignore this warning message:
# Warning message:
#   In system(paste(CXX, ARGS), ignore.stdout = TRUE, ignore.stderr = TRUE) :
#   'C:/rtools40/usr/mingw_/bin/g++' not found

# Fit model:
fit_pooled <- sampling(object = no_missing_mod, data = list(n = n,
                                                            x = x,
                                                            k = k,
                                                            feed_weights = feed_weights,
                                                            fp_carbon_dat = fp_carbon_dat,
                                                            fp_nitrogen_dat = fp_nitrogen_dat,
                                                            fp_phosphorus_dat = fp_phosphorus_dat,
                                                            fp_land_dat = fp_land_dat,
                                                            fp_water_dat = fp_water_dat),
                       iter = 10000, cores = 4,
                       control = list(adapt_delta = 0.99))
print(fit_pooled)

feeds <- c("soy", "crops", "fmfo", "animal")
feed_key <- data.frame(carbon_footprint = paste("carbon_footprint[", feeds, "]", sep = ""),
                       nitrogen_footprint = paste("nitrogen_footprint[", feeds, "]", sep = ""),
                       phosphorus_footprint = paste("phosphorus_footprint[", feeds, "]", sep = ""),
                       land_footprint = paste("land_footprint[", feeds, "]", sep = ""),
                       water_footprint = paste("water_footprint[", feeds, "]", sep = ""))

fit_pooled_clean <- fit_pooled
names(fit_pooled_clean)[grep(names(fit_pooled_clean), pattern = "feed_carbon_footprint\\[")] <- feed_key$carbon_footprint
names(fit_pooled_clean)[grep(names(fit_pooled_clean), pattern = "feed_nitrogen_footprint\\[")] <- feed_key$nitrogen_footprint
names(fit_pooled_clean)[grep(names(fit_pooled_clean), pattern = "feed_phosphorus_footprint\\[")] <- feed_key$phosphorus_footprint
names(fit_pooled_clean)[grep(names(fit_pooled_clean), pattern = "feed_land_footprint\\[")] <- feed_key$land_footprint
names(fit_pooled_clean)[grep(names(fit_pooled_clean), pattern = "feed_water_footprint\\[")] <- feed_key$water_footprint

distribution_pooled <- as.matrix(fit_pooled_clean)

# FIX IT - replace parameter names and add plots for other parameters
p_footprint <- mcmc_areas(distribution_pooled,
                          pars = vars(contains("carbon_footprint")),
                          prob = 0.8,
                          prob_outer = 0.9,
                          area_method = "scaled height",
                          point_est = "median") + 
  ggtitle("Oncorhynchus mykiss full feed footprint model", "with 80% credible intervals") +
  xlim(0, 10) +
  plot_theme 

p_footprint
ggsave(filename = file.path(outdir, "bayes-example_trout_carbon-feed-footprint.png"), width = 11, height = 8.5)


p_footprint <- mcmc_areas(distribution_pooled,
                          pars = vars(contains("nitrogen_footprint")),
                          prob = 0.8,
                          prob_outer = 0.9,
                          area_method = "scaled height",
                          point_est = "median") + 
  ggtitle("Oncorhynchus mykiss full feed footprint model", "with 80% credible intervals") +
  xlim(0, 0.01) +
  plot_theme 

p_footprint
ggsave(filename = file.path(outdir, "bayes-example_trout_nitrogen-feed-footprint.png"), width = 11, height = 8.5)

p_footprint <- mcmc_areas(distribution_pooled,
                          pars = vars(contains("phosphorus_footprint")),
                          prob = 0.8,
                          prob_outer = 0.9,
                          area_method = "scaled height",
                          point_est = "median") + 
  ggtitle("Oncorhynchus mykiss full feed footprint model", "with 80% credible intervals") +
  xlim(0, 0.001) +
  plot_theme 

p_footprint
ggsave(filename = file.path(outdir, "bayes-example_trout_phosphorus-feed-footprint.png"), width = 11, height = 8.5)

p_footprint <- mcmc_areas(distribution_pooled,
                          pars = vars(contains("land_footprint")),
                          prob = 0.8,
                          prob_outer = 0.9,
                          area_method = "scaled height",
                          point_est = "median") + 
  ggtitle("Oncorhynchus mykiss full feed footprint model", "with 80% credible intervals") +
  xlim(0, 5) +
  plot_theme 

p_footprint
ggsave(filename = file.path(outdir, "bayes-example_trout_land-feed-footprint.png"), width = 11, height = 8.5)

p_footprint <- mcmc_areas(distribution_pooled,
                          pars = vars(contains("water_footprint")),
                          prob = 0.8,
                          prob_outer = 0.9,
                          area_method = "scaled height",
                          point_est = "median") + 
  ggtitle("Oncorhynchus mykiss full feed footprint model", "with 80% credible intervals") +
  xlim(0, 0.2) +
  plot_theme 

p_footprint
ggsave(filename = file.path(outdir, "bayes-example_trout_water-feed-footprint.png"), width = 11, height = 8.5)
