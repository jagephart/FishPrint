# Bayesian regression to estimate variables based on taxa-group, intensity, and production system

rm(list=ls())
library(tidyverse)
library(rstan)
library(data.table)
library(countrycode) # part of clean.lca
library(bayesplot) # for mcmc_areas_ridges
library(brms)

# Mac
datadir <- "/Volumes/jgephart/BFA Environment 2/Data"
outdir <- "/Volumes/jgephart/BFA Environment 2/Outputs"
# Windows
# datadir <- "K:/BFA Environment 2/Data"
# outdir <- "K:BFA Environment 2/Outputs"

lca_dat <- read.csv(file.path(datadir, "LCA_compiled_20201109.csv"), fileEncoding="UTF-8-BOM") #fileEncoding needed when reading in file from windows computer (suppresses BOM hidden characters)
source("Functions.R")

# Clean LCA data
lca_dat_clean <- clean.lca(LCA_data = lca_dat)

# Rebuild FAO fish production from zip file
fishstat_dat <- rebuild_fish("/Volumes/jgephart/FishStatR/Data/Production-Global/ZippedFiles/GlobalProduction_2019.1.0.zip")

# Classify species into taxa groupings
# Use ISSCAAP grouping (in FAO data) to help with classification
lca_dat_clean <- add_taxa_group(lca_dat_clean, fishstat_dat)

# Output taxa groupings and sample sizes:
#data.frame(table(lca_dat_clean$taxa_group_name))
#write.csv(lca_dat_clean %>% select(taxa_group_name, clean_sci_name) %>% unique() %>% arrange(taxa_group_name), "taxa_group_update.csv")
######################################################################################################
# Model 1: Remove ALL NA's (predictors and response), and estimate feed conversion ratio for all scientific names, non-hierarchical

# Remove NAs (which also removes bivalves) and 0's
lca_dat_groups <- lca_dat_clean %>%
  select(clean_sci_name, FCR, taxa = taxa_group_name, intensity = Intensity, system = Production_system_group) %>%
  filter(is.na(FCR) == FALSE) %>% 
  filter(FCR != 0) %>%
  drop_na() %>%
  arrange(clean_sci_name, intensity, system) %>%
  mutate(taxa = as.factor(taxa),
         tx = as.numeric(taxa),
         intensity = as.factor(intensity),
         its = as.numeric(intensity),
         system = as.factor(system),
         ps = as.numeric(system)) %>%
  select(clean_sci_name, FCR, taxa, tx, intensity, its, system, ps)

# Does variance increase with mean? (values should bunch up around zero because there is no room for negative FCR)
ggplot(lca_dat_groups %>%
         mutate(combo = paste(taxa, system, intensity, sep = " ")), aes(x = combo, y = FCR)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 45))

# Set data:
N = nrow(lca_dat_groups)
y = lca_dat_groups$FCR

# Creat model matrix, note: automatically creates column for intercept
X <- model.matrix(object = ~taxa + intensity + system, 
                  data = lca_dat_groups %>% select(taxa, intensity, system)) 
K = ncol(X) # number of predictors

# Center all non-intercept variables and scale by 2 standard deviations
covars.sd<-apply(X[,-1], MARGIN=2, FUN=sd)
X_scaled <- scale(X[,-1], center=TRUE, scale=2*covars.sd) 

# Recreate intercept column
intercept_col <- rep(1, nrow(X_scaled))
X_scaled <- cbind(intercept_col, X_scaled)

######################################################################################################
# Model 1.1: Using STAN

# If desired, replicate data based on clean_sample_size column:
# lca_dat_clean <- rep_data(lca_dat_clean)

stan_data <- list(N = N, y = y, K = K, X_scaled = X_scaled)
# FIX IT - is canonical inverse link more appropriate for categorical variables?

# FIX IT - reparamterize on shape and mean?
# I found it is better to use
# alpha = shape;
# beta = shape/mu;
# You can keep your definition of mu the same, but instead of using phi as a precision term, we specify a shape parameter. I believe this parameterization is what rstanarm uses.

stan_one_level <- 'data {
  int N; // number of observations
  real y[N]; // response
  int K; // number of predictors, i.e., number of columns in the design matrix X_scaled
  matrix [N, K] X_scaled; // design matrix X_scaled
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
  mu = exp(X_scaled * beta); // log-link since shape and rate params are positive
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

# Use bayesplot for posterior predictive checks:
distribution_simple <- as.matrix(fit_simple)
y_reps <- distribution_simple[,grep(pattern = "y_rep", colnames(distribution_simple))]

color_scheme_set("brightblue")
# Distributions of y vs yreps (specify number of samples to display from y_rep)
ppc_dens_overlay(y, y_reps[1:50,])
# Errors: normally distributed?
ppc_error_hist(y, y_reps[1:5,])
ppc_error_scatter_avg(y, y_reps)
# Errors vs x
for (i in 2:ncol(X_scaled)){ # skip intercept column
  col_i <- colnames(X_scaled)[i]
  p <- ppc_error_scatter_avg_vs_x(y, y_reps, X_scaled[,col_i]) + ggtitle(col_i)
  print(p)
}

# Clean up param names:
fit_simple_clean <- fit_simple

# Observation-level param names (species - intensity - system - study number combo)
study_index_key <- lca_dat_groups %>%
  select(clean_sci_name, intensity, system) %>%
  group_by(clean_sci_name, intensity, system) %>%
  # Add row number to differentiate 
  mutate(study_number = row_number()) %>%
  mutate(mu_param_name = paste("mu[", clean_sci_name, "-", intensity, "-", system, "-", study_number, "]", sep =""),
         shape_param_name = paste("shape[", clean_sci_name, "-", intensity, "-", system, "-", study_number, "]", sep =""),
         rate_param_name = paste("rate[", clean_sci_name, "-", intensity, "-", system, "-", study_number, "]", sep = ""),
         yrep_param_name = paste("y_rep[", clean_sci_name, "-", intensity, "-", system, "-", study_number, "]", sep = ""))
#write.csv(sp_index_key, file.path(outdir, "sp_info.csv"), row.names = FALSE)

# Clean up beta parameter names
beta_index_key <- data.frame(beta_name = colnames(X_scaled) %>% str_remove(pattern = "taxa|intensity|system|\\(|\\)")) %>%
  mutate(beta_param_name = paste("beta[", beta_name, "]", sep = ""))

# Replace param names
names(fit_simple_clean)[grep(names(fit_simple_clean), pattern = "beta")] <- beta_index_key$beta_param_name
names(fit_simple_clean)[grep(names(fit_simple_clean), pattern = "mu")] <- study_index_key$mu_param_name
names(fit_simple_clean)[grep(names(fit_simple_clean), pattern = "shape")] <- study_index_key$shape_param_name
names(fit_simple_clean)[grep(names(fit_simple_clean), pattern = "rate")] <- study_index_key$rate_param_name
names(fit_simple_clean)[grep(names(fit_simple_clean), pattern = "y_rep")] <- study_index_key$yrep_param_name

# FIX IT - continue parsing through outputs:

######################################################################################################
# Model 1.2: Using brms

library(brms)

# assemble dataframe for brms (no missing data)
# don't need to include intercept as part of design matrix with brms
#brms_data <- data.frame(cbind(y, X[,-1]))
brms_data <- data.frame(cbind(y, X_scaled[,-1]))

# No priors:
brms_gamma <- brm(y ~ ., data = brms_data, family = Gamma(link = "log"), seed = "11729", cores = 4)
# equivalent to: 
brms_gamma <- brm(y ~ 1 + taxasalmon.char + taxatilapia + taxatrout + intensitySemi.intensive + systemopen + systemsemi.open,
                  data = brms_data, family = Gamma(link = "log"), seed = "11729")
# DEFAULT PRIORS in brm are crap

get_prior(y ~ ., data = brms_data, family = Gamma(link = "log"))

# Set brm priors to be the same as those used by rstanarm (see next section)
brms_gamma <- brm(y ~ ., data = brms_data, family = Gamma(link = "log"), seed = "11729", cores = 4,
                  set_prior("normal(0,5)", class = "b"), set_prior("normal(0,2.5", class = "Intercept"), set_prior("exponential(rate = 1)", class = "shape"))

# equivalent:
brms_gamma <- brm(y ~ 1 + taxasalmon.char + taxatilapia + taxatrout + intensitySemi.intensive + systemopen + systemsemi.open,
                  data = brms_data, family = Gamma(link = "log"), seed = "11729", cores = 4,
                  set_prior("normal(0,5)", class = "b"), set_prior("normal(0,2.5", class = "Intercept"), set_prior("exponential(rate = 1)", class = "shape"))


# Posterior predictive checks
pp_check(brms_gamma, nsamples = 50) # density plots
pp_check(brms_gamma, type = "error_hist", nsamples = 5)
pp_check(brms_gamma, type = "scatter_avg", nsamples = 100)
pp_check(brms_gamma, type = "stat_2d")
pp_check(brms_gamma, type = "rootogram") 
# error message is fixed in the development version of brms: https://discourse.mc-stan.org/t/error-using-pp-check-using-truncated-response-variable-in-brms/7555
pp_check(brms_gamma, type = "loo_pit")
pp_check(brms_gamma, type = "xyz") # Gives an overview of all valid types
######################################################################################################
# Model 1.3: Use rstanarm for guidance on priors
# Note: "adjusted" priors are the priors implemented by rstanarm to help stabilize computation
# See: https://cran.r-project.org/web/packages/rstanarm/vignettes/priors.html

rstanarm_data <- data.frame(cbind(y, X_scaled[,-1]))

# Check the priors used by rstanarm for gamma regression model
rstanarm_gamma <- stan_glm(y ~ ., data = rstanarm_data, family = Gamma(link = "log"), seed = "11729", cores = 4)

pp_check(rstanarm_gamma, plotfun = "stat", stat = "mean")
pp_check(rstanarm_gamma, plotfun = "dens_overlay")

# Can extract the stan code with the following:
stancode <- rstan::get_stancode(rstanarm_gamma$stanfit)
cat(stancode)

# Get the (specified or default) priors
prior_summary(rstanarm_gamma)

# Check the priors used by rstanarm for the gaussian portion of the model
# Remove constant variables - taxafreshwater.crustacean, marine.shrimp and tuna
rstanarm_gaus <- stan_glm(intensitySemi.intensive ~ 1 + taxasalmon.char + taxatilapia + taxatrout +
                            systemopen + systemsemi.open, data = rstanarm_data, family = gaussian(), seed = "11729", cores = 4)
prior_summary(rstanarm_gaus)
######################################################################################################
# Model 2: Model FCR, including NAs in the predictors, and predict missing FCR values

# clean up taxa_group_name
lca_categories <- lca_dat_clean %>%
  select(clean_sci_name, FCR, taxa_group_name, intensity = Intensity, system = Production_system_group) %>%
  mutate(taxa = case_when(taxa_group_name == "Cods, hakes, haddocks" ~ "cod",
                          taxa_group_name == "Common carp" ~ "com_carp",
                          taxa_group_name == "Crabs, sea-spiders" ~ "crab",
                          taxa_group_name == "Freshwater crustaceans" ~ "fresh_crust",
                          taxa_group_name == "Milkfish" ~ "milkfish",
                          taxa_group_name == "Miscellaneous diadromous fishes" ~ "misc_diad",
                          taxa_group_name == "Miscellaneous freshwater fishes" ~ "misc_fresh",
                          taxa_group_name == "Miscellaneous marine fishes" ~ "misc_marine",
                          taxa_group_name == "Mussels" ~ "mussel",
                          taxa_group_name == "Other carps, barbels and cyprinids" ~ "oth_carp",
                          taxa_group_name == "Salmon" ~ "salmon",
                          taxa_group_name == "Shrimps, prawns" ~ "shrimp",
                          taxa_group_name == "Tilapias and other cichlids" ~ "tilapia",
                          taxa_group_name == "Trout" ~ "trout",
                          taxa_group_name == "Tunas, bonitos, billfishes" ~ "tuna",
                          TRUE ~ "unassigned")) %>%
  select(FCR, clean_sci_name, taxa, intensity, system) %>%
  arrange(clean_sci_name, taxa, intensity, system) %>%
  filter(taxa %in% c("mussel")==FALSE) # Remove taxa that don't belong in FCR/feed analysis mussels


######################################################################################################

# Model 2.1
# Keep NA's in FCR (response variable), but ONLY if they have a complete set of predictors (no NAs)
# remove FCR == 0 (species that aren't fed)
lca_complete_predictors <- lca_categories %>%
  filter(FCR != 0 | is.na(FCR))  %>% # Have to explicitly include is.na(FCR) otherwise NA's get dropped by FCR != 0
  filter(is.na(intensity)==FALSE & is.na(system)==FALSE) %>% # complete predictors - i.e., both intensity AND system are non-NA
  arrange(clean_sci_name, intensity, system)

# Set data:
y = lca_complete_predictors$FCR

# Create model matrix, but keep the NA's
# First change default options for handling missing data
options(na.action='na.pass')
X <- model.matrix(object = ~taxa + intensity + system, 
                  data = lca_complete_predictors %>% select(taxa, intensity, system)) 
# Return option back to the default
options(na.action='na.omit')

# Center all non-intercept variables and scale by 2 standard deviations (ignoring NAs)
covars.sd<-apply(X[,-1], MARGIN=2, FUN=sd, na.rm=TRUE)
# First change default options for handling missing data
options(na.action='na.pass')
X_scaled <- scale(X[,-1], center=TRUE, scale=2*covars.sd) 
# Return option back to the default
options(na.action='na.omit')

# Note: don't need to create an intercept column for brms
brms_data <- data.frame(cbind(y, X_scaled))

# Use mi() to explicitly include NAs among the predictors
# Which predictors have missing data:
# X_where_na <- apply(brms_data, MARGIN = 2, is.na)
# colSums(X_where_na)
# including mi() on the right hand side of the formula means you want to model NA's in the predictors, i.e., not drop them which means you have to specify how they should be modeled to (next section)
# y_brms <- brmsformula(y | mi() ~ 1 + taxacom_carp + taxacrab + taxafresh_crust + taxamilkfish + taxamisc_diad + taxamisc_fresh + taxamisc_marine + 
#                         taxamussel + taxaoth_carp + taxasalmon + taxashrimp + taxatilapia + taxatrout + taxatuna + 
#                         mi(intensityintensive) + mi(intensitysemi) + mi(systemopen) + mi(systemsemi), family = Gamma("log"))

y_brms <- brmsformula(y | mi() ~ ., family = Gamma("log"))

all_priors <- c(set_prior("normal(0,5)", class = "b"), 
                set_prior("normal(0,2.5)", class = "Intercept"), 
                set_prior("exponential(1)", class = "shape"))

# Model converges after increasing the adapt_delta and iterations from default values
fit_complete_predictors <- brm(y_brms, data = brms_data,
                   prior = all_priors, cores = 4, seed = "11729", control = list(adapt_delta = 0.99), iter = 10000)

summary(fit_complete_predictors) # Number of observations = 103 - i.e., default is to drop NA's
prior_summary(fit_complete_predictors)

plot(conditional_effects(fit_complete_predictors), ask = FALSE)

# Posterior predictive checks
pp_check(fit_complete_predictors, nsamples = 50) # density plots
pp_check(fit_complete_predictors, type = "error_hist", nsamples = 5)
pp_check(fit_complete_predictors, type = "scatter_avg", nsamples = 100)
pp_check(fit_complete_predictors, type = "stat_2d")
pp_check(fit_complete_predictors, type = "rootogram") 
# error message is fixed in the development version of brms: https://discourse.mc-stan.org/t/error-using-pp-check-using-truncated-response-variable-in-brms/7555
pp_check(fit_complete_predictors, type = "loo_pit")
pp_check(fit_complete_predictors, type = "xyz") # Gives an overview of all valid

# Get the predicted responses:
predict(fit_complete_predictors)

# Use the fit_complete_predictors model to predict FCRs with incomplete predictors (setting NA to 0)

lca_incomplete_predictors <- lca_categories %>%
  filter(FCR != 0 | is.na(FCR))  %>% # Have to explicitly include is.na(FCR) otherwise NA's get dropped by FCR != 0
  filter(is.na(intensity)==TRUE | is.na(system)==TRUE) %>% # incomplete predictors - i.e., either intensity OR system are NA
  arrange(clean_sci_name, intensity, system)


######################################################################################################
# Model 2.2: Model FCR's with NA while imputing NA's in the predictors

lca_with_na <- lca_categories %>%
  filter(FCR != 0 | is.na(FCR))  %>% # Have to explicitly include is.na(FCR) otherwise NA's get dropped by FCR != 0
  arrange(clean_sci_name, taxa, intensity, system)

# Set data:
y = lca_with_na$FCR

# Create model matrix, but keep the NA's
# First change default options for handling missing data
options(na.action='na.pass')
X <- model.matrix(object = ~taxa + intensity + system, 
                  data = lca_with_na %>% select(taxa, intensity, system)) 
# Return option back to the default
options(na.action='na.omit')

# Center all non-intercept variables and scale by 2 standard deviations (ignoring NAs)
covars.sd<-apply(X[,-1], MARGIN=2, FUN=sd, na.rm=TRUE)
# First change default options for handling missing data
options(na.action='na.pass')
X_scaled <- scale(X[,-1], center=TRUE, scale=2*covars.sd) 
# Return option back to the default
options(na.action='na.omit')

# Note: don't need to create an intercept column for brms
brms_data <- data.frame(cbind(y, X_scaled))

######################################################################################################
# Model 2.2a: model the missing predictors with just an intercept (give it the mean value across all values): https://cran.r-project.org/web/packages/brms/vignettes/brms_missings.html

# Which predictors have missing data:
X_where_na <- apply(brms_data, MARGIN = 2, is.na)
colSums(X_where_na)

y_brms <- brmsformula(y | mi() ~ 1 + taxacom_carp + taxacrab + taxafresh_crust + taxamilkfish + taxamisc_diad + taxamisc_fresh + 
                        taxamisc_marine + taxaoth_carp + taxasalmon + taxashrimp + taxatilapia + taxatrout + taxatuna + 
                        mi(intensityintensive) + mi(intensitysemi) + mi(systemopen) + mi(systemsemi), family = Gamma("log"))

# Check the priors used by rstanarm for the gaussian portion of the model
# Note: this automatically drops missing data
# library(rstanarm)
# rstanarm_gaus <- stan_glm(systemsemi ~ 1, data = brms_data, family = gaussian(), seed = "11729", cores = 4)
# summary(rstanarm_gaus)
# prior_summary(rstanarm_gaus)

# Notes: Same rstanarm priors for all variables: intensityintensive, intensitysemi, systemopen, systemsemi:
# intercept ~ normal(0, 1.2)
# sigma ~ exponential(2)

intensity_intensive_mi <- brmsformula(intensityintensive | mi() ~ 1,
                                      family = gaussian())

intensity_semi_mi <- brmsformula(intensitysemi | mi() ~ 1,
                                 family = gaussian())

system_open_mi  <- brmsformula(systemopen | mi() ~ 1,
                               family = gaussian())

system_semi_mi <- brmsformula(systemsemi | mi() ~ 1,
                              family = gaussian())

# Use "resp = <response_variable>" to specify different priors for different response variables
all_priors <- c(set_prior("normal(0,5)", class = "b", resp = "y"), # priors for y response variables
                set_prior("normal(0,2.5)", class = "Intercept", resp = "y"), 
                set_prior("exponential(1)", class = "shape", resp = "y"),
                set_prior("normal(0,1.2)", class = "Intercept", resp = c("intensityintensive", "intensitysemi", "systemopen", "systemsemi")),
                set_prior("exponential(2)", class = "sigma", resp = c("intensityintensive", "intensitysemi", "systemopen", "systemsemi")))

# Model converges after increasing the adapt_delta and iterations from default values
# Rule of thumb: bulk and tail effective sample sizes should be 100 x number of chains (i.e., at least 400)
# increasing max_treedepth is more about efficiency (instead of validity)
# See: https://mc-stan.org/misc/warnings.html
fit_with_na <- brm(y_brms + intensity_intensive_mi + intensity_semi_mi + system_open_mi + system_semi_mi + set_rescor(FALSE), data = brms_data,
                   prior = all_priors, cores = 4, seed = "11729", iter = 20000, control = list(adapt_delta = 0.99, max_treedepth = 15))
# Bulk ESS for y_intensityintensive still < 400; try increasing iterto 50000?
stancode(fit_with_na)

# Get posterior predictions of responses and missing data
pp_with_na <- predict(fit_with_na)

######################################################################################################
# Model 2.2b: Model the missing predictors with just an intercept (same as 2.2a) but first combine rare taxa into an "other taxa" category

lca_with_na <- lca_categories %>%
  filter(FCR != 0 | is.na(FCR))  %>% # Have to explicitly include is.na(FCR) otherwise NA's get dropped by FCR != 0
  arrange(clean_sci_name, taxa, intensity, system)

table(lca_with_na$taxa)

lca_simplified_taxa <- lca_with_na %>%
  mutate(taxa = if_else(taxa %in% c("cod", "com_carp", "crab", "fresh_crust", "milkfish", "oth_carp", "tuna"), true = "oth_taxa", false = taxa)) # use cutoff of N=5

table(lca_simplified_taxa$taxa)

# Set data:
y = lca_simplified_taxa$FCR

# Create model matrix, but keep the NA's
# First change default options for handling missing data
options(na.action='na.pass')
X <- model.matrix(object = ~taxa + intensity + system, 
                  data = lca_simplified_taxa %>% select(taxa, intensity, system)) 
# Return option back to the default
options(na.action='na.omit')

# Center all non-intercept variables and scale by 2 standard deviations (ignoring NAs)
covars.sd<-apply(X[,-1], MARGIN=2, FUN=sd, na.rm=TRUE)
# First change default options for handling missing data
options(na.action='na.pass')
X_scaled <- scale(X[,-1], center=TRUE, scale=2*covars.sd) 
# Return option back to the default
options(na.action='na.omit')

# Note: don't need to create an intercept column for brms
brms_data <- data.frame(cbind(y, X_scaled))

# Which predictors have missing data:
X_where_na <- apply(brms_data, MARGIN = 2, is.na)
colSums(X_where_na)

y_brms <- brmsformula(y | mi() ~ 1 + taxamisc_fresh + 
                        taxamisc_marine + taxaoth_taxa + taxasalmon + taxashrimp + taxatilapia + taxatrout + 
                        mi(intensityintensive) + mi(intensitysemi) + mi(systemopen) + mi(systemsemi), family = Gamma("log"))

# Check the priors used by rstanarm for the gaussian portion of the model
# Note: this automatically drops missing data
# library(rstanarm)
# rstanarm_gaus <- stan_glm(systemsemi ~ 1, data = brms_data, family = gaussian(), seed = "11729", cores = 4)
# summary(rstanarm_gaus)
# prior_summary(rstanarm_gaus)

# Notes: Same rstanarm priors for all variables: intensityintensive, intensitysemi, systemopen, systemsemi:
# intercept ~ normal(0, 1.2)
# sigma ~ exponential(2)

intensity_intensive_mi <- brmsformula(intensityintensive | mi() ~ 1,
                                      family = gaussian())

intensity_semi_mi <- brmsformula(intensitysemi | mi() ~ 1,
                                 family = gaussian())

system_open_mi  <- brmsformula(systemopen | mi() ~ 1,
                               family = gaussian())

system_semi_mi <- brmsformula(systemsemi | mi() ~ 1,
                              family = gaussian())

# Use "resp = <response_variable>" to specify different priors for different response variables
all_priors <- c(set_prior("normal(0,5)", class = "b", resp = "y"), # priors for y response variables
                set_prior("normal(0,2.5)", class = "Intercept", resp = "y"), 
                set_prior("exponential(1)", class = "shape", resp = "y"),
                set_prior("normal(0,1.2)", class = "Intercept", resp = c("intensityintensive", "intensitysemi", "systemopen", "systemsemi")),
                set_prior("exponential(2)", class = "sigma", resp = c("intensityintensive", "intensitysemi", "systemopen", "systemsemi")))

# Model converges after increasing the adapt_delta and iterations from default values
# Rule of thumb: bulk and tail effective sample sizes should be 100 x number of chains (i.e., at least 400)
# increasing max_treedepth is more about efficiency (instead of validity)
# See: https://mc-stan.org/misc/warnings.html
fit_simplified_taxa <- brm(y_brms + intensity_intensive_mi + intensity_semi_mi + system_open_mi + system_semi_mi + set_rescor(FALSE), data = brms_data,
                   prior = all_priors, cores = 4, seed = "11729", iter = 20000, control = list(adapt_delta = 0.99))
