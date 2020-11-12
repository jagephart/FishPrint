# Bayesian regression to estimate variables based on taxa-group, intensity, and production system

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
# Model 1.3: Using rstanarm

library(rstanarm)

rstanarm_data <- data.frame(cbind(y, X_scaled[,-1]))
rstanarm_gamma <- stan_glm(y ~ ., data = rstanarm_data, family = Gamma(link = "log"), seed = "11729", cores = 4)

pp_check(rstanarm_gamma, plotfun = "stat", stat = "mean")
pp_check(rstanarm_gamma, plotfun = "dens_overlay")

# Can extract the stan code with the following:
stancode <- rstan::get_stancode(rstanarm_gamma$stanfit)
cat(stancode)

# Get the (specified or default) priors
prior_summary(rstanarm_gamma)

# Note: adjusted priors are the priors implemented by rstanarm to help stabilize computation
# See: https://cran.r-project.org/web/packages/rstanarm/vignettes/priors.html

######################################################################################################
# Model 2: Include NAs in the predictors

# But still remove NA's in FCR (response variable)
lca_dat_with_na <- lca_dat_clean %>%
  select(clean_sci_name, FCR, taxa = taxa_group_name, intensity = Intensity, system = Production_system_group) %>%
  filter(is.na(FCR) == FALSE) %>% 
  filter(FCR != 0) %>%
  # drop_na() %>%
  arrange(clean_sci_name, intensity, system) %>%
  mutate(taxa = as.factor(taxa),
         tx = as.numeric(taxa),
         intensity = as.factor(intensity),
         its = as.numeric(intensity),
         system = as.factor(system),
         ps = as.numeric(system)) %>%
  select(clean_sci_name, FCR, taxa, tx, intensity, its, system, ps)

# Set data:
N = nrow(lca_dat_with_na)
y = lca_dat_with_na$FCR
K = ncol(X) # number of predictors

# Create model matrix, but keep the NA's
# First change default options for handling missing data
options(na.action='na.pass')
X <- model.matrix(object = ~taxa + intensity + system, 
                  data = lca_dat_with_na %>% select(taxa, intensity, system)) 
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

# Two approaches to Missing data imputation with brms: https://cran.r-project.org/web/packages/brms/vignettes/brms_missings.html

# Option 1: Imputation before model fitting (more computationally intensive, but might be OK for a simple model)
# FIX IT - try option 1

# Option 2: Imputation while model fitting:

# Which predictors have missing data:
X_where_na <- apply(X_scaled, MARGIN = 2, is.na)
colSums(X_where_na)

# For reference, here's the no missing data formula with priors
brms_gamma <- brm(y ~ 1 + taxasalmon.char + taxatilapia + taxatrout + intensitySemi.intensive + systemopen + systemsemi.open,
                  data = brms_data, family = Gamma(link = "log"), seed = "11729", cores = 4,
                  set_prior("normal(0,5)", class = "b"), set_prior("normal(0,2.5", class = "Intercept"), set_prior("exponential(rate = 1)", class = "shape"))

# Set up brms model formula that explains where there is missing data
missing_dat_form <- brmsformula(y ~ 1 + taxasalmon.char + taxatilapia + taxatrout + mi(intensitySemi.intensive) + mi(systemopen) + mi(systemsemi.open)) +
  brmsformula(intensitySemi.intensive | mi() ~ .) +
  brmsformula(systemopen | mi() ~ .) +
  brmsformula(systemsemi.open | mi() ~ .) +
  set_rescor(FALSE) # FALSE - do not model the residual correlations between the response variables
# FIX IT try set_rescor(TRUE)

# FIX IT: LEFT OFF HERE - model not running with missing data
brms_gamma_with_na <- brm(missing_dat_form, data = brms_data, family = Gamma(link = "log"), seed = "11729", cores = 4,
                          set_prior("normal(0,5)", class = "b"), set_prior("normal(0,2.5", class = "Intercept"), set_prior("exponential(rate = 1)", class = "shape"))



# NEXT: try keeping data where FCR is NA