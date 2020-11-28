# Bayesian dirichlet regression to estimate variables based on taxa-group, intensity, and production system

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

# select relevant data columns and arrange by categorical info
# Need FCR to identify which species aren't fed (FCR == 0)
lca_categories <- lca_dat_clean %>%
  select(FCR, contains("new"), clean_sci_name, taxa, intensity = Intensity, system = Production_system_group) %>%
  arrange(clean_sci_name, taxa, intensity, system) %>%
  filter(taxa %in% c("mussel")==FALSE) # Remove taxa that don't belong in FCR/feed analysis mussels

######################################################################################################
# Model 1: No Missing data
# Remove FCR == 0 (species that aren't fed)
lca_no_na <- lca_categories %>%
  filter(FCR != 0 | is.na(FCR))  %>% # Have to explicitly include is.na(FCR) otherwise NA's get dropped by FCR != 0
  filter(is.na(intensity)==FALSE & is.na(system)==FALSE) %>% # complete predictors - i.e., both intensity AND system are non-NA
  filter(is.na(feed_soy_new)==FALSE) %>%
  rename(feed_soy = feed_soy_new,
         feed_crops = feed_crops_new,
         feed_fmfo = feed_fmfo_new,
         feed_animal = feed_animal_new)

# Set data for model:

# Create model matrix for taxa info, then center and scale
options(na.action='na.pass') # First change default options for handling missing data
X_taxa <- model.matrix(object = ~ 1 + taxa, 
                  data = lca_no_na %>% select(taxa)) 
options(na.action='na.omit') # Return option back to the default


taxa_sd <- apply(X_taxa[,-1], MARGIN=2, FUN=sd, na.rm=TRUE) # Center all non-intercept variables and scale by 2 standard deviations (ignoring NAs)
options(na.action='na.pass') # First change default options for handling missing data
X_taxa_scaled <- scale(X_taxa[,-1], center=TRUE, scale=2*taxa_sd)
options(na.action='na.omit') # Return option back to the default

# Format intensity and system as ordinal variables, then center and scale
X_ordinal <- lca_no_na %>%
  mutate(intensity = factor(intensity, levels = c("extensive", "semi", "intensive"))) %>% # set order of factors (low = extensive, high = intensive)
  mutate(system = factor(system, levels = c("open", "semi", "closed"))) %>% # set order of factors (low = open, high = closed)
  mutate(intensity = as.numeric(intensity)) %>%
  mutate(system = as.numeric(system)) %>%
  select(intensity, system) %>%
  as.matrix()
ordinal_sd<-apply(X_ordinal, MARGIN=2, FUN=sd, na.rm=TRUE) # Center all non-intercept variables and scale by 2 standard deviations (ignoring NAs)

options(na.action='na.pass') # First change default options for handling missing data
X_ordinal_scaled <- scale(X_ordinal, center=TRUE, scale=2*ordinal_sd)
options(na.action='na.omit') # Return option back to the default

# Create dataframe for brms and rename feed variables
brms_data <- data.frame(cbind(lca_no_na %>% select(contains("feed")), X_taxa_scaled, X_ordinal_scaled)) 

# Response variable must be a matrix, create function bind since cbind within the brm function is reserved for specifying multivariate models
bind <- function(...) cbind(...)

# When BOTH intensity and system are included, model does not converge, 
# CHOOSE ONE:
# For TAXA + INTENSITY 
y_brms <- brmsformula(bind(feed_soy, feed_crops, feed_fmfo, feed_animal) ~ taxamisc_marine + taxasalmon + taxatilapia + taxatrout + intensity, family = dirichlet())

# For TAXA + SYSTEM
y_brms <- brmsformula(bind(feed_soy, feed_crops, feed_fmfo, feed_animal) ~ taxamisc_marine + taxasalmon + taxatilapia + taxatrout + system, family = dirichlet())

# Model converges after increasing the adapt_delta and iterations from default values
fit_no_na <- brm(y_brms, data = brms_data,
                 cores = 4, seed = "11729", iter = 10000)

summary(fit_no_na) 

######################################################################################################
# Use model to predict NA feeds for studies 

# CHOOSE ONE: with complete taxa + INTENSITY DATA
lca_complete_predictors <- lca_categories %>%
  filter(FCR != 0 | is.na(FCR))  %>% # Have to explicitly include is.na(FCR) otherwise NA's get dropped by FCR != 0
  filter(is.na(intensity)==FALSE) %>% # complete predictors - i.e., both intensity AND system are non-NA
  filter(is.na(feed_soy_new))

# CHOOSE ONE: with complete taxa + SYSTEM DATA
lca_complete_predictors <- lca_categories %>%
  filter(FCR != 0 | is.na(FCR))  %>% # Have to explicitly include is.na(FCR) otherwise NA's get dropped by FCR != 0
  filter(is.na(system)==FALSE) %>% # complete predictors - i.e., both intensity AND system are non-NA
  filter(is.na(feed_soy_new))

# PROBLEM: lca_complete predictors has more taxa than originally model:
taxa_not_modeled <- setdiff(unique(lca_complete_predictors$taxa), unique(lca_no_na$taxa)) # these taxa were never modeled so they can't be predicted below

# DROP THESE FOR NOW:
lca_complete_predictors <- lca_complete_predictors %>%
  filter(taxa %in% taxa_not_modeled == FALSE)

# If lca_complete predictors has fewer taxa than originally model, need to convert to factor and expand/assign levels manually - having trouble automating this
# levels(lca_complete_predictors$taxa) <- list(misc_diad = "misc_diad", misc_marine = "misc_marine", salmon = "salmon", tilapia = "tilapia", trout = "trout")

# Create NEW taxa model matrix for the studies whose feeds need to be predicted
# Taxa categories:
options(na.action='na.pass') # First change default options for handling missing data
X_taxa_new <- model.matrix(object = ~ 1 + taxa, 
                       data = lca_complete_predictors %>% select(taxa)) 
options(na.action='na.omit') # Return option back to the default

# Center and Scale: BUT now center by the mean of the original modeled dataset above AND scale by the same 2*SD calculated from the original, modeled dataset above
options(na.action='na.pass') # First change default options for handling missing data
X_taxa_new_scaled <- scale(X_taxa_new[,-1], center=apply(X_taxa[,-1], MARGIN = 2, FUN = mean), scale=2*taxa_sd)
options(na.action='na.omit') # Return option back to the default

# System and Intensity variables:
# Format intensity and system as ordinal variables, then center and scale
X_ordinal_new <- lca_complete_predictors %>%
  mutate(intensity = factor(intensity, levels = c("extensive", "semi", "intensive"))) %>% # set order of factors (low = extensive, high = intensive)
  mutate(system = factor(system, levels = c("open", "semi", "closed"))) %>% # set order of factors (low = open, high = closed)
  mutate(intensity = as.numeric(intensity)) %>%
  mutate(system = as.numeric(system)) %>%
  select(intensity, system) %>%
  as.matrix()

# Center and Scale: BUT now center by the mean of the original modeled dataset above AND scale by the same 2*SD calculated from the original, modeled dataset above
options(na.action='na.pass') # First change default options for handling missing data
X_ordinal_new_scaled <- scale(X_ordinal_new, center=apply(X_ordinal, MARGIN = 2, FUN = mean), scale=2*ordinal_sd)
options(na.action='na.omit') # Return option back to the default

# Create dataframe for brms and rename feed variables
brms_new_data <- data.frame(cbind(X_taxa_new_scaled, X_ordinal_new_scaled)) 

# Make predictions
#predicted_feed_dat <- predict(fit_no_na, newdata = brms_new_data)
# Use tidybayes instead:
predicted_feed_dat <- add_predicted_draws(newdata = brms_new_data, model = fit_no_na)
