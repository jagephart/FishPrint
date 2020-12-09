# Calculate (non-feed associated) greenhouse gas footprint
# (Electricity use * country-specific GHG of electricity) + (Diesel * GHG of diesel) + (Petrol * GHG of petrol) + (Natural gas * GHG of natural gas)

# Reminder for new data: entries with data for some columns should have no blanks in other columns (these should be filled in as zeroes)

# Step 0: Run process_data_for_analysis.R, then clear environment other than:
rm(list=ls()[!(ls() %in% c("lca_dat_clean_groups", "datadir", "outdir"))])

# FIX IT - Decide: adjust zeroes upwards? or drop zeroes? See first block of code in each model


# Get model-specific data:
# SELECT STUDY ID COLUMN - use this for rejoining outputs from multiple regression models back together
# Select relevant data columns and arrange by categorical info
ghg_model_dat_categories <- lca_dat_clean_groups %>%
  select(study_id, electric = Electricity_kwh, diesel = Diesel_L, petrol = Petrol_L, natgas = NaturalGas_L, clean_sci_name, taxa, intensity = Intensity, system = Production_system_group) %>%
  arrange(clean_sci_name, taxa, intensity, system)

######################################################################################################
# Step 1: Model electricity as taxa + intensity + system

electric_no_na <- ghg_model_dat_categories %>%
  filter(is.na(electric)==FALSE)  %>% # Drop NAs (Keep zeroes)
  mutate(electric = if_else(electric == 0, true = min(electric[electric!=0]), false = electric)) %>% # Not modeling the zeroes, option 1: adjust these to the minimum value
  #filter(electric != 0) %>% # Not modeling the zeroes, option 2: drop all zeroes
  filter(is.na(intensity)==FALSE & is.na(system)==FALSE) %>% # complete predictors - i.e., both intensity AND system are non-NA
  select(study_id, electric, clean_sci_name, taxa, intensity, system)

# Create model matrix for taxa info, then center and scale
X_taxa <- model.matrix(object = ~ 1 + taxa, 
                       data = electric_no_na %>% select(taxa)) 

taxa_sd <- apply(X_taxa[,-1], MARGIN=2, FUN=sd, na.rm=TRUE) # Center all non-intercept variables and scale by 2 standard deviations (ignoring NAs)
X_taxa_scaled <- scale(X_taxa[,-1], center=TRUE, scale=2*taxa_sd)

# Format intensity and system as ordinal variables, then center and scale
X_ordinal <- electric_no_na %>%
  mutate(intensity = factor(intensity, levels = c("extensive", "semi", "intensive"))) %>% # set order of factors (low = extensive, high = intensive)
  mutate(system = factor(system, levels = c("open", "semi", "closed"))) %>% # set order of factors (low = open, high = closed)
  mutate(intensity = as.numeric(intensity)) %>%
  mutate(system = as.numeric(system)) %>%
  select(intensity, system) %>%
  as.matrix()
ordinal_sd<-apply(X_ordinal, MARGIN=2, FUN=sd, na.rm=TRUE) # Center all non-intercept variables and scale by 2 standard deviations (ignoring NAs)
X_ordinal_scaled <- scale(X_ordinal, center=TRUE, scale=2*ordinal_sd)

# Create dataframe for brms and rename feed variables
electric_brms_data <- data.frame(y = electric_no_na$electric, X_taxa_scaled, X_ordinal_scaled)

names(electric_brms_data)

# Set model formula
electric_brms <- brmsformula(y ~ 1 + ., family = Gamma("log"))
# Equivalent to:
# electric_brms <- brmsformula(y ~ 1 + taxatilapia + taxatrout +
#                            intensity + system, family = Gamma("log"))


# Use "resp = <response_variable>" to specify different priors for different response variables
all_priors <- c(set_prior("normal(0,5)", class = "b"), # priors for y response variables
                set_prior("normal(0,2.5)", class = "Intercept"), 
                set_prior("exponential(1)", class = "shape"))

# Model converges after increasing the adapt_delta and iterations from default values
# Rule of thumb: bulk and tail effective sample sizes should be 100 x number of chains (i.e., at least 400)
# increasing max_treedepth is more about efficiency (instead of validity)
# See: https://mc-stan.org/misc/warnings.html
fit_electric_no_na <- brm(electric_brms, data = electric_brms_data,
                       prior = all_priors, cores = 4, seed = "11729", iter = 20000, control = list(adapt_delta = 0.99))

# Get stan code
#stancode(fit_electric_no_na)

######################################################################################################
# Use model to predict NAs for studies with complete set of predictors
# Both intensity AND system are non-NA
electric_complete_predictors <- ghg_model_dat_categories %>%
  filter(is.na(electric)) %>% 
  filter(is.na(intensity)==FALSE & is.na(system)== FALSE) 

# PROBLEM: lca_complete predictors has more taxa than originally model:
taxa_not_modeled <- setdiff(unique(electric_complete_predictors$taxa), unique(electric_no_na$taxa)) # these taxa were never modeled so they can't be predicted below

# DROP THESE FOR NOW:
electric_complete_predictors <- electric_complete_predictors %>%
  filter(taxa %in% taxa_not_modeled == FALSE)

# Now check the other way, which taxa were in the original model but not a part of the data that needs to be predicted:
setdiff(unique(electric_no_na$taxa), unique(electric_complete_predictors$taxa))

# If original model has taxa that are not part of electric_complete_predictors, 
# Use list of unique taxa in original model and use this to expand/assign levels manually - having trouble automating this
# Include all levels here, but remember the first level won't show up in design matrix - instead, it's part of the "contrasts"
# sort(unique(electric_no_na$taxa))
# electric_complete_predictors <- electric_complete_predictors %>%
#   mutate(taxa = as.factor(taxa))
# levels(electric_complete_predictors$taxa) <- list(fresh_crust = "fresh_crust", misc_fresh = "misc_fresh", misc_marine = "misc_marine", oth_carp = "oth_carp", tilapia = "tilapia")

# Create NEW taxa model matrix for the studies to be predicted
# Taxa categories:
X_taxa_new <- model.matrix(object = ~ 1 + taxa, 
                           data = electric_complete_predictors %>% select(taxa)) 

# Center and Scale: BUT now center by the mean of the original modeled dataset above AND scale by the same 2*SD calculated from the original, modeled dataset above
X_taxa_new_scaled <- scale(X_taxa_new[,-1], center=apply(X_taxa[,-1], MARGIN = 2, FUN = mean), scale=2*taxa_sd)

# System and Intensity variables:
# Format intensity and system as ordinal variables, then center and scale
X_ordinal_new <- electric_complete_predictors %>%
  mutate(intensity = factor(intensity, levels = c("extensive", "semi", "intensive"))) %>% # set order of factors (low = extensive, high = intensive)
  mutate(system = factor(system, levels = c("open", "semi", "closed"))) %>% # set order of factors (low = open, high = closed)
  mutate(intensity = as.numeric(intensity)) %>%
  mutate(system = as.numeric(system)) %>%
  select(intensity, system) %>%
  as.matrix()

# Center and Scale: BUT now center by the mean of the original modeled dataset above AND scale by the same 2*SD calculated from the original, modeled dataset above
X_ordinal_new_scaled <- scale(X_ordinal_new, center=apply(X_ordinal, MARGIN = 2, FUN = mean), scale=2*ordinal_sd)

# Create dataframe for brms and rename feed variables
brms_new_electric_dat <- data.frame(cbind(X_taxa_new_scaled, X_ordinal_new_scaled)) 

# Make predictions
#predicted_electric_dat <- predict(fit_no_na, newdata = brms_new_electric_data)
# Use tidybayes instead:
predicted_electric_dat <- add_predicted_draws(newdata = brms_new_electric_dat, model = fit_electric_no_na)

# Get point and interval estimates from predicted data
# Select just the prediction columns
# Join these with the modeled data (electric_complete_predictors) to get metadata on taxa/intensity/syste,
electric_dat_intervals <- predicted_electric_dat %>%
  median_qi(.value = .prediction) %>% # Rename prediction to value
  ungroup() %>%
  select(contains("."))

# .row is equivalent to the row number in the modeled dataset (electric_complete_predictors) - create a join column for this
electric_metadat<- electric_complete_predictors %>%
  select(study_id, clean_sci_name, taxa, intensity, system) %>%
  mutate(.row = row_number())

electric_predictions <- electric_dat_intervals %>%
  left_join(electric_metadat, by = ".row") %>%
  rename(electric = .value)

######################################################################################################
# Model electricity with intensity + system (no taxa)

# Create dataframe for brms and rename feed variables
electric_brms_data_no_taxa <- data.frame(y = electric_no_na$electric, X_ordinal_scaled)

names(electric_brms_data_no_taxa)

# Set model formula
electric_brms_no_taxa <- brmsformula(y ~ 1 + ., family = Gamma("log"))

# Use "resp = <response_variable>" to specify different priors for different response variables
all_priors <- c(set_prior("normal(0,5)", class = "b"), # priors for y response variables
                set_prior("normal(0,2.5)", class = "Intercept"), 
                set_prior("exponential(1)", class = "shape"))

# Model converges after increasing the adapt_delta and iterations from default values
# Rule of thumb: bulk and tail effective sample sizes should be 100 x number of chains (i.e., at least 400)
# increasing max_treedepth is more about efficiency (instead of validity)
# See: https://mc-stan.org/misc/warnings.html
fit_electric_no_taxa <- brm(electric_brms_no_taxa, data = electric_brms_data_no_taxa,
                         prior = all_priors, cores = 4, seed = "11729", iter = 20000, control = list(adapt_delta = 0.99))

# Get stan code
#stancode(fit_electric_no_taxa)

######################################################################################################
# Use intensity + system model to predict electricity for taxa that were not part of the previous model
# Intensity must be non-NA
electric_complete_predictors <- ghg_model_dat_categories %>%
  filter(is.na(intensity)==FALSE & is.na(system)==FALSE) %>%
  filter(is.na(electric)) # Now, filter to just the NAs

taxa_not_modeled <- setdiff(unique(electric_complete_predictors$taxa), unique(electric_no_na$taxa)) # these taxa were never modeled so they can't be predicted below

# Only keep taxa that were not modeled previously
electric_complete_predictors_no_taxa <- electric_complete_predictors %>%
  filter(taxa %in% taxa_not_modeled) 

# Format intensity and system as ordinal variable, then center and scale
X_ordinal_no_taxa <- electric_complete_predictors_no_taxa %>%
  mutate(intensity = factor(intensity, levels = c("extensive", "semi", "intensive"))) %>% # set order of factors (low = extensive, high = intensive)
  mutate(intensity = as.numeric(intensity)) %>%
  mutate(system = factor(system, levels = c("open", "semi", "closed"))) %>% # set order of factors (low = open, high = closed)
  mutate(system = as.numeric(system)) %>%
  select(intensity, system) %>%
  as.matrix()

# Center and Scale: BUT now center by the mean of the original modeled dataset above AND scale by the same 2*SD calculated from the original, modeled dataset above
X_ordinal_scaled_no_taxa <- scale(X_ordinal_no_taxa, center=apply(X_ordinal, MARGIN = 2, FUN = mean), scale=2*ordinal_sd)

# Create dataframe for brms
brms_electric_dat_no_taxa <- data.frame(X_ordinal_scaled_no_taxa)

# Make predictions
#predicted_electric_dat <- predict(fit_no_na, newdata = brms_new_electric_data)
# Use tidybayes instead:
predicted_electric_dat_no_taxa <- add_predicted_draws(newdata = brms_electric_dat_no_taxa, model = fit_electric_no_taxa)

# Get point and interval estimates from predicted data
# Select just the prediction columns
# Join these with the modeled data (electric_complete_predictors_no_taxa) to get metadata on taxa/intensity/syste,
electric_dat_intervals_no_taxa <- predicted_electric_dat_no_taxa %>%
  median_qi(.value = .prediction) %>% # Rename prediction to value
  ungroup() %>%
  select(contains("."))

# .row is equivalent to the row number in the modeled dataset (electric_complete_predictors_no_taxa) - create a join column for this
electric_metadat_no_taxa <- electric_complete_predictors_no_taxa %>%
  select(study_id, clean_sci_name, taxa, intensity, system) %>%
  mutate(.row = row_number())

electric_predictions_no_taxa <- electric_dat_intervals_no_taxa %>%
  left_join(electric_metadat_no_taxa, by = ".row") %>%
  rename(electric = .value)

# Bind electric_no_na (data), electric_predictions, and electric_predictions_no_taxa
full_electric_dat <- electric_predictions %>%
  bind_rows(electric_no_na) %>%
  bind_rows(electric_predictions_no_taxa) %>%
  mutate(data_type = if_else(is.na(.point), true = "data", false = "prediction")) %>%
  arrange(taxa, intensity, system, clean_sci_name) %>%
  rownames_to_column() # Arrange by taxa first, then create dummy column for plotting 

# NEXT: use plot_brms_gamma_regression to produce figures for SI

######################################################################################################
# Step 2: Model diesel

diesel_no_na <- ghg_model_dat_categories %>%
  filter(is.na(diesel)==FALSE)  %>% # Drop NAs (Keep zeroes)
  mutate(diesel = if_else(diesel == 0, true = min(diesel[diesel!=0]), false = diesel)) %>% # Not modeling the zeroes, option 1: adjust these to the minimum value
  #filter(diesel != 0) %>% # Not modeling the zeroes, option 2: drop all zeroes
  filter(is.na(intensity)==FALSE & is.na(system)==FALSE) %>% # complete predictors - i.e., both intensity AND system are non-NA
  select(study_id, diesel, clean_sci_name, taxa, intensity, system)

# Create model matrix for taxa info, then center and scale
X_taxa <- model.matrix(object = ~ 1 + taxa, 
                       data = diesel_no_na %>% select(taxa)) 

taxa_sd <- apply(X_taxa[,-1], MARGIN=2, FUN=sd, na.rm=TRUE) # Center all non-intercept variables and scale by 2 standard deviations (ignoring NAs)
X_taxa_scaled <- scale(X_taxa[,-1], center=TRUE, scale=2*taxa_sd)

# Format intensity and system as ordinal variables, then center and scale
X_ordinal <- diesel_no_na %>%
  mutate(intensity = factor(intensity, levels = c("extensive", "semi", "intensive"))) %>% # set order of factors (low = extensive, high = intensive)
  mutate(system = factor(system, levels = c("open", "semi", "closed"))) %>% # set order of factors (low = open, high = closed)
  mutate(intensity = as.numeric(intensity)) %>%
  mutate(system = as.numeric(system)) %>%
  select(intensity, system) %>%
  as.matrix()
ordinal_sd<-apply(X_ordinal, MARGIN=2, FUN=sd, na.rm=TRUE) # Center all non-intercept variables and scale by 2 standard deviations (ignoring NAs)
X_ordinal_scaled <- scale(X_ordinal, center=TRUE, scale=2*ordinal_sd)

# Create dataframe for brms and rename feed variables
diesel_brms_data <- data.frame(y = diesel_no_na$diesel, X_taxa_scaled, X_ordinal_scaled)

names(diesel_brms_data)

# Set model formula
diesel_brms <- brmsformula(y ~ 1 + ., family = Gamma("log"))
# Equivalent to:
# diesel_brms <- brmsformula(y ~ 1 + taxatilapia + taxatrout +
#                            intensity + system, family = Gamma("log"))


# Use "resp = <response_variable>" to specify different priors for different response variables
all_priors <- c(set_prior("normal(0,5)", class = "b"), # priors for y response variables
                set_prior("normal(0,2.5)", class = "Intercept"), 
                set_prior("exponential(1)", class = "shape"))

# Model converges after increasing the adapt_delta and iterations from default values
# Rule of thumb: bulk and tail effective sample sizes should be 100 x number of chains (i.e., at least 400)
# increasing max_treedepth is more about efficiency (instead of validity)
# See: https://mc-stan.org/misc/warnings.html
fit_diesel_no_na <- brm(diesel_brms, data = diesel_brms_data,
                      prior = all_priors, cores = 4, seed = "11729", iter = 20000, control = list(adapt_delta = 0.99))

# Get stan code
#stancode(fit_diesel_no_na)

######################################################################################################
# Use model to predict NAs for studies with complete set of predictors
# Both intensity AND system are non-NA
diesel_complete_predictors <- ghg_model_dat_categories %>%
  filter(is.na(diesel)) %>% 
  filter(is.na(intensity)==FALSE & is.na(system)== FALSE) 

# PROBLEM: lca_complete predictors has more taxa than originally model:
taxa_not_modeled <- setdiff(unique(diesel_complete_predictors$taxa), unique(diesel_no_na$taxa)) # these taxa were never modeled so they can't be predicted below

# DROP THESE FOR NOW:
diesel_complete_predictors <- diesel_complete_predictors %>%
  filter(taxa %in% taxa_not_modeled == FALSE)

# Now check the other way, which taxa were in the original model but not a part of the data that needs to be predicted:
setdiff(unique(diesel_no_na$taxa), unique(diesel_complete_predictors$taxa))

# If original model has taxa that are not part of diesel_complete_predictors, 
# Use list of unique taxa in original model and use this to expand/assign levels manually - having trouble automating this
# Include all levels here, but remember the first level won't show up in design matrix - instead, it's part of the "contrasts"
# sort(unique(diesel_no_na$taxa))
# diesel_complete_predictors <- diesel_complete_predictors %>%
#   mutate(taxa = as.factor(taxa))
# levels(diesel_complete_predictors$taxa) <- list(fresh_crust = "fresh_crust", misc_fresh = "misc_fresh", misc_marine = "misc_marine", oth_carp = "oth_carp", tilapia = "tilapia")

# Create NEW taxa model matrix for the studies to be predicted
# Taxa categories:
X_taxa_new <- model.matrix(object = ~ 1 + taxa, 
                           data = diesel_complete_predictors %>% select(taxa)) 

# Center and Scale: BUT now center by the mean of the original modeled dataset above AND scale by the same 2*SD calculated from the original, modeled dataset above
X_taxa_new_scaled <- scale(X_taxa_new[,-1], center=apply(X_taxa[,-1], MARGIN = 2, FUN = mean), scale=2*taxa_sd)

# System and Intensity variables:
# Format intensity and system as ordinal variables, then center and scale
X_ordinal_new <- diesel_complete_predictors %>%
  mutate(intensity = factor(intensity, levels = c("extensive", "semi", "intensive"))) %>% # set order of factors (low = extensive, high = intensive)
  mutate(system = factor(system, levels = c("open", "semi", "closed"))) %>% # set order of factors (low = open, high = closed)
  mutate(intensity = as.numeric(intensity)) %>%
  mutate(system = as.numeric(system)) %>%
  select(intensity, system) %>%
  as.matrix()

# Center and Scale: BUT now center by the mean of the original modeled dataset above AND scale by the same 2*SD calculated from the original, modeled dataset above
X_ordinal_new_scaled <- scale(X_ordinal_new, center=apply(X_ordinal, MARGIN = 2, FUN = mean), scale=2*ordinal_sd)

# Create dataframe for brms and rename feed variables
brms_new_diesel_dat <- data.frame(cbind(X_taxa_new_scaled, X_ordinal_new_scaled)) 

# Make predictions
#predicted_diesel_dat <- predict(fit_no_na, newdata = brms_new_diesel_data)
# Use tidybayes instead:
predicted_diesel_dat <- add_predicted_draws(newdata = brms_new_diesel_dat, model = fit_diesel_no_na)

# Get point and interval estimates from predicted data
# Select just the prediction columns
# Join these with the modeled data (diesel_complete_predictors) to get metadata on taxa/intensity/syste,
diesel_dat_intervals <- predicted_diesel_dat %>%
  median_qi(.value = .prediction) %>% # Rename prediction to value
  ungroup() %>%
  select(contains("."))

# .row is equivalent to the row number in the modeled dataset (diesel_complete_predictors) - create a join column for this
diesel_metadat<- diesel_complete_predictors %>%
  select(study_id, clean_sci_name, taxa, intensity, system) %>%
  mutate(.row = row_number())

diesel_predictions <- diesel_dat_intervals %>%
  left_join(diesel_metadat, by = ".row") %>%
  rename(diesel = .value)

######################################################################################################
# Model diesel with intensity + system (no taxa)

# Create dataframe for brms and rename feed variables
diesel_brms_data_no_taxa <- data.frame(y = diesel_no_na$diesel, X_ordinal_scaled)

names(diesel_brms_data_no_taxa)

# Set model formula
diesel_brms_no_taxa <- brmsformula(y ~ 1 + ., family = Gamma("log"))

# Use "resp = <response_variable>" to specify different priors for different response variables
all_priors <- c(set_prior("normal(0,5)", class = "b"), # priors for y response variables
                set_prior("normal(0,2.5)", class = "Intercept"), 
                set_prior("exponential(1)", class = "shape"))

# Model converges after increasing the adapt_delta and iterations from default values
# Rule of thumb: bulk and tail effective sample sizes should be 100 x number of chains (i.e., at least 400)
# increasing max_treedepth is more about efficiency (instead of validity)
# See: https://mc-stan.org/misc/warnings.html
fit_diesel_no_taxa <- brm(diesel_brms_no_taxa, data = diesel_brms_data_no_taxa,
                        prior = all_priors, cores = 4, seed = "11729", iter = 20000, control = list(adapt_delta = 0.99))

# Get stan code
#stancode(fit_diesel_no_taxa)

######################################################################################################
# Use intensity + system model to predict diesel for taxa that were not part of the previous model
# Intensity must be non-NA
diesel_complete_predictors <- ghg_model_dat_categories %>%
  filter(is.na(intensity)==FALSE & is.na(system)==FALSE) %>%
  filter(is.na(diesel)) # Now, filter to just the NAs

taxa_not_modeled <- setdiff(unique(diesel_complete_predictors$taxa), unique(diesel_no_na$taxa)) # these taxa were never modeled so they can't be predicted below

# Only keep taxa that were not modeled previously
diesel_complete_predictors_no_taxa <- diesel_complete_predictors %>%
  filter(taxa %in% taxa_not_modeled) 

# Format intensity and system as ordinal variable, then center and scale
X_ordinal_no_taxa <- diesel_complete_predictors_no_taxa %>%
  mutate(intensity = factor(intensity, levels = c("extensive", "semi", "intensive"))) %>% # set order of factors (low = extensive, high = intensive)
  mutate(intensity = as.numeric(intensity)) %>%
  mutate(system = factor(system, levels = c("open", "semi", "closed"))) %>% # set order of factors (low = open, high = closed)
  mutate(system = as.numeric(system)) %>%
  select(intensity, system) %>%
  as.matrix()

# Center and Scale: BUT now center by the mean of the original modeled dataset above AND scale by the same 2*SD calculated from the original, modeled dataset above
X_ordinal_scaled_no_taxa <- scale(X_ordinal_no_taxa, center=apply(X_ordinal, MARGIN = 2, FUN = mean), scale=2*ordinal_sd)

# Create dataframe for brms
brms_diesel_dat_no_taxa <- data.frame(X_ordinal_scaled_no_taxa)

# Make predictions
#predicted_diesel_dat <- predict(fit_no_na, newdata = brms_new_diesel_data)
# Use tidybayes instead:
predicted_diesel_dat_no_taxa <- add_predicted_draws(newdata = brms_diesel_dat_no_taxa, model = fit_diesel_no_taxa)

# Get point and interval estimates from predicted data
# Select just the prediction columns
# Join these with the modeled data (diesel_complete_predictors_no_taxa) to get metadata on taxa/intensity/syste,
diesel_dat_intervals_no_taxa <- predicted_diesel_dat_no_taxa %>%
  median_qi(.value = .prediction) %>% # Rename prediction to value
  ungroup() %>%
  select(contains("."))

# .row is equivalent to the row number in the modeled dataset (diesel_complete_predictors_no_taxa) - create a join column for this
diesel_metadat_no_taxa <- diesel_complete_predictors_no_taxa %>%
  select(study_id, clean_sci_name, taxa, intensity, system) %>%
  mutate(.row = row_number())

diesel_predictions_no_taxa <- diesel_dat_intervals_no_taxa %>%
  left_join(diesel_metadat_no_taxa, by = ".row") %>%
  rename(diesel = .value)

# Bind diesel_no_na (data), diesel_predictions, and diesel_predictions_no_taxa
full_diesel_dat <- diesel_predictions %>%
  bind_rows(diesel_no_na) %>%
  bind_rows(diesel_predictions_no_taxa) %>%
  mutate(data_type = if_else(is.na(.point), true = "data", false = "prediction")) %>%
  arrange(taxa, intensity, system, clean_sci_name) %>%
  rownames_to_column() # Arrange by taxa first, then create dummy column for plotting 

# NEXT: use plot_brms_gamma_regression to produce figures for SI

######################################################################################################
# Step 3: Model petrol

petrol_no_na <- ghg_model_dat_categories %>%
  filter(is.na(petrol)==FALSE)  %>% # Drop NAs (Keep zeroes)
  mutate(petrol = if_else(petrol == 0, true = min(petrol[petrol!=0]), false = petrol)) %>% # Not modeling the zeroes, option 1: adjust these to the minimum value
  #filter(petrol != 0) %>% # Not modeling the zeroes, option 2: drop zeroes
  filter(is.na(intensity)==FALSE & is.na(system)==FALSE) %>% # complete predictors - i.e., both intensity AND system are non-NA
  select(study_id, petrol, clean_sci_name, taxa, intensity, system)

# Create model matrix for taxa info, then center and scale
X_taxa <- model.matrix(object = ~ 1 + taxa, 
                       data = petrol_no_na %>% select(taxa)) 

taxa_sd <- apply(X_taxa[,-1], MARGIN=2, FUN=sd, na.rm=TRUE) # Center all non-intercept variables and scale by 2 standard deviations (ignoring NAs)
X_taxa_scaled <- scale(X_taxa[,-1], center=TRUE, scale=2*taxa_sd)

# FIX IT - if only one column (e.g., after dropping zeroes, only one taxa column remaining in model matrix) no need to use apply function accross columns:
# Also, need to restore column name to X_taxa_scaled since vector doesn't retain this
#taxa_sd <- sd(X_taxa[,-1], na.rm = TRUE)
#X_taxa_scaled <- scale(X_taxa[,-1], center=TRUE, scale=2*taxa_sd)
#colnames(X_taxa_scaled) <- colnames(X_taxa)[2]

# Format intensity and system as ordinal variables, then center and scale
X_ordinal <- petrol_no_na %>%
  mutate(intensity = factor(intensity, levels = c("extensive", "semi", "intensive"))) %>% # set order of factors (low = extensive, high = intensive)
  mutate(system = factor(system, levels = c("open", "semi", "closed"))) %>% # set order of factors (low = open, high = closed)
  mutate(intensity = as.numeric(intensity)) %>%
  mutate(system = as.numeric(system)) %>%
  select(intensity, system) %>%
  as.matrix()
ordinal_sd<-apply(X_ordinal, MARGIN=2, FUN=sd, na.rm=TRUE) # Center all non-intercept variables and scale by 2 standard deviations (ignoring NAs)
X_ordinal_scaled <- scale(X_ordinal, center=TRUE, scale=2*ordinal_sd)

# Create dataframe for brms and rename feed variables
petrol_brms_data <- data.frame(y = petrol_no_na$petrol, X_taxa_scaled, X_ordinal_scaled)

names(petrol_brms_data)

# Set model formula
petrol_brms <- brmsformula(y ~ 1 + ., family = Gamma("log"))
# Equivalent to:
# petrol_brms <- brmsformula(y ~ 1 + taxatilapia + taxatrout +
#                            intensity + system, family = Gamma("log"))


# Use "resp = <response_variable>" to specify different priors for different response variables
all_priors <- c(set_prior("normal(0,5)", class = "b"), # priors for y response variables
                set_prior("normal(0,2.5)", class = "Intercept"), 
                set_prior("exponential(1)", class = "shape"))

# Model converges after increasing the adapt_delta and iterations from default values
# Rule of thumb: bulk and tail effective sample sizes should be 100 x number of chains (i.e., at least 400)
# increasing max_treedepth is more about efficiency (instead of validity)
# See: https://mc-stan.org/misc/warnings.html
fit_petrol_no_na <- brm(petrol_brms, data = petrol_brms_data,
                        prior = all_priors, cores = 4, seed = "11729", iter = 20000, control = list(adapt_delta = 0.99))

# Get stan code
#stancode(fit_petrol_no_na)

######################################################################################################
# Use model to predict NAs for studies with complete set of predictors
# Both intensity AND system are non-NA
petrol_complete_predictors <- ghg_model_dat_categories %>%
  filter(is.na(petrol)) %>% 
  filter(is.na(intensity)==FALSE & is.na(system)== FALSE) 

# PROBLEM: lca_complete predictors has more taxa than originally model:
taxa_not_modeled <- setdiff(unique(petrol_complete_predictors$taxa), unique(petrol_no_na$taxa)) # these taxa were never modeled so they can't be predicted below

# DROP THESE FOR NOW:
petrol_complete_predictors <- petrol_complete_predictors %>%
  filter(taxa %in% taxa_not_modeled == FALSE)

# Now check the other way, which taxa were in the original model but not a part of the data that needs to be predicted:
setdiff(unique(petrol_no_na$taxa), unique(petrol_complete_predictors$taxa))

# If original model has taxa that are not part of petrol_complete_predictors, 
# Use list of unique taxa in original model and use this to expand/assign levels manually - having trouble automating this
# Include all levels here, but remember the first level won't show up in design matrix - instead, it's part of the "contrasts"
# sort(unique(petrol_no_na$taxa))
# petrol_complete_predictors <- petrol_complete_predictors %>%
#   mutate(taxa = as.factor(taxa))
# levels(petrol_complete_predictors$taxa) <- list(fresh_crust = "fresh_crust", misc_fresh = "misc_fresh", misc_marine = "misc_marine", oth_carp = "oth_carp", tilapia = "tilapia")

# Create NEW taxa model matrix for the studies to be predicted
# Taxa categories:
X_taxa_new <- model.matrix(object = ~ 1 + taxa, 
                           data = petrol_complete_predictors %>% select(taxa)) 

# Center and Scale: BUT now center by the mean of the original modeled dataset above AND scale by the same 2*SD calculated from the original, modeled dataset above
X_taxa_new_scaled <- scale(X_taxa_new[,-1], center=apply(X_taxa[,-1], MARGIN = 2, FUN = mean), scale=2*taxa_sd)

# FIX IT - if X_taxa is only one column (e.g., if dropping zeroes, only one taxa column in model matrix) no need to use apply to get mean across columns
# Also need to restore column name
#X_taxa_new_scaled <- scale(X_taxa_new[,-1], center=mean(X_taxa[,-1]), scale=2*taxa_sd)
#colnames(X_taxa_new_scaled) <- colnames(X_taxa_new)[2]

# System and Intensity variables:
# Format intensity and system as ordinal variables, then center and scale
X_ordinal_new <- petrol_complete_predictors %>%
  mutate(intensity = factor(intensity, levels = c("extensive", "semi", "intensive"))) %>% # set order of factors (low = extensive, high = intensive)
  mutate(system = factor(system, levels = c("open", "semi", "closed"))) %>% # set order of factors (low = open, high = closed)
  mutate(intensity = as.numeric(intensity)) %>%
  mutate(system = as.numeric(system)) %>%
  select(intensity, system) %>%
  as.matrix()

# Center and Scale: BUT now center by the mean of the original modeled dataset above AND scale by the same 2*SD calculated from the original, modeled dataset above
X_ordinal_new_scaled <- scale(X_ordinal_new, center=apply(X_ordinal, MARGIN = 2, FUN = mean), scale=2*ordinal_sd)

# Create dataframe for brms and rename feed variables
brms_new_petrol_dat <- data.frame(cbind(X_taxa_new_scaled, X_ordinal_new_scaled)) 

# Make predictions
#predicted_petrol_dat <- predict(fit_no_na, newdata = brms_new_petrol_data)
# Use tidybayes instead:
predicted_petrol_dat <- add_predicted_draws(newdata = brms_new_petrol_dat, model = fit_petrol_no_na)

# Get point and interval estimates from predicted data
# Select just the prediction columns
# Join these with the modeled data (petrol_complete_predictors) to get metadata on taxa/intensity/syste,
petrol_dat_intervals <- predicted_petrol_dat %>%
  median_qi(.value = .prediction) %>% # Rename prediction to value
  ungroup() %>%
  select(contains("."))

# .row is equivalent to the row number in the modeled dataset (petrol_complete_predictors) - create a join column for this
petrol_metadat<- petrol_complete_predictors %>%
  select(study_id, clean_sci_name, taxa, intensity, system) %>%
  mutate(.row = row_number())

petrol_predictions <- petrol_dat_intervals %>%
  left_join(petrol_metadat, by = ".row") %>%
  rename(petrol = .value)

######################################################################################################
# Model petrol with intensity + system (no taxa)

# Create dataframe for brms and rename feed variables
petrol_brms_data_no_taxa <- data.frame(y = petrol_no_na$petrol, X_ordinal_scaled)

names(petrol_brms_data_no_taxa)

# Set model formula
petrol_brms_no_taxa <- brmsformula(y ~ 1 + ., family = Gamma("log"))

# Use "resp = <response_variable>" to specify different priors for different response variables
all_priors <- c(set_prior("normal(0,5)", class = "b"), # priors for y response variables
                set_prior("normal(0,2.5)", class = "Intercept"), 
                set_prior("exponential(1)", class = "shape"))

# Model converges after increasing the adapt_delta and iterations from default values
# Rule of thumb: bulk and tail effective sample sizes should be 100 x number of chains (i.e., at least 400)
# increasing max_treedepth is more about efficiency (instead of validity)
# See: https://mc-stan.org/misc/warnings.html
fit_petrol_no_taxa <- brm(petrol_brms_no_taxa, data = petrol_brms_data_no_taxa,
                          prior = all_priors, cores = 4, seed = "11729", iter = 20000, control = list(adapt_delta = 0.99))

# Get stan code
#stancode(fit_petrol_no_taxa)

######################################################################################################
# Use intensity + system model to predict petrol for taxa that were not part of the previous model
# Intensity must be non-NA
petrol_complete_predictors <- ghg_model_dat_categories %>%
  filter(is.na(intensity)==FALSE & is.na(system)==FALSE) %>%
  filter(is.na(petrol)) # Now, filter to just the NAs

taxa_not_modeled <- setdiff(unique(petrol_complete_predictors$taxa), unique(petrol_no_na$taxa)) # these taxa were never modeled so they can't be predicted below

# Only keep taxa that were not modeled previously
petrol_complete_predictors_no_taxa <- petrol_complete_predictors %>%
  filter(taxa %in% taxa_not_modeled) 

# Format intensity and system as ordinal variable, then center and scale
X_ordinal_no_taxa <- petrol_complete_predictors_no_taxa %>%
  mutate(intensity = factor(intensity, levels = c("extensive", "semi", "intensive"))) %>% # set order of factors (low = extensive, high = intensive)
  mutate(intensity = as.numeric(intensity)) %>%
  mutate(system = factor(system, levels = c("open", "semi", "closed"))) %>% # set order of factors (low = open, high = closed)
  mutate(system = as.numeric(system)) %>%
  select(intensity, system) %>%
  as.matrix()

# Center and Scale: BUT now center by the mean of the original modeled dataset above AND scale by the same 2*SD calculated from the original, modeled dataset above
X_ordinal_scaled_no_taxa <- scale(X_ordinal_no_taxa, center=apply(X_ordinal, MARGIN = 2, FUN = mean), scale=2*ordinal_sd)

# Create dataframe for brms
brms_petrol_dat_no_taxa <- data.frame(X_ordinal_scaled_no_taxa)

# Make predictions
#predicted_petrol_dat <- predict(fit_no_na, newdata = brms_new_petrol_data)
# Use tidybayes instead:
predicted_petrol_dat_no_taxa <- add_predicted_draws(newdata = brms_petrol_dat_no_taxa, model = fit_petrol_no_taxa)

# Get point and interval estimates from predicted data
# Select just the prediction columns
# Join these with the modeled data (petrol_complete_predictors_no_taxa) to get metadata on taxa/intensity/syste,
petrol_dat_intervals_no_taxa <- predicted_petrol_dat_no_taxa %>%
  median_qi(.value = .prediction) %>% # Rename prediction to value
  ungroup() %>%
  select(contains("."))

# .row is equivalent to the row number in the modeled dataset (petrol_complete_predictors_no_taxa) - create a join column for this
petrol_metadat_no_taxa <- petrol_complete_predictors_no_taxa %>%
  select(study_id, clean_sci_name, taxa, intensity, system) %>%
  mutate(.row = row_number())

petrol_predictions_no_taxa <- petrol_dat_intervals_no_taxa %>%
  left_join(petrol_metadat_no_taxa, by = ".row") %>%
  rename(petrol = .value)

# Bind petrol_no_na (data), petrol_predictions, and petrol_predictions_no_taxa
full_petrol_dat <- petrol_predictions %>%
  bind_rows(petrol_no_na) %>%
  bind_rows(petrol_predictions_no_taxa) %>%
  mutate(data_type = if_else(is.na(.point), true = "data", false = "prediction")) %>%
  arrange(taxa, intensity, system, clean_sci_name) %>%
  rownames_to_column() # Arrange by taxa first, then create dummy column for plotting 

# NEXT: use plot_brms_gamma_regression to produce figures for SI


######################################################################################################
# Step 4: Model natural gas

natgas_no_na <- ghg_model_dat_categories %>%
  filter(is.na(natgas)==FALSE)  %>% # Drop NAs (Keep zeroes)
  mutate(natgas = if_else(natgas == 0, true = min(natgas[natgas!=0]), false = natgas)) %>% # Not modeling the zeroes, option 1: adjust these to the minimum value
  #filter(natgas != 0) %>% # Not modeling the zeroes option 2: drop zeroes
  filter(is.na(intensity)==FALSE & is.na(system)==FALSE) %>% # complete predictors - i.e., both intensity AND system are non-NA
  select(study_id, natgas, clean_sci_name, taxa, intensity, system)

# FIX IT - option 2 (dropping zeroes) is not a viable option, after dropping zeroes, there's no varaition in taxa (all salmon), intensity (all intensive) or system (all open)

# Create model matrix for taxa info, then center and scale
X_taxa <- model.matrix(object = ~ 1 + taxa, 
                       data = natgas_no_na %>% select(taxa)) 

taxa_sd <- apply(X_taxa[,-1], MARGIN=2, FUN=sd, na.rm=TRUE) # Center all non-intercept variables and scale by 2 standard deviations (ignoring NAs)
X_taxa_scaled <- scale(X_taxa[,-1], center=TRUE, scale=2*taxa_sd)

# Format intensity and system as ordinal variables, then center and scale
X_ordinal <- natgas_no_na %>%
  mutate(intensity = factor(intensity, levels = c("extensive", "semi", "intensive"))) %>% # set order of factors (low = extensive, high = intensive)
  mutate(system = factor(system, levels = c("open", "semi", "closed"))) %>% # set order of factors (low = open, high = closed)
  mutate(intensity = as.numeric(intensity)) %>%
  mutate(system = as.numeric(system)) %>%
  select(intensity, system) %>%
  as.matrix()
ordinal_sd<-apply(X_ordinal, MARGIN=2, FUN=sd, na.rm=TRUE) # Center all non-intercept variables and scale by 2 standard deviations (ignoring NAs)
X_ordinal_scaled <- scale(X_ordinal, center=TRUE, scale=2*ordinal_sd)

# Create dataframe for brms and rename feed variables
natgas_brms_data <- data.frame(y = natgas_no_na$natgas, X_taxa_scaled, X_ordinal_scaled)

names(natgas_brms_data)

# Set model formula
natgas_brms <- brmsformula(y ~ 1 + ., family = Gamma("log"))
# Equivalent to:
# natgas_brms <- brmsformula(y ~ 1 + taxatilapia + taxatrout +
#                            intensity + system, family = Gamma("log"))


# Use "resp = <response_variable>" to specify different priors for different response variables
all_priors <- c(set_prior("normal(0,5)", class = "b"), # priors for y response variables
                set_prior("normal(0,2.5)", class = "Intercept"), 
                set_prior("exponential(1)", class = "shape"))

# Model converges after increasing the adapt_delta and iterations from default values
# Rule of thumb: bulk and tail effective sample sizes should be 100 x number of chains (i.e., at least 400)
# increasing max_treedepth is more about efficiency (instead of validity)
# See: https://mc-stan.org/misc/warnings.html
fit_natgas_no_na <- brm(natgas_brms, data = natgas_brms_data,
                        prior = all_priors, cores = 4, seed = "11729", iter = 20000, control = list(adapt_delta = 0.99))

# Get stan code
#stancode(fit_natgas_no_na)

######################################################################################################
# Use model to predict NAs for studies with complete set of predictors
# Both intensity AND system are non-NA
natgas_complete_predictors <- ghg_model_dat_categories %>%
  filter(is.na(natgas)) %>% 
  filter(is.na(intensity)==FALSE & is.na(system)== FALSE) 

# PROBLEM: lca_complete predictors has more taxa than originally model:
taxa_not_modeled <- setdiff(unique(natgas_complete_predictors$taxa), unique(natgas_no_na$taxa)) # these taxa were never modeled so they can't be predicted below

# DROP THESE FOR NOW:
natgas_complete_predictors <- natgas_complete_predictors %>%
  filter(taxa %in% taxa_not_modeled == FALSE)

# Now check the other way, which taxa were in the original model but not a part of the data that needs to be predicted:
setdiff(unique(natgas_no_na$taxa), unique(natgas_complete_predictors$taxa))

# If original model has taxa that are not part of natgas_complete_predictors, 
# Use list of unique taxa in original model and use this to expand/assign levels manually - having trouble automating this
# Include all levels here, but remember the first level won't show up in design matrix - instead, it's part of the "contrasts"
# sort(unique(natgas_no_na$taxa))
# natgas_complete_predictors <- natgas_complete_predictors %>%
#   mutate(taxa = as.factor(taxa))
# levels(natgas_complete_predictors$taxa) <- list(fresh_crust = "fresh_crust", misc_fresh = "misc_fresh", misc_marine = "misc_marine", oth_carp = "oth_carp", tilapia = "tilapia")

# Create NEW taxa model matrix for the studies to be predicted
# Taxa categories:
X_taxa_new <- model.matrix(object = ~ 1 + taxa, 
                           data = natgas_complete_predictors %>% select(taxa)) 

# Center and Scale: BUT now center by the mean of the original modeled dataset above AND scale by the same 2*SD calculated from the original, modeled dataset above
X_taxa_new_scaled <- scale(X_taxa_new[,-1], center=apply(X_taxa[,-1], MARGIN = 2, FUN = mean), scale=2*taxa_sd)

# System and Intensity variables:
# Format intensity and system as ordinal variables, then center and scale
X_ordinal_new <- natgas_complete_predictors %>%
  mutate(intensity = factor(intensity, levels = c("extensive", "semi", "intensive"))) %>% # set order of factors (low = extensive, high = intensive)
  mutate(system = factor(system, levels = c("open", "semi", "closed"))) %>% # set order of factors (low = open, high = closed)
  mutate(intensity = as.numeric(intensity)) %>%
  mutate(system = as.numeric(system)) %>%
  select(intensity, system) %>%
  as.matrix()

# Center and Scale: BUT now center by the mean of the original modeled dataset above AND scale by the same 2*SD calculated from the original, modeled dataset above
X_ordinal_new_scaled <- scale(X_ordinal_new, center=apply(X_ordinal, MARGIN = 2, FUN = mean), scale=2*ordinal_sd)

# Create dataframe for brms and rename feed variables
brms_new_natgas_dat <- data.frame(cbind(X_taxa_new_scaled, X_ordinal_new_scaled)) 

# Make predictions
#predicted_natgas_dat <- predict(fit_no_na, newdata = brms_new_natgas_data)
# Use tidybayes instead:
predicted_natgas_dat <- add_predicted_draws(newdata = brms_new_natgas_dat, model = fit_natgas_no_na)

# Get point and interval estimates from predicted data
# Select just the prediction columns
# Join these with the modeled data (natgas_complete_predictors) to get metadata on taxa/intensity/syste,
natgas_dat_intervals <- predicted_natgas_dat %>%
  median_qi(.value = .prediction) %>% # Rename prediction to value
  ungroup() %>%
  select(contains("."))

# .row is equivalent to the row number in the modeled dataset (natgas_complete_predictors) - create a join column for this
natgas_metadat<- natgas_complete_predictors %>%
  select(study_id, clean_sci_name, taxa, intensity, system) %>%
  mutate(.row = row_number())

natgas_predictions <- natgas_dat_intervals %>%
  left_join(natgas_metadat, by = ".row") %>%
  rename(natgas = .value)

######################################################################################################
# Model natgas with intensity + system (no taxa)

# Create dataframe for brms and rename feed variables
natgas_brms_data_no_taxa <- data.frame(y = natgas_no_na$natgas, X_ordinal_scaled)

names(natgas_brms_data_no_taxa)

# Set model formula
natgas_brms_no_taxa <- brmsformula(y ~ 1 + ., family = Gamma("log"))

# Use "resp = <response_variable>" to specify different priors for different response variables
all_priors <- c(set_prior("normal(0,5)", class = "b"), # priors for y response variables
                set_prior("normal(0,2.5)", class = "Intercept"), 
                set_prior("exponential(1)", class = "shape"))

# Model converges after increasing the adapt_delta and iterations from default values
# Rule of thumb: bulk and tail effective sample sizes should be 100 x number of chains (i.e., at least 400)
# increasing max_treedepth is more about efficiency (instead of validity)
# See: https://mc-stan.org/misc/warnings.html
fit_natgas_no_taxa <- brm(natgas_brms_no_taxa, data = natgas_brms_data_no_taxa,
                          prior = all_priors, cores = 4, seed = "11729", iter = 20000, control = list(adapt_delta = 0.99))

# Get stan code
#stancode(fit_natgas_no_taxa)

######################################################################################################
# Use intensity + system model to predict natgas for taxa that were not part of the previous model
# Intensity must be non-NA
natgas_complete_predictors <- ghg_model_dat_categories %>%
  filter(is.na(intensity)==FALSE & is.na(system)==FALSE) %>%
  filter(is.na(natgas)) # Now, filter to just the NAs

taxa_not_modeled <- setdiff(unique(natgas_complete_predictors$taxa), unique(natgas_no_na$taxa)) # these taxa were never modeled so they can't be predicted below

# Only keep taxa that were not modeled previously
natgas_complete_predictors_no_taxa <- natgas_complete_predictors %>%
  filter(taxa %in% taxa_not_modeled) 

# Format intensity and system as ordinal variable, then center and scale
X_ordinal_no_taxa <- natgas_complete_predictors_no_taxa %>%
  mutate(intensity = factor(intensity, levels = c("extensive", "semi", "intensive"))) %>% # set order of factors (low = extensive, high = intensive)
  mutate(intensity = as.numeric(intensity)) %>%
  mutate(system = factor(system, levels = c("open", "semi", "closed"))) %>% # set order of factors (low = open, high = closed)
  mutate(system = as.numeric(system)) %>%
  select(intensity, system) %>%
  as.matrix()

# Center and Scale: BUT now center by the mean of the original modeled dataset above AND scale by the same 2*SD calculated from the original, modeled dataset above
X_ordinal_scaled_no_taxa <- scale(X_ordinal_no_taxa, center=apply(X_ordinal, MARGIN = 2, FUN = mean), scale=2*ordinal_sd)

# Create dataframe for brms
brms_natgas_dat_no_taxa <- data.frame(X_ordinal_scaled_no_taxa)

# Make predictions
#predicted_natgas_dat <- predict(fit_no_na, newdata = brms_new_natgas_data)
# Use tidybayes instead:
predicted_natgas_dat_no_taxa <- add_predicted_draws(newdata = brms_natgas_dat_no_taxa, model = fit_natgas_no_taxa)

# Get point and interval estimates from predicted data
# Select just the prediction columns
# Join these with the modeled data (natgas_complete_predictors_no_taxa) to get metadata on taxa/intensity/syste,
natgas_dat_intervals_no_taxa <- predicted_natgas_dat_no_taxa %>%
  median_qi(.value = .prediction) %>% # Rename prediction to value
  ungroup() %>%
  select(contains("."))

# .row is equivalent to the row number in the modeled dataset (natgas_complete_predictors_no_taxa) - create a join column for this
natgas_metadat_no_taxa <- natgas_complete_predictors_no_taxa %>%
  select(study_id, clean_sci_name, taxa, intensity, system) %>%
  mutate(.row = row_number())

natgas_predictions_no_taxa <- natgas_dat_intervals_no_taxa %>%
  left_join(natgas_metadat_no_taxa, by = ".row") %>%
  rename(natgas = .value)

# Bind natgas_no_na (data), natgas_predictions, and natgas_predictions_no_taxa
full_natgas_dat <- natgas_predictions %>%
  bind_rows(natgas_no_na) %>%
  bind_rows(natgas_predictions_no_taxa) %>%
  mutate(data_type = if_else(is.na(.point), true = "data", false = "prediction")) %>%
  arrange(taxa, intensity, system, clean_sci_name) %>%
  rownames_to_column() # Arrange by taxa first, then create dummy column for plotting 

# NEXT: use plot_brms_gamma_regression to produce figures for SI
