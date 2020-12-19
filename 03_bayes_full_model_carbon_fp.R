# Impute on-farm CARBON impacts then combine with off-farm (feed) CARBON impacts 

# Calculate on-farm (non-feed associated) carbon footprint as:
# (Electricity use * country-specific GHG of electricity) + (Diesel * GHG of diesel) + (Petrol * GHG of petrol) + (Natural gas * GHG of natural gas)

# DECISION: adjusting zeroes to be very small amount for all energy datasets
# Helps with modeling the NAs
# Have to do this anyway for the final aggregation gamma model (must be positive)

# Get model-specific data:
# SELECT STUDY ID COLUMN - use this for rejoining outputs from multiple regression models back together
# Select relevant data columns and arrange by categorical info
# Select iso3c so this can be joined with country-specific GHG emissions for electricity
ghg_model_dat_categories <- lca_dat_clean_groups %>%
  select(study_id, Country, iso3c, electric = Electricity_kwh, diesel = Diesel_L, petrol = Petrol_L, natgas = NaturalGas_L, clean_sci_name, taxa, intensity = Intensity, system = Production_system_group) %>%
  arrange(clean_sci_name, taxa, intensity, system) %>%
  mutate(electric = if_else(electric == 0, true = 0.01, false = electric),
         diesel = if_else(diesel == 0, true = 0.01, false = diesel),
         petrol = if_else(petrol == 0, true = 0.01, false = petrol),
         natgas = if_else(natgas == 0, true = 0.01, false = natgas))

# Clear all memory except for final stan model:
rm(list=ls()[!(ls() %in% c("datadir", "outdir", 
                           "lca_dat_clean_groups","ghg_model_dat_categories", 
                           "electric_fp_dat", "other_energy_fp_dat"))])

######################################################################################################
# Step 1: Model electricity as taxa + intensity + system

electric_no_na <- ghg_model_dat_categories %>%
  filter(is.na(electric)==FALSE)  %>% # Drop NAs (Keep zeroes)
  #mutate(electric = if_else(electric == 0, true = min(electric[electric!=0]), false = electric)) %>% # Not modeling the zeroes, option 1: adjust these to the minimum value
  filter(electric != 0) %>% # Not modeling the zeroes, option 2: drop all zeroes
  filter(is.na(intensity)==FALSE & is.na(system)==FALSE) %>% # complete predictors - i.e., both intensity AND system are non-NA
  select(study_id, Country, iso3c, electric, clean_sci_name, taxa, intensity, system)

# Create model matrix for taxa info, then center and scale
X_taxa <- model.matrix(object = ~ 1 + taxa, 
                       data = electric_no_na %>% select(taxa)) 

taxa_sd <- apply(X_taxa[,-1], MARGIN=2, FUN=sd, na.rm=TRUE) # Center all non-intercept variables and scale by 2 standard deviations (ignoring NAs)
#taxa_sd <- sd(X_taxa[,-1], na.rm = TRUE)
X_taxa_scaled <- scale(X_taxa[,-1], center=TRUE, scale=2*taxa_sd)

# Format intensity and system as ordinal variables, then center and scale
X_ordinal <- electric_no_na %>%
  mutate(intensity = factor(intensity, levels = c("Extensive", "Imp. extensive", "Semi-intensive", "Intensive"))) %>% # set order of factors (low = extensive, high = intensive)
  mutate(system = factor(system, levels = c("On- and off-bottom", "Cages & pens", "Ponds", "Recirculating and tanks"))) %>% # set order of factors (low = open, high = closed) 
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
                          prior = all_priors, cores = 4, seed = "11729", iter = 5000, control = list(adapt_delta = 0.99))

# Get stan code
#stancode(fit_electric_no_na)

######################################################################################################
# Use model to predict NAs for studies with complete set of predictors
# Both intensity AND system are non-NA
electric_complete_predictors <- ghg_model_dat_categories %>%
  filter(is.na(electric)) %>% 
  filter(is.na(intensity)==FALSE & is.na(system)==FALSE) 

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
sort(unique(electric_no_na$taxa))
electric_complete_predictors <- electric_complete_predictors %>%
  mutate(taxa = as.factor(taxa))
levels(electric_complete_predictors$taxa) <- list(bivalves = "bivalves",
                                                  catfish = "catfish",
                                                  hypoph_carp = "hypoph_carp",
                                                  milkfish = "milkfish",
                                                  misc_diad = "misc_diad",
                                                  misc_marine = "misc_marine",
                                                  oth_carp = "oth_carp",
                                                  plants = "plants",
                                                  salmon = "salmon",
                                                  shrimp = "shrimp",
                                                  tilapia = "tilapia",
                                                  trout = "trout")

# Create NEW taxa model matrix for the studies to be predicted
# Taxa categories:
X_taxa_new <- model.matrix(object = ~ 1 + taxa, 
                           data = electric_complete_predictors %>% select(taxa)) 

# Center and Scale: BUT now center by the mean of the original modeled dataset above AND scale by the same 2*SD calculated from the original, modeled dataset above
X_taxa_new_scaled <- scale(X_taxa_new[,-1], center=apply(X_taxa[,-1], MARGIN = 2, FUN = mean), scale=2*taxa_sd)

# System and Intensity variables:
# Format intensity and system as ordinal variables, then center and scale
X_ordinal_new <- electric_complete_predictors %>%
  mutate(intensity = factor(intensity, levels = c("Extensive", "Imp. extensive", "Semi-intensive", "Intensive"))) %>% # set order of factors (low = extensive, high = intensive)
  mutate(system = factor(system, levels = c("On- and off-bottom", "Cages & pens", "Ponds", "Recirculating and tanks"))) %>% # set order of factors (low = open, high = closed)
  mutate(intensity = as.numeric(intensity)) %>%
  mutate(system = as.numeric(system)) %>%
  select(intensity, system) %>%
  #select(system) %>%
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
  select(study_id, Country, iso3c, clean_sci_name, taxa, intensity, system) %>%
  mutate(.row = row_number())

electric_predictions <- electric_dat_intervals %>%
  left_join(electric_metadat, by = ".row") %>%
  rename(electric = .value)

######################################################################################################
# The following section does not apply to electricity (all taxa predicted):
# Model electricity with intensity + system (no taxa)

# # Create dataframe for brms and rename feed variables
# electric_brms_data_no_taxa <- data.frame(y = electric_no_na$electric, X_ordinal_scaled)
# 
# names(electric_brms_data_no_taxa)
# 
# # Set model formula
# electric_brms_no_taxa <- brmsformula(y ~ 1 + ., family = Gamma("log"))
# 
# # Use "resp = <response_variable>" to specify different priors for different response variables
# all_priors <- c(set_prior("normal(0,5)", class = "b"), # priors for y response variables
#                 set_prior("normal(0,2.5)", class = "Intercept"), 
#                 set_prior("exponential(1)", class = "shape"))
# 
# # Model converges after increasing the adapt_delta and iterations from default values
# # Rule of thumb: bulk and tail effective sample sizes should be 100 x number of chains (i.e., at least 400)
# # increasing max_treedepth is more about efficiency (instead of validity)
# # See: https://mc-stan.org/misc/warnings.html
# fit_electric_no_taxa <- brm(electric_brms_no_taxa, data = electric_brms_data_no_taxa,
#                          prior = all_priors, cores = 4, seed = "11729", iter = 5000, control = list(adapt_delta = 0.99))

# Get stan code
#stancode(fit_electric_no_taxa)

######################################################################################################

# Use intensity + system model to predict electricity for taxa that were not part of the previous model
# Intensity must be non-NA
# electric_complete_predictors <- ghg_model_dat_categories %>%
#   filter(is.na(intensity)==FALSE & is.na(system)==FALSE) %>%
#   filter(is.na(electric)) # Now, filter to just the NAs
# 
# taxa_not_modeled <- setdiff(unique(electric_complete_predictors$taxa), unique(electric_no_na$taxa)) # these taxa were never modeled so they can't be predicted below
# 
# # Only keep taxa that were not modeled previously
# electric_complete_predictors_no_taxa <- electric_complete_predictors %>%
#   filter(taxa %in% taxa_not_modeled) 
# 
# # Format intensity and system as ordinal variable, then center and scale
# X_ordinal_no_taxa <- electric_complete_predictors_no_taxa %>%
#   mutate(intensity = factor(intensity, levels = c("Extensive", "Imp. extensive", "Semi-intensive", "Intensive"))) %>% # set order of factors (low = extensive, high = intensive)
#   mutate(intensity = as.numeric(intensity)) %>%
#   mutate(system = factor(system, levels = c("On- and off-bottom", "Cages & pens", "Ponds", "Recirculating and tanks"))) %>% # set order of factors (low = open, high = closed)
#   mutate(system = as.numeric(system)) %>%
#   select(intensity, system) %>%
#   #select(system) %>%
#   as.matrix()
# 
# # Center and Scale: BUT now center by the mean of the original modeled dataset above AND scale by the same 2*SD calculated from the original, modeled dataset above
# X_ordinal_scaled_no_taxa <- scale(X_ordinal_no_taxa, center=apply(X_ordinal, MARGIN = 2, FUN = mean), scale=2*ordinal_sd)
# 
# # Create dataframe for brms
# brms_electric_dat_no_taxa <- data.frame(X_ordinal_scaled_no_taxa)
# 
# # Make predictions
# #predicted_electric_dat <- predict(fit_no_na, newdata = brms_new_electric_data)
# # Use tidybayes instead:
# predicted_electric_dat_no_taxa <- add_predicted_draws(newdata = brms_electric_dat_no_taxa, model = fit_electric_no_taxa)
# 
# # Get point and interval estimates from predicted data
# # Select just the prediction columns
# # Join these with the modeled data (electric_complete_predictors_no_taxa) to get metadata on taxa/intensity/syste,
# electric_dat_intervals_no_taxa <- predicted_electric_dat_no_taxa %>%
#   median_qi(.value = .prediction) %>% # Rename prediction to value
#   ungroup() %>%
#   select(contains("."))
# 
# # .row is equivalent to the row number in the modeled dataset (electric_complete_predictors_no_taxa) - create a join column for this
# electric_metadat_no_taxa <- electric_complete_predictors_no_taxa %>%
#   select(study_id, Country, iso3c, clean_sci_name, taxa, intensity, system) %>%
#   mutate(.row = row_number())
# 
# electric_predictions_no_taxa <- electric_dat_intervals_no_taxa %>%
#   left_join(electric_metadat_no_taxa, by = ".row") %>%
#   rename(electric = .value)

######################################################################################################
# Bind electric_no_na (data), electric_predictions, and electric_predictions_no_taxa

full_electric_dat <- electric_predictions %>%
  bind_rows(electric_no_na) %>%
  #bind_rows(electric_predictions_no_taxa) %>%
  mutate(data_type = if_else(is.na(.point), true = "data", false = "prediction")) %>%
  arrange(taxa, intensity, system, clean_sci_name) %>%
  rownames_to_column() # Arrange by taxa first, then create dummy column for plotting 


# Quick check: PLOT DATA + PREDICTIONS
ggplot(full_electric_dat, aes(x = electric, y = taxa)) + geom_boxplot(outlier.shape = NA) + geom_jitter(aes(color = data_type))

## FIX IT - temporary fix - filter out predictions with high uncertainty based on highest and lowest data values
## shouldn't be predicting out of the range of observations
min_dat <- full_electric_dat %>%
  filter(data_type == "data") %>%
  pull(electric) %>%
  min()

max_dat <- full_electric_dat %>%
  filter(data_type == "data") %>%
  pull(electric) %>%
  max()

full_electric_dat <- full_electric_dat %>%
  filter(electric <= max_dat & electric >= min_dat)

# NEXT: use plot_brms_gamma_regression to produce figures for SI

# Clear all memory except for final stan model:
rm(list=ls()[!(ls() %in% c("datadir", "outdir", 
                           "lca_dat_clean_groups", "ghg_model_dat_categories", 
                           "electric_fp_dat", "other_energy_fp_dat", 
                           "full_electric_dat"))])

######################################################################################################
# Step 2: Model diesel

diesel_no_na <- ghg_model_dat_categories %>%
  filter(is.na(diesel)==FALSE)  %>% # Drop NAs (Keep zeroes)
  #mutate(diesel = if_else(diesel == 0, true = 0.01, false = diesel)) %>% # Not modeling the zeroes, option 1: adjust these to the minimum value
  filter(diesel != 0) %>% # Not modeling the zeroes, option 2: drop all zeroes
  filter(is.na(intensity)==FALSE & is.na(system)==FALSE) %>% # complete predictors - i.e., both intensity AND system are non-NA
  select(study_id, Country, iso3c, diesel, clean_sci_name, taxa, intensity, system)

# Create model matrix for taxa info, then center and scale
X_taxa <- model.matrix(object = ~ 1 + taxa, 
                       data = diesel_no_na %>% select(taxa)) 

taxa_sd <- apply(X_taxa[,-1], MARGIN=2, FUN=sd, na.rm=TRUE) # Center all non-intercept variables and scale by 2 standard deviations (ignoring NAs)
X_taxa_scaled <- scale(X_taxa[,-1], center=TRUE, scale=2*taxa_sd)

# Format intensity and system as ordinal variables, then center and scale
X_ordinal <- diesel_no_na %>%
  mutate(intensity = factor(intensity, levels = c("Extensive", "Imp. extensive", "Semi-intensive", "Intensive"))) %>% # set order of factors (low = extensive, high = intensive)
  mutate(system = factor(system, levels = c("On- and off-bottom", "Cages & pens", "Ponds", "Recirculating and tanks"))) %>% # set order of factors (low = open, high = closed)
  mutate(intensity = as.numeric(intensity)) %>%
  mutate(system = as.numeric(system)) %>%
  select(intensity, system) %>%
  #select(system) %>%
  as.matrix()
ordinal_sd<-apply(X_ordinal, MARGIN=2, FUN=sd, na.rm=TRUE) # Center all non-intercept variables and scale by 2 standard deviations (ignoring NAs)
X_ordinal_scaled <- scale(X_ordinal, center=TRUE, scale=2*ordinal_sd)

# Create dataframe for brms and rename feed variables
diesel_brms_data <- data.frame(y = diesel_no_na$diesel, X_taxa_scaled, X_ordinal_scaled)

names(diesel_brms_data)

# Set model formula
diesel_brms <- brmsformula(y ~ 1 + ., family = Gamma("log"))

# Use "resp = <response_variable>" to specify different priors for different response variables
all_priors <- c(set_prior("normal(0,5)", class = "b"), # priors for y response variables
                set_prior("normal(0,2.5)", class = "Intercept"), 
                set_prior("exponential(1)", class = "shape"))

# Model converges after increasing the adapt_delta and iterations from default values
# Rule of thumb: bulk and tail effective sample sizes should be 100 x number of chains (i.e., at least 400)
# increasing max_treedepth is more about efficiency (instead of validity)
# See: https://mc-stan.org/misc/warnings.html
fit_diesel_no_na <- brm(diesel_brms, data = diesel_brms_data,
                        prior = all_priors, cores = 4, seed = "11729", iter = 5000, control = list(adapt_delta = 0.99))

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
sort(unique(diesel_no_na$taxa))
diesel_complete_predictors <- diesel_complete_predictors %>%
  mutate(taxa = as.factor(taxa))

levels(diesel_complete_predictors$taxa) <- list(bivalves = "bivalves",
                                                catfish = "catfish",
                                                hypoph_carp = "hypoph_carp",
                                                milkfish = "milkfish",
                                                misc_diad = "misc_diad",
                                                misc_marine = "misc_marine",
                                                oth_carp = "oth_carp",
                                                plants = "plants",
                                                salmon = "salmon",
                                                shrimp = "shrimp",
                                                tilapia = "tilapia",
                                                trout = "trout")

# Create NEW taxa model matrix for the studies to be predicted
# Taxa categories:
X_taxa_new <- model.matrix(object = ~ 1 + taxa, 
                           data = diesel_complete_predictors %>% select(taxa)) 

# Center and Scale: BUT now center by the mean of the original modeled dataset above AND scale by the same 2*SD calculated from the original, modeled dataset above
X_taxa_new_scaled <- scale(X_taxa_new[,-1], center=apply(X_taxa[,-1], MARGIN = 2, FUN = mean), scale=2*taxa_sd)

# System and Intensity variables:
# Format intensity and system as ordinal variables, then center and scale
X_ordinal_new <- diesel_complete_predictors %>%
  mutate(intensity = factor(intensity, levels = c("Extensive", "Imp. extensive", "Semi-intensive", "Intensive"))) %>% # set order of factors (low = extensive, high = intensive)
  mutate(system = factor(system, levels = c("On- and off-bottom", "Cages & pens", "Ponds", "Recirculating and tanks"))) %>% # set order of factors (low = open, high = closed)
  mutate(intensity = as.numeric(intensity)) %>%
  mutate(system = as.numeric(system)) %>%
  select(intensity, system) %>%
  #select(system) %>%
  as.matrix()

# Center and Scale: BUT now center by the mean of the original modeled dataset above AND scale by the same 2*SD calculated from the original, modeled dataset above
X_ordinal_new_scaled <- scale(X_ordinal_new, center=apply(X_ordinal, MARGIN = 2, FUN = mean), scale=2*ordinal_sd)

# Create dataframe for brms and rename feed variables
brms_new_diesel_dat <- data.frame(cbind(X_taxa_new_scaled, X_ordinal_new_scaled)) 

# Make predictions
#predicted_diesel_dat <- predict(fit_no_na, newdata = brms_new_diesel_data)
# Use tidybayes instead:
predicted_diesel_dat <- add_predicted_draws(newdata = brms_new_diesel_dat, model = fit_diesel_no_na)
#predicted_diesel_dat <- add_fitted_draws(newdata = brms_new_diesel_dat, model = fit_diesel_no_na)

# Get point and interval estimates from predicted data
# Select just the prediction columns
# Join these with the modeled data (diesel_complete_predictors) to get metadata on taxa/intensity/syste,
diesel_dat_intervals <- predicted_diesel_dat %>%
  median_qi(.value = .prediction) %>% # Rename prediction to value
  ungroup() %>%
  select(contains("."))

# .row is equivalent to the row number in the modeled dataset (diesel_complete_predictors) - create a join column for this
diesel_metadat<- diesel_complete_predictors %>%
  select(study_id, Country, iso3c, clean_sci_name, taxa, intensity, system) %>%
  mutate(.row = row_number())

diesel_predictions <- diesel_dat_intervals %>%
  left_join(diesel_metadat, by = ".row") %>%
  rename(diesel = .value)

######################################################################################################
# Model diesel with intensity + system (no taxa) 
# This section doesn't apply to diesel (all taxa have data)

# Create dataframe for brms and rename feed variables
# diesel_brms_data_no_taxa <- data.frame(y = diesel_no_na$diesel, X_ordinal_scaled)
# 
# names(diesel_brms_data_no_taxa)
# 
# # Set model formula
# diesel_brms_no_taxa <- brmsformula(y ~ 1 + ., family = Gamma("log"))
# 
# # Use "resp = <response_variable>" to specify different priors for different response variables
# all_priors <- c(set_prior("normal(0,5)", class = "b"), # priors for y response variables
#                 set_prior("normal(0,2.5)", class = "Intercept"),
#                 set_prior("exponential(1)", class = "shape"))
# 
# # Model converges after increasing the adapt_delta and iterations from default values
# # Rule of thumb: bulk and tail effective sample sizes should be 100 x number of chains (i.e., at least 400)
# # increasing max_treedepth is more about efficiency (instead of validity)
# # See: https://mc-stan.org/misc/warnings.html
# fit_diesel_no_taxa <- brm(diesel_brms_no_taxa, data = diesel_brms_data_no_taxa,
#                         prior = all_priors, cores = 4, seed = "11729", iter = 5000, control = list(adapt_delta = 0.99))
# 
# # Get stan code
# #stancode(fit_diesel_no_taxa)
# 
# ######################################################################################################
# # Use intensity + system model to predict diesel for taxa that were not part of the previous model
# # Intensity must be non-NA
# diesel_complete_predictors <- ghg_model_dat_categories %>%
#   filter(is.na(intensity)==FALSE & is.na(system)==FALSE) %>%
#   filter(is.na(diesel)) # Now, filter to just the NAs
# 
# taxa_not_modeled <- setdiff(unique(diesel_complete_predictors$taxa), unique(diesel_no_na$taxa)) # these taxa were never modeled so they can't be predicted below
# 
# # Only keep taxa that were not modeled previously
# diesel_complete_predictors_no_taxa <- diesel_complete_predictors %>%
#   filter(taxa %in% taxa_not_modeled)
# 
# # Format intensity and system as ordinal variable, then center and scale
# X_ordinal_no_taxa <- diesel_complete_predictors_no_taxa %>%
#   mutate(intensity = factor(intensity, levels = c("Extensive", "Imp. extensive", "Semi-intensive", "Intensive"))) %>% # set order of factors (low = extensive, high = intensive)
#   mutate(intensity = as.numeric(intensity)) %>%
#   mutate(system = factor(system, levels = c("On- and off-bottom", "Cages & pens", "Ponds", "Recirculating and tanks"))) %>% # set order of factors (low = open, high = closed)
#   mutate(system = as.numeric(system)) %>%
#   select(intensity, system) %>%
#   #select(system) %>%
#   as.matrix()
# 
# # Center and Scale: BUT now center by the mean of the original modeled dataset above AND scale by the same 2*SD calculated from the original, modeled dataset above
# X_ordinal_scaled_no_taxa <- scale(X_ordinal_no_taxa, center=apply(X_ordinal, MARGIN = 2, FUN = mean), scale=2*ordinal_sd)
# 
# # Create dataframe for brms
# brms_diesel_dat_no_taxa <- data.frame(X_ordinal_scaled_no_taxa)
# 
# # Make predictions
# #predicted_diesel_dat <- predict(fit_no_na, newdata = brms_new_diesel_data)
# # Use tidybayes instead:
# predicted_diesel_dat_no_taxa <- add_predicted_draws(newdata = brms_diesel_dat_no_taxa, model = fit_diesel_no_taxa)
# 
# # Get point and interval estimates from predicted data
# # Select just the prediction columns
# # Join these with the modeled data (diesel_complete_predictors_no_taxa) to get metadata on taxa/intensity/syste,
# diesel_dat_intervals_no_taxa <- predicted_diesel_dat_no_taxa %>%
#   median_qi(.value = .prediction) %>% # Rename prediction to value
#   ungroup() %>%
#   select(contains("."))
# 
# # .row is equivalent to the row number in the modeled dataset (diesel_complete_predictors_no_taxa) - create a join column for this
# diesel_metadat_no_taxa <- diesel_complete_predictors_no_taxa %>%
#   select(study_id, Country, iso3c, clean_sci_name, taxa, intensity, system) %>%
#   mutate(.row = row_number())
# 
# diesel_predictions_no_taxa <- diesel_dat_intervals_no_taxa %>%
#   left_join(diesel_metadat_no_taxa, by = ".row") %>%
#   rename(diesel = .value)

#######################################################################################################
# Bind diesel_no_na (data), diesel_predictions, and diesel_predictions_no_taxa
full_diesel_dat <- diesel_predictions %>%
  bind_rows(diesel_no_na) %>%
  #bind_rows(diesel_predictions_no_taxa) %>%
  mutate(data_type = if_else(is.na(.point), true = "data", false = "prediction")) %>%
  arrange(taxa, intensity, system, clean_sci_name) %>%
  rownames_to_column() # Arrange by taxa first, then create dummy column for plotting 

# Quick check: PLOT DATA + PREDICTIONS
ggplot(full_diesel_dat, aes(x = diesel, y = taxa)) + geom_boxplot(outlier.shape = NA) + geom_jitter(aes(color = data_type))

## FIX IT - temporary fix - filter out predictions with high uncertainty based on highest and lowest data values
## shouldn't be predicting out of the range of observations
min_dat <- full_diesel_dat %>%
  filter(data_type == "data") %>%
  pull(diesel) %>%
  min()

max_dat <- full_diesel_dat %>%
  filter(data_type == "data") %>%
  pull(diesel) %>%
  max()

full_diesel_dat <- full_diesel_dat %>%
  filter(diesel <= max_dat & diesel >= min_dat)

# Check again: PLOT DATA + PREDICTIONS
ggplot(full_diesel_dat, aes(x = diesel, y = taxa)) + geom_boxplot(outlier.shape = NA) + geom_jitter(aes(color = data_type))

# NEXT: use plot_brms_gamma_regression to produce figures for SI

# Clear all memory except for final stan model:
rm(list=ls()[!(ls() %in% c("datadir", "outdir", 
                           "lca_dat_clean_groups","ghg_model_dat_categories", 
                           "electric_fp_dat", "other_energy_fp_dat", 
                           "full_electric_dat", "full_diesel_dat"))])

######################################################################################################
# Step 3: Model petrol

petrol_no_na <- ghg_model_dat_categories %>%
  filter(is.na(petrol)==FALSE)  %>% # Drop NAs (Keep zeroes)
  #mutate(petrol = if_else(petrol == 0, true = min(petrol[petrol!=0]), false = petrol)) %>% # Not modeling the zeroes, option 1: adjust these to the minimum value
  filter(petrol != 0) %>% # Not modeling the zeroes, option 2: drop zeroes
  filter(is.na(intensity)==FALSE & is.na(system)==FALSE) %>% # complete predictors - i.e., both intensity AND system are non-NA
  select(study_id, Country, iso3c, petrol, clean_sci_name, taxa, intensity, system)

# Create model matrix for taxa info, then center and scale
X_taxa <- model.matrix(object = ~ 1 + taxa, 
                       data = petrol_no_na %>% select(taxa)) 

taxa_sd <- apply(X_taxa[,-1], MARGIN=2, FUN=sd, na.rm=TRUE) # Center all non-intercept variables and scale by 2 standard deviations (ignoring NAs)
X_taxa_scaled <- scale(X_taxa[,-1], center=TRUE, scale=2*taxa_sd)

# Format intensity and system as ordinal variables, then center and scale
X_ordinal <- petrol_no_na %>%
  mutate(intensity = factor(intensity, levels = c("Extensive", "Imp. extensive", "Semi-intensive", "Intensive"))) %>% # set order of factors (low = extensive, high = intensive)
  mutate(system = factor(system, levels = c("On- and off-bottom", "Cages & pens", "Ponds", "Recirculating and tanks"))) %>% # set order of factors (low = open, high = closed)
  mutate(intensity = as.numeric(intensity)) %>%
  mutate(system = as.numeric(system)) %>%
  select(intensity, system) %>%
  #select(system) %>%
  as.matrix()
ordinal_sd<-apply(X_ordinal, MARGIN=2, FUN=sd, na.rm=TRUE) # Center all non-intercept variables and scale by 2 standard deviations (ignoring NAs)
X_ordinal_scaled <- scale(X_ordinal, center=TRUE, scale=2*ordinal_sd)

# Create dataframe for brms and rename feed variables
petrol_brms_data <- data.frame(y = petrol_no_na$petrol, X_taxa_scaled, X_ordinal_scaled)

names(petrol_brms_data)

# Set model formula
petrol_brms <- brmsformula(y ~ 1 + ., family = Gamma("log"))

# Use "resp = <response_variable>" to specify different priors for different response variables
all_priors <- c(set_prior("normal(0,5)", class = "b"), # priors for y response variables
                set_prior("normal(0,2.5)", class = "Intercept"), 
                set_prior("exponential(1)", class = "shape"))

# Model converges after increasing the adapt_delta and iterations from default values
# Rule of thumb: bulk and tail effective sample sizes should be 100 x number of chains (i.e., at least 400)
# increasing max_treedepth is more about efficiency (instead of validity)
# See: https://mc-stan.org/misc/warnings.html
fit_petrol_no_na <- brm(petrol_brms, data = petrol_brms_data,
                        prior = all_priors, cores = 4, seed = "11729", iter = 5000, control = list(adapt_delta = 0.99))

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
sort(unique(petrol_no_na$taxa))
petrol_complete_predictors <- petrol_complete_predictors %>%
  mutate(taxa = as.factor(taxa))
levels(petrol_complete_predictors$taxa) <- list(bivalves = "bivalves",
                                                catfish = "catfish",
                                                hypoph_carp = "hypoph_carp",
                                                milkfish = "milkfish",
                                                misc_diad = "misc_diad",
                                                misc_marine = "misc_marine",
                                                oth_carp = "oth_carp",
                                                plants = "plants",
                                                salmon = "salmon",
                                                shrimp = "shrimp",
                                                tilapia = "tilapia",
                                                trout = "trout")

# Create NEW taxa model matrix for the studies to be predicted
# Taxa categories:
X_taxa_new <- model.matrix(object = ~ 1 + taxa, 
                           data = petrol_complete_predictors %>% select(taxa)) 

# Center and Scale: BUT now center by the mean of the original modeled dataset above AND scale by the same 2*SD calculated from the original, modeled dataset above
X_taxa_new_scaled <- scale(X_taxa_new[,-1], center=apply(X_taxa[,-1], MARGIN = 2, FUN = mean), scale=2*taxa_sd)

# Format intensity and system as ordinal variable, then center and scale
# System and Intensity variables:
# Format intensity and system as ordinal variables, then center and scale
X_ordinal_new <- petrol_complete_predictors %>%
  mutate(intensity = factor(intensity, levels = c("Extensive", "Imp. extensive", "Semi-intensive", "Intensive"))) %>% # set order of factors (low = extensive, high = intensive)
  mutate(system = factor(system, levels = c("On- and off-bottom", "Cages & pens", "Ponds", "Recirculating and tanks"))) %>% # set order of factors (low = open, high = closed)
  mutate(intensity = as.numeric(intensity)) %>%
  mutate(system = as.numeric(system)) %>%
  select(intensity, system) %>%
  #select(system) %>%
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
  select(study_id, Country, iso3c, clean_sci_name, taxa, intensity, system) %>%
  mutate(.row = row_number())

petrol_predictions <- petrol_dat_intervals %>%
  left_join(petrol_metadat, by = ".row") %>%
  rename(petrol = .value)

######################################################################################################
# Model petrol with intensity + system (no taxa)
# FOLLOWING SECTION NOT APPLICABLE

# Create dataframe for brms and rename feed variables
# petrol_brms_data_no_taxa <- data.frame(y = petrol_no_na$petrol, X_ordinal_scaled)
# 
# names(petrol_brms_data_no_taxa)
# 
# # Set model formula
# petrol_brms_no_taxa <- brmsformula(y ~ 1 + ., family = Gamma("log"))
# 
# # Use "resp = <response_variable>" to specify different priors for different response variables
# all_priors <- c(set_prior("normal(0,5)", class = "b"), # priors for y response variables
#                 set_prior("normal(0,2.5)", class = "Intercept"), 
#                 set_prior("exponential(1)", class = "shape"))
# 
# # Model converges after increasing the adapt_delta and iterations from default values
# # Rule of thumb: bulk and tail effective sample sizes should be 100 x number of chains (i.e., at least 400)
# # increasing max_treedepth is more about efficiency (instead of validity)
# # See: https://mc-stan.org/misc/warnings.html
# fit_petrol_no_taxa <- brm(petrol_brms_no_taxa, data = petrol_brms_data_no_taxa,
#                           prior = all_priors, cores = 4, seed = "11729", iter = 5000, control = list(adapt_delta = 0.99))
# 
# # Get stan code
# #stancode(fit_petrol_no_taxa)
# 
# ######################################################################################################
# # Use intensity + system model to predict petrol for taxa that were not part of the previous model
# # Intensity must be non-NA
# petrol_complete_predictors <- ghg_model_dat_categories %>%
#   filter(is.na(intensity)==FALSE & is.na(system)==FALSE) %>%
#   filter(is.na(petrol)) # Now, filter to just the NAs
# 
# taxa_not_modeled <- setdiff(unique(petrol_complete_predictors$taxa), unique(petrol_no_na$taxa)) # these taxa were never modeled so they can't be predicted below
# 
# # Only keep taxa that were not modeled previously
# petrol_complete_predictors_no_taxa <- petrol_complete_predictors %>%
#   filter(taxa %in% taxa_not_modeled) 
# 
# # Format intensity and system as ordinal variable, then center and scale
# X_ordinal_no_taxa <- petrol_complete_predictors_no_taxa %>%
#   mutate(intensity = factor(intensity, levels = c("Extensive", "Imp. extensive", "Semi-intensive", "Intensive"))) %>% # set order of factors (low = extensive, high = intensive)
#   mutate(intensity = as.numeric(intensity)) %>%
#   mutate(system = factor(system, levels = c("On- and off-bottom", "Cages & pens", "Ponds", "Recirculating and tanks"))) %>% # set order of factors (low = open, high = closed)
#   mutate(system = as.numeric(system)) %>%
#   select(intensity, system) %>%
#   #select(system) %>%
#   as.matrix()
# 
# # Center and Scale: BUT now center by the mean of the original modeled dataset above AND scale by the same 2*SD calculated from the original, modeled dataset above
# X_ordinal_scaled_no_taxa <- scale(X_ordinal_no_taxa, center=apply(X_ordinal, MARGIN = 2, FUN = mean), scale=2*ordinal_sd)
# 
# # Create dataframe for brms
# brms_petrol_dat_no_taxa <- data.frame(X_ordinal_scaled_no_taxa)
# 
# # Make predictions
# #predicted_petrol_dat <- predict(fit_no_na, newdata = brms_new_petrol_data)
# # Use tidybayes instead:
# predicted_petrol_dat_no_taxa <- add_predicted_draws(newdata = brms_petrol_dat_no_taxa, model = fit_petrol_no_taxa)
# 
# # Get point and interval estimates from predicted data
# # Select just the prediction columns
# # Join these with the modeled data (petrol_complete_predictors_no_taxa) to get metadata on taxa/intensity/syste,
# petrol_dat_intervals_no_taxa <- predicted_petrol_dat_no_taxa %>%
#   median_qi(.value = .prediction) %>% # Rename prediction to value
#   ungroup() %>%
#   select(contains("."))
# 
# # .row is equivalent to the row number in the modeled dataset (petrol_complete_predictors_no_taxa) - create a join column for this
# petrol_metadat_no_taxa <- petrol_complete_predictors_no_taxa %>%
#   select(study_id, Country, iso3c, clean_sci_name, taxa, intensity, system) %>%
#   mutate(.row = row_number())
# 
# petrol_predictions_no_taxa <- petrol_dat_intervals_no_taxa %>%
#   left_join(petrol_metadat_no_taxa, by = ".row") %>%
#   rename(petrol = .value)

######################################################################################################
# Bind petrol_no_na (data), petrol_predictions, and petrol_predictions_no_taxa
full_petrol_dat <- petrol_predictions %>%
  bind_rows(petrol_no_na) %>%
  #bind_rows(petrol_predictions_no_taxa) %>%
  mutate(data_type = if_else(is.na(.point), true = "data", false = "prediction")) %>%
  arrange(taxa, intensity, system, clean_sci_name) %>%
  rownames_to_column() # Arrange by taxa first, then create dummy column for plotting 

# Quick Check: PLOT DATA + PREDICTIONS
ggplot(full_petrol_dat, aes(x = petrol, y = taxa)) + geom_boxplot(outlier.shape = NA) + geom_jitter(aes(color = data_type))

# NEXT: use plot_brms_gamma_regression to produce figures for SI

# Clear all memory except for final stan model:
rm(list=ls()[!(ls() %in% c("datadir", "outdir", 
                           "lca_dat_clean_groups","ghg_model_dat_categories", 
                           "electric_fp_dat", "other_energy_fp_dat", 
                           "full_electric_dat", "full_diesel_dat", "full_petrol_dat"))])
######################################################################################################
# Step 4: Model natural gas

natgas_no_na <- ghg_model_dat_categories %>%
  filter(is.na(natgas)==FALSE)  %>% # Drop NAs (Keep zeroes)
  #mutate(natgas = if_else(natgas == 0, true = min(natgas[natgas!=0]), false = natgas)) %>% # Not modeling the zeroes, option 1: adjust these to the minimum value
  filter(natgas != 0) %>% # Not modeling the zeroes option 2: drop zeroes
  filter(is.na(intensity)==FALSE & is.na(system)==FALSE) %>% # complete predictors - i.e., both intensity AND system are non-NA
  select(study_id, Country, iso3c, natgas, clean_sci_name, taxa, intensity, system)

# Create model matrix for taxa info, then center and scale
X_taxa <- model.matrix(object = ~ 1 + taxa, 
                       data = natgas_no_na %>% select(taxa)) 

# If only one column, no need to use the apply() function
taxa_sd <- apply(X_taxa[,-1], MARGIN=2, FUN=sd, na.rm=TRUE) # Center all non-intercept variables and scale by 2 standard deviations (ignoring NAs)
#taxa_sd <- sd(X_taxa[,-1], na.rm = TRUE)
X_taxa_scaled <- scale(X_taxa[,-1], center=TRUE, scale=2*taxa_sd)

# For natural gas, only select predictors system since there's no variation in intensity
# Format intensity and system as ordinal variables, then center and scale
X_ordinal <- natgas_no_na %>%
  mutate(intensity = factor(intensity, levels = c("Extensive", "Imp. extensive", "Semi-intensive", "Intensive"))) %>% # set order of factors (low = extensive, high = intensive)
  mutate(system = factor(system, levels = c("On- and off-bottom", "Cages & pens", "Ponds", "Recirculating and tanks"))) %>% # set order of factors (low = open, high = closed)
  mutate(intensity = as.numeric(intensity)) %>%
  mutate(system = as.numeric(system)) %>%
  select(intensity, system) %>%
  #select(system) %>%
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
                        prior = all_priors, cores = 4, seed = "11729", iter = 5000, control = list(adapt_delta = 0.99))

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
sort(unique(natgas_no_na$taxa))
natgas_complete_predictors <- natgas_complete_predictors %>%
  mutate(taxa = as.factor(taxa))
levels(natgas_complete_predictors$taxa) <- list(bivalves = "bivalves",
                                                catfish = "catfish",
                                                hypoph_carp = "hypoph_carp",
                                                milkfish = "milkfish",
                                                misc_diad = "misc_diad",
                                                misc_marine = "misc_marine",
                                                oth_carp = "oth_carp",
                                                plants = "plants",
                                                salmon = "salmon",
                                                shrimp = "shrimp",
                                                tilapia = "tilapia",
                                                trout = "trout")


# Create NEW taxa model matrix for the studies to be predicted
# Taxa categories:
X_taxa_new <- model.matrix(object = ~ 1 + taxa, 
                           data = natgas_complete_predictors %>% select(taxa)) 

# There's only one taxa column so don't need to use function apply()
# Center and Scale: BUT now center by the mean of the original modeled dataset above AND scale by the same 2*SD calculated from the original, modeled dataset above
X_taxa_new_scaled <- scale(X_taxa_new[,-1], center=apply(X_taxa[,-1], MARGIN = 2, FUN = mean), scale=2*taxa_sd)

# System and Intensity variables:
# Format intensity and system as ordinal variables, then center and scale
X_ordinal_new <- natgas_complete_predictors %>%
  mutate(intensity = factor(intensity, levels = c("Extensive", "Imp. extensive", "Semi-intensive", "Intensive"))) %>% # set order of factors (low = extensive, high = intensive)
  mutate(system = factor(system, levels = c("On- and off-bottom", "Cages & pens", "Ponds", "Recirculating and tanks"))) %>% # set order of factors (low = open, high = closed)
  mutate(intensity = as.numeric(intensity)) %>%
  mutate(system = as.numeric(system)) %>%
  select(intensity, system) %>%
  #select(system) %>%
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
  select(study_id, Country, iso3c, clean_sci_name, taxa, intensity, system) %>%
  mutate(.row = row_number())

natgas_predictions <- natgas_dat_intervals %>%
  left_join(natgas_metadat, by = ".row") %>%
  rename(natgas = .value)

######################################################################################################
# Model natgas with intensity + system (no taxa)
# FOLLOWING SECTION NOT APPLICABLE

# Create dataframe for brms and rename feed variables
# natgas_brms_data_no_taxa <- data.frame(y = natgas_no_na$natgas, X_ordinal_scaled)
# 
# names(natgas_brms_data_no_taxa)
# 
# # Set model formula
# natgas_brms_no_taxa <- brmsformula(y ~ 1 + ., family = Gamma("log"))
# 
# # Use "resp = <response_variable>" to specify different priors for different response variables
# all_priors <- c(set_prior("normal(0,5)", class = "b"), # priors for y response variables
#                 set_prior("normal(0,2.5)", class = "Intercept"), 
#                 set_prior("exponential(1)", class = "shape"))
# 
# # Model converges after increasing the adapt_delta and iterations from default values
# # Rule of thumb: bulk and tail effective sample sizes should be 100 x number of chains (i.e., at least 400)
# # increasing max_treedepth is more about efficiency (instead of validity)
# # See: https://mc-stan.org/misc/warnings.html
# fit_natgas_no_taxa <- brm(natgas_brms_no_taxa, data = natgas_brms_data_no_taxa,
#                           prior = all_priors, cores = 4, seed = "11729", iter = 5000, control = list(adapt_delta = 0.99))
# 
# # Get stan code
# #stancode(fit_natgas_no_taxa)
# 
# ######################################################################################################
# # Use intensity + system model to predict natgas for taxa that were not part of the previous model
# # Intensity must be non-NA
# natgas_complete_predictors <- ghg_model_dat_categories %>%
#   filter(is.na(intensity)==FALSE & is.na(system)==FALSE) %>%
#   filter(is.na(natgas)) # Now, filter to just the NAs
# 
# taxa_not_modeled <- setdiff(unique(natgas_complete_predictors$taxa), unique(natgas_no_na$taxa)) # these taxa were never modeled so they can't be predicted below
# 
# # Only keep taxa that were not modeled previously
# natgas_complete_predictors_no_taxa <- natgas_complete_predictors %>%
#   filter(taxa %in% taxa_not_modeled) 
# 
# # REMOVE INTENSITY - If predictor model did not include intensity or system, make sure not to include them here
# # Format intensity and system as ordinal variable, then center and scale
# X_ordinal_no_taxa <- natgas_complete_predictors_no_taxa %>%
#   mutate(intensity = factor(intensity, levels = c("Extensive", "Imp. extensive", "Semi-intensive", "Intensive"))) %>% # set order of factors (low = extensive, high = intensive)
#   mutate(intensity = as.numeric(intensity)) %>%
#   mutate(system = factor(system, levels = c("On- and off-bottom", "Cages & pens", "Ponds", "Recirculating and tanks"))) %>% # set order of factors (low = open, high = closed)
#   mutate(system = as.numeric(system)) %>%
#   #select(intensity, system) %>%
#   select(system) %>%
#   as.matrix()
# 
# # Center and Scale: BUT now center by the mean of the original modeled dataset above AND scale by the same 2*SD calculated from the original, modeled dataset above
# X_ordinal_scaled_no_taxa <- scale(X_ordinal_no_taxa, center=apply(X_ordinal, MARGIN = 2, FUN = mean), scale=2*ordinal_sd)
# 
# # Create dataframe for brms
# brms_natgas_dat_no_taxa <- data.frame(X_ordinal_scaled_no_taxa)
# 
# # Make predictions
# #predicted_natgas_dat <- predict(fit_no_na, newdata = brms_new_natgas_data)
# # Use tidybayes instead:
# predicted_natgas_dat_no_taxa <- add_predicted_draws(newdata = brms_natgas_dat_no_taxa, model = fit_natgas_no_taxa)
# 
# # Get point and interval estimates from predicted data
# # Select just the prediction columns
# # Join these with the modeled data (natgas_complete_predictors_no_taxa) to get metadata on taxa/intensity/syste,
# natgas_dat_intervals_no_taxa <- predicted_natgas_dat_no_taxa %>%
#   median_qi(.value = .prediction) %>% # Rename prediction to value
#   ungroup() %>%
#   select(contains("."))
# 
# # .row is equivalent to the row number in the modeled dataset (natgas_complete_predictors_no_taxa) - create a join column for this
# natgas_metadat_no_taxa <- natgas_complete_predictors_no_taxa %>%
#   select(study_id, Country, iso3c, clean_sci_name, taxa, intensity, system) %>%
#   mutate(.row = row_number())
# 
# natgas_predictions_no_taxa <- natgas_dat_intervals_no_taxa %>%
#   left_join(natgas_metadat_no_taxa, by = ".row") %>%
#   rename(natgas = .value)

######################################################################################################
# Bind natgas_no_na (data), natgas_predictions, and natgas_predictions_no_taxa
full_natgas_dat <- natgas_predictions %>%
  bind_rows(natgas_no_na) %>%
  #bind_rows(natgas_predictions_no_taxa) %>%
  mutate(data_type = if_else(is.na(.point), true = "data", false = "prediction")) %>%
  arrange(taxa, intensity, system, clean_sci_name) %>%
  rownames_to_column() # Arrange by taxa first, then create dummy column for plotting 

# PLOT DATA + PREDICTIONS
ggplot(full_natgas_dat, aes(x = natgas, y = taxa)) + geom_boxplot(outlier.shape = NA) + geom_jitter(aes(color = data_type))

# NEXT: use plot_brms_gamma_regression to produce figures for SI

# Clear all memory except for final stan model:
rm(list=ls()[!(ls() %in% c("datadir", "outdir", 
                           "lca_dat_clean_groups","ghg_model_dat_categories", 
                           "electric_fp_dat", "other_energy_fp_dat", 
                           "full_electric_dat", "full_diesel_dat", "full_petrol_dat", "full_natgas_dat"))])



#save.image(file.path(outdir, paste(Sys.Date(), "_on-farm-carbon-all-data-prior-to-aggregation.RData", sep = "")))
# datadir <- "/Volumes/jgephart/BFA Environment 2/Data"
# outdir <- "/Volumes/jgephart/BFA Environment 2/Outputs"
# load(file.path(outdir, "2020-12-16_on-farm-ghg-all-data-prior-to-aggregation.RData"))

######################################################################################################
# RESTARTING POINT

# Load on-farm carbon (GHG) variables
# load(file.path(outdir, "2020-12-16_on-farm-ghg-all-data-prior-to-aggregation.RData"))

# Load off-farm (feed-associated) variables
load(file.path(outdir, "2020-12-17_feed-impacts-all-data-prior-to-aggregation.RData"))
######################################################################################################

# STEP 3 - aggregate up to taxa level and estimate total feed footprint
# Set data for model

# FIX IT - change clean.lca in Functions.R so no less than 0.01
# Adjust feed proportions to be no less than 0.01
# TEMPORARY FIX:
# Normalize the FINAL feed proportion values to be greater than 0 and no less than 0.01
full_feed_dat <- full_feed_dat %>%
  mutate(feed_proportion = if_else(feed_proportion < 0.01, true = 0.0105, false = feed_proportion))

# Merge fcr and feed datasets by study_id, remove interval info, only use medians for rest of analysis
# Format feed_dat
# Standardize number of significant digits (round to digist = 3)
# Re-scale so rowsums equals 1
feed_dat_merge <- full_feed_dat %>%
  select(study_id, clean_sci_name, taxa, intensity, system, feed_data_type = data_type, .category, contains("feed")) %>%
  pivot_wider(names_from = .category, values_from = feed_proportion) %>%
  mutate(feed_soy = round(feed_soy, digits = 4),
         feed_crops = round(feed_crops, digits = 4),
         feed_fmfo = round(feed_fmfo, digits = 4),
         feed_animal = round(feed_animal, digits = 4)) %>%
  mutate(new_scale = feed_soy + feed_crops + feed_fmfo + feed_animal) %>%
  mutate(feed_soy = feed_soy / new_scale,
         feed_crops = feed_crops / new_scale,
         feed_fmfo = feed_fmfo / new_scale,
         feed_animal = feed_animal / new_scale)



# Format FCR dat
# Rename data_type to keep track of both feed vs fcr data types
fcr_dat_merge <- full_fcr_dat %>%
  select(study_id, clean_sci_name, taxa, intensity, system, fcr_data_type = data_type, fcr)

# MERGE
feed_footprint_dat <- feed_dat_merge %>%
  full_join(fcr_dat_merge, by = intersect(names(feed_dat_merge), names(fcr_dat_merge))) %>%
  drop_na() %>%  # Make sure to drop na before creating sci and tx indices (otherwise some indices drop out)
  arrange(clean_sci_name, taxa) %>%
  mutate(clean_sci_name = as.factor(clean_sci_name),
         sci = as.numeric(clean_sci_name),
         taxa = as.factor(taxa),
         tx = as.numeric(taxa)) 

# Set data
# overall model:
N = nrow(feed_footprint_dat)
N_SCI <- length(unique(feed_footprint_dat$sci))
n_to_sci <- feed_footprint_dat$sci
N_TX <- length(unique(feed_footprint_dat$tx))
sci_to_tx = feed_footprint_dat %>%
  select(sci, tx) %>%
  unique() %>%
  pull(tx)


# for FCR model:
x <- feed_footprint_dat$fcr

# for Feed proportion model:
K = 4
feed_weights <- feed_footprint_dat %>%
  select(feed_soy, feed_crops, feed_fmfo, feed_animal) %>%
  as.matrix()

# Get counts per sci name and counts per taxa group (also included as data in the model):
sci_kappa <- feed_footprint_dat %>% 
  group_by(sci) %>% 
  summarise(n_obs = n()) %>%
  ungroup() %>%
  arrange(sci) %>%
  pull(n_obs)

tx_kappa <- feed_footprint_dat %>% 
  group_by(tx) %>% 
  summarise(n_obs = n()) %>%
  ungroup() %>%
  arrange(tx) %>%
  pull(n_obs)

# Get footprint data
# DECISION: Choose allocation method
# IMPORTANT - multiply all values by 1000 to convert to kg CO2 per tonne (currently in kg CO2 per kg)

set_allocation <- "Mass"

fp_dat <- read.csv(file.path(datadir, "20201217_weighted_feed_fp.csv")) %>%
  filter(Allocation == set_allocation) %>%
  mutate(ave_stressor_per_tonne = ave_stressor * 1000)

# ORDER OF feed_weights data vector: soy, crops, fmfo, animal
head(feed_weights)
# ORDER of footprint data vector should match this:
set_fp_order <- c("Soy", "Crop", "Fishery", "Livestock")

fp_c_dat <- fp_dat %>%
  filter(Impact.category == "Global warming potential") %>% 
  arrange(match(Input.type, set_fp_order)) %>% # Match index and arrange by custom order
  select(ave_stressor_per_tonne) %>%
  as.matrix() %>%
  c()

fp_n_dat <- fp_dat %>%
  filter(Impact.category == "Marine eutrophication") %>% 
  arrange(match(Input.type, set_fp_order)) %>% # Match index and arrange by custom order
  select(ave_stressor_per_tonne) %>%
  as.matrix() %>%
  c()

fp_p_dat <- fp_dat %>%
  filter(Impact.category == "Freshwater eutrophication") %>% 
  arrange(match(Input.type, set_fp_order)) %>% # Match index and arrange by custom order
  select(ave_stressor_per_tonne) %>%
  as.matrix() %>%
  c()

fp_land_dat <- fp_dat %>%
  filter(Impact.category == "Land use") %>% 
  arrange(match(Input.type, set_fp_order)) %>% # Match index and arrange by custom order
  select(ave_stressor_per_tonne) %>%
  as.matrix() %>%
  c()

fp_water_dat <- fp_dat %>%
  filter(Impact.category == "Water consumption") %>% 
  arrange(match(Input.type, set_fp_order)) %>% # Match index and arrange by custom order
  select(ave_stressor_per_tonne) %>%
  as.matrix() %>%
  c()

# Get sci-level weightings for generating taxa-level quantities:
# IMPORTANT arrange by clean_sci_name so that the order matches data
sci_prod_weights <- read.csv(file.path(datadir, "aqua_prod_weightings.csv")) %>%
  arrange(clean_sci_name)

# Drop sci-names not found in the data
sci_prod_weights <- sci_prod_weights %>%
  filter(clean_sci_name %in% feed_footprint_dat$clean_sci_name)

# Check that the order of sci names in both weights and data are the same (sum = 0)
sum(sci_prod_weights$clean_sci_name != unique(feed_footprint_dat$clean_sci_name))
sci_w <- sci_prod_weights$weighting


# Need the following for generating ragged array in STAN - indexing vector of sci_mu's based on their taxa identity (each slice will have a different length)
where_tx <- order(sci_to_tx) # i.e., give the positions in sci_to_tx in order of their taxa level (where in sci_to_tx is taxa level 1, followed by taxa level 2, etc)
# How many sci are in each taxa - need this to declare length of each sci_mu vector
n_sci_in_tx <- feed_footprint_dat %>%
  select(sci, tx) %>%
  unique() %>%
  group_by(tx) %>%
  mutate(n_sci_in_tx = n()) %>%
  ungroup() %>%
  arrange(tx) %>%
  select(tx, n_sci_in_tx) %>%
  unique() %>%
  pull(n_sci_in_tx)
slice_where_tx <- cumsum(n_sci_in_tx) # These are the breaks in where_tx corresponding to each taxa level - need this to split up where_tx
slice_where_tx <- c(0, slice_where_tx)

stan_data <- list(N = N,
                  N_SCI = N_SCI, 
                  n_to_sci = n_to_sci,
                  N_TX = N_TX,
                  sci_to_tx = sci_to_tx,
                  x = x,
                  K = K,
                  feed_weights = feed_weights,
                  sci_kappa = sci_kappa,
                  tx_kappa = tx_kappa,
                  fp_c_dat = fp_c_dat,
                  fp_n_dat = fp_n_dat,
                  fp_p_dat = fp_p_dat,
                  fp_land_dat = fp_land_dat,
                  fp_water_dat = fp_water_dat, 
                  sci_w = sci_w,
                  where_tx = where_tx,
                  n_sci_in_tx = n_sci_in_tx,
                  slice_where_tx = slice_where_tx)

# Estimate foot print for all scientific names and taxa groups (removed the "all-seafood" level for simplicity)
# GAMMA distribution hierarchical model
# stan_no_na <- 'data {
#   // data for gamma model for FCR
#   int<lower=0> N;  // number of observations
#   vector<lower=0>[N] x; // data
#   int N_TX; // number of taxa groups
#   int N_SCI; // number of scientific names
#   int n_to_sci[N]; // sciname index
#   int sci_to_tx[N_SCI]; // taxa-group indices
# 
#   // data for dirichlet model for feed
#   int K; // number of feed types
#   simplex[K] feed_weights[N]; // array of observed feed weights simplexes
#   int sci_kappa[N_SCI]; // number of observations per sci-name
#   int tx_kappa[N_TX]; // number of observations per taxa group
# 
#   // constants for generated quantities
#   vector[K] fp_c_dat;
#   vector[K] fp_n_dat;
#   vector[K] fp_p_dat;
#   vector[K] fp_land_dat;
#   vector[K] fp_water_dat;
# }
# parameters {
#   // FCR model:
#   vector<lower=0>[N_TX] tx_mu;
#   vector<lower=0>[N_SCI] sci_mu;
#   // only need sigmas if defining shape and rate with mu and sigma
#   //real<lower=0> tx_sigma;
#   //real<lower=0> sci_sigma;
#   // if using variance instead of st dev
#   real<lower=0> tx_sigma_sq;
#   real<lower=0> sci_sigma_sq;
# 
#   // Feed proportion model:
#   simplex[K] sci_phi[N_SCI];
#   simplex[K] sci_theta[N_SCI]; // vectors of estimated sci-level feed weight simplexes
#   simplex[K] tx_phi[N_TX];
#   simplex[K] tx_theta[N_TX];
# 
#   // Params for the dirichlet priors:
#   // real<lower=0> sigma_1;
#   // real<lower=0> sigma_2;
# }
# transformed parameters {
#   // define transofrmed params for gamma model for FCRs
#   vector<lower=0>[N_SCI] sci_shape;
#   vector<lower=0>[N_SCI] sci_rate;
#   vector<lower=0>[N_TX] tx_shape;
#   vector<lower=0>[N_TX] tx_rate;
# 
#   // define params for dirichlet model for feed proportions
#   vector<lower=0>[K] sci_alpha[N_SCI];
#   vector<lower=0>[K] tx_alpha[N_TX];
# 
#   // gamma model reparameterization
#   // option 1: reparamaterize gamma to get mu and sigma; defining these here instead of the model section allows us to see these parameters in the output
#   // taxa group level
#   //for (n_tx in 1:N_TX){
#   //  tx_shape[n_tx] = square(tx_mu[n_tx]) ./ square(tx_sigma);
#   //  tx_rate[n_tx] = tx_mu[n_tx] ./ square(tx_sigma);
#   //}
# 
#   // sci level
#   //for (n_sci in 1:N_SCI){
#   //  sci_shape[n_sci] = square(sci_mu[n_sci]) ./ square(sci_sigma);
#   //  sci_rate[n_sci] = sci_mu[n_sci] ./ square(sci_sigma);
#   //}
# 
#   // option 2: reparameterize gamma to get just mu
#   //for (n_tx in 1:N_TX){
#   //  tx_rate[n_tx] = tx_shape[n_tx] ./ tx_mu[n_tx];
#   //}
#   // sci level
#   //for (n_sci in 1:N_SCI){
#   //  sci_rate[n_sci] = sci_shape[n_sci] ./ sci_mu[n_sci];
#   //}
#   
#   // option 3: reparameterize shape and rate as inverse(va) and inverse(va)/mu
#   for (n_tx in 1:N_TX){
#     tx_shape[n_tx] = 1 / tx_sigma_sq;
#     tx_rate[n_tx] = 1 / (tx_sigma_sq * tx_mu[n_tx]);
#   }
#   for (n_sci in 1:N_SCI){
#     sci_shape[n_sci] = 1 / sci_sigma_sq;
#     sci_rate[n_sci] = 1 / (sci_sigma_sq * sci_mu[n_sci]);
#   }
#   
#   
#   // dirichlet model reparameterization
#   // reparameterize alphas as a vector of means (phi) and counts (kappas)
#   // theta is expected value of mean feed weights
#   // kappa is strength of the prior measured in number of prior observations (minus K)
#   for (n_tx in 1:N_TX) {
#     tx_alpha[n_tx] = tx_kappa[n_tx] * tx_theta[n_tx];
#   }
# 
#   for (n_sci in 1:N_SCI) {
#     sci_alpha[n_sci] = sci_kappa[n_sci] * sci_theta[n_sci];
#   }
# }
# model {
#   // define priors for gamma model for FCRs
#   // Put priors on mu and sigma (instead of shape and rate) since this is more intuitive:
#   tx_mu ~ uniform(0, 100); // note: uniform(0,100) for all of these doesnt help much with convergence
#   sci_mu ~ uniform(0, 100);
#   //tx_sigma ~ uniform(0, 100); // only need sigmas if calculating shape and rate with mu and sigma
#   //sci_sigma ~ uniform(0, 100);
# 
#   // define priors for dirichlet model for feed proportions
#   // sci_phi defined as sci_phi[n_to_sci][K]
# 
#   // option 1: define feed proportion priors as lower upper bounds
#   // sci_phi[2][1] ~ uniform(0.1, 0.2); // hypothetical lower and upper bounds for Oncorhynchus mykiss soy
# 
#   // option 2: define feed proportions as means (need to define sigmas in parameters block: real<lower=0> sigma_1, sigma_2 etc;)
#   // sci_phi[2][1] ~ normal(0.13, sigma_1); // mean for Oncorhynhchus mykiss soy feed
# 
#   // likelihood
#   // gamma model sci-name and taxa-level
#   for (n in 1:N){
#     x[n] ~ gamma(sci_shape[n_to_sci[n]], sci_rate[n_to_sci[n]]);
#   }
# 
#   for (n_sci in 1:N_SCI){
#     sci_mu[n_sci] ~ gamma(tx_shape[sci_to_tx[n_sci]], tx_rate[sci_to_tx[n_sci]]);
#   }
# 
#   // dirichlet model
#   // use data to estimate sci-level dirichlet shape param (alphas)
#   for (n in 1:N) {
#     feed_weights[n] ~ dirichlet(to_vector(sci_alpha[n_to_sci[n]]));
#   }
#   // use sci-level thetas to estimate taxa-level dirichlet shape param (alphas)
#   for (n_sci in 1:N_SCI) {
#     sci_theta[n_sci] ~ dirichlet(to_vector(tx_alpha[sci_to_tx[n_sci]]));
#   }
# }
# generated quantities {
#   // Carbon
#   vector[N_TX] tx_c_fp;
#   vector[K] tx_feed_c_fp[N_TX]; // array of k feeds and N_TX taxa groups
#   vector[N_TX] tx_sum_feed_c_fp;
#   vector[N_SCI] sci_c_fp;
#   vector[K] sci_feed_c_fp[N_SCI];
#   vector[N_SCI] sci_sum_feed_c_fp;
#   // Nitrogen
# 
#   // Phosphorus
# 
#   // Land
# 
#   // Water
# 
#   // Calculations
#   for (n_tx in 1:N_TX) {
#     tx_feed_c_fp[n_tx] = fp_c_dat .* tx_theta[n_tx];
#     tx_sum_feed_c_fp[n_tx] = sum(tx_feed_c_fp[n_tx]);
#     tx_c_fp[n_tx] = tx_mu[n_tx] * tx_sum_feed_c_fp[n_tx];
#   }
#   for (n_sci in 1:N_SCI) {
#     sci_feed_c_fp[n_sci] = fp_c_dat .* sci_theta[n_sci];
#     sci_sum_feed_c_fp[n_sci] = sum(sci_feed_c_fp[n_sci]);
#     sci_c_fp[n_sci] = sci_mu[n_sci] * sci_sum_feed_c_fp[n_sci];
#   }
# }'


# NORMAL distribution hierarchical model
stan_no_na <- 'data {
  // data for normal model for FCR
  int<lower=0> N;  // number of observations
  vector<lower=0>[N] x; // fcr data
  int N_TX; // number of taxa groups
  int N_SCI; // number of scientific names
  int n_to_sci[N]; // sciname index
  int sci_to_tx[N_SCI]; // taxa-group indices

  // data for dirichlet model for feed
  int K; // number of feed types
  simplex[K] feed_weights[N]; // array of observed feed weights simplexes
  int sci_kappa[N_SCI]; // number of observations per sci-name
  int tx_kappa[N_TX]; // number of observations per taxa group

  // constants for generated quantities
  vector[K] fp_c_dat;
  //vector[K] fp_n_dat;
  //vector[K] fp_p_dat;
  //vector[K] fp_land_dat;
  //vector[K] fp_water_dat;
  
  // data for slicing vectors for calculating weighted means
  vector<lower=0>[N_SCI] sci_w; // sci-level production weights
  int where_tx[N_SCI]; // order sci_to_tx by taxa-level
  int n_sci_in_tx[N_TX]; // number of sci in each taxa-level, ordered by taxa-level
  int slice_where_tx[N_TX + 1]; // breaks in where_tx by taxa-level
}
parameters {
  // FCR model:
  vector<lower=0>[N_TX] tx_mu;
  vector<lower=0>[N_SCI] sci_mu;
  vector<lower=0>[N_TX] tx_sigma;
  vector<lower=0>[N_SCI] sci_sigma;
  //real<lower=0> tx_sigma;
  //real<lower=0> sci_sigma;
  

  // Feed proportion model:
  simplex[K] sci_phi[N_SCI];
  simplex[K] sci_theta[N_SCI]; // vectors of estimated sci-level feed weight simplexes
  simplex[K] tx_phi[N_TX];
  simplex[K] tx_theta[N_TX];
}
transformed parameters {
  // define params for dirichlet model for feed proportions
  vector<lower=0>[K] sci_alpha[N_SCI];
  vector<lower=0>[K] tx_alpha[N_TX];

  // dirichlet model reparameterization
  // reparameterize alphas as a vector of means (phi) and counts (kappas)
  // theta is expected value of mean feed weights
  // kappa is strength of the prior measured in number of prior observations (minus K)
  for (n_tx in 1:N_TX) {
    tx_alpha[n_tx] = tx_kappa[n_tx] * tx_theta[n_tx];
  }

  for (n_sci in 1:N_SCI) {
    sci_alpha[n_sci] = sci_kappa[n_sci] * sci_theta[n_sci];
  }
}
model {
  // example priors for gamma model for FCRs
  // Put priors on mu and sigma (instead of shape and rate) since this is more intuitive:
  //tx_mu ~ uniform(0, 100);
  //tx_sigma ~ uniform(0, 100);

  // example priors for dirichlet model for feed proportions
  // sci_phi defined as sci_phi[n_to_sci][K]
  // sci_phi[2][1] ~ normal(0.13, 5); // mean for Oncorhynhchus mykiss soy feed

  // likelihood
  // normal model sci-name and taxa-level
  for (n in 1:N){
    x[n] ~ normal(sci_mu[n_to_sci[n]], sci_sigma[n_to_sci[n]]);
  }

  for (n_sci in 1:N_SCI){
    sci_mu[n_sci] ~ normal(tx_mu[sci_to_tx[n_sci]], tx_sigma[sci_to_tx[n_sci]]);
  }

  // dirichlet model
  // use data to estimate sci-level dirichlet shape param (alphas)
  for (n in 1:N) {
    feed_weights[n] ~ dirichlet(to_vector(sci_alpha[n_to_sci[n]]));
  }
  // use sci-level thetas to estimate taxa-level dirichlet shape param (alphas)
  for (n_sci in 1:N_SCI) {
    sci_theta[n_sci] ~ dirichlet(to_vector(tx_alpha[sci_to_tx[n_sci]]));
  }
}
generated quantities {
  // Carbon
  vector[N_TX] tx_c_fp;
  vector[N_SCI] sci_c_fp;

  // Nitrogen
  //vector[N_TX] tx_n_fp;
  //vector[N_SCI] sci_n_fp;

  // Phosphorus
  //vector[N_TX] tx_p_fp;
  //vector[N_SCI] sci_p_fp;

  // Land
  //vector[N_TX] tx_land_fp;
  //vector[N_SCI] sci_land_fp;

  // Water
  //vector[N_TX] tx_water_fp;
  //vector[N_SCI] sci_water_fp;
  
  // Declare vectors for weightings
// Carbon
vector[N_SCI] sci_mu_w_c; // weighted means
vector[N_TX] tx_out_c; // weighted means pooled to taxa level

// Nitrogen
//vector[N_SCI] sci_mu_w_n; // weighted means
//vector[N_TX] tx_out_n; // weighted means pooled to taxa level

// Phosphorus
//vector[N_SCI] sci_mu_w_p; // weighted means
//vector[N_TX] tx_out_p; // weighted means pooled to taxa level

// Land
//vector[N_SCI] sci_mu_w_land; // weighted means
//vector[N_TX] tx_out_land; // weighted means pooled to taxa level

// Water
//vector[N_SCI] sci_mu_w_water; // weighted means
//vector[N_TX] tx_out_water; // weighted means pooled to taxa level


  // Calculations
  for (n_tx in 1:N_TX) {
    tx_c_fp[n_tx] = tx_mu[n_tx] * sum(fp_c_dat .* tx_theta[n_tx]);
    //tx_n_fp[n_tx] = tx_mu[n_tx] * sum(fp_n_dat .* tx_theta[n_tx]);
    //tx_p_fp[n_tx] = tx_mu[n_tx] * sum(fp_p_dat .* tx_theta[n_tx]);
    //tx_land_fp[n_tx] = tx_mu[n_tx] * sum(fp_land_dat .* tx_theta[n_tx]);
    //tx_water_fp[n_tx] = tx_mu[n_tx] * sum(fp_water_dat .* tx_theta[n_tx]);
  }
  for (n_sci in 1:N_SCI) {
    sci_c_fp[n_sci] = sci_mu[n_sci] * sum(fp_c_dat .* sci_theta[n_sci]);
    //sci_n_fp[n_sci] = sci_mu[n_sci] * sum(fp_n_dat .* sci_theta[n_sci]);
    //sci_p_fp[n_sci] = sci_mu[n_sci] * sum(fp_p_dat .* sci_theta[n_sci]);
    //sci_land_fp[n_sci] = sci_mu[n_sci] * sum(fp_land_dat .* sci_theta[n_sci]);
    //sci_water_fp[n_sci] = sci_mu[n_sci] * sum(fp_water_dat .* sci_theta[n_sci]);
  }

// Apply weightings
sci_mu_w_c = sci_c_fp .* sci_w; // weighted sci-level means
//sci_mu_w_n = sci_n_fp .* sci_w; // weighted sci-level means
//sci_mu_w_p = sci_p_fp .* sci_w; // weighted sci-level means
//sci_mu_w_land = sci_land_fp .* sci_w; // weighted sci-level means
//sci_mu_w_water = sci_water_fp .* sci_w; // weighted sci-level means


for (n_tx in 1:N_TX){
  vector[n_sci_in_tx[n_tx]] sci_mu_w_vec_c; // declare vector of sci_mu in taxa-level n_tx
  //vector[n_sci_in_tx[n_tx]] sci_mu_w_vec_n; // declare vector of sci_mu in taxa-level n_tx
  //vector[n_sci_in_tx[n_tx]] sci_mu_w_vec_p; // declare vector of sci_mu in taxa-level n_tx
  //vector[n_sci_in_tx[n_tx]] sci_mu_w_vec_land; // declare vector of sci_mu in taxa-level n_tx
  //vector[n_sci_in_tx[n_tx]] sci_mu_w_vec_water; // declare vector of sci_mu in taxa-level n_tx
  
  // Carbon
  sci_mu_w_vec_c = sci_mu_w_c[where_tx[slice_where_tx[n_tx]+1:slice_where_tx[n_tx+1]]]; // get all the sci_mu in taxa-level n_tx
  tx_out_c[n_tx] = sum(sci_mu_w_vec_c); // sum sci_mu_w_vec to pool to taxa-level
  
  // Nitrogen
  //sci_mu_w_vec_n = sci_mu_w_n[where_tx[slice_where_tx[n_tx]+1:slice_where_tx[n_tx+1]]]; // get all the sci_mu in taxa-level n_tx
  //tx_out_n[n_tx] = sum(sci_mu_w_vec_n); // sum sci_mu_w_vec to pool to taxa-level
  
  // Phosphorus
  //sci_mu_w_vec_p = sci_mu_w_p[where_tx[slice_where_tx[n_tx]+1:slice_where_tx[n_tx+1]]]; // get all the sci_mu in taxa-level n_tx
  //tx_out_p[n_tx] = sum(sci_mu_w_vec_p); // sum sci_mu_w_vec to pool to taxa-level
  
  // Land
  //sci_mu_w_vec_land = sci_mu_w_land[where_tx[slice_where_tx[n_tx]+1:slice_where_tx[n_tx+1]]]; // get all the sci_mu in taxa-level n_tx
  //tx_out_land[n_tx] = sum(sci_mu_w_vec_land); // sum sci_mu_w_vec to pool to taxa-level
  
  // Water
  //sci_mu_w_vec_water = sci_mu_w_water[where_tx[slice_where_tx[n_tx]+1:slice_where_tx[n_tx+1]]]; // get all the sci_mu in taxa-level n_tx
  //tx_out_water[n_tx] = sum(sci_mu_w_vec_c); // sum sci_mu_w_vec to pool to taxa-level
}

}'


no_na_mod <- stan_model(model_code = stan_no_na)
# Note: For Windows, apparently OK to ignore this warning message:
# Warning message:
#   In system(paste(CXX, ARGS), ignore.stdout = TRUE, ignore.stderr = TRUE) :
#   'C:/rtools40/usr/mingw_/bin/g++' not found

# Fit model:
# Set seed while testing
fit_no_na <- sampling(object = no_na_mod, data = stan_data, cores = 4, seed = "11729", iter = 2500, control = list(adapt_delta = 0.99))
#fit_no_na <- sampling(object = no_na_mod, data = stan_data, cores = 4, iter = 5000, control = list(adapt_delta = 0.99))
summary(fit_no_na)$summary

#launch_shinystan(fit_no_na)
######################################################################################################
# RESTARTING POINT
# rm(list=ls()[!(ls() %in% c("datadir", "outdir", 
#                            "lca_dat_clean_groups", "feed_model_dat_categories",
#                            "full_feed_dat", "full_fcr_dat", "feed_footprint_dat", "fit_no_na"))])
#save.image(file.path(outdir, paste(Sys.Date(), "_feed-impacts_", set_allocation, "-allocation_all-data-prior-to-plotting.RData", sep = "")))

# datadir <- "/Volumes/jgephart/BFA Environment 2/Data"
# outdir <- "/Volumes/jgephart/BFA Environment 2/Outputs"
# load(file.path(outdir, "2020-12-17_feed-impacts_Mass-allocation_all-data-prior-to-plotting.RData"))
# set_allocation <- "Mass"
######################################################################################################
# PLOT RESULTS

# SET THEME
sci_plot_theme <- theme(title = element_text(size = 18),
                    axis.title.x = element_text(size = 16),
                    axis.text=element_text(size=14, color = "black"))
tx_plot_theme <- theme(title = element_text(size = 20),
                        axis.title.x = element_text(size = 20),
                        axis.text=element_text(size=20, color = "black"))

# Key for naming sci and taxa levels
# Get full taxa group names back
index_key <- feed_footprint_dat %>%
  group_by(clean_sci_name) %>%
  mutate(n_obs = n()) %>%
  ungroup() %>%
  select(clean_sci_name, sci, taxa, tx, n_obs) %>%
  unique() %>%
  arrange(taxa) %>%
  mutate(taxa = as.character(taxa),
         full_taxa_name = case_when(taxa == "hypoph_carp" ~ "bighead/silverhead carp",
                                    taxa == "misc_marine" ~ "misc marine fishes",
                                    taxa == "misc_fresh" ~ "misc freshwater fishes",
                                    taxa == "misc_diad" ~ "misc diadromous fishes",
                                    taxa == "fresh_crust" ~ "freshwater crustaceans",
                                    TRUE ~ taxa),
         taxa = as.factor(taxa),
         full_taxa_name = as.factor(full_taxa_name))

# Use tidybayes + ggdist for finer control of aes mapping (instead of bayesplots) 
get_variables(fit_no_na)

# WEIGHTED OUTPUTS:
# Sci-level feed footprints as point intervals:
# Carbon
fit_no_na %>%
  spread_draws(sci_mu_w_c[sci]) %>%
  median_qi(.width = 0.95) %>%
  left_join(index_key, by = "sci") %>% # Join with index key to get sci and taxa names
  mutate(clean_sci_name = fct_reorder(clean_sci_name, as.character(full_taxa_name))) %>%
  ggplot(aes(y = clean_sci_name, x = sci_mu_w_c, xmin = .lower, xmax = .upper, color = full_taxa_name)) +
  geom_pointinterval() +
  theme_classic() + 
  sci_plot_theme + 
  labs(x = "kg CO2-eq per tonne", y = "", title = "Carbon", color = "taxa group")
ggsave(filename = file.path(outdir, paste("plot_feed-impact_carbon_", set_allocation, "-allocation_sci-level-WEIGHTED.png", sep = "")), width = 11, height = 8.5)

# Nitrogen
fit_no_na %>%
  spread_draws(sci_mu_w_n[sci]) %>%
  median_qi(.width = 0.95) %>%
  left_join(index_key, by = "sci") %>% # Join with index key to get sci and taxa names
  mutate(clean_sci_name = fct_reorder(clean_sci_name, as.character(full_taxa_name))) %>%
  ggplot(aes(y = clean_sci_name, x = sci_mu_w_n, xmin = .lower, xmax = .upper, color = full_taxa_name)) +
  geom_pointinterval() +
  theme_classic() + 
  sci_plot_theme + 
  labs(x = "kg N-eq per tonne", y = "", title = "Nitrogen", color = "taxa group")
ggsave(filename = file.path(outdir, paste("plot_feed-impact_nitrogen_", set_allocation, "-allocation_sci-level-WEIGHTED.png", sep = "")), width = 11, height = 8.5)

# Phosphorus
fit_no_na %>%
  spread_draws(sci_mu_w_p[sci]) %>%
  median_qi(.width = 0.95) %>%
  left_join(index_key, by = "sci") %>% # Join with index key to get sci and taxa names
  mutate(clean_sci_name = fct_reorder(clean_sci_name, as.character(full_taxa_name))) %>%
  ggplot(aes(y = clean_sci_name, x = sci_mu_w_p, xmin = .lower, xmax = .upper, color = full_taxa_name)) +
  geom_pointinterval() +
  theme_classic() + 
  sci_plot_theme + 
  labs(x = "kg P-eq per tonne", y = "", title = "Phosphorus", color = "taxa group")
ggsave(filename = file.path(outdir, paste("plot_feed-impact_phosphorus_", set_allocation, "-allocation_sci-level-WEIGHTED.png", sep = "")), width = 11, height = 8.5)

# Land
fit_no_na %>%
  spread_draws(sci_mu_w_land[sci]) %>%
  median_qi(.width = 0.95) %>%
  left_join(index_key, by = "sci") %>% # Join with index key to get sci and taxa names
  mutate(clean_sci_name = fct_reorder(clean_sci_name, as.character(full_taxa_name))) %>%
  ggplot(aes(y = clean_sci_name, x = sci_mu_w_land, xmin = .lower, xmax = .upper, color = full_taxa_name)) +
  geom_pointinterval() +
  theme_classic() + 
  sci_plot_theme + 
  labs(x = bquote('m'^2~'a per tonne'), y = "", title = "Land", color = "taxa group")
ggsave(filename = file.path(outdir, paste("plot_feed-impact_land_", set_allocation, "-allocation_sci-level-WEIGHTED.png", sep = "")), width = 11, height = 8.5)

# Water
fit_no_na %>%
  spread_draws(sci_mu_w_water[sci]) %>%
  median_qi(.width = 0.95) %>%
  left_join(index_key, by = "sci") %>% # Join with index key to get sci and taxa names
  mutate(clean_sci_name = fct_reorder(clean_sci_name, as.character(full_taxa_name))) %>%
  ggplot(aes(y = clean_sci_name, x = sci_mu_w_water, xmin = .lower, xmax = .upper, color = full_taxa_name)) +
  geom_pointinterval() +
  theme_classic() + 
  sci_plot_theme + 
  labs(x = bquote('m'^3~'per tonne'), y = "", title = "Water", color = "taxa group")
ggsave(filename = file.path(outdir, paste("plot_feed-impact_water_", set_allocation, "-allocation_sci-level-WEIGHTED.png", sep = "")), width = 11, height = 8.5)

# Taxa-level feed footprints as point intervals:
# Carbon
fit_no_na %>%
  spread_draws(tx_out_c[tx]) %>%
  median_qi(.width = 0.95) %>%
  left_join(index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  ggplot(aes(y = full_taxa_name, x = tx_out_c, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() +
  #stat_halfeye(aes(slab_fill = full_taxa_name)) +
  theme_classic() + 
  tx_plot_theme + 
  theme(legend.position = "none") +
  labs(x = "kg CO2-eq per tonne", y = "", title = "Carbon")
ggsave(filename = file.path(outdir, paste("plot_feed-impact_carbon_", set_allocation, "-allocation_taxa-level-WEIGHTED.png", sep = "")), width = 11, height = 8.5)

# Nitrogen
fit_no_na %>%
  spread_draws(tx_out_n[tx]) %>%
  median_qi(.width = 0.95) %>%
  left_join(index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  ggplot(aes(y = full_taxa_name, x = tx_out_n, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() +
  #stat_halfeye(aes(slab_fill = full_taxa_name)) +
  theme_classic() + 
  tx_plot_theme + 
  theme(legend.position = "none") +
  labs(x = "kg N-eq per tonne", y = "", title = "Nitrogen")
ggsave(filename = file.path(outdir, paste("plot_feed-impact_nitrogen_", set_allocation, "-allocation_taxa-level-WEIGHTED.png", sep = "")), width = 11, height = 8.5)


# Phosphorus
fit_no_na %>%
  spread_draws(tx_out_p[tx]) %>%
  median_qi(.width = 0.95) %>%
  left_join(index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  ggplot(aes(y = full_taxa_name, x = tx_out_p, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() +
  #stat_halfeye(aes(slab_fill = full_taxa_name)) +
  theme_classic() + 
  tx_plot_theme + 
  theme(legend.position = "none") +
  labs(x = "kg P-eq per tonne", y = "", title = "Phosphorus")
ggsave(filename = file.path(outdir, paste("plot_feed-impact_phosphorus_", set_allocation, "-allocation_taxa-level-WEIGHTED.png", sep = "")), width = 11, height = 8.5)


# Land
fit_no_na %>%
  spread_draws(tx_out_land[tx]) %>%
  median_qi(.width = 0.95) %>%
  left_join(index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  ggplot(aes(y = full_taxa_name, x = tx_out_land, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() +
  #stat_halfeye(aes(slab_fill = full_taxa_name)) +
  theme_classic() + 
  tx_plot_theme + 
  theme(legend.position = "none") +
  labs(x = bquote('m'^2~'a per tonne'), y = "", title = "Land")
ggsave(filename = file.path(outdir, paste("plot_feed-impact_land_", set_allocation, "-allocation_taxa-level-WEIGHTED.png", sep = "")), width = 11, height = 8.5)

# Water
fit_no_na %>%
  spread_draws(tx_out_water[tx]) %>%
  median_qi(.width = 0.95) %>%
  left_join(index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  ggplot(aes(y = full_taxa_name, x = tx_out_water, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() +
  #stat_halfeye(aes(slab_fill = full_taxa_name)) +
  theme_classic() + 
  tx_plot_theme + 
  theme(legend.position = "none") +
  labs(x = bquote('m'^3~'per tonne'), y = "", title = "Water")
ggsave(filename = file.path(outdir, paste("plot_feed-impact_water_", set_allocation, "-allocation_taxa-level-WEIGHTED.png", sep = "")), width = 11, height = 8.5)



# UNWEIGHTED OUTPUTS:
# Sci-level feed footprints as point intervals:
# Carbon
fit_no_na %>%
  spread_draws(sci_c_fp[sci]) %>%
  median_qi(.width = 0.95) %>%
  left_join(index_key, by = "sci") %>% # Join with index key to get sci and taxa names
  mutate(clean_sci_name = fct_reorder(clean_sci_name, as.character(full_taxa_name))) %>%
  ggplot(aes(y = clean_sci_name, x = sci_c_fp, xmin = .lower, xmax = .upper, color = full_taxa_name)) +
  geom_pointinterval() +
  theme_classic() + 
  sci_plot_theme + 
  labs(x = "kg CO2-eq per tonne", y = "", title = "Carbon", color = "taxa group")
ggsave(filename = file.path(outdir, paste("plot_feed-impact_carbon_", set_allocation, "-allocation_sci-level.png", sep = "")), width = 11, height = 8.5)

# Nitrogen
fit_no_na %>%
  spread_draws(sci_n_fp[sci]) %>%
  median_qi(.width = 0.95) %>%
  left_join(index_key, by = "sci") %>% # Join with index key to get sci and taxa names
  mutate(clean_sci_name = fct_reorder(clean_sci_name, as.character(full_taxa_name))) %>%
  ggplot(aes(y = clean_sci_name, x = sci_n_fp, xmin = .lower, xmax = .upper, color = full_taxa_name)) +
  geom_pointinterval() +
  theme_classic() + 
  sci_plot_theme + 
  labs(x = "kg N-eq per tonne", y = "", title = "Nitrogen", color = "taxa group")
ggsave(filename = file.path(outdir, paste("plot_feed-impact_nitrogen_", set_allocation, "-allocation_sci-level.png", sep = "")), width = 11, height = 8.5)

# Phosphorus
fit_no_na %>%
  spread_draws(sci_p_fp[sci]) %>%
  median_qi(.width = 0.95) %>%
  left_join(index_key, by = "sci") %>% # Join with index key to get sci and taxa names
  mutate(clean_sci_name = fct_reorder(clean_sci_name, as.character(full_taxa_name))) %>%
  ggplot(aes(y = clean_sci_name, x = sci_p_fp, xmin = .lower, xmax = .upper, color = full_taxa_name)) +
  geom_pointinterval() +
  theme_classic() + 
  sci_plot_theme + 
  labs(x = "kg P-eq per tonne", y = "", title = "Phosphorus", color = "taxa group")
ggsave(filename = file.path(outdir, paste("plot_feed-impact_phosphorus_", set_allocation, "-allocation_sci-level.png", sep = "")), width = 11, height = 8.5)

# Land
fit_no_na %>%
  spread_draws(sci_land_fp[sci]) %>%
  median_qi(.width = 0.95) %>%
  left_join(index_key, by = "sci") %>% # Join with index key to get sci and taxa names
  mutate(clean_sci_name = fct_reorder(clean_sci_name, as.character(full_taxa_name))) %>%
  ggplot(aes(y = clean_sci_name, x = sci_land_fp, xmin = .lower, xmax = .upper, color = full_taxa_name)) +
  geom_pointinterval() +
  theme_classic() + 
  sci_plot_theme + 
  labs(x = bquote('m'^2~'a per tonne'), y = "", title = "Land", color = "taxa group")
ggsave(filename = file.path(outdir, paste("plot_feed-impact_land_", set_allocation, "-allocation_sci-level.png", sep = "")), width = 11, height = 8.5)

# Water
fit_no_na %>%
  spread_draws(sci_water_fp[sci]) %>%
  median_qi(.width = 0.95) %>%
  left_join(index_key, by = "sci") %>% # Join with index key to get sci and taxa names
  mutate(clean_sci_name = fct_reorder(clean_sci_name, as.character(full_taxa_name))) %>%
  ggplot(aes(y = clean_sci_name, x = sci_water_fp, xmin = .lower, xmax = .upper, color = full_taxa_name)) +
  geom_pointinterval() +
  theme_classic() + 
  sci_plot_theme + 
  labs(x = bquote('m'^3~'per tonne'), y = "", title = "Water", color = "taxa group")
ggsave(filename = file.path(outdir, paste("plot_feed-impact_water_", set_allocation, "-allocation_sci-level.png", sep = "")), width = 11, height = 8.5)

# Taxa-level feed footprints as point intervals:
# Carbon
fit_no_na %>%
  spread_draws(tx_c_fp[tx]) %>%
  median_qi(.width = 0.95) %>%
  left_join(index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  ggplot(aes(y = full_taxa_name, x = tx_c_fp, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() +
  #stat_halfeye(aes(slab_fill = full_taxa_name)) +
  theme_classic() + 
  tx_plot_theme + 
  theme(legend.position = "none") +
  labs(x = "kg CO2-eq per tonne", y = "", title = "Carbon")
ggsave(filename = file.path(outdir, paste("plot_feed-impact_carbon_", set_allocation, "-allocation_taxa-level.png", sep = "")), width = 11, height = 8.5)

# Nitrogen
fit_no_na %>%
  spread_draws(tx_n_fp[tx]) %>%
  median_qi(.width = 0.95) %>%
  left_join(index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  ggplot(aes(y = full_taxa_name, x = tx_n_fp, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() +
  #stat_halfeye(aes(slab_fill = full_taxa_name)) +
  theme_classic() + 
  tx_plot_theme + 
  theme(legend.position = "none") +
  labs(x = "kg N-eq per tonne", y = "", title = "Nitrogen")
ggsave(filename = file.path(outdir, paste("plot_feed-impact_nitrogen_", set_allocation, "-allocation_taxa-level.png", sep = "")), width = 11, height = 8.5)


# Phosphorus
fit_no_na %>%
  spread_draws(tx_p_fp[tx]) %>%
  median_qi(.width = 0.95) %>%
  left_join(index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  ggplot(aes(y = full_taxa_name, x = tx_p_fp, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() +
  #stat_halfeye(aes(slab_fill = full_taxa_name)) +
  theme_classic() + 
  tx_plot_theme + 
  theme(legend.position = "none") +
  labs(x = "kg P-eq per tonne", y = "", title = "Phosphorus")
ggsave(filename = file.path(outdir, paste("plot_feed-impact_phosphorus_", set_allocation, "-allocation_taxa-level.png", sep = "")), width = 11, height = 8.5)


# Land
fit_no_na %>%
  spread_draws(tx_land_fp[tx]) %>%
  median_qi(.width = 0.95) %>%
  left_join(index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  ggplot(aes(y = full_taxa_name, x = tx_land_fp, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() +
  #stat_halfeye(aes(slab_fill = full_taxa_name)) +
  theme_classic() + 
  tx_plot_theme + 
  theme(legend.position = "none") +
  labs(x = bquote('m'^2~'a per tonne'), y = "", title = "Land")
ggsave(filename = file.path(outdir, paste("plot_feed-impact_land_", set_allocation, "-allocation_taxa-level.png", sep = "")), width = 11, height = 8.5)

# Water
fit_no_na %>%
  spread_draws(tx_water_fp[tx]) %>%
  median_qi(.width = 0.95) %>%
  left_join(index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  ggplot(aes(y = full_taxa_name, x = tx_water_fp, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() +
  #stat_halfeye(aes(slab_fill = full_taxa_name)) +
  theme_classic() + 
  tx_plot_theme + 
  theme(legend.position = "none") +
  labs(x = bquote('m'^3~'per tonne'), y = "", title = "Water")
ggsave(filename = file.path(outdir, paste("plot_feed-impact_water_", set_allocation, "-allocation_taxa-level.png", sep = "")), width = 11, height = 8.5)





# OPTION 2: Taxa-level plots with color themes:
# If we want to mimic bayesplot color schemes, can get hexadecimal colors and input manually to stat_halfeye aesthetics
color_scheme_get("blue")
color_scheme_get("green")





