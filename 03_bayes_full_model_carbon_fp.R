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

######################################################################################################
# RESTARTING POINT
#datadir <- "/Volumes/jgephart/BFA Environment 2/Data"
#outdir <- "/Volumes/jgephart/BFA Environment 2/Outputs"

# Load on-farm carbon (GHG) variables
# load(file.path(outdir, "2020-12-16_on-farm-carbon-all-data-prior-to-aggregation.RData"))

######################################################################################################
# STEP 3 - aggregate up to taxa level and estimate total feed footprint
# Set data for model

# Load off-farm (feed-associated) variables
load(file.path(outdir, "2020-12-20_off-farm-all-impacts-all-data-prior-to-aggregation.RData"))

# Organize off-farm feed data (FCR and feed proportions)
# FIX IT - change clean.lca in Functions.R so no less than 0.01
# Adjust feed proportions to be no less than 0.01
# TEMPORARY FIX:
# Normalize the FINAL feed proportion values to be greater than 0 and no less than 0.01
full_feed_dat <- full_feed_dat %>%
  mutate(feed_proportion = if_else(feed_proportion < 0.01, true = 0.0105, false = feed_proportion))

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

# Organize on-farm carbon data
# Merge electric, diesel, petrol, and natgas datasets
electric_dat_merge <- full_electric_dat %>%
  select(study_id, Country, iso3c, clean_sci_name, taxa, intensity, system, electric_data_type = data_type, electric)

diesel_dat_merge <- full_diesel_dat %>%
  select(study_id, Country, iso3c, clean_sci_name, taxa, intensity, system, diesel_data_type = data_type, diesel)

petrol_dat_merge <- full_petrol_dat %>%
  select(study_id, Country, iso3c, clean_sci_name, taxa, intensity, system, petrol_data_type = data_type, petrol)

natgas_dat_merge <- full_natgas_dat %>%
  select(study_id, Country, iso3c, clean_sci_name, taxa, intensity, system, natgas_data_type = data_type, natgas)

# Merge and format on-farm data
# Calculate n_in_sci and n_in_taxa in case we want to remove small sample sizes
ghg_footprint_merge <- electric_dat_merge %>%
  full_join(diesel_dat_merge, by = intersect(names(electric_dat_merge), names(diesel_dat_merge))) %>%
  full_join(petrol_dat_merge, by = intersect(names(.), names(petrol_dat_merge))) %>%
  full_join(natgas_dat_merge, by = intersect(names(.), names(natgas_dat_merge))) %>%
  drop_na() %>% # Make sure to drop na before creating sci and tx indices (otherwise some indices drop out)
  group_by(clean_sci_name) %>%
  mutate(n_in_sci = n()) %>%
  ungroup() %>%
  group_by(taxa) %>%
  mutate(n_in_taxa = n()) %>%
  ungroup() 

# MERGE on and off-farm data
carbon_footprint_dat_raw <- feed_dat_merge %>%
  full_join(fcr_dat_merge, by = intersect(names(feed_dat_merge), names(fcr_dat_merge))) %>%
  full_join(ghg_footprint_merge, by = intersect(names(ghg_footprint_merge), names(.))) %>%
  drop_na() %>%  # Make sure to drop na before creating sci and tx indices (otherwise some indices drop out)
  arrange(clean_sci_name, taxa) %>%
  mutate(clean_sci_name = as.factor(clean_sci_name),
         sci = as.numeric(clean_sci_name),
         taxa = as.factor(taxa),
         tx = as.numeric(taxa)) 

# OUTPUT data used for analysis:
write.csv(carbon_footprint_dat_raw, file.path(outdir, "dat_for_bayesian-carbon.csv"), row.names = FALSE)

# Multiply energy use by their GHG footprint - in STAN model generated quantities will be just the sum of each energy source's footprint
# (Electricity use * country-specific GHG of electricity) + (Diesel * GHG of diesel) + (Petrol * GHG of petrol) + (Natural gas * GHG of natural gas)

# Get farm-associated carbon footprints data
# Add iso3c
electric_fp_dat <- read.csv(file.path(datadir, "electricity_GWP.csv")) %>%
  mutate(iso3c = countrycode(Country, origin = "country.name", destination = "iso3c"))
other_energy_fp_dat <- read.csv(file.path(datadir, "energy_carriers_impact_factors.csv"))

# units for diesel and petrol constants are in kg CO2-eq / L
diesel_fp <- other_energy_fp_dat %>% filter(Impact.category == "Global warming potential" & Input == "Diesel") %>% pull(Value)
petrol_fp <- other_energy_fp_dat %>% filter(Impact.category == "Global warming potential" & Input == "Petrol") %>% pull(Value)

# NOTE: natural gas is in kg CO2-eq / m3, need to convert to kg CO2-eq / L - multiply by 0.001 m3 / L
natgas_fp <- other_energy_fp_dat %>% filter(Impact.category == "Global warming potential" & Input == "Natural gas") %>% pull(Value) * 0.001

# FIX IT - join with full lca_dat_clean dataset - and incorporate non-fed taxa (plants, bivalves) into model output (these would only have on-farm impacts)
# FIX IT - check studies with no country data, might be able to fix some of these manually
# Calculate GHG footprint for each energy source
carbon_footprint_dat <- carbon_footprint_dat_raw %>%
  # Calculate electriciy GHG footprint
  left_join(electric_fp_dat %>% select(-Country), by = "iso3c") %>% 
  mutate(GWP_perkWh_kgCO2eq = if_else(is.na(GWP_perkWh_kgCO2eq), true = mean(GWP_perkWh_kgCO2eq, na.rm = TRUE), false = GWP_perkWh_kgCO2eq)) %>% # for studies with no country data, just use the average across countries
  mutate(electric_ghg = electric * GWP_perkWh_kgCO2eq) %>%
  # Calculate diesel GHG footprint
  mutate(diesel_ghg = diesel * diesel_fp) %>%
  # Calculate petrol GHG footprint
  mutate(petrol_ghg = petrol * petrol_fp) %>%
  # Calculate natural gas GHG footprint
  mutate(natgas_ghg = natgas * natgas_fp) %>%
  # Calculate sum total of GHG footprint
  mutate(total_ghg = electric_ghg + diesel_ghg + petrol_ghg + natgas_ghg)

######################################################################################################
# END all cleaning: carbon_footprint_dat is the data used for Bayesian

# Set Data, indices, constants, weights:

# VARIABLE-SPECIFIC DATA:

# For on_farm ghg
farm_c <- carbon_footprint_dat$total_ghg
# For FCR model:
fcr <- carbon_footprint_dat$fcr
# For Feed proportion model:
K = 4
feed_weights <- carbon_footprint_dat %>%
  select(feed_soy, feed_crops, feed_fmfo, feed_animal) %>%
  as.matrix()
# Also needed for feed proportion dirichlet: counts per sci name and counts per taxa group (also included as data in the model):
sci_kappa <- carbon_footprint_dat %>% 
  group_by(sci) %>% 
  summarise(n_obs = n()) %>%
  ungroup() %>%
  arrange(sci) %>%
  pull(n_obs)
tx_kappa <- carbon_footprint_dat %>% 
  group_by(tx) %>% 
  summarise(n_obs = n()) %>%
  ungroup() %>%
  arrange(tx) %>%
  pull(n_obs)

# INDICES:
N = nrow(carbon_footprint_dat)
N_SCI <- length(unique(carbon_footprint_dat$sci))
n_to_sci <- carbon_footprint_dat$sci
n_to_tx <- carbon_footprint_dat$tx
N_TX <- length(unique(carbon_footprint_dat$tx))
sci_to_tx <- carbon_footprint_dat %>%
  select(sci, tx) %>%
  unique() %>%
  pull(tx)

# FEED IMPACT CONSTANTS:
# Choose allocation method
# Choose CARBON for Impact.category
# IMPORTANT - multiply all values by 1000 to convert to kg CO2 per tonne (currently in kg CO2 per kg)
impact <- "Global warming potential" # i.e., Carbon impact
set_allocation <- "Mass"

fp_dat <- read.csv(file.path(datadir, "20201217_weighted_feed_fp.csv")) %>%
  filter(Allocation == set_allocation) %>%
  mutate(ave_stressor_per_tonne = ave_stressor * 1000)

# ORDER OF feed_weights data vector: soy, crops, fmfo, animal
head(feed_weights)
# ORDER of footprint data vector should match this:
set_fp_order <- c("Soy", "Crop", "Fishery", "Livestock")

fp_constant <- fp_dat %>%
  filter(Impact.category == impact) %>% 
  arrange(match(Input.type, set_fp_order)) %>% # Match index and arrange by custom order
  select(ave_stressor_per_tonne) %>%
  as.matrix() %>%
  c()

# WEIGHTS:
# Get sci-level weightings for generating taxa-level quantities:
# IMPORTANT arrange by clean_sci_name so that the order matches data
sci_prod_weights <- read.csv(file.path(datadir, "aqua_prod_weightings.csv")) %>%
  arrange(clean_sci_name)
# Drop sci-names not found in the data
sci_prod_weights <- sci_prod_weights %>%
  filter(clean_sci_name %in% carbon_footprint_dat$clean_sci_name)
# Check that the order of sci names in both weights and data are the same (sum = 0)
sum(sci_prod_weights$clean_sci_name != unique(carbon_footprint_dat$clean_sci_name))
sci_w <- sci_prod_weights$weighting

# Need the following for generating ragged array in STAN - indexing vector of sci_mu's based on their taxa identity (each slice will have a different length)
where_tx <- order(sci_to_tx) # i.e., give the positions in sci_to_tx in order of their taxa level (where in sci_to_tx is taxa level 1, followed by taxa level 2, etc)
# How many sci are in each taxa - need this to declare length of each sci_mu vector
n_sci_in_tx <- carbon_footprint_dat %>%
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

# Set data for stan:
stan_data <- list(N = N,
                  N_SCI = N_SCI, 
                  n_to_sci = n_to_sci,
                  N_TX = N_TX,
                  sci_to_tx = sci_to_tx,
                  fcr = fcr,
                  K = K,
                  feed_weights = feed_weights,
                  farm_c = farm_c,
                  sci_kappa = sci_kappa,
                  tx_kappa = tx_kappa,
                  fp_constant = fp_constant,
                  sci_w = sci_w,
                  where_tx = where_tx,
                  n_sci_in_tx = n_sci_in_tx,
                  slice_where_tx = slice_where_tx)

# NORMAL distribution hierarchical model
stan_no_na <- 'data {
  // indices
  int<lower=0> N;  // number of observations
  int N_TX; // number of taxa groups
  int N_SCI; // number of scientific names
  int n_to_sci[N]; // sciname index
  int sci_to_tx[N_SCI]; // taxa-group indices

  // data for feed footrpint
  vector<lower=0>[N] fcr; // fcr data
  int K; // number of feed types
  simplex[K] feed_weights[N]; // array of observed feed weights simplexes
  int sci_kappa[N_SCI]; // number of observations per sci-name
  int tx_kappa[N_TX]; // number of observations per taxa group
  
  // data for on-farm footrpint
  vector<lower=0>[N] farm_c; // data

  // constants to apply to feed footrpint
  //vector[K] fp_constant;
  
  // data for slicing vectors for calculating weighted means
  vector<lower=0>[N_SCI] sci_w; // sci-level production weights
  int where_tx[N_SCI]; // order sci_to_tx by taxa-level
  int n_sci_in_tx[N_TX]; // number of sci in each taxa-level, ordered by taxa-level
  int slice_where_tx[N_TX + 1]; // breaks in where_tx by taxa-level
}
parameters {
  // FCR model
  vector<lower=0>[N_TX] tx_mu_fcr;
  vector<lower=0>[N_SCI] sci_mu_fcr;
  real<lower=0> tx_sigma_fcr;
  real<lower=0> sci_sigma_fcr;
  //vector<lower=0>[N_TX] tx_sigma_fcr;
  //vector<lower=0>[N_SCI] sci_sigma_fcr;
  
  // Feed proportion model:
  simplex[K] sci_theta[N_SCI]; // vectors of estimated sci-level feed weight simplexes
  simplex[K] tx_theta[N_TX];
  
  // On farm model
  vector<lower=0>[N_TX] tx_mu_farm;
  vector<lower=0>[N_SCI] sci_mu_farm;
  real<lower=0> tx_sigma_farm;
  real<lower=0> sci_sigma_farm;
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

  // weak priors on sigma
  tx_sigma_fcr ~ cauchy(0, 1);
  sci_sigma_fcr ~ cauchy(0, 1);
  tx_sigma_farm ~ cauchy(0, 10);
  sci_sigma_farm ~ cauchy(0, 10);

  // likelihood
  // normal model sci-name and taxa-level for FCR
  fcr ~ normal(sci_mu_fcr[n_to_sci], sci_sigma_fcr);
  sci_mu_fcr ~ normal(tx_mu_fcr[sci_to_tx], tx_sigma_fcr);

  // dirichlet model for feed proportions
  // use data to estimate sci-level dirichlet shape param (alphas)
  // LOOPING feed proportions:
  for (n in 1:N){
    feed_weights[n] ~ dirichlet(to_vector(sci_alpha[n_to_sci[n]]));
  }
  // use sci-level thetas to estimate taxa-level dirichlet shape param (alphas)
  for (n_sci in 1:N_SCI){
    sci_theta[n_sci] ~ dirichlet(to_vector(tx_alpha[sci_to_tx[n_sci]]));
  }
  
  // normal model for sci and taxa-level on-farm footrpint
  farm_c ~ normal(sci_mu_farm[n_to_sci], sci_sigma_farm);
  sci_mu_farm ~ normal(tx_mu_farm[sci_to_tx], tx_sigma_farm);
}
generated quantities {
  // Declare vectors for weightings
  vector[N_TX] tx_feed_fp;
  vector[N_SCI] sci_feed_fp;
  vector[N_SCI] sci_total_fp; // unweighted sci-level
  vector[N_TX] tx_total_fp; //
  vector[N_SCI] sci_total_fp_w; // weighted sci-level
  vector[N_TX] tx_total_fp_w; // weighted tx-level

  // Feed footrpint calculations
  for (n_tx in 1:N_TX) {
    tx_feed_fp[n_tx] = tx_mu_feed[n_tx] * sum(fp_constant .* tx_theta[n_tx]);
  }
  for (n_sci in 1:N_SCI) {
    sci_feed_fp[n_sci] = sci_mu_feed[n_sci] * sum(fp_constant .* sci_theta[n_sci]);
  }

  // Sum off farm (feed) and on farm footprints (UNWEIGHTED sci and tx-level outputs)
  for (n_sci in 1:N_SCI) {
    sci_total_fp[n_sci] = sci_feed_fp[n_sci] + sci_mu_farm[n_sci];
  }
  for (n_tx in 1:N_TX) {
    tx_total_fp[n_tx] = tx_feed_fp[n_tx] + tx_mu_farm[n_tx];
  }

  // Apply weightings
  //sci_total_fp_w = sci_total_fp .* sci_w; // WEIGHTED sci-level outputs

  //for (n_tx in 1:N_TX){
  //  vector[n_sci_in_tx[n_tx]] sci_mu_w_vec; // declare vector of sci_mu in taxa-level n_tx
  //  sci_mu_w_vec = sci_total_fp_w[where_tx[slice_where_tx[n_tx]+1:slice_where_tx[n_tx+1]]]; // get all the sci_mu in taxa-level n_tx
  //  tx_total_fp_w[n_tx] = sum(sci_mu_w_vec); // sum sci_mu_w_vec to get WEIGHTED tx-level outputs
  //}
}'

no_na_mod <- stan_model(model_code = stan_no_na)
# Note: For Windows, apparently OK to ignore this warning message:
# Warning message:
#   In system(paste(CXX, ARGS), ignore.stdout = TRUE, ignore.stderr = TRUE) :
#   'C:/rtools40/usr/mingw_/bin/g++' not found

# Fit model:
# Set seed while testing
fit_no_na <- sampling(object = no_na_mod, 
                      data = stan_data, 
                      cores = 4, 
                      seed = "11729", 
                      iter = 2500, 
                      control = list(adapt_delta = 0.99, max_treedepth = 15))
#fit_no_na <- sampling(object = no_na_mod, data = stan_data, cores = 4, iter = 5000, control = list(adapt_delta = 0.99))
summary(fit_no_na)$summary

#launch_shinystan(fit_no_na)
######################################################################################################
# RESTARTING POINT
# FIX IT - which objects to clear before saving?
# rm(list=ls()[!(ls() %in% c("datadir", "outdir", 
#                            "lca_dat_clean_groups", "feed_model_dat_categories",
#                            "full_feed_dat", "full_fcr_dat", "feed_footprint_dat", "fit_no_na"))])
#save.image(file.path(outdir, paste(Sys.Date(), "_full-model_", impact, "_", set_allocation, "-allocation_all-data-prior-to-plotting.RData", sep = "")))

# datadir <- "/Volumes/jgephart/BFA Environment 2/Data"
# outdir <- "/Volumes/jgephart/BFA Environment 2/Outputs"
# load(file.path(outdir, "<file-name>.RData"))
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
index_key <- carbon_footprint_dat %>%
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

# Key for naming [sci, feed] and [taxa, feed] levels
# Order of feeds: soy, crops, fmfo, animal
tx_feed_key <- carbon_footprint_dat %>%
  select(contains(c("taxa", "tx", "soy", "crops", "fmfo", "animal"))) %>%
  pivot_longer(cols = contains(c("soy", "crops", "fmfo", "animal")), names_to = "feed") %>%
  select(-value) %>%
  unique() %>%
  mutate(feed_index = case_when(str_detect(feed, "soy") ~ 1,
                                str_detect(feed, "crops") ~ 2,
                                str_detect(feed, "fmfo") ~ 3,
                                str_detect(feed, "animal") ~ 4)) %>%
  # Clean feed names
  mutate(feed = gsub(feed, pattern = "feed_", replacement = "")) %>%
  mutate(index = paste("[", tx, ",", feed_index, "]", sep = "")) %>%
  mutate(tx_theta_name = paste("theta[", taxa, ", ", feed, "]", sep = "")) %>%
  mutate(tx_alpha_name = paste("alpha[", taxa, ", ", feed, "]", sep = "")) %>%
  # IMPORTANT before replaceing param names: ARRANGE BY FEED, THEN TAXA NAME TO MATCH HOW NAMES ARE ARRANGED IN STANFIT OBJECT
  arrange(feed_index, tx)

sci_feed_key <- carbon_footprint_dat %>%
  select(contains(c("clean_sci_name", "taxa", "sci", "soy", "crops", "fmfo", "animal"))) %>%
  pivot_longer(cols = contains(c("soy", "crops", "fmfo", "animal")), names_to = "feed") %>%
  select(-value) %>%
  unique() %>%
  mutate(feed_index = case_when(str_detect(feed, "soy") ~ 1,
                                str_detect(feed, "crops") ~ 2,
                                str_detect(feed, "fmfo") ~ 3,
                                str_detect(feed, "animal") ~ 4)) %>%
  # Clean feed names
  mutate(feed = gsub(feed, pattern = "feed_", replacement = "")) %>%
  mutate(index = paste("[", sci, ",", feed_index, "]", sep = "")) %>%
  mutate(sci_theta_name = paste("theta[", clean_sci_name, ", ", feed, "]", sep = "")) %>%
  mutate(sci_alpha_name = paste("alpha[", clean_sci_name, ", ", feed, "]", sep = "")) %>%
  # IMPORTANT before replaceing param names: ARRANGE BY FEED, THEN TAXA NAME TO MATCH HOW NAMES ARE ARRANGED IN STANFIT OBJECT
  arrange(feed_index, sci)


# Use tidybayes + ggdist for finer control of aes mapping (instead of bayesplots) 
get_variables(fit_no_na)

######################################################################################################
# PLOT final outputs (total on + off farm impacts)

######################################################################################################
# PLOT other intermediate-level calculations
######################################################################################################
# WEIGHTED OUTPUTS:
# Sci-level feed footprints as point intervals:
# Carbon
fit_no_na %>%
  spread_draws(sci_total_fp_w[sci]) %>%
  median_qi(.width = 0.95) %>%
  left_join(index_key, by = "sci") %>% # Join with index key to get sci and taxa names
  mutate(clean_sci_name = fct_reorder(clean_sci_name, as.character(full_taxa_name))) %>%
  ggplot(aes(y = clean_sci_name, x = sci_total_fp_w, xmin = .lower, xmax = .upper, color = full_taxa_name)) +
  geom_pointinterval() +
  theme_classic() + 
  sci_plot_theme + 
  labs(x = "kg CO2-eq per tonne", y = "", title = "Carbon", color = "taxa group")
ggsave(filename = file.path(outdir, paste("plot_feed-impact_carbon_", set_allocation, "-allocation_sci-level-WEIGHTED.png", sep = "")), width = 11, height = 8.5)

# Taxa-level feed footprints as point intervals:
# Carbon
fit_no_na %>%
  spread_draws(tx_total_fp_w[tx]) %>%
  median_qi(.width = 0.95) %>%
  left_join(index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  ggplot(aes(y = full_taxa_name, x = tx_total_fp_w, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() +
  #stat_halfeye(aes(slab_fill = full_taxa_name)) +
  theme_classic() + 
  tx_plot_theme + 
  theme(legend.position = "none") +
  labs(x = "kg CO2-eq per tonne", y = "", title = "Carbon")
ggsave(filename = file.path(outdir, paste("plot_feed-impact_carbon_", set_allocation, "-allocation_taxa-level-WEIGHTED.png", sep = "")), width = 11, height = 8.5)


# UNWEIGHTED OUTPUTS:
# Sci-level feed footprints as point intervals:
# Carbon
fit_no_na %>%
  spread_draws(sci_feed_fp[sci]) %>%
  median_qi(.width = 0.95) %>%
  left_join(index_key, by = "sci") %>% # Join with index key to get sci and taxa names
  mutate(clean_sci_name = fct_reorder(clean_sci_name, as.character(full_taxa_name))) %>%
  ggplot(aes(y = clean_sci_name, x = sci_feed_fp, xmin = .lower, xmax = .upper, color = full_taxa_name)) +
  geom_pointinterval() +
  theme_classic() + 
  sci_plot_theme + 
  labs(x = "kg CO2-eq per tonne", y = "", title = "Carbon", color = "taxa group")
ggsave(filename = file.path(outdir, paste("plot_feed-impact_carbon_", set_allocation, "-allocation_sci-level.png", sep = "")), width = 11, height = 8.5)

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

# FIX IT - fct_reorder not working for feed proportions?
# Sci-level theta (feed proportions)
# Make separate plot for each feed component
feed_component <- c("soy", "crops", "fmfo", "animal")
for (i in 1:length(feed_component)){
  plot_dat <- fit_no_na %>%
    spread_draws(sci_theta[sci, feed_index]) %>%
    median_qi(.width = 0.95) %>%
    left_join(sci_feed_key, by = c("sci" = "sci", "feed_index" = "feed_index")) %>% # Join with index key to get sci and taxa names
    filter(feed == feed_component[i]) %>%
    mutate(clean_sci_name = fct_reorder(clean_sci_name, as.character(taxa)))
    #mutate(clean_sci_name = fct_reorder(clean_sci_name, sci_theta))
    p <- ggplot(plot_dat, aes(y = clean_sci_name, x = sci_theta, xmin = .lower, xmax = .upper, color = taxa)) +
      geom_pointinterval() +
      #stat_halfeye(aes(slab_fill = full_taxa_name)) +
      coord_cartesian(xlim = c(0, 1)) +
      theme_classic() + 
      sci_plot_theme + 
      labs(x = "Feed Proportion", y = "", title = feed_component[i], color = "taxa group")
    print(p)
    ggsave(filename = file.path(outdir, paste("plot_feed-proportion_", feed_component[i], "_sci-level.png", sep = "")), width = 11, height = 8.5)
}

# Taxa-level theta (feed proportions)
# Make separate plot for each feed component
feed_component <- c("soy", "crops", "fmfo", "animal")
for (i in 1:length(feed_component)){
  plot_dat <- fit_no_na %>%
    spread_draws(tx_theta[tx, feed_index]) %>%
    median_qi(.width = 0.95) %>%
    left_join(tx_feed_key, by = c("tx" = "tx", "feed_index" = "feed_index")) %>% # Join with index key to get sci and taxa names
    filter(feed == feed_component[i])
  p <- ggplot(plot_dat, aes(y = taxa, x = tx_theta, xmin = .lower, xmax = .upper)) +
    geom_pointinterval() +
    #stat_halfeye(aes(slab_fill = full_taxa_name)) +
    coord_cartesian(xlim = c(0, 1)) +
    theme_classic() + 
    tx_plot_theme + 
    theme(legend.position = "none") +
    labs(x = "Feed Proportion", y = "", title = feed_component[i])
  print(p)
  ggsave(filename = file.path(outdir, paste("plot_feed-proportion_", feed_component[i], "_taxa-level.png", sep = "")), width = 11, height = 8.5)
}

######################################################################################################

# Taxa-level FCR
fit_no_na %>%
  spread_draws(tx_mu_fcr[tx]) %>%
  median_qi(.width = 0.95) %>%
  left_join(index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  ggplot(aes(y = full_taxa_name, x = tx_mu_fcr, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() +
  #stat_halfeye(aes(slab_fill = full_taxa_name)) +
  theme_classic() + 
  tx_plot_theme + 
  theme(legend.position = "none") +
  labs(x = "", y = "", title = "FCR")
ggsave(filename = file.path(outdir, paste("plot_fcr_taxa-level.png", sep = "")), width = 11, height = 8.5)




# OPTION 2: Taxa-level plots with color themes:
# If we want to mimic bayesplot color schemes, can get hexadecimal colors and input manually to stat_halfeye aesthetics
color_scheme_get("blue")
color_scheme_get("green")





