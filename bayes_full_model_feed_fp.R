# Estiamte on-farm, feed-associated footprint
# Reminder for new data: entries with data for some columns should have no blanks in other columns (these should be filled in as zeroes)

# First run: process_data_for_analysis.R, then clear environment other than:
rm(list=ls()[!(ls() %in% c("lca_dat_clean_groups", "datadir", "outdir"))])

# Get model-specific data:
# Remove FCR == 0 (species that aren't fed)
# SELECT STUDY ID COLUMN - use this for rejoining outputs from multiple regression models back together
# Select relevant data columns and arrange by categorical info
feed_model_dat_categories <- lca_dat_clean_groups %>%
  select(study_id, fcr = FCR, contains("new"), clean_sci_name, taxa, intensity = Intensity, system = Production_system_group) %>%
  arrange(clean_sci_name, taxa, intensity, system) %>%
  filter(taxa %in% c("bivalves", "plants")==FALSE) # Remove taxa that don't belong in FCR/feed analysis - bivalves

fcr_no_na <- feed_model_dat_categories %>%
  #filter(fcr != 0 | is.na(fcr))  %>% # If we want to retain NAs, have to explicitly include is.na(fcr) otherwise NA's get dropped by fcr != 0
  filter(fcr != 0) %>% # This also automatically drops NAs
  filter(is.na(intensity)==FALSE & is.na(system)==FALSE) %>% # complete predictors - i.e., both intensity AND system are non-NA
  select(study_id, fcr, clean_sci_name, taxa, intensity, system)

######################################################################################################
# Step 1: Model FCR as taxa + intensity + system
# ie, remove all NAs and model these data, then use the model to predict FCR for those data with a complete set of predictors

# Create model matrix for taxa info, then center and scale
options(na.action='na.pass') # First change default options for handling missing data
X_taxa <- model.matrix(object = ~ 1 + taxa, 
                       data = fcr_no_na %>% select(taxa)) 
options(na.action='na.omit') # Return option back to the default

taxa_sd <- apply(X_taxa[,-1], MARGIN=2, FUN=sd, na.rm=TRUE) # Center all non-intercept variables and scale by 2 standard deviations (ignoring NAs)
options(na.action='na.pass') # First change default options for handling missing data
X_taxa_scaled <- scale(X_taxa[,-1], center=TRUE, scale=2*taxa_sd)
options(na.action='na.omit') # Return option back to the default

# Format intensity and system as ordinal variables, then center and scale
X_ordinal <- fcr_no_na %>%
  mutate(intensity = factor(intensity, levels = c("Extensive", "Imp. extensive", "Semi-intensive", "Intensive"))) %>% # set order of factors (low = extensive, high = intensive)
  mutate(system = factor(system, levels = c("On- and off-bottom", "Cages & pens", "Ponds", "Recirculating and tanks"))) %>% # set order of factors (low = open, high = closed)
  mutate(intensity = as.numeric(intensity)) %>%
  mutate(system = as.numeric(system)) %>%
  select(intensity, system) %>%
  as.matrix()
ordinal_sd<-apply(X_ordinal, MARGIN=2, FUN=sd, na.rm=TRUE) # Center all non-intercept variables and scale by 2 standard deviations (ignoring NAs)

options(na.action='na.pass') # First change default options for handling missing data
X_ordinal_scaled <- scale(X_ordinal, center=TRUE, scale=2*ordinal_sd)
options(na.action='na.omit') # Return option back to the default

# Create dataframe for brms and rename feed variables
fcr_brms_data <- data.frame(y = fcr_no_na$fcr, X_taxa_scaled, X_ordinal_scaled)

# Set model formula
fcr_brms <- brmsformula(y ~ 1 + ., family = Gamma("log"))

# Use "resp = <response_variable>" to specify different priors for different response variables
all_priors <- c(set_prior("normal(0,5)", class = "b"), # priors for y response variables
                set_prior("normal(0,2.5)", class = "Intercept"), 
                set_prior("exponential(1)", class = "shape"))

# Model converges after increasing the adapt_delta and iterations from default values
# Rule of thumb: bulk and tail effective sample sizes should be 100 x number of chains (i.e., at least 400)
# increasing max_treedepth is more about efficiency (instead of validity)
# See: https://mc-stan.org/misc/warnings.html
fit_fcr_no_na <- brm(fcr_brms, data = fcr_brms_data,
                       prior = all_priors, cores = 4, seed = "11729", iter = 5000, control = list(adapt_delta = 0.99))

# Get stan code
#stancode(fit_fcr_no_na)

######################################################################################################
# Use fcr model to predict NA fcrs for studies with complete set of predictors
# Both intensity AND system are non-NA
fcr_complete_predictors <- feed_model_dat_categories %>%
  filter(fcr != 0 | is.na(fcr))  %>% # Have to explicitly include is.na(fcr) otherwise NA's get dropped by fcr != 0
  filter(is.na(intensity)==FALSE & is.na(system)==FALSE) %>%
  filter(is.na(fcr))

# PROBLEM: lca_complete predictors has more taxa than originally model:
taxa_not_modeled <- setdiff(unique(fcr_complete_predictors$taxa), unique(fcr_no_na$taxa)) # these taxa were never modeled so they can't be predicted below

# DROP THESE FOR NOW:
fcr_complete_predictors <- fcr_complete_predictors %>%
  filter(taxa %in% taxa_not_modeled == FALSE)

# Now check the other way, which taxa were in the original model but not a part of the data that needs to be predicted:
setdiff(unique(fcr_no_na$taxa), unique(fcr_complete_predictors$taxa))

# If original model has taxa that are not part of fcr_complete_predictors, need to convert to factor and expand/assign levels manually - having trouble automating this
# See list of unique taxa in fcr_no_na - remember the first level is part of the "contrasts" in design matrix
sort(unique(fcr_no_na$taxa))
fcr_complete_predictors <- fcr_complete_predictors %>%
  mutate(taxa = as.factor(taxa))
levels(fcr_complete_predictors$taxa) <- list(catfish = "catfish",
                                             hypoph_carp = "hypoph_carp",
                                             milkfish = "milkfish",
                                             misc_diad = "misc_diad", 
                                             misc_marine = "misc_marine",
                                             oth_carp = "oth_carp",
                                             salmon = "salmon", 
                                             shrimp = "shrimp", 
                                             tilapia = "tilapia", 
                                             trout = "trout")

# Create NEW taxa model matrix for the studies whose feeds need to be predicted
# Taxa categories:
options(na.action='na.pass') # First change default options for handling missing data
X_taxa_new <- model.matrix(object = ~ 1 + taxa, 
                           data = fcr_complete_predictors %>% select(taxa)) 
options(na.action='na.omit') # Return option back to the default

# Center and Scale: BUT now center by the mean of the original modeled dataset above AND scale by the same 2*SD calculated from the original, modeled dataset above
options(na.action='na.pass') # First change default options for handling missing data
X_taxa_new_scaled <- scale(X_taxa_new[,-1], center=apply(X_taxa[,-1], MARGIN = 2, FUN = mean), scale=2*taxa_sd)
options(na.action='na.omit') # Return option back to the default

# System and Intensity variables:
# Format intensity and system as ordinal variables, then center and scale
X_ordinal_new <- fcr_complete_predictors %>%
  mutate(intensity = factor(intensity, levels = c("Extensive", "Imp. extensive", "Semi-intensive", "Intensive"))) %>% # set order of factors (low = extensive, high = intensive)
  mutate(system = factor(system, levels = c("On- and off-bottom", "Cages & pens", "Ponds", "Recirculating and tanks"))) %>% # set order of factors (low = open, high = closed)
  mutate(intensity = as.numeric(intensity)) %>%
  mutate(system = as.numeric(system)) %>%
  select(intensity, system) %>%
  as.matrix()

# Center and Scale: BUT now center by the mean of the original modeled dataset above AND scale by the same 2*SD calculated from the original, modeled dataset above
options(na.action='na.pass') # First change default options for handling missing data
X_ordinal_new_scaled <- scale(X_ordinal_new, center=apply(X_ordinal, MARGIN = 2, FUN = mean), scale=2*ordinal_sd)
options(na.action='na.omit') # Return option back to the default

# Create dataframe for brms and rename feed variables
brms_new_fcr_data <- data.frame(cbind(X_taxa_new_scaled, X_ordinal_new_scaled)) 

# Make predictions
#predicted_fcr_dat <- predict(fit_no_na, newdata = brms_new_fcr_data)
# Use tidybayes instead:
predicted_fcr_dat <- add_predicted_draws(newdata = brms_new_fcr_data, model = fit_fcr_no_na)

# Get point and interval estimates from predicted data
# Select just the prediction columns
# Join these with the modeled data (fcr_complete_predictors) to get metadata on taxa/intensity/syste,
fcr_dat_intervals <- predicted_fcr_dat %>%
  median_qi(.value = .prediction) %>% # Rename prediction to value
  ungroup() %>%
  select(contains("."))

# .row is equivalent to the row number in the modeled dataset (fcr_complete_predictors) - create a join column for this
fcr_metadat<- fcr_complete_predictors %>%
  select(study_id, clean_sci_name, taxa, intensity, system) %>%
  mutate(.row = row_number())

fcr_predictions <- fcr_dat_intervals %>%
  left_join(fcr_metadat, by = ".row") %>%
  rename(fcr = .value)

######################################################################################################
# Model fcr with intensity + system (no taxa)
# FOLLOWING SECTION DOES NOT APPLY

# Create dataframe for brms
# fcr_brms_data_no_taxa <- data.frame(y = fcr_no_na$fcr, X_ordinal_scaled)
# 
# names(fcr_brms_data_no_taxa)
# 
# # Set model formula
# fcr_brms_no_taxa <- brmsformula(y ~ 1 + ., family = Gamma("log"))
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
# fit_fcr_no_taxa <- brm(fcr_brms_no_taxa, data = fcr_brms_data_no_taxa,
#                          prior = all_priors, cores = 4, seed = "11729", iter = 5000, control = list(adapt_delta = 0.99))
# 
# # Get stan code
# #stancode(fit_fcr_no_taxa)
# 
# ######################################################################################################
# # Use model of intensity + system to predict fcr for taxa that were not part of the previous model
# # Intensity must be non-NA
# fcr_complete_predictors <- feed_model_dat_categories %>%
#   filter(fcr != 0 | is.na(fcr))  %>% # First drop the zeroes; Have to explicitly include is.na(fcr) otherwise NA's get dropped by fcr != 0
#   filter(is.na(intensity)==FALSE & is.na(system)==FALSE) %>%
#   filter(is.na(fcr)) # Now, filter to just the NAs
# 
# # PROBLEM: lca_complete predictors has more taxa than originally model:
# taxa_not_modeled <- setdiff(unique(fcr_complete_predictors$taxa), unique(fcr_no_na$taxa)) # these taxa were never modeled so they can't be predicted below
# 
# # Only keep taxa that were not modeled previously
# fcr_complete_predictors_no_taxa <- fcr_complete_predictors %>%
#   filter(taxa %in% taxa_not_modeled) 
# 
# # Format intensity as ordinal variable, then center and scale
# X_ordinal_no_taxa <- fcr_complete_predictors_no_taxa %>%
#   mutate(intensity = factor(intensity, levels = c("Extensive", "Imp. extensive", "Semi-intensive", "Intensive"))) %>% # set order of factors (low = extensive, high = intensive)
#   mutate(intensity = as.numeric(intensity))  %>%
#   mutate(system = factor(system, levels = c("On- and off-bottom", "Cages & pens", "Ponds", "Recirculating and tanks"))) %>% # set order of factors (low = open, high = closed)
#   mutate(system = as.numeric(system)) %>%
#   select(intensity, system) %>%
#   as.matrix()
# 
# # Center and Scale: BUT now center by the mean of the original modeled dataset above AND scale by the same 2*SD calculated from the original, modeled dataset above
# X_ordinal_scaled_no_taxa <- scale(X_ordinal_no_taxa, center=apply(X_ordinal, MARGIN = 2, FUN = mean), scale=2*ordinal_sd)
# 
# # Create dataframe for brms
# brms_fcr_dat_no_taxa <- data.frame(X_ordinal_scaled_no_taxa)
# 
# # Make predictions
# #predicted_fcr_dat <- predict(fit_no_na, newdata = brms_new_fcr_data)
# # Use tidybayes instead:
# predicted_fcr_dat_no_taxa <- add_predicted_draws(newdata = brms_fcr_dat_no_taxa, model = fit_fcr_no_taxa)
# 
# # Get point and interval estimates from predicted data
# # Select just the prediction columns
# # Join these with the modeled data (fcr_complete_predictors_no_taxa) to get metadata on taxa/intensity/syste,
# fcr_dat_intervals_no_taxa <- predicted_fcr_dat_no_taxa %>%
#   median_qi(.value = .prediction) %>% # Rename prediction to value
#   ungroup() %>%
#   select(contains("."))
# 
# # .row is equivalent to the row number in the modeled dataset (fcr_complete_predictors_no_taxa) - create a join column for this
# fcr_metadat_no_taxa <- fcr_complete_predictors_no_taxa %>%
#   select(study_id, clean_sci_name, taxa, intensity, system) %>%
#   mutate(.row = row_number())
# 
# fcr_predictions_no_taxa <- fcr_dat_intervals_no_taxa %>%
#   left_join(fcr_metadat_no_taxa, by = ".row") %>%
#   rename(fcr = .value)

######################################################################################################
# Bind fcr_no_na (data), fcr_predictions, and fcr_predictions_no_taxa
full_fcr_dat <- fcr_predictions %>%
  bind_rows(fcr_no_na) %>%
  #bind_rows(fcr_predictions_no_taxa) %>%
  mutate(data_type = if_else(is.na(.point), true = "data", false = "prediction")) %>%
  arrange(taxa, intensity, system, clean_sci_name) %>%
  rownames_to_column() # Arrange by taxa first, then create dummy column for plotting 


# Quick check: PLOT DATA + PREDICTIONS
ggplot(full_fcr_dat, aes(x = fcr, y = taxa)) + geom_boxplot(outlier.shape = NA) + geom_jitter(aes(color = data_type))

# NEXT: use plot_brms_outputs_gamma_regression to produce figures for SI

######################################################################################################
# Step 2: Model feed proportions
# Brms dirichlet regression doesn't support missing data, so strategy here is to remove all NAs and model these data
# Then use the model to predict missing feed proportion data for those data that have a complete set of predictors

# Remove fcr == 0 (species that aren't fed)
feed_no_na <- feed_model_dat_categories %>%
  filter(fcr != 0 | is.na(fcr))  %>% # Have to explicitly include is.na(fcr) otherwise NA's get dropped by fcr != 0
  filter(is.na(intensity)==FALSE & is.na(system)==FALSE) %>% # complete predictors - i.e., both intensity AND system are non-NA
  filter(is.na(feed_soy_new)==FALSE) %>%
  rename(feed_soy = feed_soy_new,
         feed_crops = feed_crops_new,
         feed_fmfo = feed_fmfo_new,
         feed_animal = feed_animal_new) %>%
  select(study_id, clean_sci_name, taxa, intensity, system, contains("feed"))

# Set data for model:

# Create model matrix for taxa info, then center and scale
options(na.action='na.pass') # First change default options for handling missing data
X_taxa <- model.matrix(object = ~ 1 + taxa, 
                       data = feed_no_na %>% select(taxa)) 
options(na.action='na.omit') # Return option back to the default

taxa_sd <- apply(X_taxa[,-1], MARGIN=2, FUN=sd, na.rm=TRUE) # Center all non-intercept variables and scale by 2 standard deviations (ignoring NAs)
options(na.action='na.pass') # First change default options for handling missing data
X_taxa_scaled <- scale(X_taxa[,-1], center=TRUE, scale=2*taxa_sd)
options(na.action='na.omit') # Return option back to the default

# Format intensity and system as ordinal variables, then center and scale
X_ordinal <- feed_no_na %>%
  mutate(intensity = factor(intensity, levels = c("Extensive", "Imp. extensive", "Semi-intensive", "Intensive"))) %>% # set order of factors (low = extensive, high = intensive)
  mutate(system = factor(system, levels = c("On- and off-bottom", "Cages & pens", "Ponds", "Recirculating and tanks"))) %>% # set order of factors (low = open, high = closed)
  mutate(intensity = as.numeric(intensity)) %>%
  mutate(system = as.numeric(system)) %>%
  select(intensity, system) %>%
  as.matrix()
ordinal_sd<-apply(X_ordinal, MARGIN=2, FUN=sd, na.rm=TRUE) # Center all non-intercept variables and scale by 2 standard deviations (ignoring NAs)

options(na.action='na.pass') # First change default options for handling missing data
X_ordinal_scaled <- scale(X_ordinal, center=TRUE, scale=2*ordinal_sd)
options(na.action='na.omit') # Return option back to the default

# Create dataframe for brms and rename feed variables
feed_brms_data <- data.frame(cbind(feed_no_na %>% select(contains("feed")), X_taxa_scaled, X_ordinal_scaled)) 

# Response variable must be a matrix, create function bind since cbind within the brm function is reserved for specifying multivariate models
bind <- function(...) cbind(...)

feed_brms <- brmsformula(bind(feed_soy, feed_crops, feed_fmfo, feed_animal) ~ 1 + ., family = dirichlet())

# Model converges after increasing the adapt_delta and iterations from default values
fit_feed_no_na <- brm(feed_brms, data = feed_brms_data,
                 cores = 4, seed = "11729", iter = 5000)

summary(fit_feed_no_na) 

######################################################################################################
# Use feed proportion model to predict NA feeds for studies with complete set of predictors

# Both intensity AND system are non-NA
feed_complete_predictors <- feed_model_dat_categories %>%
  filter(fcr != 0 | is.na(fcr))  %>% # Have to explicitly include is.na(fcr) otherwise NA's get dropped by fcr != 0
  filter(is.na(intensity)==FALSE & is.na(system)== FALSE) %>%
  filter(is.na(feed_soy_new))

# PROBLEM: lca_complete predictors has more taxa than originally model:
taxa_not_modeled <- setdiff(unique(feed_complete_predictors$taxa), unique(feed_no_na$taxa)) # these taxa were never modeled so they can't be predicted below

# DROP THESE FOR NOW:
feed_complete_predictors <- feed_complete_predictors %>%
  filter(taxa %in% taxa_not_modeled == FALSE)

# Now check the other way, which taxa were in the original model but not a part of the data that needs to be predicted:
setdiff(unique(feed_no_na$taxa), unique(feed_complete_predictors$taxa))

# If original model has taxa that are not part of feed_complete_predictors, need to convert to factor and expand/assign levels manually - having trouble automating this
# See list of unique taxa in feed_no_na - remember the first level is part of the "contrasts" in design matrix
sort(unique(feed_no_na$taxa))
feed_complete_predictors <- feed_complete_predictors %>%
  mutate(taxa = as.factor(taxa))
levels(feed_complete_predictors$taxa) <- list(catfish = "catfish",
                                             hypoph_carp = "hypoph_carp",
                                             milkfish = "milkfish",
                                             misc_diad = "misc_diad", 
                                             misc_marine = "misc_marine",
                                             oth_carp = "oth_carp",
                                             salmon = "salmon", 
                                             shrimp = "shrimp", 
                                             tilapia = "tilapia", 
                                             trout = "trout")

# Create NEW taxa model matrix for the studies whose feeds need to be predicted
# Taxa categories:
options(na.action='na.pass') # First change default options for handling missing data
X_taxa_new <- model.matrix(object = ~ 1 + taxa, 
                           data = feed_complete_predictors %>% select(taxa)) 
options(na.action='na.omit') # Return option back to the default

# Center and Scale: BUT now center by the mean of the original modeled dataset above AND scale by the same 2*SD calculated from the original, modeled dataset above
options(na.action='na.pass') # First change default options for handling missing data
X_taxa_new_scaled <- scale(X_taxa_new[,-1], center=apply(X_taxa[,-1], MARGIN = 2, FUN = mean), scale=2*taxa_sd)
options(na.action='na.omit') # Return option back to the default

# System and Intensity variables:
# Format intensity and system as ordinal variables, then center and scale
X_ordinal_new <- feed_complete_predictors %>%
  mutate(intensity = factor(intensity, levels = c("Extensive", "Imp. extensive", "Semi-intensive", "Intensive"))) %>% # set order of factors (low = extensive, high = intensive)
  mutate(system = factor(system, levels = c("On- and off-bottom", "Cages & pens", "Ponds", "Recirculating and tanks"))) %>% # set order of factors (low = open, high = closed)
  mutate(intensity = as.numeric(intensity)) %>%
  mutate(system = as.numeric(system)) %>%
  select(intensity, system) %>%
  as.matrix()

# Center and Scale: BUT now center by the mean of the original modeled dataset above AND scale by the same 2*SD calculated from the original, modeled dataset above
options(na.action='na.pass') # First change default options for handling missing data
X_ordinal_new_scaled <- scale(X_ordinal_new, center=apply(X_ordinal, MARGIN = 2, FUN = mean), scale=2*ordinal_sd)
options(na.action='na.omit') # Return option back to the default

# Create dataframe for brms and rename feed variables
brms_new_feed_data <- data.frame(cbind(X_taxa_new_scaled, X_ordinal_new_scaled)) 

# Make predictions
#predicted_feed_dat <- predict(fit_no_na, newdata = brms_new_feed_data)
# Use tidybayes instead:
predicted_feed_dat <- add_predicted_draws(newdata = brms_new_feed_data, model = fit_feed_no_na)

# Get point and interval estimates from predicted data
# Select just the prediction columns
# Join these with the original lca data (feed_complete_predictors) to get metadata on taxa/intensity/syste,
feed_dat_intervals <- predicted_feed_dat %>%
  median_qi(.value = .prediction) %>% # Rename prediction to value
  ungroup() %>%
  select(contains("."))

# .row is equivalent to the row number in the original dataset (feed_complete_predictors) - create a join column for this
feed_metadat<- feed_complete_predictors %>%
  select(study_id, clean_sci_name, taxa, intensity, system) %>%
  mutate(.row = row_number())

feed_predictions <- feed_dat_intervals %>%
  left_join(feed_metadat, by = ".row") %>%
  rename(feed_proportion = .value)

# ######################################################################################################
# # Step 2: Model fcr with intensity + system (no taxa)
# # FOLLOWING SECTION IS NOT APPLIACABLE
# 
# # Create dataframe for brms
# feed_brms_data_no_taxa <- data.frame(cbind(feed_no_na %>% select(contains("feed")), X_ordinal_scaled)) 
# 
# names(feed_brms_data_no_taxa)
# 
# # Set model formula
# feed_brms_no_taxa <- brmsformula(bind(feed_soy, feed_crops, feed_fmfo, feed_animal) ~ 1 + ., family = dirichlet())
# 
# # Rule of thumb: bulk and tail effective sample sizes should be 100 x number of chains (i.e., at least 400)
# # increasing max_treedepth is more about efficiency (instead of validity)
# # See: https://mc-stan.org/misc/warnings.html
# fit_feed_no_taxa <- brm(feed_brms_no_taxa, data = feed_brms_data_no_taxa,
#                         cores = 4, seed = "11729", iter = 5000)
# 
# 
# # Get stan code
# #stancode(fit_feed_no_taxa)
# 
# ######################################################################################################
# # Use model of intensity + system to predict feed for taxa that were not part of the previous model
# # Intensity must be non-NA
# feed_complete_predictors <- feed_model_dat_categories %>%
#   filter(fcr != 0 | is.na(fcr))  %>% # First drop the zeroes; Have to explicitly include is.na(fcr) otherwise NA's get dropped by fcr != 0
#   filter(is.na(intensity)==FALSE & is.na(system)==FALSE) %>%
#   filter(is.na(fcr)) # Now, filter to just the NAs
# 
# # PROBLEM: lca_complete predictors has more taxa than originally model:
# taxa_not_modeled <- setdiff(unique(feed_complete_predictors$taxa), unique(feed_no_na$taxa)) # these taxa were never modeled so they can't be predicted below
# 
# # Only keep taxa that were not modeled previously
# feed_complete_predictors_no_taxa <- feed_complete_predictors %>%
#   filter(taxa %in% taxa_not_modeled) 
# 
# # Format intensity as ordinal variable, then center and scale
# X_ordinal_no_taxa <- feed_complete_predictors_no_taxa %>%
#   mutate(intensity = factor(intensity, levels = c("Extensive", "Imp. extensive", "Semi-intensive", "Intensive"))) %>% # set order of factors (low = extensive, high = intensive)
#   mutate(intensity = as.numeric(intensity))  %>%
#   mutate(system = factor(system, levels = c("On- and off-bottom", "Cages & pens", "Ponds", "Recirculating and tanks"))) %>% # set order of factors (low = open, high = closed)
#   mutate(system = as.numeric(system)) %>%
#   select(intensity, system) %>%
#   as.matrix()
# 
# # Center and Scale: BUT now center by the mean of the original modeled dataset above AND scale by the same 2*SD calculated from the original, modeled dataset above
# X_ordinal_scaled_no_taxa <- scale(X_ordinal_no_taxa, center=apply(X_ordinal, MARGIN = 2, FUN = mean), scale=2*ordinal_sd)
# 
# # Create dataframe for brms
# brms_feed_dat_no_taxa <- data.frame(X_ordinal_scaled_no_taxa)
# 
# # Make predictions
# #predicted_feed_dat <- predict(fit_no_na, newdata = brms_new_feed_data)
# # Use tidybayes instead:
# predicted_feed_dat_no_taxa <- add_predicted_draws(newdata = brms_feed_dat_no_taxa, model = fit_feed_no_taxa)
# 
# # Get point and interval estimates from predicted data
# # Select just the prediction columns
# # Join these with the modeled data (feed_complete_predictors_no_taxa) to get metadata on taxa/intensity/syste,
# feed_dat_intervals_no_taxa <- predicted_feed_dat_no_taxa %>%
#   median_qi(.value = .prediction) %>% # Rename prediction to value
#   ungroup() %>%
#   select(contains("."))
# 
# # .row is equivalent to the row number in the modeled dataset (feed_complete_predictors_no_taxa) - create a join column for this
# feed_metadat_no_taxa <- feed_complete_predictors_no_taxa %>%
#   select(study_id, clean_sci_name, taxa, intensity, system) %>%
#   mutate(.row = row_number())
# 
# feed_predictions_no_taxa <- feed_dat_intervals_no_taxa %>%
#   left_join(feed_metadat_no_taxa, by = ".row") %>%
#   rename(feed_proportion = .value)


######################################################################################################
# Reformat data so it can row-bind with predictions
feed_no_na_long <- feed_no_na %>%
  pivot_longer(cols = c("feed_soy", "feed_crops", "feed_fmfo", "feed_animal"), names_to = ".category", values_to = "feed_proportion")

# Bind feed_no_na (data), feed_predictions, and feed_predictions_no_taxa
full_feed_dat <- feed_predictions %>%
  bind_rows(feed_no_na_long) %>%
  #bind_rows(feed_predictions_no_taxa) %>%
  mutate(data_type = if_else(is.na(.point), true = "data", false = "prediction")) %>%
  arrange(taxa, intensity, system, clean_sci_name) %>%
  rownames_to_column() # Arrange by taxa first, then create dummy column for plotting 

# Check data + predictions
feed_vars <- c("soy", "crops", "fmfo", "animal")
for (i in 1:length(feed_vars)) {
  p <- ggplot(data = full_feed_dat %>% 
                filter(str_detect(.category, feed_vars[i])), aes(y = taxa, x = feed_proportion)) +
    geom_boxplot(outlier.shape = NA) +
    #geom_violin(aes(color = taxa), scale = "width") +
    geom_jitter(aes(color = data_type), size = 3) +
    theme_classic() +
    plot_theme +
    labs(title = paste("Boxplots of ", feed_vars[i], " feed proportions", sep = ""),
         x = "",
         y = "")  +
    theme(axis.text.x = element_text(hjust = 1))
  print(p)
}

# NEXT: use plot_brms_outputs_dirichlet_regression to produce figures for SI
######################################################################################################
# RESTARTING POINT
rm(list=ls()[!(ls() %in% c("datadir", "outdir", 
                           "lca_dat_clean_groups", "feed_model_dat_categories",
                           "full_feed_dat", "full_fcr_dat"))])

#save.image(file.path(outdir, paste(Sys.Date(), "_feed-impacts-all-data-prior-to-aggregation.RData", sep = "")))
# datadir <- "/Volumes/jgephart/BFA Environment 2/Data"
# outdir <- "/Volumes/jgephart/BFA Environment 2/Outputs"
# load(file.path(outdir, "<name of file>.RData"))

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

set_allocation <- "Gross energy content"

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
                  fp_water_dat = fp_water_dat)

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
  // data for gamma model for FCR
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
  vector[K] fp_n_dat;
  vector[K] fp_p_dat;
  vector[K] fp_land_dat;
  vector[K] fp_water_dat;
}
parameters {
  // FCR model:
  vector<lower=0>[N_TX] tx_mu;
  real<lower=0> tx_sigma;
  vector<lower=0>[N_SCI] sci_mu;
  real<lower=0> sci_sigma;

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
    x[n] ~ normal(sci_mu[n_to_sci[n]], sci_sigma);
  }

  for (n_sci in 1:N_SCI){
    sci_mu[n_sci] ~ normal(tx_mu[sci_to_tx[n_sci]], tx_sigma);
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
  vector[N_TX] tx_n_fp;
  vector[N_SCI] sci_n_fp;

  // Phosphorus
  vector[N_TX] tx_p_fp;
  vector[N_SCI] sci_p_fp;

  // Land
  vector[N_TX] tx_land_fp;
  vector[N_SCI] sci_land_fp;

  // Water
  vector[N_TX] tx_water_fp;
  vector[N_SCI] sci_water_fp;

  // Calculations
  for (n_tx in 1:N_TX) {
    tx_c_fp[n_tx] = tx_mu[n_tx] * sum(fp_c_dat .* tx_theta[n_tx]);
    tx_n_fp[n_tx] = tx_mu[n_tx] * sum(fp_n_dat .* tx_theta[n_tx]);
    tx_p_fp[n_tx] = tx_mu[n_tx] * sum(fp_p_dat .* tx_theta[n_tx]);
    tx_land_fp[n_tx] = tx_mu[n_tx] * sum(fp_land_dat .* tx_theta[n_tx]);
    tx_water_fp[n_tx] = tx_mu[n_tx] * sum(fp_water_dat .* tx_theta[n_tx]);
  }
  for (n_sci in 1:N_SCI) {
    sci_c_fp[n_sci] = sci_mu[n_sci] * sum(fp_c_dat .* sci_theta[n_sci]);
    sci_n_fp[n_sci] = sci_mu[n_sci] * sum(fp_n_dat .* sci_theta[n_sci]);
    sci_p_fp[n_sci] = sci_mu[n_sci] * sum(fp_p_dat .* sci_theta[n_sci]);
    sci_land_fp[n_sci] = sci_mu[n_sci] * sum(fp_land_dat .* sci_theta[n_sci]);
    sci_water_fp[n_sci] = sci_mu[n_sci] * sum(fp_water_dat .* sci_theta[n_sci]);
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





