# Impute all missing values for FCR and feed proportions (will need these data for each of the total on-farm + off-farm impacts)

# First run: process_data_for_analysis.R, to make lca_dat_clean_groups
lca_dat_clean_group <- read.csv(lca_dat_clean_groups, file.path(datadir, "lca_clean_with_groups.csv"), row.names = FALSE)

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
# OUTPUT
rm(list=ls()[!(ls() %in% c("datadir", "outdir", 
                           "lca_dat_clean_groups", "feed_model_dat_categories",
                           "full_feed_dat", "full_fcr_dat"))])
# When data is finalized, write to CSV
#write.csv()
#save.image(file.path(outdir, paste(Sys.Date(), "_feed-impacts-all-data-prior-to-aggregation.RData", sep = "")))

######################################################################################################