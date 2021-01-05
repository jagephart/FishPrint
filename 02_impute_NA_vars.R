# Impute all missing values for:
# FCR
# Feed proportions
# All energy inputs: Electricity, Diesel, Petrol, Natural Gas
# Yield

# First run: 01_process_data_for_analysis.R, to make lca_dat_clean_groups
rm(list=ls()[!(ls() %in% c("lca_dat_clean_groups"))])

# Or just read in:
#lca_dat_clean_groups <- read.csv(lca_dat_clean_groups, file.path(datadir, "lca_clean_with_groups.csv"), row.names = FALSE)

######################################################################################################
# Section 1: Create feed_model_dat_categories for modeling FCR and feed proportions

# Remove species that aren't fed: bivalves, plants
# SELECT STUDY ID COLUMN - use this for rejoining outputs from multiple regression models back together
# Select relevant data columns and arrange by categorical info
feed_model_dat_categories <- lca_dat_clean_groups %>%
  select(study_id, fcr = FCR, contains("new"), clean_sci_name, taxa, intensity = Intensity, system = Production_system_group) %>%
  arrange(clean_sci_name, taxa, intensity, system) %>%
  filter(taxa %in% c("bivalves", "plants")==FALSE) # Remove taxa that don't belong in FCR/feed analysis - bivalves

######################################################################################################
# Step 1: Model FCR as taxa + intensity + system
# ie, remove all NAs and model these data, then use the model to predict FCR for those data with a complete set of predictors

fcr_no_na <- feed_model_dat_categories %>%
  #filter(fcr != 0 | is.na(fcr))  %>% # If we want to retain NAs, have to explicitly include is.na(fcr) otherwise NA's get dropped by fcr != 0
  filter(fcr != 0) %>% # This also automatically drops NAs
  filter(is.na(intensity)==FALSE & is.na(system)==FALSE) %>% # complete predictors - i.e., both intensity AND system are non-NA
  select(study_id, fcr, clean_sci_name, taxa, intensity, system)

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
# No need to create model of just inetnsity + system because all taxa were predicted
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

# Generate summaries for supplementary information:
source("Functions.R") 
plot_for_si(name_of_fit = "fit_fcr_no_na", name_of_data = "full_fcr_dat", name_of_var = "fcr")

rm(list=ls()[!(ls() %in% c("lca_dat_clean_groups", "feed_model_dat_categories",
                           "full_fcr_dat"))])

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

######################################################################################################
# No need to create model of just inetnsity + system because all taxa were predicted
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
    labs(title = paste("Boxplots of ", feed_vars[i], " feed proportions", sep = ""),
         x = "",
         y = "")  +
    theme(axis.text.x = element_text(hjust = 1))
  print(p)
}

# Generate summaries for supplementary information:
source("Functions.R") 
plot_for_si(name_of_fit = "fit_feed_no_na", name_of_data = "full_feed_dat", name_of_var = "feed_proportion", regression_type = "dirichlet")

rm(list=ls()[!(ls() %in% c("lca_dat_clean_groups", "feed_model_dat_categories",
                           "full_fcr_dat",
                           "full_feed_dat"))])
######################################################################################################
# Section 2: Create ghg_model_dat_categories for modeling all energy inputs

# DECISION: adjusting zeroes to be very small amount for all energy datasets
# Need to do this for gamma regression

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
# No need to create model of just inetnsity + system because all taxa were predicted
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
# Outlier in misc_diad looks like it's being pulled to a much lower electricity value because system == "Cages & pens"
# See plot for just the Cages and pens studies:
ggplot(full_electric_dat %>% filter(system == "Cages & pens"), aes(x = electric, y = taxa)) + geom_boxplot(outlier.shape = NA) + geom_jitter(aes(color = data_type))

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

# Generate summaries for supplementary information:
source("Functions.R") 
plot_for_si(name_of_fit = "fit_electric_no_na", name_of_data = "full_electric_dat", name_of_var = "electric")

rm(list=ls()[!(ls() %in% c("lca_dat_clean_groups", 
                           "feed_model_dat_categories",
                           "full_fcr_dat",
                           "full_feed_dat",
                           "ghg_model_dat_categories", 
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
# No need to create model of just inetnsity + system because all taxa were predicted
######################################################################################################

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

# Generate summaries for supplementary information
source("Functions.R") 
plot_for_si(name_of_fit = "fit_diesel_no_na", name_of_data = "full_diesel_dat", name_of_var = "diesel")

rm(list=ls()[!(ls() %in% c("lca_dat_clean_groups", 
                           "feed_model_dat_categories",
                           "full_fcr_dat",
                           "full_feed_dat",
                           "ghg_model_dat_categories", 
                           "full_electric_dat",
                           "full_diesel_dat"))])

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
# No need to create model of just inetnsity + system because all taxa were predicted
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

# Generate summaries for supplementary information
source("Functions.R") 
plot_for_si(name_of_fit = "fit_petrol_no_na", name_of_data = "full_petrol_dat", name_of_var = "petrol")

rm(list=ls()[!(ls() %in% c("lca_dat_clean_groups", 
                           "feed_model_dat_categories",
                           "full_fcr_dat",
                           "full_feed_dat",
                           "ghg_model_dat_categories", 
                           "full_electric_dat",
                           "full_diesel_dat",
                           "full_petrol_dat"))])

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
# No need to create model of just inetnsity + system because all taxa were predicted
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

# Generate summaries for supplementary information
source("Functions.R") 
plot_for_si(name_of_fit = "fit_natgas_no_na", name_of_data = "full_natgas_dat", name_of_var = "natgas")

rm(list=ls()[!(ls() %in% c("lca_dat_clean_groups", 
                           "feed_model_dat_categories",
                           "full_fcr_dat",
                           "full_feed_dat",
                           "ghg_model_dat_categories", 
                           "full_electric_dat",
                           "full_diesel_dat",
                           "full_petrol_dat",
                           "full_natgas_dat"))])

######################################################################################################
# Section 3: Create land_model_dat_categories for modeling yield (i.e., land)

# Get model-specific data: Only need ponds and recirculating and tanks for this section
# SELECT STUDY ID COLUMN - use this for rejoining outputs from multiple regression models back together
# Select relevant data columns and arrange by categorical info
# Clean harvest data (do all data cleaning like this up top) or move this to Functions.R
land_model_dat_categories <- lca_dat_clean_groups %>%
  select(study_id, yield = Yield_m2_per_t, clean_sci_name, taxa, intensity = Intensity, system = Production_system_group) %>%
  filter(system %in% c("Ponds", "Recirculating and tanks")) %>%
  arrange(clean_sci_name, taxa, intensity, system)

######################################################################################################

# Step 1: Model yield as taxa + intensity + system
# Remove all NAs and model these data, then use the model to predict yield for those data with a complete set of predictors
yield_no_na <- land_model_dat_categories %>%
  filter(yield != 0)  %>% # This also automatically drops yield == NA
  filter(is.na(intensity)==FALSE & is.na(system)==FALSE) %>% # complete predictors - i.e., both intensity AND system are non-NA
  select(study_id, yield, clean_sci_name, taxa, intensity, system)

# See which ones have missing predictors: land_model_dat_categories %>% filter(is.na(yield)==FALSE & (is.na(intensity) | is.na(system)))

# Set data for model:
# Create model matrix for taxa info, then center and scale
X_taxa <- model.matrix(object = ~ 1 + taxa, 
                       data = yield_no_na %>% select(taxa)) 

taxa_sd <- apply(X_taxa[,-1], MARGIN=2, FUN=sd, na.rm=TRUE) # Center all non-intercept variables and scale by 2 standard deviations (ignoring NAs)
X_taxa_scaled <- scale(X_taxa[,-1], center=TRUE, scale=2*taxa_sd)

# Format intensity and system as ordinal variable, then center and scale
X_ordinal <- yield_no_na %>%
  mutate(intensity = factor(intensity, levels = c("Extensive", "Imp. extensive", "Semi-intensive", "Intensive"))) %>% # set order of factors (low = extensive, high = intensive)
  mutate(system = factor(system, levels = c("On- and off-bottom", "Cages & pens", "Ponds", "Recirculating and tanks"))) %>% # set order of factors (low = open, high = closed) 
  mutate(intensity = as.numeric(intensity)) %>%
  mutate(system = as.numeric(system)) %>%
  select(intensity, system) %>%
  as.matrix()
ordinal_sd<-apply(X_ordinal, MARGIN=2, FUN=sd, na.rm=TRUE) # Center all non-intercept variables and scale by 2 standard deviations (ignoring NAs)
X_ordinal_scaled <- scale(X_ordinal, center=TRUE, scale=2*ordinal_sd)

# Create dataframe for brms
yield_brms_data <- data.frame(y = yield_no_na$yield, X_taxa_scaled, X_ordinal_scaled)

names(yield_brms_data)

# Set model formula
yield_brms <- brmsformula(y ~ 1 + ., family = Gamma("log"))

# Use "resp = <response_variable>" to specify different priors for different response variables
all_priors <- c(set_prior("normal(0,5)", class = "b"), # priors for y response variables
                set_prior("normal(0,2.5)", class = "Intercept"), 
                set_prior("exponential(1)", class = "shape"))

# Model converges after increasing the adapt_delta and iterations from default values
# Rule of thumb: bulk and tail effective sample sizes should be 100 x number of chains (i.e., at least 400)
# increasing max_treedepth is more about efficiency (instead of validity)
# See: https://mc-stan.org/misc/warnings.html
fit_yield_no_na <- brm(yield_brms, data = yield_brms_data,
                       prior = all_priors, cores = 4, seed = "11729", iter = 5000, control = list(adapt_delta = 0.99))

# Get stan code
#stancode(fit_yield_no_na)

######################################################################################################
# Use model to predict NAs for studies with complete set of predictors
# Intensity and system must be non-NA
yield_complete_predictors <- land_model_dat_categories %>%
  filter(is.na(intensity)==FALSE & is.na(system)==FALSE) %>%
  filter(is.na(yield)) # Now, filter to just the NAs

taxa_not_modeled <- setdiff(unique(yield_complete_predictors$taxa), unique(yield_no_na$taxa)) # these taxa were never modeled so they can't be predicted below

# DROP THESE FOR NOW:
yield_complete_predictors <- yield_complete_predictors %>%
  filter(taxa %in% taxa_not_modeled == FALSE)

# Now check the other way, which taxa were in the original model but not a part of the data that needs to be predicted:
setdiff(unique(yield_no_na$taxa), unique(yield_complete_predictors$taxa))

# If original model has taxa that are not part of yield_complete_predictors, 
# Use list of unique taxa in original model and use this to expand/assign levels manually - having trouble automating this
# Include all levels here, but remember the first level won't show up in design matrix - instead, it's part of the "contrasts"
sort(unique(yield_no_na$taxa))
yield_complete_predictors <- yield_complete_predictors %>%
  mutate(taxa = as.factor(taxa))
levels(yield_complete_predictors$taxa) <- list(catfish = "catfish",
                                               hypoph_carp = "hypoph_carp",
                                               milkfish = "milkfish",
                                               misc_marine = "misc_marine",
                                               oth_carp = "oth_carp", 
                                               salmon = "salmon", 
                                               shrimp = "shrimp",
                                               tilapia = "tilapia",
                                               trout = "trout")

# Create NEW taxa model matrix for the studies whose feeds need to be predicted
# Taxa categories:
X_taxa_new <- model.matrix(object = ~ 1 + taxa, 
                           data = yield_complete_predictors %>% select(taxa)) 

# Center and Scale: BUT now center by the mean of the original modeled dataset above AND scale by the same 2*SD calculated from the original, modeled dataset above
X_taxa_new_scaled <- scale(X_taxa_new[,-1], center=apply(X_taxa[,-1], MARGIN = 2, FUN = mean), scale=2*taxa_sd)

# Still only using intensity (no system)
# Format intensity as ordinal variable, then center and scale
X_ordinal_new <- yield_complete_predictors %>%
  mutate(intensity = factor(intensity, levels = c("Extensive", "Imp. extensive", "Semi-intensive", "Intensive"))) %>% # set order of factors (low = extensive, high = intensive)
  mutate(system = factor(system, levels = c("On- and off-bottom", "Cages & pens", "Ponds", "Recirculating and tanks"))) %>% # set order of factors (low = open, high = closed) 
  mutate(intensity = as.numeric(intensity)) %>%
  mutate(system = as.numeric(system)) %>%
  select(intensity, system) %>%
  as.matrix()

# Center and Scale: BUT now center by the mean of the original modeled dataset above AND scale by the same 2*SD calculated from the original, modeled dataset above
X_ordinal_new_scaled <- scale(X_ordinal_new, center=apply(X_ordinal, MARGIN = 2, FUN = mean), scale=2*ordinal_sd)

# Create dataframe for brms
brms_new_yield_dat <- data.frame(cbind(X_taxa_new_scaled, X_ordinal_new_scaled)) 

# Make predictions
#predicted_yield_dat <- predict(fit_no_na, newdata = brms_new_yield_data)
# Use tidybayes instead:
predicted_yield_dat <- add_predicted_draws(newdata = brms_new_yield_dat, model = fit_yield_no_na)

# Get point and interval estimates from predicted data
# Select just the prediction columns
# Join these with the modeled data (yield_complete_predictors) to get metadata on taxa/intensity/syste,
yield_dat_intervals <- predicted_yield_dat %>%
  median_qi(.value = .prediction) %>% # Rename prediction to value
  ungroup() %>%
  select(contains("."))

# .row is equivalent to the row number in the modeled dataset (yield_complete_predictors) - create a join column for this
yield_metadat<- yield_complete_predictors %>%
  select(study_id, clean_sci_name, taxa, intensity, system) %>%
  mutate(.row = row_number())

yield_predictions <- yield_dat_intervals %>%
  left_join(yield_metadat, by = ".row") %>%
  rename(yield = .value)

######################################################################################################
# Step 2: Model yield with intensity + system (no taxa)

# Create dataframe for brms
yield_brms_data_no_taxa <- data.frame(y = yield_no_na$yield, X_ordinal_scaled)

names(yield_brms_data_no_taxa)

# Set model formula
yield_brms_no_taxa <- brmsformula(y ~ 1 + ., family = Gamma("log"))

# Use "resp = <response_variable>" to specify different priors for different response variables
all_priors <- c(set_prior("normal(0,5)", class = "b"), # priors for y response variables
                set_prior("normal(0,2.5)", class = "Intercept"), 
                set_prior("exponential(1)", class = "shape"))

# Model converges after increasing the adapt_delta and iterations from default values
# Rule of thumb: bulk and tail effective sample sizes should be 100 x number of chains (i.e., at least 400)
# increasing max_treedepth is more about efficiency (instead of validity)
# See: https://mc-stan.org/misc/warnings.html
fit_yield_no_taxa <- brm(yield_brms_no_taxa, data = yield_brms_data_no_taxa,
                         prior = all_priors, cores = 4, seed = "11729", iter = 5000, control = list(adapt_delta = 0.99))

# Get stan code
#stancode(fit_yield_no_taxa)

######################################################################################################
# Use model of intensity + system to predict yield for taxa that were not part of the previous model
# Intensity must be non-NA
yield_complete_predictors <- land_model_dat_categories %>%
  filter(is.na(intensity)==FALSE & is.na(system)==FALSE) %>%
  filter(is.na(yield)) # Now, filter to just the NAs

taxa_not_modeled <- setdiff(unique(yield_complete_predictors$taxa), unique(yield_no_na$taxa)) # these taxa were never modeled so they can't be predicted below

# Only keep taxa that were not modeled previously
yield_complete_predictors_no_taxa <- yield_complete_predictors %>%
  filter(taxa %in% taxa_not_modeled) 

# Format intensity as ordinal variable, then center and scale
X_ordinal_no_taxa <- yield_complete_predictors_no_taxa %>%
  mutate(intensity = factor(intensity, levels = c("Extensive", "Imp. extensive", "Semi-intensive", "Intensive"))) %>% # set order of factors (low = extensive, high = intensive)
  mutate(system = factor(system, levels = c("On- and off-bottom", "Cages & pens", "Ponds", "Recirculating and tanks"))) %>% # set order of factors (low = open, high = closed) 
  mutate(intensity = as.numeric(intensity)) %>%
  mutate(system = as.numeric(system)) %>%
  select(intensity, system) %>%
  as.matrix()

# Center and Scale: BUT now center by the mean of the original modeled dataset above AND scale by the same 2*SD calculated from the original, modeled dataset above
X_ordinal_scaled_no_taxa <- scale(X_ordinal_no_taxa, center=apply(X_ordinal, MARGIN = 2, FUN = mean), scale=2*ordinal_sd)

# Create dataframe for brms
brms_yield_dat_no_taxa <- data.frame(X_ordinal_scaled_no_taxa)

# Make predictions
#predicted_yield_dat <- predict(fit_no_na, newdata = brms_new_yield_data)
# Use tidybayes instead:
predicted_yield_dat_no_taxa <- add_predicted_draws(newdata = brms_yield_dat_no_taxa, model = fit_yield_no_taxa)

# Get point and interval estimates from predicted data
# Select just the prediction columns
# Join these with the modeled data (yield_complete_predictors_no_taxa) to get metadata on taxa/intensity/syste,
yield_dat_intervals_no_taxa <- predicted_yield_dat_no_taxa %>%
  median_qi(.value = .prediction) %>% # Rename prediction to value
  ungroup() %>%
  select(contains("."))

# .row is equivalent to the row number in the modeled dataset (yield_complete_predictors_no_taxa) - create a join column for this
yield_metadat_no_taxa <- yield_complete_predictors_no_taxa %>%
  select(study_id, clean_sci_name, taxa, intensity, system) %>%
  mutate(.row = row_number())

yield_predictions_no_taxa <- yield_dat_intervals_no_taxa %>%
  left_join(yield_metadat_no_taxa, by = ".row") %>%
  rename(yield = .value)

# Bind yield_no_na (data), yield_predictions, and yield_predictions_no_taxa
full_yield_dat <- yield_predictions %>%
  bind_rows(yield_no_na) %>%
  bind_rows(yield_predictions_no_taxa) %>%
  mutate(data_type = if_else(is.na(.point), true = "data", false = "prediction")) %>%
  arrange(taxa, intensity, system, clean_sci_name) %>%
  rownames_to_column() # Arrange by taxa first, then create dummy column for plotting 

# Quick check: PLOT DATA + PREDICTIONS
ggplot(full_yield_dat, aes(x = yield, y = taxa)) + geom_boxplot(outlier.shape = NA) + geom_jitter(aes(color = data_type))

# Generate summaries for supplementary information
source("Functions.R") 
plot_for_si(name_of_fit = "fit_yield_no_na", name_of_data = "full_yield_dat", name_of_var = "yield")

rm(list=ls()[!(ls() %in% c("lca_dat_clean_groups", 
                           "feed_model_dat_categories",
                           "full_fcr_dat",
                           "full_feed_dat",
                           "ghg_model_dat_categories", 
                           "full_electric_dat",
                           "full_diesel_dat",
                           "full_petrol_dat",
                           "full_natgas_dat",
                           "land_model_dat_categories",
                           "full_yield_dat"))])

######################################################################################################
# Create master data frame of all variables before beginning bayesian hierarchical analyses

# Rename columns in lca_dat_clean_groups for merging
lca_dat_clean_groups_merge <- lca_dat_clean_groups %>%
  rename(intensity = Intensity, system = Production_system_group,
         fcr = FCR,
         feed_soy = feed_soy_new,
         feed_crops = feed_crops_new,
         feed_fmfo = feed_fmfo_new,
         feed_animal = feed_animal_new)

######################################################################################################
# Get back data that have incomplete system/intensity predictors (i.e., were not part of the regressions above and merge with full_feed_dat)

# FIX IT - change clean.lca in Functions.R so no less than 0.01
# Adjust feed proportions to be no less than 0.01
# TEMPORARY FIX:
# Normalize the FINAL feed proportion values to be greater than 0 and no less than 0.01
full_feed_dat <- full_feed_dat %>%
  mutate(feed_proportion = if_else(feed_proportion < 0.01, true = 0.0105, false = feed_proportion))

# Feed
feed_dat_incomplete_predictors <- lca_dat_clean_groups_merge %>%
  filter(study_id %in% setdiff(lca_dat_clean_groups_merge$study_id, full_feed_dat$study_id)) %>%
  select(study_id, clean_sci_name, taxa, intensity, system,  contains("feed")) %>%
  mutate(feed_data_type = "data") %>%
  # Do the same temporary fix for the original LCA data
  mutate(feed_soy = if_else(feed_soy < 0.01, true = 0.0105, false = feed_soy),
         feed_crops = if_else(feed_crops < 0.01, true = 0.0105, false = feed_crops),
         feed_fmfo = if_else(feed_fmfo < 0.01, true = 0.0105, false = feed_fmfo),
         feed_animal = if_else(feed_animal < 0.01, true = 0.0105, false = feed_animal)) %>%
  mutate(feed_soy = round(feed_soy, digits = 4),
         feed_crops = round(feed_crops, digits = 4),
         feed_fmfo = round(feed_fmfo, digits = 4),
         feed_animal = round(feed_animal, digits = 4)) %>%
  mutate(new_scale = feed_soy + feed_crops + feed_fmfo + feed_animal) %>%
  mutate(feed_soy = feed_soy / new_scale,
         feed_crops = feed_crops / new_scale,
         feed_fmfo = feed_fmfo / new_scale,
         feed_animal = feed_animal / new_scale) %>%
  select(-c(new_scale, Feed_type)) 

# FCR
fcr_dat_incomplete_predictors <- lca_dat_clean_groups_merge %>%
  filter(study_id %in% setdiff(lca_dat_clean_groups_merge$study_id, full_fcr_dat$study_id)) %>%
  select(study_id, clean_sci_name, taxa, intensity, system, fcr) %>%
  mutate(fcr_data_type = "data") 

# Electricity
electric_dat_incomplete_predictors <- lca_dat_clean_groups_merge %>%
  filter(study_id %in% setdiff(lca_dat_clean_groups_merge$study_id, full_electric_dat$study_id)) %>%
  select(study_id, Country, iso3c, clean_sci_name, taxa, intensity, system, Electricity_kwh) %>%
  mutate(electric_data_type = "data") 

# Diesel
diesel_dat_incomplete_predictors <- lca_dat_clean_groups_merge %>%
  filter(study_id %in% setdiff(lca_dat_clean_groups_merge$study_id, full_diesel_dat$study_id)) %>%
  select(study_id, Country, iso3c, clean_sci_name, taxa, intensity, system, Diesel_L) %>%
  mutate(diesel_data_type = "data") 

# Petrol
petrol_dat_incomplete_predictors <- lca_dat_clean_groups_merge %>%
  filter(study_id %in% setdiff(lca_dat_clean_groups_merge$study_id, full_petrol_dat$study_id)) %>%
  select(study_id, Country, iso3c, clean_sci_name, taxa, intensity, system, Petrol_L) %>%
  mutate(petrol_data_type = "data") 

# Natural Gas
natgas_dat_incomplete_predictors <- lca_dat_clean_groups_merge %>%
  filter(study_id %in% setdiff(lca_dat_clean_groups_merge$study_id, full_natgas_dat$study_id)) %>%
  select(study_id, Country, iso3c, clean_sci_name, taxa, intensity, system, NaturalGas_L) %>%
  mutate(natgas_data_type = "data") 

# Yield/Land
yield_dat_incomplete_predictors <- lca_dat_clean_groups_merge %>%
  filter(study_id %in% setdiff(lca_dat_clean_groups_merge$study_id, full_yield_dat$study_id)) %>%
  select(study_id, clean_sci_name, taxa, intensity, system, Yield_m2_per_t) %>%
  mutate(yield_data_type = "data") 

######################################################################################################
# Join dat with incomplete predictors with regression dat (full_dat)

# Feed
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
         feed_animal = feed_animal / new_scale) %>%
  select(-new_scale) %>%
  full_join(feed_dat_incomplete_predictors, by = intersect(names(.), names(feed_dat_incomplete_predictors)))

# FCR
fcr_dat_merge <- full_fcr_dat %>%
  select(study_id, clean_sci_name, taxa, intensity, system, fcr_data_type = data_type, fcr) %>%
  full_join(fcr_dat_incomplete_predictors, by = intersect(names(.), names(fcr_dat_incomplete_predictors)))

# Electricity - rename to get units back
electric_dat_merge <- full_electric_dat %>%
  select(study_id, Country, iso3c, clean_sci_name, taxa, intensity, system, electric_data_type = data_type, Electricity_kwh = electric) %>%
  full_join(electric_dat_incomplete_predictors, by = intersect(names(.), names(electric_dat_incomplete_predictors)))

# Diesel - rename to get units back
diesel_dat_merge <- full_diesel_dat %>%
  select(study_id, Country, iso3c, clean_sci_name, taxa, intensity, system, diesel_data_type = data_type, Diesel_L = diesel) %>%
  full_join(diesel_dat_incomplete_predictors, by = intersect(names(.), names(diesel_dat_incomplete_predictors)))

# Petrol - rename to get units back
petrol_dat_merge <- full_petrol_dat %>%
  select(study_id, Country, iso3c, clean_sci_name, taxa, intensity, system, petrol_data_type = data_type, Petrol_L = petrol) %>%
  full_join(petrol_dat_incomplete_predictors, by = intersect(names(.), names(petrol_dat_incomplete_predictors)))

# Natural gas - rename to get units back
natgas_dat_merge <- full_natgas_dat %>%
  select(study_id, Country, iso3c, clean_sci_name, taxa, intensity, system, natgas_data_type = data_type, NaturalGas_L = natgas) %>%
  full_join(natgas_dat_incomplete_predictors, by = intersect(names(.), names(natgas_dat_incomplete_predictors)))

# Yield- rename to get units back
yield_dat_merge <- full_yield_dat %>%
  select(study_id, clean_sci_name, taxa, intensity, system, yield_data_type = data_type, Yield_m2_per_t = yield) %>%
  full_join(yield_dat_incomplete_predictors, by = intersect(names(.), names(yield_dat_incomplete_predictors)))

######################################################################################################
# MERGE:
lca_dat_imputed <- feed_dat_merge %>%
  full_join(fcr_dat_merge, by = intersect(names(feed_dat_merge), names(fcr_dat_merge))) %>%
  full_join(electric_dat_merge, by = intersect(names(.), names(electric_dat_merge))) %>%
  full_join(diesel_dat_merge, by = intersect(names(.), names(diesel_dat_merge))) %>%
  full_join(petrol_dat_merge, by = intersect(names(.), names(petrol_dat_merge))) %>%
  full_join(natgas_dat_merge, by = intersect(names(.), names(natgas_dat_merge))) %>%
  full_join(yield_dat_merge, by = intersect(names(.), names(yield_dat_merge))) %>%
  # BIVALVES AND PLANTS SHOULD GET AN FCR OF 0 - this way, NA is specific to MISSING DATA
  mutate(fcr = if_else(is.na(fcr) & taxa %in% c("bivalves", "plants"), true = 0, false = fcr)) %>%
  # FEEDS SHOULD BE 0 (not NA) WHENEVER FCR IS 0
  mutate(feed_soy = if_else(fcr==0, true = 0, false = feed_soy),
         feed_crops = if_else(fcr==0, true = 0, false = feed_crops),
         feed_fmfo = if_else(fcr==0, true = 0, false = feed_fmfo),
         feed_animal = if_else(fcr==0, true = 0, false = feed_animal))

datadir <- "/Volumes/jgephart/BFA Environment 2/Data"
write.csv(lca_dat_imputed, file.path(datadir, "lca-dat-imputed-vars.csv"), row.names = FALSE)
