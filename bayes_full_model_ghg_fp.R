# Calculate (non-feed associated) greenhouse gas footprint
# (Electricity use * country-specific GHG of electricity) + (Diesel * GHG of diesel) + (Petrol * GHG of petrol) + (Natural gas * GHG of natural gas)

# Reminder for new data: entries with data for some columns should have no blanks in other columns (these should be filled in as zeroes)

# Step 0: Run process_data_for_analysis.R, then clear environment other than:
rm(list=ls()[!(ls() %in% c("lca_dat_clean_groups", "datadir", "outdir"))])

# Get model-specific data:
# SELECT STUDY ID COLUMN - use this for rejoining outputs from multiple regression models back together
# Select relevant data columns and arrange by categorical info
ghg_model_dat_categories <- lca_dat_clean_groups %>%
  select(study_id, electricity = Electricity_kwh, diesel = Diesel_L, petrol = Petrol_L, natural_gas = NaturalGas_L, clean_sci_name, taxa, intensity = Intensity, system = Production_system_group) %>%
  arrange(clean_sci_name, taxa, intensity, system)

######################################################################################################
# Step 1: Model electricity

elec_no_na <- ghg_model_dat_categories %>%
  filter(is.na(electricity)==FALSE)  %>% # Drop NAs (Keep zeroes)
  mutate(electricity = if_else(electricity == 0, true = min(electricity[electricity!=0]), false = electricity)) %>% # Not modeling the zeroes, adjust these to the minimum value
  filter(is.na(intensity)==FALSE & is.na(system)==FALSE) %>% # complete predictors - i.e., both intensity AND system are non-NA
  select(study_id, electricity, clean_sci_name, taxa, intensity, system)

# Create model matrix for taxa info, then center and scale
X_taxa <- model.matrix(object = ~ 1 + taxa, 
                       data = elec_no_na %>% select(taxa)) 

taxa_sd <- apply(X_taxa[,-1], MARGIN=2, FUN=sd, na.rm=TRUE) # Center all non-intercept variables and scale by 2 standard deviations (ignoring NAs)
X_taxa_scaled <- scale(X_taxa[,-1], center=TRUE, scale=2*taxa_sd)

# Format intensity and system as ordinal variables, then center and scale
X_ordinal <- elec_no_na %>%
  mutate(intensity = factor(intensity, levels = c("extensive", "semi", "intensive"))) %>% # set order of factors (low = extensive, high = intensive)
  mutate(system = factor(system, levels = c("open", "semi", "closed"))) %>% # set order of factors (low = open, high = closed)
  mutate(intensity = as.numeric(intensity)) %>%
  mutate(system = as.numeric(system)) %>%
  select(intensity, system) %>%
  as.matrix()
ordinal_sd<-apply(X_ordinal, MARGIN=2, FUN=sd, na.rm=TRUE) # Center all non-intercept variables and scale by 2 standard deviations (ignoring NAs)
X_ordinal_scaled <- scale(X_ordinal, center=TRUE, scale=2*ordinal_sd)

# Create dataframe for brms and rename feed variables
elec_brms_data <- data.frame(y = elec_no_na$electricity, X_taxa_scaled, X_ordinal_scaled)

names(elec_brms_data)

# Set model formula
# elec_brms <- brmsformula(y ~ 1 + ., family = Gamma("log"))
# Equivalent to:
# elec_brms <- brmsformula(y ~ 1 + taxatilapia + taxatrout +
#                            intensity + system, family = Gamma("log"))


# Use "resp = <response_variable>" to specify different priors for different response variables
all_priors <- c(set_prior("normal(0,5)", class = "b"), # priors for y response variables
                set_prior("normal(0,2.5)", class = "Intercept"), 
                set_prior("exponential(1)", class = "shape"))

# Model converges after increasing the adapt_delta and iterations from default values
# Rule of thumb: bulk and tail effective sample sizes should be 100 x number of chains (i.e., at least 400)
# increasing max_treedepth is more about efficiency (instead of validity)
# See: https://mc-stan.org/misc/warnings.html
fit_elec_no_na <- brm(elec_brms, data = elec_brms_data,
                       prior = all_priors, cores = 4, seed = "11729", iter = 20000, control = list(adapt_delta = 0.99))

# Get stan code
stancode(fit_elec_no_na)

######################################################################################################
# Use model to predict NAs for studies with complete set of predictors
# Both intensity AND system are non-NA
elec_complete_predictors <- ghg_model_dat_categories %>%
  filter(is.na(electricity)) %>% 
  filter(is.na(intensity)==FALSE & is.na(system)== FALSE) 

# PROBLEM: lca_complete predictors has more taxa than originally model:
taxa_not_modeled <- setdiff(unique(elec_complete_predictors$taxa), unique(elec_no_na$taxa)) # these taxa were never modeled so they can't be predicted below

# DROP THESE FOR NOW:
elec_complete_predictors <- elec_complete_predictors %>%
  filter(taxa %in% taxa_not_modeled == FALSE)

# Now check the other way, which taxa were in the original model but not a part of the data that needs to be predicted:
setdiff(unique(elec_no_na$taxa), unique(elec_complete_predictors$taxa))

# If original model has taxa that are not part of elec_complete_predictors, 
# Use list of unique taxa in original model and use this to expand/assign levels manually - having trouble automating this
# Include all levels here, but remember the first level won't show up in design matrix - instead, it's part of the "contrasts"
# sort(unique(elec_no_na$taxa))
# elec_complete_predictors <- elec_complete_predictors %>%
#   mutate(taxa = as.factor(taxa))
# levels(elec_complete_predictors$taxa) <- list(fresh_crust = "fresh_crust", misc_fresh = "misc_fresh", misc_marine = "misc_marine", oth_carp = "oth_carp", tilapia = "tilapia")

# Create NEW taxa model matrix for the studies to be predicted
# Taxa categories:
X_taxa_new <- model.matrix(object = ~ 1 + taxa, 
                           data = elec_complete_predictors %>% select(taxa)) 

# Center and Scale: BUT now center by the mean of the original modeled dataset above AND scale by the same 2*SD calculated from the original, modeled dataset above
X_taxa_new_scaled <- scale(X_taxa_new[,-1], center=apply(X_taxa[,-1], MARGIN = 2, FUN = mean), scale=2*taxa_sd)

# System and Intensity variables:
# Format intensity and system as ordinal variables, then center and scale
X_ordinal_new <- elec_complete_predictors %>%
  mutate(intensity = factor(intensity, levels = c("extensive", "semi", "intensive"))) %>% # set order of factors (low = extensive, high = intensive)
  mutate(system = factor(system, levels = c("open", "semi", "closed"))) %>% # set order of factors (low = open, high = closed)
  mutate(intensity = as.numeric(intensity)) %>%
  mutate(system = as.numeric(system)) %>%
  select(intensity, system) %>%
  as.matrix()

# Center and Scale: BUT now center by the mean of the original modeled dataset above AND scale by the same 2*SD calculated from the original, modeled dataset above
X_ordinal_new_scaled <- scale(X_ordinal_new, center=apply(X_ordinal, MARGIN = 2, FUN = mean), scale=2*ordinal_sd)

# Create dataframe for brms and rename feed variables
brms_new_elec_dat <- data.frame(cbind(X_taxa_new_scaled, X_ordinal_new_scaled)) 

# Make predictions
#predicted_elec_dat <- predict(fit_no_na, newdata = brms_new_elec_data)
# Use tidybayes instead:
predicted_elec_dat <- add_predicted_draws(newdata = brms_new_elec_dat, model = fit_elec_no_na)

# Get point and interval estimates from predicted data
# Select just the prediction columns
# Join these with the modeled data (elec_complete_predictors) to get metadata on taxa/intensity/syste,
elec_dat_intervals <- predicted_elec_dat %>%
  median_qi(.value = .prediction) %>% # Rename prediction to value
  ungroup() %>%
  select(contains("."))

# .row is equivalent to the row number in the modeled dataset (elec_complete_predictors) - create a join column for this
elec_metadat<- elec_complete_predictors %>%
  select(study_id, clean_sci_name, taxa, intensity, system) %>%
  mutate(.row = row_number())

elec_predictions <- elec_dat_intervals %>%
  left_join(elec_metadat, by = ".row") %>%
  rename(electricity = .value)

# Bind
full_elec_dat <- elec_predictions %>%
  bind_rows(elec_no_na) %>%
  mutate(data_type = if_else(is.na(.point), true = "data", false = "prediction")) %>%
  arrange(taxa, intensity, system, clean_sci_name) %>%
  rownames_to_column() # Arrange by taxa first, then create dummy column for plotting 

# NEXT: use plot_brms_gamma_regression to produce figures for SI

######################################################################################################
# Step 2: Model diesel

######################################################################################################
# Step 3: Model petrol

######################################################################################################
# Step 4: Model natural gas