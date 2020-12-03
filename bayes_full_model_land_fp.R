# Use Bayes gamma regression to estimate harvest and yield
# Calculate land footprint as harvest / yield

# Step 0: If needed, run process_data_for_analysis.R, then clear environment other than:
rm(list=ls()[!(ls() %in% c("lca_dat_clean_groups", "datadir", "outdir"))])

# MOVE THIS TO process_data_for_analysis.R:
# Standardize yield columns as tonnes per Ha
# 1000 kg in 1 tonne
# 10,000 sq meters in 1 hectare

# Check that no entries report both units:
# lca_dat_clean %>%
#   filter(is.na(Yield_t_per_Ha)==FALSE & is.na(Yield_kg_per_m3)==FALSE)

# FIX IT: volume to area not possible without assumptions
# lca_dat_clean %>%
#   mutate(Yield_t_per_Ha == if_else(is.na(Yield_kg_per_m3)==FALSE, true = ))

# Get model-specific data
# ie, remove all NAs and model these data, then use the model to predict yield for those data with a complete set of predictors
# CREATE STUDY ID COLUMN - use this for rejoining outputs from multiple regression models back together
# Select relevant data columns and arrange by categorical info

# FIX IT - ADD "HARVEST" column to this:

land_model_dat_categories <- lca_dat_clean_groups %>%
  select(yield = Yield_t_per_Ha, clean_sci_name, taxa, intensity = Intensity, system = Production_system_group) %>%
  arrange(clean_sci_name, taxa, intensity, system) %>%
  mutate(study_id = row_number())

######################################################################################################
# Step 1: Model yield

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

# Format intensity and system as ordinal variables, then center and scale
X_ordinal <- yield_no_na %>%
  mutate(intensity = factor(intensity, levels = c("extensive", "semi", "intensive"))) %>% # set order of factors (low = extensive, high = intensive)
  mutate(system = factor(system, levels = c("open", "semi", "closed"))) %>% # set order of factors (low = open, high = closed)
  mutate(intensity = as.numeric(intensity)) %>%
  mutate(system = as.numeric(system)) %>%
  select(intensity, system) %>%
  as.matrix()
ordinal_sd<-apply(X_ordinal, MARGIN=2, FUN=sd, na.rm=TRUE) # Center all non-intercept variables and scale by 2 standard deviations (ignoring NAs)
X_ordinal_scaled <- scale(X_ordinal, center=TRUE, scale=2*ordinal_sd)

# Create dataframe for brms and rename feed variables
gamma_brms_data <- data.frame(y = yield_no_na$yield, X_taxa_scaled, X_ordinal_scaled)

names(gamma_brms_data)

# Remove NA columns (predictors with no variation)
gamma_brms_data <- gamma_brms_data[ ,colSums(is.na(gamma_brms_data)) == 0]

# Set model formula
# FIX IT - for now removing system since there's no variation in available data
yield_brms <- brmsformula(y ~ 1 + taxamisc_fresh + taxamisc_marine + taxaoth_carp + taxatilapia + 
                        intensity, family = Gamma("log"))

# Versus:
#yield_brms_shortcut <- brmsformula(y ~ 1 + ., family = Gamma("log"))


# Use "resp = <response_variable>" to specify different priors for different response variables
all_priors <- c(set_prior("normal(0,5)", class = "b"), # priors for y response variables
                set_prior("normal(0,2.5)", class = "Intercept"), 
                set_prior("exponential(1)", class = "shape"))

# Model converges after increasing the adapt_delta and iterations from default values
# Rule of thumb: bulk and tail effective sample sizes should be 100 x number of chains (i.e., at least 400)
# increasing max_treedepth is more about efficiency (instead of validity)
# See: https://mc-stan.org/misc/warnings.html
fit_yield_no_na <- brm(yield_brms, data = gamma_brms_data,
                     prior = all_priors, cores = 4, seed = "11729", iter = 20000, control = list(adapt_delta = 0.99))

# Versus:
# fit_yield_no_na_shortcut <- brm(yield_brms, data = gamma_brms_data,
#                        prior = all_priors, cores = 4, seed = "11729", iter = 20000, control = list(adapt_delta = 0.99))

# Get stan code
stancode(fit_yield_no_na)

######################################################################################################
# Use model to predict NAs for studies with complete set of predictors
# Both intensity AND system are non-NA
land_complete_predictors <- land_model_dat_categories %>%
  filter(yield != 0 | is.na(yield))  %>% # Have to explicitly include is.na(yield) otherwise NA's get dropped by yield != 0
  filter(is.na(intensity)==FALSE & is.na(system)== FALSE) %>%
  filter(is.na(yield)) # filter to just the NAs

# PROBLEM: lca_complete predictors has more taxa than originally model:
taxa_not_modeled <- setdiff(unique(land_complete_predictors$taxa), unique(yield_no_na$taxa)) # these taxa were never modeled so they can't be predicted below

# DROP THESE FOR NOW:
land_complete_predictors <- land_complete_predictors %>%
  filter(taxa %in% taxa_not_modeled == FALSE)

# Now check the other way, which taxa were in the original model but not a part of the data that needs to be predicted:
setdiff(unique(yield_no_na$taxa), unique(land_complete_predictors$taxa))

# If original model has taxa that are not part of land_complete_predictors, 
# Use list of unique taxa in original model and use this to expand/assign levels manually - having trouble automating this
# Include all levels here, but remember the first level won't show up in design matrix - instead, it's part of the "contrasts"
sort(unique(yield_no_na$taxa))
land_complete_predictors <- land_complete_predictors %>%
  mutate(taxa = as.factor(taxa))
levels(land_complete_predictors$taxa) <- list(fresh_crust = "fresh_crust", misc_fresh = "misc_fresh", misc_marine = "misc_marine", oth_carp = "oth_carp", tilapia = "tilapia")

# Create NEW taxa model matrix for the studies whose feeds need to be predicted
# Taxa categories:
X_taxa_new <- model.matrix(object = ~ 1 + taxa, 
                           data = land_complete_predictors %>% select(taxa)) 

# Center and Scale: BUT now center by the mean of the original modeled dataset above AND scale by the same 2*SD calculated from the original, modeled dataset above
X_taxa_new_scaled <- scale(X_taxa_new[,-1], center=apply(X_taxa[,-1], MARGIN = 2, FUN = mean), scale=2*taxa_sd)

# System and Intensity variables:
# Format intensity and system as ordinal variables, then center and scale
X_ordinal_new <- land_complete_predictors %>%
  mutate(intensity = factor(intensity, levels = c("extensive", "semi", "intensive"))) %>% # set order of factors (low = extensive, high = intensive)
  mutate(system = factor(system, levels = c("open", "semi", "closed"))) %>% # set order of factors (low = open, high = closed)
  mutate(intensity = as.numeric(intensity)) %>%
  mutate(system = as.numeric(system)) %>%
  select(intensity, system) %>%
  as.matrix()

# Center and Scale: BUT now center by the mean of the original modeled dataset above AND scale by the same 2*SD calculated from the original, modeled dataset above
X_ordinal_new_scaled <- scale(X_ordinal_new, center=apply(X_ordinal, MARGIN = 2, FUN = mean), scale=2*ordinal_sd)

# Create dataframe for brms and rename feed variables
brms_new_yield_dat <- data.frame(cbind(X_taxa_new_scaled, X_ordinal_new_scaled)) 

# FIX IT - system was not part of the original model (no variation in available data) so removing this column for now:
brms_new_yield_dat <- brms_new_yield_dat %>%
  select(-system)

# Make predictions
#predicted_yield_dat <- predict(fit_no_na, newdata = brms_new_yield_data)
# Use tidybayes instead:
predicted_yield_dat <- add_predicted_draws(newdata = brms_new_yield_dat, model = fit_yield_no_na)

# Get point and interval estimates from predicted data
# Select just the prediction columns
# Join these with the modeled data (land_complete_predictors) to get metadata on taxa/intensity/syste,
yield_dat_intervals <- predicted_yield_dat %>%
  median_qi(.value = .prediction) %>% # Rename prediction to value
  ungroup() %>%
  select(contains("."))

# .row is equivalent to the row number in the modeled dataset (land_complete_predictors) - create a join column for this
yield_metadat<- land_complete_predictors %>%
  select(study_id, clean_sci_name, taxa, intensity, system) %>%
  mutate(.row = row_number())

yield_predictions <- yield_dat_intervals %>%
  left_join(yield_metadat, by = ".row") %>%
  rename(yield = .value)

# Bind
full_yield_dat <- yield_predictions %>%
  bind_rows(yield_no_na) %>%
  mutate(data_type = if_else(is.na(.point), true = "data", false = "prediction")) %>%
  arrange(taxa, intensity, system, clean_sci_name) %>%
  rownames_to_column() # Arrange by taxa first, then create dummy column for plotting 

# NEXT: use plot_brms_gamma_regression to produce figures for SI

######################################################################################################
# Step 3: Model harvest

harvest_no_na <- lca_model_dat_categories

