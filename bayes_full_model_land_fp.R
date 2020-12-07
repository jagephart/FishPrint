# Estimate land footprint from harvest and yield
# Land = harvest / yield

# Step 0: Run process_data_for_analysis.R, then clear environment other than:
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

# Get model-specific data:
# SELECT STUDY ID COLUMN - use this for rejoining outputs from multiple regression models back together
# Select relevant data columns and arrange by categorical info
# Clean harvest data (do all data cleaning like this up top) or move this to Functions.R
land_model_dat_categories <- lca_dat_clean_groups %>%
  select(study_id, yield = Yield_t_per_Ha, harvest_raw = Product, clean_sci_name, taxa, intensity = Intensity, system = Production_system_group) %>%
  arrange(clean_sci_name, taxa, intensity, system) %>%
  # First remove any harvest_raw entries that don't make sense
  mutate(harvest_raw = case_when(harvest_raw == "t tonne fresh" ~ "",
                                 harvest_raw != "" ~ harvest_raw)) %>% # insert NAs for missing data (case_when inserts NAs when no cases match)
  # convert "tonne" and "tonnes" to "t"
  mutate(harvest_raw = str_replace(harvest_raw, pattern = "\\btonne\\b|\\btonnes\\b", replacement = "t")) %>% # \\b indicate word boundaries
  # Separate tonnage and kg from state/presentation description
  separate(col = harvest_raw, into = c("weight", "state_presentation"), sep = "(?<=\\bt\\b)|(?<=\\bkg\\b)", remove = FALSE) %>%
  # Find entries that have to be dealt with manually
  #filter(is.na(weight)==FALSE & is.na(state_presentation))
  mutate(weight = case_when(str_detect(harvest_raw, "= 460 g") ~ "0.46 kg",
                            TRUE ~ weight)) %>%
  mutate(state_presentation = case_when(str_detect(harvest_raw, "= 460 g") ~ "live weight",
                                        TRUE ~ state_presentation)) %>%
  # Now separate column weight into weight and weight_unit
  separate(weight, into = c("weight", "weight_unit"), sep = " ", remove = TRUE) %>%
  # Create standardized column "weight" (in units tonnes)
  mutate(weight = as.numeric(weight, na.rm = TRUE)) %>%
  mutate(weight = case_when(weight_unit == "t" ~ weight,
                            weight_unit == "kg" ~ weight / 1000,
                            TRUE ~ NA_real_)) %>%
  # FIX IT - Create standardized unit "harvest" by multiplying weight by conversion factor?? does it matter?
  #mutate(harvest = case_when(str_detect(state_presentation, pattern = "fresh|live" ~ weight)))
  # For now, just use weight:
  rename(harvest = weight)

# Get footprint data
source("Functions.R")
fp_dat <- read.csv(file.path(datadir, "Feed_impact_factors_20201203.csv"))
fp_clean <- clean.feedFP(fp_dat)

######################################################################################################
# Step 1: Model yield
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
yield_brms_data <- data.frame(y = yield_no_na$yield, X_taxa_scaled, X_ordinal_scaled)

names(yield_brms_data)

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
fit_yield_no_na <- brm(yield_brms, data = yield_brms_data,
                     prior = all_priors, cores = 4, seed = "11729", iter = 20000, control = list(adapt_delta = 0.99))

# Get stan code
stancode(fit_yield_no_na)

######################################################################################################
# Use model to predict NAs for studies with complete set of predictors
# Both intensity AND system are non-NA
yield_complete_predictors <- land_model_dat_categories %>%
  filter(yield != 0 | is.na(yield))  %>% # First drop the zeroes; Have to explicitly include is.na(yield) otherwise NA's get dropped by yield != 0
  filter(is.na(intensity)==FALSE & is.na(system)== FALSE) %>%
  filter(is.na(yield)) # Now, filter to just the NAs

# PROBLEM: lca_complete predictors has more taxa than originally model:
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
levels(yield_complete_predictors$taxa) <- list(fresh_crust = "fresh_crust", misc_fresh = "misc_fresh", misc_marine = "misc_marine", oth_carp = "oth_carp", tilapia = "tilapia")

# Create NEW taxa model matrix for the studies whose feeds need to be predicted
# Taxa categories:
X_taxa_new <- model.matrix(object = ~ 1 + taxa, 
                           data = yield_complete_predictors %>% select(taxa)) 

# Center and Scale: BUT now center by the mean of the original modeled dataset above AND scale by the same 2*SD calculated from the original, modeled dataset above
X_taxa_new_scaled <- scale(X_taxa_new[,-1], center=apply(X_taxa[,-1], MARGIN = 2, FUN = mean), scale=2*taxa_sd)

# System and Intensity variables:
# Format intensity and system as ordinal variables, then center and scale
X_ordinal_new <- yield_complete_predictors %>%
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

# Bind
full_yield_dat <- yield_predictions %>%
  bind_rows(yield_no_na) %>%
  mutate(data_type = if_else(is.na(.point), true = "data", false = "prediction")) %>%
  arrange(taxa, intensity, system, clean_sci_name) %>%
  rownames_to_column() # Arrange by taxa first, then create dummy column for plotting 

# NEXT: use plot_brms_gamma_regression to produce figures for SI

######################################################################################################
# Step 2: Model harvest

harvest_no_na <- land_model_dat_categories %>%
  # Now filter to complete data (response and predictors)
  filter(is.na(intensity)==FALSE & is.na(system)== FALSE) %>%
  filter(is.na(harvest)==FALSE) %>% # Now, filter to just the NAs
  select(study_id, harvest, clean_sci_name, taxa, intensity, system)

# Set data for model:
# Create model matrix for taxa info, then center and scale
X_taxa <- model.matrix(object = ~ 1 + taxa, 
                       data = harvest_no_na %>% select(taxa)) 

taxa_sd <- apply(X_taxa[,-1], MARGIN=2, FUN=sd, na.rm=TRUE) # Center all non-intercept variables and scale by 2 standard deviations (ignoring NAs)
X_taxa_scaled <- scale(X_taxa[,-1], center=TRUE, scale=2*taxa_sd)

# Format intensity and system as ordinal variables, then center and scale
X_ordinal <- harvest_no_na %>%
  mutate(intensity = factor(intensity, levels = c("extensive", "semi", "intensive"))) %>% # set order of factors (low = extensive, high = intensive)
  mutate(system = factor(system, levels = c("open", "semi", "closed"))) %>% # set order of factors (low = open, high = closed)
  mutate(intensity = as.numeric(intensity)) %>%
  mutate(system = as.numeric(system)) %>%
  select(intensity, system) %>%
  as.matrix()
ordinal_sd<-apply(X_ordinal, MARGIN=2, FUN=sd, na.rm=TRUE) # Center all non-intercept variables and scale by 2 standard deviations (ignoring NAs)
X_ordinal_scaled <- scale(X_ordinal, center=TRUE, scale=2*ordinal_sd)

# Create dataframe for brms and rename feed variables
harvest_brms_data <- data.frame(y = harvest_no_na$harvest, X_taxa_scaled, X_ordinal_scaled)

names(harvest_brms_data)

# Set model formula
harvest_brms <- brmsformula(y ~ 1 + ., family = Gamma("log"))


# Use "resp = <response_variable>" to specify different priors for different response variables
all_priors <- c(set_prior("normal(0,5)", class = "b"), # priors for y response variables
                set_prior("normal(0,2.5)", class = "Intercept"), 
                set_prior("exponential(1)", class = "shape"))

# Model converges after increasing the adapt_delta and iterations from default values
# Rule of thumb: bulk and tail effective sample sizes should be 100 x number of chains (i.e., at least 400)
# increasing max_treedepth is more about efficiency (instead of validity)
# See: https://mc-stan.org/misc/warnings.html
fit_harvest_no_na <- brm(harvest_brms, data = harvest_brms_data,
                       prior = all_priors, cores = 4, seed = "11729", iter = 20000, control = list(adapt_delta = 0.99))

# Get stan code
stancode(fit_harvest_no_na)

######################################################################################################
# Use model to predict NAs for studies with complete set of predictors
# Both intensity AND system are non-NA
harvest_complete_predictors <- land_model_dat_categories %>%
  filter(is.na(intensity)==FALSE & is.na(system)== FALSE) %>%
  filter(is.na(harvest)) # Now, filter to just the NAs

# PROBLEM: lca_complete predictors has more taxa than originally model:
taxa_not_modeled <- setdiff(unique(harvest_complete_predictors$taxa), unique(harvest_no_na$taxa)) # these taxa were never modeled so they can't be predicted below

# DROP THESE FOR NOW:
harvest_complete_predictors <- harvest_complete_predictors %>%
  filter(taxa %in% taxa_not_modeled == FALSE)

# Now check the other way, which taxa were in the original model but not a part of the data that needs to be predicted:
setdiff(unique(harvest_no_na$taxa), unique(harvest_complete_predictors$taxa))

# If original model has taxa that are not part of harvest_complete_predictors, 
# Use sorted list of unique taxa in original model and use this to expand/assign levels manually - having trouble automating this
# Include all levels here, but remember the first level won't show up in design matrix - instead, it's part of the "contrasts"
sort(unique(harvest_no_na$taxa))
harvest_complete_predictors <- harvest_complete_predictors %>%
  mutate(taxa = as.factor(taxa))
levels(harvest_complete_predictors$taxa) <- list(com_carp = "com_carp", milkfish = "milkfish", misc_diad = "misc_diad", misc_fresh = "misc_fresh", misc_marine = "misc_marine", 
                                                 mussel = "mussel", oth_carp = "oth_carp", salmon = "salmon", shrimp = "shrimp", tilapia = "tilapia", trout = "trout")

# Create NEW taxa model matrix for the studies whose feeds need to be predicted
# Taxa categories:
X_taxa_new <- model.matrix(object = ~ 1 + taxa, 
                           data = harvest_complete_predictors %>% select(taxa)) 

# Center and Scale: BUT now center by the mean of the original modeled dataset above AND scale by the same 2*SD calculated from the original, modeled dataset above
X_taxa_new_scaled <- scale(X_taxa_new[,-1], center=apply(X_taxa[,-1], MARGIN = 2, FUN = mean), scale=2*taxa_sd)

# System and Intensity variables:
# Format intensity and system as ordinal variables, then center and scale
X_ordinal_new <- harvest_complete_predictors %>%
  mutate(intensity = factor(intensity, levels = c("extensive", "semi", "intensive"))) %>% # set order of factors (low = extensive, high = intensive)
  mutate(system = factor(system, levels = c("open", "semi", "closed"))) %>% # set order of factors (low = open, high = closed)
  mutate(intensity = as.numeric(intensity)) %>%
  mutate(system = as.numeric(system)) %>%
  select(intensity, system) %>%
  as.matrix()

# Center and Scale: BUT now center by the mean of the original modeled dataset above AND scale by the same 2*SD calculated from the original, modeled dataset above
X_ordinal_new_scaled <- scale(X_ordinal_new, center=apply(X_ordinal, MARGIN = 2, FUN = mean), scale=2*ordinal_sd)

# Create dataframe for brms and rename feed variables
brms_new_harvest_dat <- data.frame(cbind(X_taxa_new_scaled, X_ordinal_new_scaled)) 

# Make predictions
#predicted_harvest_dat <- predict(fit_no_na, newdata = brms_new_harvest_data)
# Use tidybayes instead:
predicted_harvest_dat <- add_predicted_draws(newdata = brms_new_harvest_dat, model = fit_harvest_no_na)

# Get point and interval estimates from predicted data
# Select just the prediction columns
# Join these with the modeled data (harvest_complete_predictors) to get metadata on taxa/intensity/syste,
harvest_dat_intervals <- predicted_harvest_dat %>%
  median_qi(.value = .prediction) %>% # Rename prediction to value
  ungroup() %>%
  select(contains("."))

# .row is equivalent to the row number in the modeled dataset (harvest_complete_predictors) - create a join column for this
harvest_metadat<- harvest_complete_predictors %>%
  select(study_id, clean_sci_name, taxa, intensity, system) %>%
  mutate(.row = row_number())

harvest_predictions <- harvest_dat_intervals %>%
  left_join(harvest_metadat, by = ".row") %>%
  rename(harvest = .value)

# Bind
full_harvest_dat <- harvest_predictions %>%
  bind_rows(harvest_no_na) %>%
  mutate(data_type = if_else(is.na(.point), true = "data", false = "prediction")) %>%
  arrange(taxa, intensity, system, clean_sci_name) %>%
  rownames_to_column() # Arrange by taxa first, then create dummy column for plotting 

# NEXT: use plot_brms_gamma_regression to produce figures for SI

######################################################################################################
# STEP 3 - aggregate up to taxa level and estimate total feed footprint
# Set data for model

# Merge datasets by study_id, remove interval info, only use medians for rest of analysis
# Format dat for merge
# Rename data_type to keep track of different datasets
yield_dat_merge <- full_yield_dat %>%
  select(study_id, clean_sci_name, taxa, intensity, system, yield_data_type = data_type, yield)

harvest_dat_merge <- full_harvest_dat %>%
  select(study_id, clean_sci_name, taxa, intensity, system, harvest_data_type = data_type, harvest)


# Merge as full_join and remove NAs
# FIX IT - check dim of full_yield_dat and 
land_footprint_dat <- yield_dat_merge %>%
  full_join(harvest_dat_merge, by = intersect(names(yield_dat_merge), names(harvest_dat_merge))) %>%
  drop_na() %>%
  arrange(clean_sci_name, taxa) %>%
  mutate(clean_sci_name = as.factor(clean_sci_name),
       sci = as.numeric(clean_sci_name),
       taxa = as.factor(taxa),
       tx = as.numeric(taxa))

# Set data
N = nrow(land_footprint_dat)
N_SCI <- length(unique(land_footprint_dat$sci))
n_to_sci <- land_footprint_dat$sci
N_TX <- length(unique(land_footprint_dat$tx))
sci_to_tx = land_footprint_dat %>%
  select(sci, tx) %>%
  unique() %>%
  pull(tx)
harvest <- land_footprint_dat$harvest
yield <- land_footprint_dat$yield

stan_data <- list(N = N,
                  N_SCI = N_SCI, 
                  n_to_sci = n_to_sci,
                  N_TX = N_TX,
                  sci_to_tx = sci_to_tx,
                  harvest = harvest,
                  yield = yield)

# Estimate foot print for all scientific names and taxa groups (removed the "all-seafood" level for simplicity)
# GAMMA distribution hierarchical model
stan_no_na <- 'data {
  // data for gamma model for FCR
  int<lower=0> N;  // number of observations
  vector<lower=0>[N] harvest; // data
  vector<lower=0>[N] yield; // data
  int N_TX; // number of taxa groups
  int N_SCI; // number of scientific names
  int n_to_sci[N]; // sciname index
  int sci_to_tx[N_SCI]; // taxa-group indices
}
parameters {
  vector<lower=0>[N_TX] tx_mu;
  vector<lower=0>[N_SCI] sci_mu;
  // only need sigmas if defining shape and rate with mu and sigma
  real<lower=0> tx_sigma;
  real<lower=0> sci_sigma;
  // if using variance instead of st dev
  //real<lower=0> tx_sigma_sq;
  //real<lower=0> sci_sigma_sq;
}
transformed parameters {
  // define transofrmed params for gamma model for FCRs
  vector<lower=0>[N_SCI] sci_shape;
  vector<lower=0>[N_SCI] sci_rate;
  vector<lower=0>[N_TX] tx_shape;
  vector<lower=0>[N_TX] tx_rate;

  // gamma model reparameterization
  // option 1: reparamaterize gamma to get mu and sigma; defining these here instead of the model section allows us to see these parameters in the output
  // taxa group level
  for (n_tx in 1:N_TX){
    tx_shape[n_tx] = square(tx_mu[n_tx]) ./ square(tx_sigma);
    tx_rate[n_tx] = tx_mu[n_tx] ./ square(tx_sigma);
  }

  // sci level
  for (n_sci in 1:N_SCI){
    sci_shape[n_sci] = square(sci_mu[n_sci]) ./ square(sci_sigma);
    sci_rate[n_sci] = sci_mu[n_sci] ./ square(sci_sigma);
  }

  // option 2: reparameterize gamma to get just mu
  //for (n_tx in 1:N_TX){
  //  tx_rate[n_tx] = tx_shape[n_tx] ./ tx_mu[n_tx];
  //}
  // sci level
  //for (n_sci in 1:N_SCI){
  //  sci_rate[n_sci] = sci_shape[n_sci] ./ sci_mu[n_sci];
  //}
  
  // option 3: reparameterize shape and rate as inverse(va) and inverse(va)/mu
  //for (n_tx in 1:N_TX){
  //  tx_shape[n_tx] = 1 / tx_sigma_sq;
  //  tx_rate[n_tx] = 1 / (tx_sigma_sq * tx_mu[n_tx]);
  //}
  //for (n_sci in 1:N_SCI){
  //  sci_shape[n_sci] = 1 / sci_sigma_sq;
  //  sci_rate[n_sci] = 1 / (sci_sigma_sq * sci_mu[n_sci]);
  //}
}
model {
  // define priors for gamma model
  // Put priors on mu and sigma (instead of shape and rate) since this is more intuitive:
  tx_mu ~ uniform(0, 100); // note: uniform(0,100) for all of these doesnt help much with convergence
  sci_mu ~ uniform(0, 100);
  //tx_sigma ~ uniform(0, 100); // only need sigmas if calculating shape and rate with mu and sigma
  //sci_sigma ~ uniform(0, 100);

  // likelihood
  // gamma model sci-name and taxa-level
  for (n in 1:N){
    harvest[n]/yield[n] ~ gamma(sci_shape[n_to_sci[n]], sci_rate[n_to_sci[n]]);
  }

  for (n_sci in 1:N_SCI){
    sci_mu[n_sci] ~ gamma(tx_shape[sci_to_tx[n_sci]], tx_rate[sci_to_tx[n_sci]]);
  }
}'


no_na_mod <- stan_model(model_code = stan_no_na)
# Note: For Windows, apparently OK to ignore this warning message:
# Warning message:
#   In system(paste(CXX, ARGS), ignore.stdout = TRUE, ignore.stderr = TRUE) :
#   'C:/rtools40/usr/mingw_/bin/g++' not found

# Fit model:
# Set seed while testing
#fit_no_na <- sampling(object = no_na_mod, data = stan_data, cores = 4, seed = "11729", iter = 10000, control = list(adapt_delta = 0.99))
fit_no_na <- sampling(object = no_na_mod, data = stan_data, cores = 4, iter = 10000, control = list(adapt_delta = 0.999))
summary(fit_no_na)

launch_shinystan(fit_no_na)
