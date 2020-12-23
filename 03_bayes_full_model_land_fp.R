# Impute on-farm LAND impacts then combine with off-farm (feed) LAND impacts 

# Calculate on-farm (non-feed associated) land impact as: 1 / yield

##################################### REMEMBER TO ONLY DO THIS FOR SYSTEM == PONDS or RECIRCULATING and TANKS

# Step 0: Run process_data_for_analysis.R, then clear environment other than:
rm(list=ls()[!(ls() %in% c("lca_dat_clean_groups", "datadir", "outdir"))])

# MOVE THIS TO process_data_for_analysis.R:
# Standardize yield columns as tonnes per Ha
# 1000 kg in 1 tonne
# 10,000 sq meters in 1 hectare

# Check that no entries report both units:
# lca_dat_clean %>%
#   filter(is.na(Yield_t_per_Ha)==FALSE & is.na(Yield_kg_per_m3)==FALSE)

# Get model-specific data:
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
                     prior = all_priors, cores = 4, seed = "11729", iter = 20000, control = list(adapt_delta = 0.99))

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
                       prior = all_priors, cores = 4, seed = "11729", iter = 20000, control = list(adapt_delta = 0.99))

# Get stan code
#stancode(fit_yield_no_taxa)

######################################################################################################
# Use model of intensity + system to predict yield for taxa that were not part of the previous model
# Intensity must be non-NA
yield_complete_predictors <- land_model_dat_categories %>%
  filter(is.na(intensity)==FALSE & is.na(system)==FALSE) %>%
  filter(is.na(yield)) # Now, filter to just the NAs

# PROBLEM: lca_complete predictors has more taxa than originally model:
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

# NEXT: use plot_brms_gamma_regression to produce figures for SI


# Clear all memory except for final stan model:
rm(list=ls()[!(ls() %in% c("datadir", "outdir", 
                           "lca_dat_clean_groups",
                           "land_model_dat_categories",
                           "full_yield_dat"))])

#save.image(file.path(outdir, paste(Sys.Date(), "_on-farm-land-all-data-prior-to-aggregation.RData", sep = "")))

######################################################################################################
# RESTARTING POINT
# Load on-farm land variables

# datadir <- "/Volumes/jgephart/BFA Environment 2/Data"
# outdir <- "/Volumes/jgephart/BFA Environment 2/Outputs"
# load(file.path(outdir, "2020-12-21_on-farm-land-all-data-prior-to-aggregation.RData"))

######################################################################################################
# Estimate taxa-level footprint

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

# Format on-farm land data
# Calculate n_in_sci and n_in_taxa in case we want to remove small sample sizes
land_footprint_merge <- full_yield_dat %>%
  select(study_id, clean_sci_name, taxa, intensity, system, data_type, yield) %>%
  drop_na() %>% # Make sure to drop na before creating sci and tx indices (otherwise some indices drop out)
  group_by(clean_sci_name) %>%  
  mutate(n_in_sci = n()) %>%
  ungroup() %>%
  group_by(taxa) %>%
  mutate(n_in_taxa = n()) %>%
  ungroup() 

# MERGE on and off-farm data
# Then MERGE with complete lca_data_clean_groups to get all studies excluded from on-farm land analysis
land_footprint_dat <- feed_dat_merge %>%
  full_join(fcr_dat_merge, by = intersect(names(feed_dat_merge), names(fcr_dat_merge))) %>%
  full_join(land_footprint_merge, by = intersect(names(land_footprint_merge), names(.))) %>%
  # FOR HIERARCHICAL LAND MODEL, NA's are actually zeroes (all Cages & Pens left out of LAND REGRESSION MODEL should be assigned to zero here)
  # But adjust zeroes upward to allow gamma model to work
  mutate(yield = if_else(is.na(yield), true = 0.1, false = yield)) %>%
  arrange(clean_sci_name, taxa) %>%
  mutate(clean_sci_name = as.factor(clean_sci_name),
         sci = as.numeric(clean_sci_name),
         taxa = as.factor(taxa),
         tx = as.numeric(taxa)) 

# OUTPUT data used for analysis:
write.csv(land_footprint_dat, file.path(outdir, "dat_for_bayesian-land.csv"), row.names = FALSE)

# FEED IMPACT CONSTANTS:
# Choose allocation method
# Choose Impact.category
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


# Set data
N = nrow(land_footprint_dat)
N_SCI <- length(unique(land_footprint_dat$sci))
n_to_sci <- land_footprint_dat$sci
n_to_tx <- land_footprint_dat$tx
N_TX <- length(unique(land_footprint_dat$tx))
sci_to_tx <- land_footprint_dat %>%
  select(sci, tx) %>%
  unique() %>%
  pull(tx)
yield <- land_footprint_dat$yield

stan_data <- list(N = N,
                  N_SCI = N_SCI, 
                  n_to_sci = n_to_sci,
                  N_TX = N_TX,
                  #n_to_tx = n_to_tx,
                  sci_to_tx = sci_to_tx,
                  yield = yield)

# Estimate foot print for all scientific names and taxa groups (removed the "all-seafood" level for simplicity)
# GAMMA distribution hierarchical model
# stan_no_na <- 'data {
#   // data for gamma model for FCR
#   int<lower=0> N;  // number of observations
#   vector<lower=0>[N] yield; // data
#   int N_TX; // number of taxa groups
#   int N_SCI; // number of scientific names
#   int n_to_sci[N]; // sciname index
#   int sci_to_tx[N_SCI]; // taxa-group indices
# }
# parameters {
#   vector<lower=0>[N_TX] tx_mu;
#   vector<lower=0>[N_SCI] sci_mu;
#   real<lower=0> tx_sigma; // only need to define sigmas if using option 1
#   real<lower=0> sci_sigma; 
#   //vector<lower=0>[N_SCI] sci_shape; // define shape here if not transforming (using option 2 below)
# }
# transformed parameters {
#   // define transofrmed params for gamma model for FCRs
#   vector<lower=0>[N_SCI] sci_shape;
#   vector<lower=0>[N_SCI] sci_rate;
#   vector<lower=0>[N_TX] tx_shape;
#   vector<lower=0>[N_TX] tx_rate;
# 
#   // reparamaterize gamma to get mu and sigma; defining these here instead of the model section allows us to see these parameters in the output
#   // option 1: mean and variance
#   for (n_tx in 1:N_TX){
#     tx_shape[n_tx] = square(tx_mu[n_tx]) / square(tx_sigma);
#     tx_rate[n_tx] = tx_mu[n_tx] / square(tx_sigma);
#   }
#   for (n_sci in 1:N_SCI){
#     sci_shape[n_sci] = square(sci_mu[n_sci]) / square(sci_sigma);
#     sci_rate[n_sci] = sci_mu[n_sci] / square(sci_sigma);
#   }
#   
#   // option 2: rate = shape / mean
#   //for (n_tx in 1:N_TX){
#   //  tx_rate[n_tx] = tx_shape[n_tx] / tx_mu[n_tx];
#   //}
#   //for (n_sci in 1:N_SCI){
#   //  sci_rate[n_sci] = sci_shape[n_sci] / sci_mu[n_sci];
#   //}
#   
# }
# model {
#   // define priors for gamma model
#   // Put priors on mu and sigma (instead of shape and rate) since this is more intuitive:
#   //tx_mu ~ uniform(0, 100); // note: uniform(0,100) for all of these doesnt help much with convergence
#   //sci_mu ~ uniform(0, 100);
#   //tx_sigma ~ uniform(0, 100); // only need sigmas if calculating shape and rate with mu and sigma
#   //sci_sigma ~ uniform(0, 100);
#   //tx_sigma ~ cauchy(.05, .01);
#   // try adding weakly informative priors on shape and/or rate parameters
#   // shape priors require target += notation: values come outputs of shape parameters in non-hierarchical model - did not help
#   //target += cauchy_lpdf(sci_shape[1] | 16, 8);
#   //target += cauchy_lpdf(sci_shape[2] | 200, 100);
#   //target += cauchy_lpdf(sci_shape[3] | 8, 4);
#   //target += cauchy_lpdf(sci_shape[4] | 0.3, 0.15);
#   //target += cauchy_lpdf(sci_shape[5] | 7, 3.5);
#   //target += cauchy_lpdf(sci_shape[6] | 0.4, 0.2);
#   //target += cauchy_lpdf(sci_shape[7] | 0.4, 0.2);
#   //target += cauchy_lpdf(sci_shape[8] | 14, 7);
#   //target += cauchy_lpdf(sci_shape[9] | 1, 0.5);
#   //target += cauchy_lpdf(sci_shape[10] | 1, 0.5);
#   
#   // FIX IT - create generated quantities section for land = 1 / yield
#   // likelihood
#   // gamma model sci-name and taxa-level
#   for (n in 1:N){
#     yield[n] ~ gamma(sci_shape[n_to_sci[n]], sci_rate[n_to_sci[n]]);
#   }
#   for (n_sci in 1:N_SCI){
#     sci_mu[n_sci] ~ gamma(tx_shape[sci_to_tx[n_sci]], tx_rate[sci_to_tx[n_sci]]);
#   }
# }'

# NORMAL distribution hierarchical model
stan_no_na <- 'data {
  // data for gamma model for FCR
  int<lower=0> N;  // number of observations
  vector<lower=0>[N] yield; // data
  int N_TX; // number of taxa groups
  int N_SCI; // number of scientific names
  int n_to_sci[N]; // sciname index
  int sci_to_tx[N_SCI]; // taxa-group indices
}
parameters {
  vector<lower=0>[N_TX] tx_mu;
  vector<lower=0>[N_SCI] sci_mu;
  //vector<lower=0>[N_TX] tx_sigma;
  //vector<lower=0>[N_SCI] sci_sigma; 
  real<lower=0> tx_sigma;
  real<lower=0> sci_sigma; 
}
model {
  // define priors for gamma model
  // Put priors on mu and sigma (instead of shape and rate) since this is more intuitive:
  //tx_mu ~ uniform(0, 100); // note: uniform(0,100) for all of these doesnt help much with convergence
  //sci_mu ~ uniform(0, 100);
  //tx_sigma ~ uniform(0, 100); 
  //sci_sigma ~ uniform(0, 100);
  //tx_sigma ~ cauchy(.05, .01);

  
  // likelihood
  // normal model sci-name and taxa-level
  for (n in 1:N){
    yield[n] ~ normal(sci_mu[n_to_sci[n]], sci_sigma);
  }

  for (n_sci in 1:N_SCI){
    sci_mu[n_sci] ~ normal(tx_mu[sci_to_tx[n_sci]], tx_sigma);
  }
}
generated quantities {
  // Land = 1 / Yield - model converges better if yield ~ normal() instead of 1/yield ~ normal()
  vector<lower=0>[N_TX] tx_land;
  vector<lower=0>[N_SCI] sci_land;
  
  // Calculations
  for (n_tx in 1:N_TX) {
    tx_land[n_tx] = 1 / tx_mu[n_tx];
  }
  for (n_sci in 1:N_SCI) {
    sci_land[n_sci] = 1 / sci_mu[n_sci];
  }
  
}'

no_na_mod <- stan_model(model_code = stan_no_na)
# Note: For Windows, apparently OK to ignore this warning message:
# Warning message:
#   In system(paste(CXX, ARGS), ignore.stdout = TRUE, ignore.stderr = TRUE) :
#   'C:/rtools40/usr/mingw_/bin/g++' not found

# Fit model:
# Set seed while testing
fit_no_na <- sampling(object = no_na_mod, data = stan_data, cores = 4, seed = "11729", iter = 10000, control = list(adapt_delta = 0.99))
#fit_no_na <- sampling(object = no_na_mod, data = stan_data, cores = 4, iter = 10000, control = list(adapt_delta = 0.99))
summary(fit_no_na)$summary

#launch_shinystan(fit_no_na)

######################################################################################################
# PLOT RESULTS
plot_theme <- theme(title = element_text(size = 18),
                    axis.title.x = element_text(size = 16),
                    axis.text=element_text(size=14, color = "black"))

# Key for naming sci and taxa levels
# Get full taxa group names back
index_key <- land_footprint_dat %>%
  select(clean_sci_name, sci, taxa, tx, n_in_sci, n_in_taxa) %>%
  unique() %>%
  mutate(taxa = as.character(taxa),
         full_taxa_name = case_when(taxa == "fresh_crust" ~ "freshwater crustaceans",
                                    taxa == "misc_diad" ~ "misc diadromous fishes",
                                    taxa == "misc_fresh" ~ "misc freshwater fishes",
                                    taxa == "misc_marine" ~ "misc marine fishes",
                                    TRUE ~ taxa),
         taxa = as.factor(taxa),
         full_taxa_name = as.factor(full_taxa_name))

# Use tidybayes + ggdist for finer control of aes mapping (instead of bayesplots) 
get_variables(fit_no_na)

# Sci-level land footprints as point intervals:
fit_no_na %>%
  spread_draws(sci_land[sci]) %>%
  median_qi(.width = 0.8) %>%
  left_join(index_key, by = "sci") %>% # Join with index key to get sci and taxa names
  mutate(clean_sci_name = fct_reorder(clean_sci_name, as.character(full_taxa_name))) %>%
  ggplot(aes(y = clean_sci_name, x = sci_land, xmin = .lower, xmax = .upper, color = full_taxa_name)) +
  geom_pointinterval() +
  theme_classic() + 
  plot_theme + 
  labs(x = "hectares per tonne", y = "", title = "Land", color = "taxa group")
ggsave(filename = file.path(outdir, "plot_land-footprint_sci-level.png"), width = 11, height = 8.5)

# Taxa-level land footprints as densities: 
# FIX IT - DENSITIES ARE SUPER STEEP (see just trout for example) - if this is true for final dataset, just plot as point-interval
fit_no_na %>%
  spread_draws(tx_land[tx]) %>%
  #median_qi(.width = 0.8) %>% # need at least 2 points to select a bandwidth automatically
  left_join(index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  #filter(tx_land > 100)
  filter(full_taxa_name == "freshwater crustaceans") %>%
  ggplot(aes(y = full_taxa_name, x = tx_land)) +
  stat_halfeye(aes(slab_fill = full_taxa_name)) +
  theme_classic() + 
  plot_theme + 
  theme(legend.position = "none") +
  labs(x = "hectares per tonne", y = "", title = "Land")
ggsave(filename = file.path(outdir, "plot_land-footprint_taxa-level.png"), width = 11, height = 8.5)


# OPTION 2: Taxa-level plots with color themes:
# If we want to mimic bayesplot color schemes, can get hexadecimal colors and input manually to stat_halfeye aesthetics
color_scheme_get("blue")
color_scheme_get("green")

