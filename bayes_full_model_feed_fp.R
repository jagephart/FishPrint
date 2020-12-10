# Estiamte on-farm, feed-associated footprint
# Reminder for new data: entries with data for some columns should have no blanks in other columns (these should be filled in as zeroes)

# First run: process_data_for_analysis.R, then clear environment other than:
rm(list=ls()[!(ls() %in% c("lca_dat_clean_groups", "datadir", "outdir"))])


######################################################################################################
# Step 1: Model FCR, imputing NAs in the predictors (with just the intercept - equivalent to setting this to the mean), and predict missing FCR values

# Remove species that are not fed (FCR == 0), format intensity and system variables as ordinal
# lca_with_na <- feed_model_dat_categories %>%
#   filter(FCR != 0 | is.na(FCR))  %>% # Have to explicitly include is.na(FCR) otherwise NA's get dropped by FCR != 0
#   arrange(clean_sci_name, taxa, intensity, system)
# 
# # Create model matrix from taxa data, then center and scale
# options(na.action='na.pass') # If needed, change default options for handling missing data 
# X_taxa <- model.matrix(object = ~taxa, 
#                        data = lca_with_na %>% select(taxa)) 
# options(na.action='na.omit') # Return option back to the default
# taxa_sd <- apply(X_taxa[,-1], MARGIN=2, FUN=sd, na.rm=TRUE) # Calculate standard deviation across all vars
# # Center all non-intercept variables and scale by 2 standard deviations (ignoring NAs)
# options(na.action='na.pass') # If needed, change default options for handling missing data
# X_taxa_scaled <- scale(X_taxa[,-1], center=TRUE, scale=2*taxa_sd)
# options(na.action='na.omit') # Return option back to the default
# 
# # Format intensity and system as ordinal variables, then center and scale
# X_ordinal <- lca_with_na %>%
#   mutate(intensity = factor(intensity, levels = c("extensive", "semi", "intensive"))) %>% # set order of factors (low = extensive, high = intensive)
#   mutate(system = factor(system, levels = c("open", "semi", "closed"))) %>% # set order of factors (low = open, high = closed)
#   mutate(intensity = as.numeric(intensity)) %>%
#   mutate(system = as.numeric(system)) %>%
#   select(intensity, system) %>%
#   as.matrix()
# ordinal_sd<-apply(X_ordinal, MARGIN=2, FUN=sd, na.rm=TRUE) # Calculate standard deviation across all vars
# # Center all non-intercept variables and scale by 2 standard deviations (ignoring NAs)
# options(na.action='na.pass') # First change default options for handling missing data
# X_ordinal_scaled <- scale(X_ordinal, center=TRUE, scale=2*ordinal_sd)
# options(na.action='na.omit') # Return option back to the default
# 
# # Create dataframe for brms and rename feed variables
# gamma_brms_data <- data.frame(y = lca_with_na$FCR, X_taxa_scaled, X_ordinal_scaled)
# 
# # Which predictors have missing data:
# X_where_na <- apply(gamma_brms_data, MARGIN = 2, is.na)
# colSums(X_where_na)
# 
# # Set model formulas
# y_brms <- brmsformula(y | mi() ~ 1 + taxacrab + taxafresh_crust + taxamilkfish + taxamisc_diad + taxamisc_fresh +
#                         taxamisc_marine + taxaoth_carp + taxasalmon + taxashrimp + taxatilapia + taxatrout + taxatuna +
#                         mi(intensity) + mi(system), family = Gamma("log"))
# 
# intensity_mi <- brmsformula(intensity | mi() ~ 1,
#                             family = gaussian())
# 
# system_mi  <- brmsformula(system | mi() ~ 1,
#                           family = gaussian())
# 
# # Use "resp = <response_variable>" to specify different priors for different response variables
# all_priors <- c(set_prior("normal(0,5)", class = "b", resp = "y"), # priors for y response variables
#                 set_prior("normal(0,2.5)", class = "Intercept", resp = "y"), 
#                 set_prior("exponential(1)", class = "shape", resp = "y"),
#                 set_prior("normal(0,1.2)", class = "Intercept", resp = c("intensity", "system")),
#                 set_prior("exponential(2)", class = "sigma", resp = c("intensity", "system")))
# 
# # Model converges after increasing the adapt_delta and iterations from default values
# # Rule of thumb: bulk and tail effective sample sizes should be 100 x number of chains (i.e., at least 400)
# # increasing max_treedepth is more about efficiency (instead of validity)
# # See: https://mc-stan.org/misc/warnings.html
# fit_fcr_with_na <- brm(y_brms + intensity_mi + system_mi + set_rescor(FALSE), data = gamma_brms_data,
#                    prior = all_priors, cores = 4, seed = "11729", iter = 20000, control = list(adapt_delta = 0.99))
# 
# # Bulk ESS for y_intensityintensive still < 400; try increasing iterto 50000?
# 
# # Get stan code
# stancode(fit_fcr_with_na)
# 
# # NEXT: open extract_brms_outputs.R to produce figures
# 
# # Missing FCRs are predicted by the model and listed as part of the outputed variables
# na_predictions <- fit_fcr_with_na %>%
#   spread_draws(Ymi_y[y_row]) %>%
#   median_qi() %>% # median values of all missing predicted responses
#   rename(FCR = Ymi_y)
# 
# #Get all the non-NAs in gamma_brms_data
# brms_non_na <- gamma_brms_data %>%
#   mutate(y_row = row_number()) %>%
#   filter(is.na(y)==FALSE) %>%
#   select(FCR = y, y_row)
# 
# # combine FCR data with predictions of missing FCRs
# dat_combine <- na_predictions %>%
#   bind_rows(brms_non_na) %>%
#   arrange(y_row) %>%
#   mutate(data_type = if_else(is.na(.point), true = "data", false = "prediction"))
# 
# # combine with original lca data to get meta data
# full_fcr_dat <- dat_combine %>%
#   left_join(lca_with_na %>% mutate(y_row = row_number()), by = "y_row")

######################################################################################################
# Step 1b: Model FCR but mirror dirichlet regression which doesn't imput missing predictors
# ie, remove all NAs and model these data, then use the model to predict FCR for those data with a complete set of predictors
# Remove FCR == 0 (species that aren't fed)

# Get model-specific data:
# SELECT STUDY ID COLUMN - use this for rejoining outputs from multiple regression models back together
# Select relevant data columns and arrange by categorical info
feed_model_dat_categories <- lca_dat_clean_groups %>%
  select(study_id, FCR, contains("new"), clean_sci_name, taxa, intensity = Intensity, system = Production_system_group) %>%
  arrange(clean_sci_name, taxa, intensity, system) %>%
  filter(taxa %in% c("mussel")==FALSE) # Remove taxa that don't belong in FCR/feed analysis - mussels

# Get footprint data
source("Functions.R")
fp_dat <- read.csv(file.path(datadir, "Feed_impact_factors_20201203.csv"))
fp_clean <- clean.feedFP(fp_dat)

fcr_no_na <- feed_model_dat_categories %>%
  #filter(FCR != 0 | is.na(FCR))  %>% # If we want to retain NAs, have to explicitly include is.na(FCR) otherwise NA's get dropped by FCR != 0
  filter(FCR != 0) %>% # This also automatically drops NAs
  filter(is.na(intensity)==FALSE & is.na(system)==FALSE) %>% # complete predictors - i.e., both intensity AND system are non-NA
  select(study_id, FCR, clean_sci_name, taxa, intensity, system)

# Set data for model:

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
  mutate(intensity = factor(intensity, levels = c("extensive", "semi", "intensive"))) %>% # set order of factors (low = extensive, high = intensive)
  mutate(system = factor(system, levels = c("open", "semi", "closed"))) %>% # set order of factors (low = open, high = closed)
  mutate(intensity = as.numeric(intensity)) %>%
  mutate(system = as.numeric(system)) %>%
  select(intensity, system) %>%
  as.matrix()
ordinal_sd<-apply(X_ordinal, MARGIN=2, FUN=sd, na.rm=TRUE) # Center all non-intercept variables and scale by 2 standard deviations (ignoring NAs)

options(na.action='na.pass') # First change default options for handling missing data
X_ordinal_scaled <- scale(X_ordinal, center=TRUE, scale=2*ordinal_sd)
options(na.action='na.omit') # Return option back to the default

# Create dataframe for brms and rename feed variables

# Create dataframe for brms and rename feed variables
fcr_brms_data <- data.frame(y = fcr_no_na$FCR, X_taxa_scaled, X_ordinal_scaled)

# Set model formula
fcr_brms <- brmsformula(y ~ 1 + taxamisc_diad + taxamisc_fresh +
                        taxamisc_marine + taxasalmon + taxashrimp + taxatilapia + taxatrout + taxatuna +
                        intensity + system, family = Gamma("log"))


# Use "resp = <response_variable>" to specify different priors for different response variables
all_priors <- c(set_prior("normal(0,5)", class = "b"), # priors for y response variables
                set_prior("normal(0,2.5)", class = "Intercept"), 
                set_prior("exponential(1)", class = "shape"))

# Model converges after increasing the adapt_delta and iterations from default values
# Rule of thumb: bulk and tail effective sample sizes should be 100 x number of chains (i.e., at least 400)
# increasing max_treedepth is more about efficiency (instead of validity)
# See: https://mc-stan.org/misc/warnings.html
fit_fcr_no_na <- brm(fcr_brms, data = fcr_brms_data,
                       prior = all_priors, cores = 4, seed = "11729", iter = 20000, control = list(adapt_delta = 0.99))

# Bulk ESS for y_intensityintensive still < 400; try increasing iterto 50000?

# Get stan code
stancode(fit_fcr_no_na)

######################################################################################################
# Use FCR model to predict NA FCRs for studies with complete set of predictors
# Both intensity AND system are non-NA
fcr_complete_predictors <- feed_model_dat_categories %>%
  filter(FCR != 0 | is.na(FCR))  %>% # Have to explicitly include is.na(FCR) otherwise NA's get dropped by FCR != 0
  filter(is.na(intensity)==FALSE & is.na(system)== FALSE) %>%
  filter(is.na(FCR))

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
levels(fcr_complete_predictors$taxa) <- list(fresh_crust = "fresh_crust", misc_diad = "misc_diad", misc_fresh = "misc_fresh", misc_marine = "misc_marine", salmon = "salmon", shrimp = "shrimp", tilapia = "tilapia", trout = "trout", tuna = "tuna")

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
  mutate(intensity = factor(intensity, levels = c("extensive", "semi", "intensive"))) %>% # set order of factors (low = extensive, high = intensive)
  mutate(system = factor(system, levels = c("open", "semi", "closed"))) %>% # set order of factors (low = open, high = closed)
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
  rename(FCR = .value)

# Bind
full_fcr_dat <- fcr_predictions %>%
  bind_rows(fcr_no_na) %>%
  mutate(data_type = if_else(is.na(.point), true = "data", false = "prediction")) %>%
  arrange(taxa, intensity, system, clean_sci_name) %>%
  rownames_to_column() # Arrange by taxa first, then create dummy column for plotting 

# NEXT: use extract_brms_gamma_no_imputation to produce figures for SI

######################################################################################################
# Step 2: Model feed proportions
# Brms dirichlet regression doesn't support missing data, so strategy here is to remove all NAs and model these data
# Then use the model to predict missing feed proportion data for those data that have a complete set of predictors

# Remove FCR == 0 (species that aren't fed)
feed_no_na <- feed_model_dat_categories %>%
  filter(FCR != 0 | is.na(FCR))  %>% # Have to explicitly include is.na(FCR) otherwise NA's get dropped by FCR != 0
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
  mutate(intensity = factor(intensity, levels = c("extensive", "semi", "intensive"))) %>% # set order of factors (low = extensive, high = intensive)
  mutate(system = factor(system, levels = c("open", "semi", "closed"))) %>% # set order of factors (low = open, high = closed)
  mutate(intensity = as.numeric(intensity)) %>%
  mutate(system = as.numeric(system)) %>%
  select(intensity, system) %>%
  as.matrix()
ordinal_sd<-apply(X_ordinal, MARGIN=2, FUN=sd, na.rm=TRUE) # Center all non-intercept variables and scale by 2 standard deviations (ignoring NAs)

options(na.action='na.pass') # First change default options for handling missing data
X_ordinal_scaled <- scale(X_ordinal, center=TRUE, scale=2*ordinal_sd)
options(na.action='na.omit') # Return option back to the default

# Create dataframe for brms and rename feed variables
feed_data <- data.frame(cbind(feed_no_na %>% select(contains("feed")), X_taxa_scaled, X_ordinal_scaled)) 

# Response variable must be a matrix, create function bind since cbind within the brm function is reserved for specifying multivariate models
bind <- function(...) cbind(...)

# CHOOSE ONE:
# Option1: For TAXA + INTENSITY + SYSTEM
feed_brms <- brmsformula(bind(feed_soy, feed_crops, feed_fmfo, feed_animal) ~ 
                        taxamisc_diad + taxamisc_fresh + taxamisc_marine + taxasalmon + taxashrimp + taxatilapia + taxatrout + taxatuna + 
                        intensity + system, family = dirichlet())
# Option2: For TAXA + INTENSITY 
# y_brms <- brmsformula(bind(feed_soy, feed_crops, feed_fmfo, feed_animal) ~ taxamisc_marine + taxasalmon + taxatilapia + taxatrout + intensity, family = dirichlet())
# Option3: For TAXA + SYSTEM
# y_brms <- brmsformula(bind(feed_soy, feed_crops, feed_fmfo, feed_animal) ~ taxamisc_marine + taxasalmon + taxatilapia + taxatrout + system, family = dirichlet())

# Model converges after increasing the adapt_delta and iterations from default values
fit_feed_no_na <- brm(feed_brms, data = feed_data,
                 cores = 4, seed = "11729", iter = 10000)

summary(fit_feed_no_na) 

######################################################################################################
# Use feed proportion model to predict NA feeds for studies with complete set of predictors

# Both intensity AND system are non-NA
feed_complete_predictors <- feed_model_dat_categories %>%
  filter(FCR != 0 | is.na(FCR))  %>% # Have to explicitly include is.na(FCR) otherwise NA's get dropped by FCR != 0
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
levels(feed_complete_predictors$taxa) <- list(fresh_crust = "fresh_crust", misc_diad = "misc_diad", misc_fresh = "misc_fresh", misc_marine = "misc_marine", salmon = "salmon", shrimp = "shrimp", tilapia = "tilapia", trout = "trout", tuna = "tuna")

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
  mutate(intensity = factor(intensity, levels = c("extensive", "semi", "intensive"))) %>% # set order of factors (low = extensive, high = intensive)
  mutate(system = factor(system, levels = c("open", "semi", "closed"))) %>% # set order of factors (low = open, high = closed)
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

# Reformat data so it can row-bind with predictions
feed_no_na_long <- feed_no_na %>%
  pivot_longer(cols = c("feed_soy", "feed_crops", "feed_fmfo", "feed_animal"), names_to = ".category", values_to = "feed_proportion")

# Bind and pivot wide again to match original data format
full_feed_dat <- feed_predictions %>%
  bind_rows(feed_no_na_long) %>%
  mutate(data_type = if_else(is.na(.point), true = "data", false = "prediction"))

# NEXT: use plot_brms_output and produce figures for SI

######################################################################################################
# BOX PLOTS OF DATA:
# Theme for ALL PLOTS (including mcmc plots)
plot_theme <- theme(title = element_text(size = 18),
                    axis.title.x = element_text(size = 16),
                    axis.text=element_text(size=14, color = "black"))

# FCR:
# plot_fcr <- full_fcr_dat %>%
#   select(clean_sci_name, FCR) 
#   #mutate(clean_sci_name = paste(clean_sci_name, row_number(), sep = ""))
ggplot(data = full_fcr_dat %>% 
         mutate(clean_sci_name = as.factor(clean_sci_name)) %>%
         mutate(clean_sci_name = fct_reorder(clean_sci_name, taxa, min)), 
       aes(y = clean_sci_name, x = FCR)) +
  #geom_boxplot(aes(color = taxa)) +
  geom_violin(aes(color = taxa), scale = "width") +
  geom_jitter(aes(color = taxa, shape = data_type), size = 3) +
  theme_classic() +
  plot_theme +
  labs(title = "Boxplots of FCRs",
       x = "",
       y = "") +
  theme(axis.text.x = element_text(hjust = 1))
ggsave(file.path(outdir, "boxplot_fcr.png"), height = 8, width = 11.5)


# feed proportion all taxa
feed_vars <- c("soy", "crops", "fmfo", "animal")
for (i in 1:length(feed_vars)) {
  p <- ggplot(data = full_feed_dat %>% 
                filter(str_detect(.category, feed_vars[i])) %>%
                mutate(clean_sci_name = as.factor(clean_sci_name)) %>%
                mutate(clean_sci_name = fct_reorder(clean_sci_name, taxa, min)), aes(y = clean_sci_name, x = feed_proportion)) +
    #geom_boxplot() +
    geom_violin(aes(color = taxa), scale = "width") +
    geom_jitter(aes(color = taxa, shape = data_type), size = 3) +
    theme_classic() +
    plot_theme +
    labs(title = paste("Boxplots of ", feed_vars[i], " feed proportions", sep = ""),
         x = "",
         y = "")  +
    theme(axis.text.x = element_text(hjust = 1))
  print(p)
  ggsave(file.path(outdir, paste("boxplot_feed-prop_", feed_vars[i], ".png", sep = "")), height = 8, width = 11.5)
}

# Feed proportions with facet_wrap
p <- ggplot(data = full_feed_dat %>% 
              #filter(str_detect(.category, feed_vars[i])) %>%
              mutate(clean_sci_name = as.factor(clean_sci_name)) %>%
              mutate(clean_sci_name = fct_reorder(clean_sci_name, taxa, min)), aes(y = clean_sci_name, x = feed_proportion)) +
  #geom_boxplot() +
  geom_violin(aes(color = taxa), scale = "width") +
  geom_jitter(aes(color = taxa, shape = data_type), size = 3) +
  theme_classic() +
  plot_theme +
  labs(title = paste("Boxplots of all feed proportions", sep = ""),
       x = "",
       y = "")  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),) +
  facet_wrap(~.category, nrow = 1)
print(p)
ggsave(file.path(outdir, "boxplot_feed-prop_all-feeds.png"), height = 8, width = 11.5)

######################################################################################################
# STEP 3 - aggregate up to taxa level and estimate total feed footprint
# Set data for model

# Merge fcr and feed datasets by study_id, remove interval info, only use medians for rest of analysis
# FIX IT - vectorize over different allocation methods of footprint data (now called impact factors)

# Format feed_dat
# Standardize number of significant digits (round to digist = 3)
# Re-scale so rowsums equals 1
feed_dat_merge <- full_feed_dat %>%
  select(study_id, clean_sci_name, taxa, intensity, system, feed_data_type = data_type, .category, contains("feed")) %>%
  pivot_wider(names_from = .category, values_from = feed_proportion) %>%
  mutate(feed_soy = round(feed_soy, digits = 3),
         feed_crops = round(feed_crops, digits = 3),
         feed_fmfo = round(feed_fmfo, digits = 3),
         feed_animal = round(feed_animal, digits = 3)) %>%
  mutate(new_scale = feed_soy + feed_crops + feed_fmfo + feed_animal) %>%
  mutate(feed_soy = feed_soy / new_scale,
         feed_crops = feed_crops / new_scale,
         feed_fmfo = feed_fmfo / new_scale,
         feed_animal = feed_animal / new_scale)

# Format FCR dat
# Rename data_type to keep track of both feed vs fcr data types
fcr_dat_merge <- full_fcr_dat %>%
  select(study_id, clean_sci_name, taxa, intensity, system, fcr_data_type = data_type, FCR)

# FIX IT - check dim of feed_dat_merge and fcr_dat_merge (which is bigger - i.e., is left_join appropriate?)
# MERGE
feed_footprint_dat <- feed_dat_merge %>%
  full_join(fcr_dat_merge, by = intersect(names(feed_dat_merge), names(fcr_dat_merge))) %>%
  drop_na() %>%
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
x <- feed_footprint_dat$FCR

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

# Data (constants) for final foot print calculation 
# FIX IT - need updated footprint data (for now just average across all allocation methods)
fp_dat <- fp_clean %>%
  filter(Category != "Energy") %>%
  group_by(FP, Category) %>%
  summarise(FP_val = mean(FP_val)) %>% 
  ungroup()

# ORDER: Animal, crop, FMFO, soy
fp_c_dat <- fp_dat %>%
  filter(FP == "Carbon") %>%
  select(-FP) %>%
  pivot_wider(names_from = Category, values_from = FP_val) %>%
  as.matrix() %>%
  c()

fp_n_dat <- fp_dat %>%
  filter(FP == "Nitrogen") %>%
  select(-FP) %>%
  pivot_wider(names_from = Category, values_from = FP_val) %>%
  as.matrix() %>%
  c()

fp_p_dat <- fp_dat %>%
  filter(FP == "Phosphorus") %>%
  select(-FP) %>%
  pivot_wider(names_from = Category, values_from = FP_val) %>%
  as.matrix() %>%
  c()

fp_land_dat <- fp_dat %>%
  filter(FP == "Land") %>%
  select(-FP) %>%
  pivot_wider(names_from = Category, values_from = FP_val) %>%
  as.matrix() %>%
  c()

fp_water_dat <- fp_dat %>%
  filter(FP == "Water") %>%
  select(-FP) %>%
  pivot_wider(names_from = Category, values_from = FP_val) %>%
  as.matrix() %>%
  c()

# Prior information
# Mean feed proportions per sci-name
sci_phi_mean <- feed_footprint_dat %>% 
  select(contains(c("soy", "crops", "fmfo", "animal", "sci", "clean_sci_name"))) %>%
  group_by(clean_sci_name, sci) %>% 
  summarise(across(contains(c("soy", "crops", "fmfo", "animal")), mean)) %>%
  ungroup() %>%
  arrange(sci)

# Mean feed proportions per taxa group
tx_phi_mean <- feed_footprint_dat %>% 
  select(contains(c("soy", "crops", "fmfo", "animal", "tx", "taxa"))) %>%
  group_by(taxa, tx) %>% 
  summarise(across(contains(c("soy", "crops", "fmfo", "animal")), mean)) %>%
  ungroup() %>%
  arrange(tx)


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
stan_no_na <- 'data {
  // data for gamma model for FCR
  int<lower=0> N;  // number of observations
  vector<lower=0>[N] x; // data
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
  vector<lower=0>[N_SCI] sci_mu;
  // only need sigmas if defining shape and rate with mu and sigma
  //real<lower=0> tx_sigma;
  //real<lower=0> sci_sigma;
  // if using variance instead of st dev
  real<lower=0> tx_sigma_sq;
  real<lower=0> sci_sigma_sq;

  // Feed proportion model:
  simplex[K] sci_phi[N_SCI];
  simplex[K] sci_theta[N_SCI]; // vectors of estimated sci-level feed weight simplexes
  simplex[K] tx_phi[N_TX];
  simplex[K] tx_theta[N_TX];

  // Params for the dirichlet priors:
  // real<lower=0> sigma_1;
  // real<lower=0> sigma_2;
}
transformed parameters {
  // define transofrmed params for gamma model for FCRs
  vector<lower=0>[N_SCI] sci_shape;
  vector<lower=0>[N_SCI] sci_rate;
  vector<lower=0>[N_TX] tx_shape;
  vector<lower=0>[N_TX] tx_rate;

  // define params for dirichlet model for feed proportions
  vector<lower=0>[K] sci_alpha[N_SCI];
  vector<lower=0>[K] tx_alpha[N_TX];

  // gamma model reparameterization
  // option 1: reparamaterize gamma to get mu and sigma; defining these here instead of the model section allows us to see these parameters in the output
  // taxa group level
  //for (n_tx in 1:N_TX){
  //  tx_shape[n_tx] = square(tx_mu[n_tx]) ./ square(tx_sigma);
  //  tx_rate[n_tx] = tx_mu[n_tx] ./ square(tx_sigma);
  //}

  // sci level
  //for (n_sci in 1:N_SCI){
  //  sci_shape[n_sci] = square(sci_mu[n_sci]) ./ square(sci_sigma);
  //  sci_rate[n_sci] = sci_mu[n_sci] ./ square(sci_sigma);
  //}

  // option 2: reparameterize gamma to get just mu
  //for (n_tx in 1:N_TX){
  //  tx_rate[n_tx] = tx_shape[n_tx] ./ tx_mu[n_tx];
  //}
  // sci level
  //for (n_sci in 1:N_SCI){
  //  sci_rate[n_sci] = sci_shape[n_sci] ./ sci_mu[n_sci];
  //}
  
  // option 3: reparameterize shape and rate as inverse(va) and inverse(va)/mu
  for (n_tx in 1:N_TX){
    tx_shape[n_tx] = 1 / tx_sigma_sq;
    tx_rate[n_tx] = 1 / (tx_sigma_sq * tx_mu[n_tx]);
  }
  for (n_sci in 1:N_SCI){
    sci_shape[n_sci] = 1 / sci_sigma_sq;
    sci_rate[n_sci] = 1 / (sci_sigma_sq * sci_mu[n_sci]);
  }
  
  
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
  // define priors for gamma model for FCRs
  // Put priors on mu and sigma (instead of shape and rate) since this is more intuitive:
  tx_mu ~ uniform(0, 100); // note: uniform(0,100) for all of these doesnt help much with convergence
  sci_mu ~ uniform(0, 100);
  //tx_sigma ~ uniform(0, 100); // only need sigmas if calculating shape and rate with mu and sigma
  //sci_sigma ~ uniform(0, 100);

  // define priors for dirichlet model for feed proportions
  // sci_phi defined as sci_phi[n_to_sci][K]

  // option 1: define feed proportion priors as lower upper bounds
  // sci_phi[2][1] ~ uniform(0.1, 0.2); // hypothetical lower and upper bounds for Oncorhynchus mykiss soy

  // option 2: define feed proportions as means (need to define sigmas in parameters block: real<lower=0> sigma_1, sigma_2 etc;)
  // sci_phi[2][1] ~ normal(0.13, sigma_1); // mean for Oncorhynhchus mykiss soy feed

  // likelihood
  // gamma model sci-name and taxa-level
  for (n in 1:N){
    x[n] ~ gamma(sci_shape[n_to_sci[n]], sci_rate[n_to_sci[n]]);
  }

  for (n_sci in 1:N_SCI){
    sci_mu[n_sci] ~ gamma(tx_shape[sci_to_tx[n_sci]], tx_rate[sci_to_tx[n_sci]]);
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
  vector[K] tx_feed_c_fp[N_TX]; // array of k feeds and N_TX taxa groups
  vector[N_TX] tx_sum_feed_c_fp;
  vector[N_SCI] sci_c_fp;
  vector[K] sci_feed_c_fp[N_SCI];
  vector[N_SCI] sci_sum_feed_c_fp;
  // Nitrogen

  // Phosphorus

  // Land

  // Water

  // Calculations
  for (n_tx in 1:N_TX) {
    tx_feed_c_fp[n_tx] = fp_c_dat .* tx_theta[n_tx];
    tx_sum_feed_c_fp[n_tx] = sum(tx_feed_c_fp[n_tx]);
    tx_c_fp[n_tx] = tx_mu[n_tx] * tx_sum_feed_c_fp[n_tx];
  }
  for (n_sci in 1:N_SCI) {
    sci_feed_c_fp[n_sci] = fp_c_dat .* sci_theta[n_sci];
    sci_sum_feed_c_fp[n_sci] = sum(sci_feed_c_fp[n_sci]);
    sci_c_fp[n_sci] = sci_mu[n_sci] * sci_sum_feed_c_fp[n_sci];
  }
}'




# NORMAL distribution hierarchical model
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
#   real<lower=0> tx_sigma;
#   vector<lower=0>[N_SCI] sci_mu;
#   real<lower=0> sci_sigma;
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
#   // define params for dirichlet model for feed proportions
#   vector<lower=0>[K] sci_alpha[N_SCI];
#   vector<lower=0>[K] tx_alpha[N_TX];
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
#   //tx_mu ~ uniform(0, 100); // note: uniform(0,100) for all of these doesnt help much with convergence
#   //tx_sigma ~ uniform(0, 100);
#   //sci_mu ~ uniform(0, 100);
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
#   // normal model sci-name and taxa-level
#   for (n in 1:N){
#     x[n] ~ normal(sci_mu[n_to_sci[n]], sci_sigma);
#   }
#   
#   for (n_sci in 1:N_SCI){
#     sci_mu[n_sci] ~ normal(tx_mu[sci_to_tx[n_sci]], tx_sigma);
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
#   vector[N_SCI] sci_c_fp;
# 
#   // Nitrogen
#   vector[N_TX] tx_n_fp;
#   vector[N_SCI] sci_n_fp;
# 
#   // Phosphorus
#   vector[N_TX] tx_p_fp;
#   vector[N_SCI] sci_p_fp;
#   
#   // Land
#   vector[N_TX] tx_land_fp;
#   vector[N_SCI] sci_land_fp;
# 
#   // Water
#   vector[N_TX] tx_water_fp;
#   vector[N_SCI] sci_water_fp;
#   
#   // Calculations
#   for (n_tx in 1:N_TX) {
#     tx_c_fp[n_tx] = tx_mu[n_tx] * sum(fp_c_dat .* tx_theta[n_tx]);
#     tx_n_fp[n_tx] = tx_mu[n_tx] * sum(fp_n_dat .* tx_theta[n_tx]);
#     tx_p_fp[n_tx] = tx_mu[n_tx] * sum(fp_p_dat .* tx_theta[n_tx]);
#     tx_land_fp[n_tx] = tx_mu[n_tx] * sum(fp_land_dat .* tx_theta[n_tx]);
#     tx_water_fp[n_tx] = tx_mu[n_tx] * sum(fp_water_dat .* tx_theta[n_tx]);
#   }
#   for (n_sci in 1:N_SCI) {
#     sci_c_fp[n_sci] = sci_mu[n_sci] * sum(fp_c_dat .* sci_theta[n_sci]);
#     sci_n_fp[n_sci] = sci_mu[n_sci] * sum(fp_n_dat .* sci_theta[n_sci]);
#     sci_p_fp[n_sci] = sci_mu[n_sci] * sum(fp_p_dat .* sci_theta[n_sci]);
#     sci_land_fp[n_sci] = sci_mu[n_sci] * sum(fp_land_dat .* sci_theta[n_sci]);
#     sci_water_fp[n_sci] = sci_mu[n_sci] * sum(fp_water_dat .* sci_theta[n_sci]);
#   }
# }'


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

######################################################################################################
# PLOT RESULTS

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
         full_taxa_name = case_when(taxa == "misc_marine" ~ "misc marine fishes",
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
  median_qi(.width = 0.8) %>%
  left_join(index_key, by = "sci") %>% # Join with index key to get sci and taxa names
  mutate(clean_sci_name = fct_reorder(clean_sci_name, as.character(full_taxa_name))) %>%
  ggplot(aes(y = clean_sci_name, x = sci_c_fp, xmin = .lower, xmax = .upper, color = full_taxa_name)) +
  geom_pointinterval() +
  theme_classic() + 
  plot_theme + 
  labs(x = "kg CO2-eq", y = "", title = "Carbon", color = "taxa group")
ggsave(filename = file.path(outdir, "plot_feed-footprint_carbon_obs-level.png"), width = 11, height = 8.5)

# Nitrogen
fit_no_na %>%
  spread_draws(sci_n_fp[sci]) %>%
  median_qi(.width = 0.8) %>%
  left_join(index_key, by = "sci") %>% # Join with index key to get sci and taxa names
  mutate(clean_sci_name = fct_reorder(clean_sci_name, as.character(full_taxa_name))) %>%
  ggplot(aes(y = clean_sci_name, x = sci_n_fp, xmin = .lower, xmax = .upper, color = full_taxa_name)) +
  geom_pointinterval() +
  theme_classic() + 
  plot_theme + 
  labs(x = "kg N-eq", y = "", title = "Nitrogen", color = "taxa group")
ggsave(filename = file.path(outdir, "plot_feed-footprint_nitrogen_obs-level.png"), width = 11, height = 8.5)

# Phosphorus
fit_no_na %>%
  spread_draws(sci_p_fp[sci]) %>%
  median_qi(.width = 0.8) %>%
  left_join(index_key, by = "sci") %>% # Join with index key to get sci and taxa names
  mutate(clean_sci_name = fct_reorder(clean_sci_name, as.character(full_taxa_name))) %>%
  ggplot(aes(y = clean_sci_name, x = sci_p_fp, xmin = .lower, xmax = .upper, color = full_taxa_name)) +
  geom_pointinterval() +
  theme_classic() + 
  plot_theme + 
  labs(x = "kg P-eq", y = "", title = "Phosphorus", color = "taxa group")
ggsave(filename = file.path(outdir, "plot_feed-footprint_phosphorus_obs-level.png"), width = 11, height = 8.5)

# Land
fit_no_na %>%
  spread_draws(sci_land_fp[sci]) %>%
  median_qi(.width = 0.8) %>%
  left_join(index_key, by = "sci") %>% # Join with index key to get sci and taxa names
  mutate(clean_sci_name = fct_reorder(clean_sci_name, as.character(full_taxa_name))) %>%
  ggplot(aes(y = clean_sci_name, x = sci_land_fp, xmin = .lower, xmax = .upper, color = full_taxa_name)) +
  geom_pointinterval() +
  theme_classic() + 
  plot_theme + 
  labs(x = bquote('m'^2~'a'), y = "", title = "Land", color = "taxa group")
ggsave(filename = file.path(outdir, "plot_feed-footprint_land_obs-level.png"), width = 11, height = 8.5)

# Water
fit_no_na %>%
  spread_draws(sci_water_fp[sci]) %>%
  median_qi(.width = 0.8) %>%
  left_join(index_key, by = "sci") %>% # Join with index key to get sci and taxa names
  mutate(clean_sci_name = fct_reorder(clean_sci_name, as.character(full_taxa_name))) %>%
  ggplot(aes(y = clean_sci_name, x = sci_water_fp, xmin = .lower, xmax = .upper, color = full_taxa_name)) +
  geom_pointinterval() +
  theme_classic() + 
  plot_theme + 
  labs(x = bquote('m'^3), y = "", title = "Water", color = "taxa group")
ggsave(filename = file.path(outdir, "plot_feed-footprint_water_obs-level.png"), width = 11, height = 8.5)

# Taxa-level feed footprints as densities:
# Carbon
fit_no_na %>%
  spread_draws(tx_c_fp[tx]) %>%
  left_join(index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  ggplot(aes(y = full_taxa_name, x = tx_c_fp)) +
  stat_halfeye(aes(slab_fill = full_taxa_name)) +
  theme_classic() + 
  plot_theme + 
  theme(legend.position = "none") +
  labs(x = "kg CO2-eq", y = "", title = "Carbon")
ggsave(filename = file.path(outdir, "plot_feed-footprint_carbon_taxa-level.png"), width = 11, height = 8.5)

# Nitrogen
fit_no_na %>%
  spread_draws(tx_n_fp[tx]) %>%
  left_join(index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  ggplot(aes(y = full_taxa_name, x = tx_n_fp)) +
  stat_halfeye(aes(slab_fill = full_taxa_name)) +
  theme_classic() + 
  plot_theme + 
  theme(legend.position = "none") +
  labs(x = "kg N-eq", y = "", title = "Nitrogen")
ggsave(filename = file.path(outdir, "plot_feed-footprint_nitrogen_taxa-level.png"), width = 11, height = 8.5)


# Phosphorus
fit_no_na %>%
  spread_draws(tx_p_fp[tx]) %>%
  left_join(index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  ggplot(aes(y = full_taxa_name, x = tx_p_fp)) +
  stat_halfeye(aes(slab_fill = full_taxa_name)) +
  theme_classic() + 
  plot_theme + 
  theme(legend.position = "none") +
  labs(x = "kg P-eq", y = "", title = "Phosphorus")
ggsave(filename = file.path(outdir, "plot_feed-footprint_phosphorus_taxa-level.png"), width = 11, height = 8.5)


# Land
fit_no_na %>%
  spread_draws(tx_land_fp[tx]) %>%
  left_join(index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  ggplot(aes(y = full_taxa_name, x = tx_land_fp)) +
  stat_halfeye(aes(slab_fill = full_taxa_name)) +
  theme_classic() + 
  plot_theme + 
  theme(legend.position = "none") +
  labs(x = bquote('m'^2~'a'), y = "", title = "Land")
ggsave(filename = file.path(outdir, "plot_feed-footprint_land_taxa-level.png"), width = 11, height = 8.5)

# Water
fit_no_na %>%
  spread_draws(tx_water_fp[tx]) %>%
  left_join(index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  ggplot(aes(y = full_taxa_name, x = tx_water_fp)) +
  stat_halfeye(aes(slab_fill = full_taxa_name)) +
  theme_classic() + 
  plot_theme + 
  theme(legend.position = "none") +
  labs(x = bquote('m'^3), y = "", title = "Water")
ggsave(filename = file.path(outdir, "plot_feed-footprint_water_taxa-level.png"), width = 11, height = 8.5)


# OPTION 2: Taxa-level plots with color themes:
# If we want to mimic bayesplot color schemes, can get hexadecimal colors and input manually to stat_halfeye aesthetics
color_scheme_get("blue")
color_scheme_get("green")










# Remaining code below was for initial testing/model building:
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################

feed_footprint_dat <- feed_dat_merge %>%
  left_join(fcr_dat_merge, by = intersect(names(feed_dat_merge), names(fcr_dat_merge))) %>%
  arrange(clean_sci_name, taxa) %>%
  # FILTER TO SMALLER DATASET
  #filter(taxa %in% c("trout", "tilapia")) %>%
  # FILTER OUTLIER DATA
  #filter(FCR < 2 & FCR > 1) %>%
  # FILTER TO ONLY INCLUDE sci-names with 3 or more observations
  #group_by(clean_sci_name) %>%
  #mutate(n_obs = n(),
  #       sci_mean = mean(FCR)) %>%
  #ungroup() %>% 
  #filter(n_obs > 2) %>%
  # Model doesn't converge once taxa level is included - Get taxa means to use for priors
  #group_by(taxa) %>%
  #mutate(taxa_mean = mean(FCR)) %>%
  #ungroup() %>%
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
x <- feed_footprint_dat$FCR

# Get two-level model with generated quantities for just the FCR (gamma model) to converge
stan_data <- list(N = N,
                  N_SCI = N_SCI, 
                  n_to_sci = n_to_sci,
                  N_TX = N_TX,
                  sci_to_tx = sci_to_tx,
                  x = x,
                  fp_c_dat = fp_c_dat,
                  fp_n_dat = fp_n_dat,
                  fp_p_dat = fp_p_dat,
                  fp_land_dat = fp_land_dat,
                  fp_water_dat = fp_water_dat)

# stan_no_na <- 'data {
#   // data for gamma model for FCR
#   int<lower=0> N;  // number of observations
#   vector<lower=0>[N] x; // data
#   //int N_TX; // number of taxa groups
#   int N_SCI; // number of scientific names
#   int n_to_sci[N]; // sciname index
#   //int sci_to_tx[N_SCI]; // taxa-group indices
# 
# }
# parameters {
#   // FCR model:
#   //real<lower=0> mu;
#   //real<lower=0> sigma;
#   //vector<lower=0>[N_TX] tx_mu;
#   //real<lower=0> tx_sigma;
#   vector<lower=0>[N_SCI] sci_mu;
#   real<lower=0> sci_sigma;
# }
# transformed parameters {
#   // define transofrmed params for gamma model for FCRs
#   //real shape;
#   //real rate;
#   vector[N_SCI] sci_shape;
#   vector[N_SCI] sci_rate;
#   //vector[N_TX] tx_shape;
#   //vector[N_TX] tx_rate;
# 
#   // gamma model reparameterization
#   // reparamaterize gamma to get mu and sigma; defining these here instead of the model section allows us to see these parameters in the output
#   // global-level
#   //shape = square(mu) / square(sigma);
#   //rate = mu / square(sigma);
#   
#   // taxa group level
#   //for (n_tx in 1:N_TX){
#   //  tx_shape[n_tx] = square(tx_mu[n_tx]) ./ square(tx_sigma);
#   //  tx_rate[n_tx] = tx_mu[n_tx] ./ square(tx_sigma);
#   //}
#   
#   // sci level
#   for (n_sci in 1:N_SCI){
#     sci_shape[n_sci] = square(sci_mu[n_sci]) ./ square(sci_sigma);
#     sci_rate[n_sci] = sci_mu[n_sci] ./ square(sci_sigma);
#   }
# 
# }
# model {
#   // define priors for gamma model for FCRs
#   // Put priors on mu and sigma (instead of shape and rate) since this is more intuitive:
#   //mu ~ uniform(0, 10);
#   //sigma ~ uniform(0, 10);
#   //tx_mu ~ uniform(0, 100);
#   sci_mu ~ uniform(0, 100);
#   
#   //tx_mu[1] ~ uniform(0, 2);
#   //tx_mu[2] ~ uniform(0, 2);
#   //tx_mu[3] ~ uniform(0, 2);
#   //tx_mu[4] ~ uniform(0, 2);
#   //tx_mu[5] ~ uniform(0, 2);
#   //tx_mu[6] ~ uniform(0, 3);
#   //tx_mu[7] ~ uniform(0, 2);
#   //tx_mu[8] ~ uniform(0, 2);
#   
# 
#   //sci_mu[1] ~ uniform(0, 5);
#   //sci_mu[2] ~ uniform(0, 5);
#   //sci_mu[3] ~ uniform(0, 5);
#   //sci_mu[4] ~ uniform(0, 5);
#   //sci_mu[5] ~ uniform(0, 5);
#   //sci_mu[6] ~ uniform(0, 5);
#   //sci_mu[7] ~ uniform(0, 5);
#   //sci_mu[8] ~ uniform(0, 5);
#   //sci_mu[9] ~ uniform(0, 5);
#   //sci_mu[10] ~ uniform(0, 5);
#   //sci_mu[11] ~ uniform(0, 5);
#   
#   //tx_sigma ~ uniform(0, 100);
#   sci_sigma ~ uniform(0, 100);
#   
#   // likelihood
#   // gamma model sci-name and taxa-level
#   for (n in 1:N){
#     x[n] ~ gamma(sci_shape[n_to_sci[n]], sci_rate[n_to_sci[n]]);
#   }
#   
#   //for (n_sci in 1:N_SCI){
#   //  sci_mu[n_sci] ~ gamma(tx_shape[sci_to_tx[n_sci]], tx_rate[sci_to_tx[n_sci]]);
#   //}
#   
#   //for (n_tx in 1:N_TX){
#   //  tx_mu[n_tx] ~ gamma(shape, rate);
#   //}
#   
# }'


# Try hierarchical model with normal distribution
stan_no_na <- 'data {
  // data for gamma model for FCR
  int<lower=0> N;  // number of observations
  vector<lower=0>[N] x; // data
  int N_TX; // number of taxa groups
  int N_SCI; // number of scientific names
  int n_to_sci[N]; // sciname index
  int sci_to_tx[N_SCI]; // taxa-group indices

}
parameters {
  // FCR model:
  //real<lower=0> mu;
  //real<lower=0> sigma;
  vector<lower=0>[N_TX] tx_mu;
  real<lower=0> tx_sigma;
  vector<lower=0>[N_SCI] sci_mu;
  real<lower=0> sci_sigma;
}
model {
  // define priors for FCRs
  // Put priors on mu and sigma (instead of shape and rate) since this is more intuitive:
  //mu ~ uniform(0, 10);
  //sigma ~ uniform(0, 10);
  //tx_mu ~ uniform(0, 100);
  //sci_mu ~ uniform(0, 100);
  //tx_sigma ~ uniform(0, 100);
  //sci_sigma ~ uniform(0, 100);

  // likelihood
  // normal model sci-name and taxa-level
  for (n in 1:N){
    x[n] ~ normal(sci_mu[n_to_sci[n]], sci_sigma);
  }

  for (n_sci in 1:N_SCI){
    sci_mu[n_sci] ~ normal(tx_mu[sci_to_tx[n_sci]], tx_sigma);
  }

  //for (n_tx in 1:N_TX){
  //  tx_mu[n_tx] ~ gamma(shape, rate);
  //}

}'


no_na_mod <- stan_model(model_code = stan_no_na)
# Note: For Windows, apparently OK to ignore this warning message:
# Warning message:
#   In system(paste(CXX, ARGS), ignore.stdout = TRUE, ignore.stderr = TRUE) :
#   'C:/rtools40/usr/mingw_/bin/g++' not found

# Fit model:
# Set seed while testing
fit_no_na <- sampling(object = no_na_mod, data = stan_data, cores = 4, seed = "11729", iter = 10000, control = list(adapt_delta = 0.99))

launch_shinystan(fit_no_na)






lca_dat_no_zeroes <- clean.lca(LCA_data = lca_dat) %>%
  #select(clean_sci_name, Feed_soy_percent, Feed_othercrops_percent, Feed_FMFO_percent, Feed_animal_percent, taxa_group_name) %>%
  # NOTE multinomial-dirchlet model requires all elements > 0 (change some to 0.001 for now?)
  mutate(feed_soy_new = if_else(Feed_soy_percent == 0, true = 0.01, false = Feed_soy_percent),
         feed_crops_new = if_else(Feed_othercrops_percent == 0, true = 0.01, false = Feed_othercrops_percent),
         feed_fmfo_new = if_else(Feed_FMFO_percent == 0, true = 0.01, false = Feed_FMFO_percent),
         feed_animal_new = if_else(Feed_animal_percent == 0, true = 0.01, false = Feed_animal_percent)) %>%
  # Renomoralize values so they sum to 1
  mutate(sum = rowSums(select(., contains("new")))) %>%
  mutate(feed_soy_new = feed_soy_new / sum,
         feed_crops_new = feed_crops_new / sum,
         feed_fmfo_new = feed_fmfo_new / sum,
         feed_animal_new = feed_animal_new / sum) 

fp_dat <- read.csv(file.path(datadir, "Feed_FP_raw.csv"))
fp_clean <- clean.feedFP(fp_dat)

######################################################################################################
# Model 2: Remove all NAs - estimate feed footprint for all sci names

# Remove NAs
lca_dat_no_na <- lca_dat_no_zeroes %>%
  filter(is.na(Feed_soy_percent)==FALSE) %>%
  filter(is.na(FCR) == FALSE) %>%
  filter(FCR != 0) %>%
  mutate(clean_sci_name = as.factor(clean_sci_name),
         sci = as.numeric(clean_sci_name),
         taxa_group_name = as.factor(taxa_group_name),
         tx = as.numeric(taxa_group_name)) %>%
  select(clean_sci_name, sci, taxa_group_name, tx, FCR, Feed_soy_percent, Feed_othercrops_percent, Feed_FMFO_percent, Feed_animal_percent) %>%
  # NOTE multinomial-dirchlet model requires all elements > 0 (change some to 0.001 for now?)
  mutate(feed_soy_new = if_else(Feed_soy_percent == 0, true = 0.01, false = Feed_soy_percent),
         feed_crops_new = if_else(Feed_othercrops_percent == 0, true = 0.01, false = Feed_othercrops_percent),
         feed_fmfo_new = if_else(Feed_FMFO_percent == 0, true = 0.01, false = Feed_FMFO_percent),
         feed_animal_new = if_else(Feed_animal_percent == 0, true = 0.01, false = Feed_animal_percent)) %>%
  # Renomoralize values so they sum to 1
  mutate(sum = rowSums(select(., contains("new")))) %>%
  mutate(feed_soy_new = feed_soy_new / sum,
         feed_crops_new = feed_crops_new / sum,
         feed_fmfo_new = feed_fmfo_new / sum,
         feed_animal_new = feed_animal_new / sum) %>%
  select(clean_sci_name, sci, taxa_group_name, tx, FCR, contains("new"))

# Try for just two scientific names:
lca_dat_no_na <- lca_dat_no_na %>%
  filter(clean_sci_name %in% c("Oncorhynchus mykiss", "Salmo salar"))

# BOX PLOTS OF DATA:
# Theme for ALL PLOTS (including mcmc plots)
plot_theme <- theme(axis.text=element_text(size=14, color = "black"))

# FCR:
plot_fcr <- lca_dat_no_na %>%
  select(clean_sci_name, FCR) %>%
  #mutate(clean_sci_name = paste(clean_sci_name, row_number(), sep = "")) %>%
  pivot_longer(cols = FCR)
ggplot(data = plot_fcr, aes(x = clean_sci_name, y = value)) +
  geom_boxplot() +
  theme_classic() +
  plot_theme +
  labs(title = "Boxplots of FCRs",
       x = "",
       y = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#ggsave(file.path(outdir, "boxplot_fcr_no_na.png"), height = 8, width = 11.5)


# feed proportion:
plot_feed_prop <- lca_dat_no_na %>%
  select(clean_sci_name, soy = feed_soy_new, crops = feed_crops_new, fmfo = feed_fmfo_new, animal = feed_animal_new) %>%
  pivot_longer(cols = soy:animal)



feed_vars <- c("soy", "crops", "fmfo", "animal")
for (i in 1:length(feed_vars)) {
  p <- ggplot(data = plot_feed_prop %>% filter(name == feed_vars[i]), aes(x = clean_sci_name, y = value)) +
    geom_boxplot() +
    theme_classic() +
    plot_theme +
    labs(title = paste("Boxplots of ", feed_vars[i], " feed proportions", sep = ""),
         x = "",
         y = "")  +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p)
  #ggsave(file.path(outdir, "boxplot_feed-prop_no_na.png"), height = 8, width = 11.5)
}

# Set data for model

# overall model:
N = nrow(lca_dat_no_na)
N_SCI <- length(unique(lca_dat_no_na$sci))
sci <- lca_dat_no_na$sci
N_TX <- length(unique(lca_dat_no_na$tx))
tx <- lca_dat_no_na$tx

# for FCR model:
x <- lca_dat_no_na$FCR

# for Feed proportion model:
K = 4
feed_weights <- lca_dat_no_na %>%
  select(contains("new")) %>%
  as.matrix()
  
# Get counts per sci name and counts per taxa group (also included as data in the model):
sci_kappa <- lca_dat_no_na %>% 
  select(contains(c("new", "sci", "obs"))) %>%
  group_by(sci) %>% 
  summarise(n_obs = n()) %>%
  ungroup() %>%
  arrange(sci) %>%
  pull(n_obs)

tx_kappa <- lca_dat_no_na %>% 
  select(contains(c("new", "tx", "obs"))) %>%
  group_by(tx) %>% 
  summarise(n_obs = n()) %>%
  ungroup() %>%
  arrange(tx) %>%
  pull(n_obs)

# data (constants) for final foot print calculation 
fp_dat <- fp_clean %>%
  filter(Category != "Energy") %>%
  select(FP, Category, FP_val) 

# ORDER: Animal, crop, FMFO, soy
fp_carbon_dat <- fp_dat %>%
  filter(FP == "Carbon") %>%
  select(-FP) %>%
  pivot_wider(names_from = Category, values_from = FP_val) %>%
  as.matrix() %>%
  c()

fp_nitrogen_dat <- fp_dat %>%
  filter(FP == "Nitrogen") %>%
  select(-FP) %>%
  pivot_wider(names_from = Category, values_from = FP_val) %>%
  as.matrix() %>%
  c()

fp_phosphorus_dat <- fp_dat %>%
  filter(FP == "Phosphorus") %>%
  select(-FP) %>%
  pivot_wider(names_from = Category, values_from = FP_val) %>%
  as.matrix() %>%
  c()

fp_land_dat <- fp_dat %>%
  filter(FP == "Land") %>%
  select(-FP) %>%
  pivot_wider(names_from = Category, values_from = FP_val) %>%
  as.matrix() %>%
  c()

fp_water_dat <- fp_dat %>%
  filter(FP == "Water") %>%
  select(-FP) %>%
  pivot_wider(names_from = Category, values_from = FP_val) %>%
  as.matrix() %>%
  c()

# Prior information
# Mean feed proportions per sci-name
sci_phi_mean <- lca_dat_no_na %>% 
  select(contains(c("new", "sci", "clean_sci_name"))) %>%
  group_by(clean_sci_name, sci) %>% 
  summarise(across(contains("new"), mean)) %>%
  ungroup() %>%
  arrange(sci)

# Mean feed proportions per taxa group
tx_phi_mean <- lca_dat_no_na %>% 
  select(contains(c("new", "tx", "taxa_group_name"))) %>%
  group_by(taxa_group_name, tx) %>% 
  summarise(across(contains("new"), mean)) %>%
  ungroup() %>%
  arrange(tx)


stan_data <- list(N = N,
                  N_SCI = N_SCI, 
                  sci = sci,
                  N_TX = N_TX,
                  tx = tx,
                  x = x,
                  K = K,
                  feed_weights = feed_weights,
                  sci_kappa = sci_kappa,
                  tx_kappa = tx_kappa,
                  fp_carbon_dat = fp_carbon_dat,
                  fp_nitrogen_dat = fp_nitrogen_dat,
                  fp_phosphorus_dat = fp_phosphorus_dat,
                  fp_land_dat = fp_land_dat,
                  fp_water_dat = fp_water_dat)

# Estimate foot print for all scientific names without NAs
stan_no_na <- 'data {
  // data for gamma model for FCR
  int<lower=0> N;  // number of observations
  vector<lower=0>[N] x; // data
  int N_TX; // number of taxa groups
  int tx[N]; // taxa group index (ordered by unique sci index)
  int N_SCI; // number of scientific names
  int sci[N]; // sciname index
  
  // data for dirichlet model for feed
  int K; // number of feed types
  simplex[K] feed_weights[N]; // array of observed feed weights simplexes
  int sci_kappa[N_SCI]; // number of observations per sci-name
  int tx_kappa[N_TX]; // number of observations per taxa group
  
  // constants for generated quantities
  vector[K] fp_carbon_dat;
  vector[K] fp_nitrogen_dat;
  vector[K] fp_phosphorus_dat;
  vector[K] fp_land_dat;
  vector[K] fp_water_dat;
}
parameters {
  // FCR model:
  real<lower=0> mu;
  real<lower=0> sigma;
  vector<lower=0>[N_TX] tx_mu;
  real<lower=0> tx_sigma;
  vector<lower=0>[N_SCI] sci_mu;
  real<lower=0> sci_sigma;
  
  // Feed proportion model:
  simplex[K] sci_phi[N_SCI];
  simplex[K] sci_theta[N_SCI]; // vectors of estimated sci-level feed weight simplexes
  simplex[K] tx_phi[N_TX];
  simplex[K] tx_theta[N_TX];
  simplex[K] phi;
  simplex[K] theta;
  
  // Params for the dirichlet priors:
  // real<lower=0> sigma_1;
  // real<lower=0> sigma_2;
}
transformed parameters {
  // define transofrmed params for gamma model for FCRs
  real shape;
  real rate; 
  vector[N_SCI] sci_shape;
  vector[N_SCI] sci_rate;
  vector[N_TX] tx_shape;
  vector[N_TX] tx_rate;
  
  // define params for dirichlet model for feed proportions
  vector<lower=0>[K] sci_alpha[N_SCI];
  vector<lower=0>[K] tx_alpha[N_TX];
  vector<lower=0>[K] alpha;
  
  // reparamaterize gamma to get mu and sigma; defining these here instead of the model section allows us to see these parameters in the output
  // global-level
  shape = square(mu) / square(sigma);
  rate = mu / square(sigma);
  // sci name and taxa group levels
  for (n in 1:N){
    tx_shape[tx[n]] = square(tx_mu[tx[n]]) ./ square(tx_sigma);
    tx_rate[tx[n]] = tx_mu[tx[n]] ./ square(tx_sigma);
    sci_shape[sci[n]] = square(sci_mu[sci[n]]) ./ square(sci_sigma);
    sci_rate[sci[n]] = sci_mu[sci[n]] ./ square(sci_sigma);
  }
  
  // reparameterize alphas as a vector of means (phi) and counts (kappas)
  // phi is expected value of theta (mean feed weights)
  // kappa is strength of the prior measured in number of prior observations (minus K)
  alpha = N * phi;
  for (n_tx in 1:N_TX) {
    tx_alpha[n_tx] = tx_kappa[n_tx] * tx_phi[n_tx];
  }    
  
  for (n_sci in 1:N_SCI) {
    sci_alpha[n_sci] = sci_kappa[n_sci] * sci_phi[n_sci];
  }
}
model {
  // define priors for gamma model for FCRs
  // Put priors on mu and sigma (instead of shape and rate) since this is more intuitive:
  mu ~ uniform(0, 10); // note: uniform(0,100) for all of these doesnt help much with convergence
  sigma ~ uniform(0, 10);
  tx_mu ~ uniform(0, 10);
  tx_sigma ~ uniform(0, 10);
  sci_mu ~ uniform(0, 10);
  sci_sigma ~ uniform(0, 10);
  
  // define priors for dirichlet model for feed proportions
  // sci_phi defined as sci_phi[sci][K]
  
  // option 1: define feed proportion priors as lower upper bounds
  // sci_phi[2][1] ~ uniform(0.1, 0.2); // hypothetical lower and upper bounds for Oncorhynchus mykiss soy 
  // sci_phi[6][1] ~ uniform(0.05, 0.2); // lower upper bounds fo Salmo salar soy
  
  // option 2: define feed proportions as means (need to define sigmas in parameters block: real<lower=0> sigma_1, sigma_2 etc;)
  // sci_phi[2][1] ~ normal(0.13, sigma_1); // mean for Oncorhynhchus mykiss soy feed 
  // sci_phi[6][1] ~ normal(0.7, sigma_2); // mean for Salmo salar soy feed
  // tx_phi[6][1] ~ normal(0.13, sigma_1); // mean for taxa-level trout soy feed 
  // tx_phi[4][1] ~ normal(0.7, sigma_2); // mean for taxa-level salmon/char soy feed
  
  // likelihood
  
  // global-level for dirichlet
  tx_mu ~ gamma(shape, rate);
  
  for (n in 1:N){
    // gamma model sci-name and taxa-level
    sci_mu[sci[n]] ~ gamma(tx_shape[tx[n]], tx_rate[tx[n]]);
    x[n] ~ gamma(sci_shape[sci[n]], sci_rate[sci[n]]);
    
    // dirichlet model sci-name and taxa-level
    tx_phi[tx[n]] ~ dirichlet(alpha);
    sci_phi[sci[n]] ~ dirichlet(to_vector(tx_alpha[tx[n]]));
    feed_weights[n] ~ dirichlet(to_vector(sci_alpha[sci[n]])); 
  }

  // dirichlet model - estimate feed weights based on estimated alphas
  // global level estimates
  theta ~ dirichlet(to_vector(alpha));
  // taxa level estimates
  for (n_tx in 1:N_TX) {
    tx_theta[n_tx] ~ dirichlet(to_vector(tx_alpha[n_tx]));
  }
  // sci-name level estimates
  for (n_sci in 1:N_SCI) {
    sci_theta[n_sci] ~ dirichlet(to_vector(sci_alpha[n_sci]));
  }
}'

# Add generated quantities later
# generated quantities {
#   // Carbon
#   real<lower=0> species_carbon_footprint;
#   vector[K] feed_carbon_footprint;
#   real total_feed_carbon_footprint;
#   // Nitrogen
#   real<lower=0> species_nitrogen_footprint;
#   vector[K] feed_nitrogen_footprint;
#   real total_feed_nitrogen_footprint;
#   // Phosphorus
#   real<lower=0> species_phosphorus_footprint;
#   vector[K] feed_phosphorus_footprint;
#   real total_feed_phosphorus_footprint;
#   // Land
#   real<lower=0> species_land_footprint;
#   vector[K] feed_land_footprint;
#   real total_feed_land_footprint;
#   // Water
#   real<lower=0> species_water_footprint;
#   vector[K] feed_water_footprint;
#   real total_feed_water_footprint;
#   
#   // Calculations
#   feed_carbon_footprint = fp_carbon_dat .* theta;
#   total_feed_carbon_footprint = sum(feed_carbon_footprint);
#   species_carbon_footprint = mu * total_feed_carbon_footprint;
# 
#   feed_nitrogen_footprint = fp_nitrogen_dat .* theta;
#   total_feed_nitrogen_footprint = sum(feed_nitrogen_footprint);
#   species_nitrogen_footprint = mu * total_feed_nitrogen_footprint;
# 
#   feed_phosphorus_footprint = fp_phosphorus_dat .* theta;
#   total_feed_phosphorus_footprint = sum(feed_phosphorus_footprint);
#   species_phosphorus_footprint = mu * total_feed_phosphorus_footprint;
#   
#   feed_land_footprint = fp_land_dat .* theta;
#   total_feed_land_footprint = sum(feed_land_footprint);
#   species_land_footprint = mu * total_feed_land_footprint;
#   
#   feed_water_footprint = fp_water_dat .* theta;
#   total_feed_water_footprint = sum(feed_water_footprint);
#   species_water_footprint = mu * total_feed_water_footprint;
# }'

no_na_mod <- stan_model(model_code = stan_no_na, verbose = TRUE)
# Note: For Windows, apparently OK to ignore this warning message:
# Warning message:
#   In system(paste(CXX, ARGS), ignore.stdout = TRUE, ignore.stderr = TRUE) :
#   'C:/rtools40/usr/mingw_/bin/g++' not found

# Fit model:
# Set seed while testing
fit_no_na <- sampling(object = no_na_mod, data = stan_data, cores = 4, seed = "11729")
                       #,iter = 10000, cores = 4,
                       #control = list(adapt_delta = 0.99))
print(fit_no_na)

# Check diagnostics: https://betanalpha.github.io/assets/case_studies/rstan_workflow.html
# check_n_eff(fit_no_na)
# check_rhat(fit_no_na)
# check_treedepth(fit_no_na) # if needed, re-run with larger treedepth
# check_div(fit_no_na) # if needed, re-run with larger adapt
check_all_diagnostics(fit_no_na)

launch_shinystan(fit_no_na)


######################################################################################################
# Model 1: Remove all NAs - estimate feed footprint for Oncorhynchus mykiss

# Remove NAs
lca_dat_no_na <- lca_dat_no_zeroes %>%
  filter(clean_sci_name == "Oncorhynchus mykiss") %>%
  filter(is.na(Feed_soy_percent)==FALSE) 


# BOX PLOTS OF DATA:

# Theme for ALL PLOTS (including mcmc plots)
plot_theme <- theme(axis.text=element_text(size=14, color = "black"))

# FCR:
plot_fcr <- lca_dat_no_na %>%
  select(clean_sci_name, FCR) %>%
  mutate(clean_sci_name = paste(clean_sci_name, row_number(), sep = "")) %>%
  pivot_longer(cols = FCR)
ggplot(data = plot_fcr, aes(x = name, y = value)) +
  geom_boxplot() +
  theme_classic() +
  plot_theme +
  labs(title = "Boxplots of FCRs for Oncorhynchus mykiss",
       x = "",
       y = "")
ggsave(file.path(outdir, "boxplot_fcr_trout.png"), height = 8, width = 11.5)


# feed proportion:
plot_feed_prop <- lca_dat_no_na %>%
  select(clean_sci_name, soy = feed_soy_new, crops = feed_crops_new, fmfo = feed_fmfo_new, animal = feed_animal_new) %>%
  mutate(clean_sci_name = paste(clean_sci_name, row_number(), sep = "")) %>%
  pivot_longer(cols = soy:animal)

ggplot(data = plot_feed_prop, aes(x = name, y = value)) +
  geom_boxplot() +
  theme_classic() +
  plot_theme +
  labs(title = "Boxplots of feed proportions for Oncorhynchus mykiss",
       x = "",
       y = "")
ggsave(file.path(outdir, "boxplot_feed-prop_trout.png"), height = 8, width = 11.5)

# Set data for model
# for FCR model:
x <- lca_dat_no_na$FCR

# for Feed proportion model:
k = 4
n = 3
feed_weights <- lca_dat_no_na %>%
  select(contains("new")) %>%
  as.matrix()

# for final foot print calculation
fp_dat <- fp_clean %>%
  filter(Category != "Energy") %>%
  select(FP, Category, FP_val) 

# ORDER: Animal, crop, FMFO, soy
fp_carbon_dat <- fp_dat %>%
  filter(FP == "Carbon") %>%
  select(-FP) %>%
  pivot_wider(names_from = Category, values_from = FP_val) %>%
  as.matrix() %>%
  c()

fp_nitrogen_dat <- fp_dat %>%
  filter(FP == "Nitrogen") %>%
  select(-FP) %>%
  pivot_wider(names_from = Category, values_from = FP_val) %>%
  as.matrix() %>%
  c()

fp_phosphorus_dat <- fp_dat %>%
  filter(FP == "Phosphorus") %>%
  select(-FP) %>%
  pivot_wider(names_from = Category, values_from = FP_val) %>%
  as.matrix() %>%
  c()

fp_land_dat <- fp_dat %>%
  filter(FP == "Land") %>%
  select(-FP) %>%
  pivot_wider(names_from = Category, values_from = FP_val) %>%
  as.matrix() %>%
  c()

fp_water_dat <- fp_dat %>%
  filter(FP == "Water") %>%
  select(-FP) %>%
  pivot_wider(names_from = Category, values_from = FP_val) %>%
  as.matrix() %>%
  c()

# note: dirichlet_rng is just a random number generator:
# rep_vector(x, m) creates a column consisting of m copies of x
# generated quantities {
#   vector[k] theta = dirichlet_rng(rep_vector(alpha, k));
# }

# Estimate feed component proportions for a single species
stan_pooled <- 'data {
  int<lower=0> n;  // number of observations
  vector[n] x; // data
  int<lower=1> k; // number of feed types
  simplex[k] feed_weights[n]; // array of feed weights simplexes
  vector[k] fp_carbon_dat;
  vector[k] fp_nitrogen_dat;
  vector[k] fp_phosphorus_dat;
  vector[k] fp_land_dat;
  vector[k] fp_water_dat;
}
parameters {
  // FCR model:
  real<lower=0> mu;
  real<lower=0> sigma;
  // Feed proportion model:
  vector<lower=0>[k] alpha;
  simplex[k] theta;
}
model {
  
  x ~ normal(mu, sigma); // note: stan interprets second param as standard deviation

  for (i in 1:n) {
    feed_weights[i] ~ dirichlet(alpha); // estimate vector of alphas based on the data of feed weights
  }
  theta ~ dirichlet(alpha); // now, estimate feed weights based on the vector of alphas
}
generated quantities {
  // Carbon
  real<lower=0> species_carbon_footprint;
  vector[k] feed_carbon_footprint;
  real total_feed_carbon_footprint;
  // Nitrogen
  real<lower=0> species_nitrogen_footprint;
  vector[k] feed_nitrogen_footprint;
  real total_feed_nitrogen_footprint;
  // Phosphorus
  real<lower=0> species_phosphorus_footprint;
  vector[k] feed_phosphorus_footprint;
  real total_feed_phosphorus_footprint;
  // Land
  real<lower=0> species_land_footprint;
  vector[k] feed_land_footprint;
  real total_feed_land_footprint;
  // Water
  real<lower=0> species_water_footprint;
  vector[k] feed_water_footprint;
  real total_feed_water_footprint;
  
  // Calculations
  feed_carbon_footprint = fp_carbon_dat .* theta;
  total_feed_carbon_footprint = sum(feed_carbon_footprint);
  species_carbon_footprint = mu * total_feed_carbon_footprint;

  feed_nitrogen_footprint = fp_nitrogen_dat .* theta;
  total_feed_nitrogen_footprint = sum(feed_nitrogen_footprint);
  species_nitrogen_footprint = mu * total_feed_nitrogen_footprint;

  feed_phosphorus_footprint = fp_phosphorus_dat .* theta;
  total_feed_phosphorus_footprint = sum(feed_phosphorus_footprint);
  species_phosphorus_footprint = mu * total_feed_phosphorus_footprint;
  
  feed_land_footprint = fp_land_dat .* theta;
  total_feed_land_footprint = sum(feed_land_footprint);
  species_land_footprint = mu * total_feed_land_footprint;
  
  feed_water_footprint = fp_water_dat .* theta;
  total_feed_water_footprint = sum(feed_water_footprint);
  species_water_footprint = mu * total_feed_water_footprint;
}'

no_missing_mod <- stan_model(model_code = stan_pooled, verbose = TRUE)
# Note: For Windows, apparently OK to ignore this warning message:
# Warning message:
#   In system(paste(CXX, ARGS), ignore.stdout = TRUE, ignore.stderr = TRUE) :
#   'C:/rtools40/usr/mingw_/bin/g++' not found

# Fit model:
fit_pooled <- sampling(object = no_missing_mod, data = list(n = n,
                                                            x = x,
                                                            k = k,
                                                            feed_weights = feed_weights,
                                                            fp_carbon_dat = fp_carbon_dat,
                                                            fp_nitrogen_dat = fp_nitrogen_dat,
                                                            fp_phosphorus_dat = fp_phosphorus_dat,
                                                            fp_land_dat = fp_land_dat,
                                                            fp_water_dat = fp_water_dat))#,
                       # iter = 10000, cores = 4,
                       # control = list(adapt_delta = 0.99))
print(fit_pooled)

launch_shinystan(fit_pooled)

feeds <- c("soy", "crops", "fmfo", "animal")
feed_key <- data.frame(carbon_footprint = paste("carbon_footprint[", feeds, "]", sep = ""),
                       nitrogen_footprint = paste("nitrogen_footprint[", feeds, "]", sep = ""),
                       phosphorus_footprint = paste("phosphorus_footprint[", feeds, "]", sep = ""),
                       land_footprint = paste("land_footprint[", feeds, "]", sep = ""),
                       water_footprint = paste("water_footprint[", feeds, "]", sep = ""))

fit_pooled_clean <- fit_pooled
names(fit_pooled_clean)[grep(names(fit_pooled_clean), pattern = "feed_carbon_footprint\\[")] <- feed_key$carbon_footprint
names(fit_pooled_clean)[grep(names(fit_pooled_clean), pattern = "feed_nitrogen_footprint\\[")] <- feed_key$nitrogen_footprint
names(fit_pooled_clean)[grep(names(fit_pooled_clean), pattern = "feed_phosphorus_footprint\\[")] <- feed_key$phosphorus_footprint
names(fit_pooled_clean)[grep(names(fit_pooled_clean), pattern = "feed_land_footprint\\[")] <- feed_key$land_footprint
names(fit_pooled_clean)[grep(names(fit_pooled_clean), pattern = "feed_water_footprint\\[")] <- feed_key$water_footprint

distribution_pooled <- as.matrix(fit_pooled_clean)

# FIX IT - replace parameter names and add plots for other parameters
p_footprint <- mcmc_areas(distribution_pooled,
                          pars = vars(contains("carbon_footprint")),
                          prob = 0.8,
                          prob_outer = 0.9,
                          area_method = "scaled height",
                          point_est = "median") + 
  ggtitle("Oncorhynchus mykiss full feed footprint model", "with 80% credible intervals") +
  xlim(0, 10) +
  plot_theme 

p_footprint
ggsave(filename = file.path(outdir, "bayes-example_trout_carbon-feed-footprint.png"), width = 11, height = 8.5)


p_footprint <- mcmc_areas(distribution_pooled,
                          pars = vars(contains("nitrogen_footprint")),
                          prob = 0.8,
                          prob_outer = 0.9,
                          area_method = "scaled height",
                          point_est = "median") + 
  ggtitle("Oncorhynchus mykiss full feed footprint model", "with 80% credible intervals") +
  xlim(0, 0.01) +
  plot_theme 

p_footprint
ggsave(filename = file.path(outdir, "bayes-example_trout_nitrogen-feed-footprint.png"), width = 11, height = 8.5)

p_footprint <- mcmc_areas(distribution_pooled,
                          pars = vars(contains("phosphorus_footprint")),
                          prob = 0.8,
                          prob_outer = 0.9,
                          area_method = "scaled height",
                          point_est = "median") + 
  ggtitle("Oncorhynchus mykiss full feed footprint model", "with 80% credible intervals") +
  xlim(0, 0.001) +
  plot_theme 

p_footprint
ggsave(filename = file.path(outdir, "bayes-example_trout_phosphorus-feed-footprint.png"), width = 11, height = 8.5)

p_footprint <- mcmc_areas(distribution_pooled,
                          pars = vars(contains("land_footprint")),
                          prob = 0.8,
                          prob_outer = 0.9,
                          area_method = "scaled height",
                          point_est = "median") + 
  ggtitle("Oncorhynchus mykiss full feed footprint model", "with 80% credible intervals") +
  xlim(0, 5) +
  plot_theme 

p_footprint
ggsave(filename = file.path(outdir, "bayes-example_trout_land-feed-footprint.png"), width = 11, height = 8.5)

p_footprint <- mcmc_areas(distribution_pooled,
                          pars = vars(contains("water_footprint")),
                          prob = 0.8,
                          prob_outer = 0.9,
                          area_method = "scaled height",
                          point_est = "median") + 
  ggtitle("Oncorhynchus mykiss full feed footprint model", "with 80% credible intervals") +
  xlim(0, 0.2) +
  plot_theme 

p_footprint
ggsave(filename = file.path(outdir, "bayes-example_trout_water-feed-footprint.png"), width = 11, height = 8.5)
