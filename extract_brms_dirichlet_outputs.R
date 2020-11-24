# Create visualizations of brms dirichlet regression outputs
# Do not clear workspace at start

# Mac
datadir <- "/Volumes/jgephart/BFA Environment 2/Data"
outdir <- "/Volumes/jgephart/BFA Environment 2/Outputs"

library(tidyverse)
library(modelr)
library(ggdist)
library(tidybayes)
library(ggplot2)
library(cowplot)
library(rstan)
library(brms)
library(ggrepel)
library(RColorBrewer)

# After running brms regressions (e.g., bayes_gamma_regression.R)

# identify model output to be used for visualizations
name_of_fit <- "fit_no_na" 
brms_output <- get(name_of_fit)

summary(brms_output)

## PP_CHECK not implemented for dirichlet regression models
# pp_check(brms_output, nsamples = 50) + 
#   ggtitle("Posterior predictive check")
# ggsave(filename = file.path(outdir, "plot_gamma-regression_post-pred-checks_density.png"), width  = 11.5, height = 8)
# pp_check(brms_output, type = "error_hist", nsamples = 5)
# pp_check(brms_output, type = "scatter_avg", nsamples = 1000)
# pp_check(brms_output, type = "stat_2d")
# pp_check(brms_output, type = "stat")

get_variables(brms_output)

# Plot coefficients (separate for each feed component)
feed_coeffs <- c("b_mufeedcrops", "b_mufeedfmfo", "b_mufeedanimal")
for (i in 1:length(feed_coeffs)){
  p <- mcmc_plot(brms_output, pars = feed_coeffs[i])
  coeff_data <- p$data %>% 
    mutate(effect = case_when(m > 0 ~ "positive",
                              m < 0 ~ "negative",
                              TRUE ~ "none"))
  p_custom <- ggplot(data = coeff_data) +
    geom_segment(data = coeff_data, mapping = aes(xend = hh, yend = parameter, x = ll, y = parameter)) +
    geom_segment(data = coeff_data, mapping = aes(xend = h, yend = parameter, x = l, y = parameter), size = 2) +
    geom_point(aes(x = m, y = parameter, color = effect), size = 3) +
    labs(x = "", y = "") +
    theme_classic() +
    theme(axis.text = element_text(size = 16))
  plot(p_custom)
  
  ggsave(filename = file.path(outdir, paste("plot_dirichlet-regression_coeffs_", feed_coeffs[i], ".png", sep = "")), height = 8.5, width = 11)
}


# Get point and interval estimates from predicted data
# Select jsut the prediction columns
# Join these with the original lca data (lca_complete_predictors) to get metadata on taxa/intensity/syste,
feed_dat_intervals <- predicted_feed_dat %>%
  median_qi(.value = .prediction) %>% # Rename prediction to value
  ungroup() %>%
  select(contains("."))

# .row is equivalent to the row number in the original dataset (lca_complete_predictors) - create a join column for this
predicted_metadat<- lca_complete_predictors %>%
  select(clean_sci_name, taxa, intensity, system) %>%
  mutate(.row = row_number())

na_predictions <- feed_dat_intervals %>%
  left_join(predicted_metadat, by = ".row")
  
# Reformat data so it can row-bind with predictions
lca_no_na_long <- lca_no_na %>%
  pivot_longer(cols = c("feed_soy", "feed_crops", "feed_fmfo", "feed_animal"), names_to = ".category", values_to = ".value")

full_dat <- na_predictions %>%
  bind_rows(lca_no_na_long)

# ONLY MODELED USING INTENSITY, but display both?
# Make separate plot of data and predictions for each taxa group
# feed_cols <- c("feed_soy", "feed_crops", "feed_fmfo", "feed_animal")
for (i in 1:length(unique(full_dat$taxa))){
  taxa_i <- unique(full_dat$taxa)[i]
  #feed_j <- feed_cols[j]
  dat_taxa_i <- full_dat %>%
    filter(taxa == taxa_i) %>%
    #filter(.category == feed_j) %>%
    mutate(y = row_number()) %>%
    replace_na(replace = list(taxa = "unknown", intensity = "unknown", system = "unknown"))
  p <- ggplot() +
    geom_pointinterval(aes(y = y, x = .value, xmin = .lower, xmax = .upper, shape = system, point_color = intensity), size = 2, data = dat_taxa_i) +
    geom_point(aes(y = y, x = .value, shape = system, color = intensity), data = dat_taxa_i, size = 3) +
    #scale_interval_shape(drop = FALSE) +
    scale_shape_discrete(drop = FALSE) +
    #coord_cartesian(xlim = c(0, 10)) +
    theme_classic() +
    labs(x = "FCR", y = "", title = taxa_i) +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          axis.text.y = element_blank()) +
    facet_wrap(~.category, ncol = 1)
  plot(p)
  file_i <- paste("plot_dirichlet-regression_missing-dat-predictions_taxa-", taxa_i, ".png", sep = "")
  ggsave(filename = file.path(outdir, file_i), width = 8, height = 11.5)
}
