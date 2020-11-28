# Create visualizations of brms gamma regression outputs
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
name_of_fit <- "fit_fcr_with_na" 
brms_output <- get(name_of_fit)

summary(brms_output)

# Default posterior predictive check is a density plot:
# Specify response variable in resp
pp_check(brms_output, resp = "y", nsamples = 50) + 
ggtitle("Posterior predictive check")
ggsave(filename = file.path(outdir, "plot_gamma-regression_post-pred-checks_density.png"), width  = 11.5, height = 8)

# Other posterior predictive checks
pp_check(brms_output, type = "error_hist", nsamples = 5, resp = "y")
pp_check(brms_output, type = "scatter_avg", nsamples = 1000, resp = "y")
ggsave(filename = file.path(outdir, "plot_gamma-regression_post-pred-checks_scatter.png"), width  = 11.5, height = 8)
pp_check(brms_output, type = "stat_2d", resp = "y")
pp_check(brms_output, type = "stat", resp = "y")

get_variables(brms_output)

# Plot coefficients
p <- mcmc_plot(brms_output, pars = c("^b_y", "^bsp_y"))

coeff_data <- p$data %>% 
  mutate(effect = case_when(m > 0 ~ "positive",
                            m < 0 ~ "negative",
                            TRUE ~ "none"))

ggplot(data = coeff_data) +
  geom_segment(data = coeff_data, mapping = aes(xend = hh, yend = parameter, x = ll, y = parameter)) +
  geom_segment(data = coeff_data, mapping = aes(xend = h, yend = parameter, x = l, y = parameter), size = 2) +
  geom_point(aes(x = m, y = parameter, color = effect), size = 3) +
  labs(x = "", y = "") +
  theme_classic() +
  theme(axis.text = element_text(size = 16))

ggsave(filename = file.path(outdir, "plot_gamma-regression_coeffs.png"), height = 8.5, width = 11)
  
# Missing FCRs are predicted by the model and listed as part of the outputed variables
na_predictions <- fit_fcr_with_na %>%
  spread_draws(Ymi_y[y_row]) %>%
  median_qi() # median values of all missing predicted responses

#Get all the non-NAs in gamma_brms_data
brms_non_na <- gamma_brms_data %>%
  mutate(y_row = row_number()) %>%
  filter(is.na(y)==FALSE) %>%
  select(Ymi_y = y, y_row)

# combine FCR data with predictions of missing FCRs
dat_combine <- na_predictions %>%
  bind_rows(brms_non_na) %>%
  arrange(y_row) %>%
  mutate(data_type = if_else(is.na(.point), true = "data", false = "prediction"))

# combine with original lca data to get meta data
full_dat <- dat_combine %>%
  left_join(lca_with_na %>% mutate(y_row = row_number()), by = "y_row")

# plot predictions of missing responses
p <- ggplot() +
  geom_pointinterval(aes(y = y_row, x = Ymi_y, xmin = .lower, xmax = .upper), size = 0.5, data = full_dat) +
  geom_point(aes(y = y_row, x = Ymi_y, color = data_type), data = full_dat) +
  coord_cartesian(xlim = c(0, 10)) +
  theme_classic() +
  labs(x = "FCR", y = "LCA study") +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 20))
plot(p)
ggsave(filename = file.path(outdir, "plot_gamma-regression_missing-dat-predictions.png"), width = 8, height = 11.5)

# Reorder by response value
p <- ggplot(data = full_dat %>% 
              mutate(y_row = as.factor(y_row)) %>%
              mutate(y_row = fct_reorder(y_row, Ymi_y)) %>%
              arrange(Ymi_y)) +
  geom_pointinterval(aes(y = y_row, x = Ymi_y, xmin = .lower, xmax = .upper), size = 0.5) +
  geom_point(aes(y = y_row, x = Ymi_y, color = data_type)) +
  coord_cartesian(xlim = c(0, 10)) +
  theme_classic() +
  labs(x = "FCR", y = "LCA study") +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 20),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())
plot(p)
ggsave(filename = file.path(outdir, "plot_gamma-regression_missing-dat-predictions-ordered.png"), width = 8, height = 11.5)

# plot predictions of missing responses with meta-data
study_dat <- c("taxa", "intensity", "system")
for (i in 1:length(study_dat)){
  study_dat_i <- study_dat[i]
  p <- ggplot(data = full_dat %>% 
                mutate(y_row = as.factor(y_row)) %>%
                mutate(y_row = fct_reorder(y_row, Ymi_y)) %>%
                arrange(Ymi_y)) +
    geom_pointinterval(aes(y = y_row, x = Ymi_y, xmin = .lower, xmax = .upper), size = 0.5) +
    geom_point(aes(y = y_row, x = Ymi_y, color = !!sym(study_dat_i))) +
    coord_cartesian(xlim = c(0, 10)) +
    theme_classic() +
    labs(x = "FCR", y = "LCA study") +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank())
  plot(p)
  file_i <- paste("plot_gamma-regression_missing-dat-predictions_", study_dat_i, ".png", sep = "")
  ggsave(filename = file.path(outdir, file_i), width = 8, height = 11.5)
}

# Make separate plot for each taxa group
for (i in 1:length(unique(full_dat$taxa))){
  taxa_i <- unique(full_dat$taxa)[i]
  dat_taxa_i <- full_dat %>%
    filter(taxa == taxa_i) %>%
    mutate(y_row = row_number()) %>%
    replace_na(replace = list(taxa = "unknown", intensity = "unknown", system = "unknown"))
  p <- ggplot() +
    geom_pointinterval(aes(y = y_row, x = Ymi_y, xmin = .lower, xmax = .upper, shape = system, point_color = intensity), size = 2, data = dat_taxa_i) +
    geom_point(aes(y = y_row, x = Ymi_y, shape = system, color = intensity), data = dat_taxa_i, size = 3) +
    #scale_interval_shape(drop = FALSE) +
    scale_shape_discrete(drop = FALSE) +
    coord_cartesian(xlim = c(0, 10)) +
    theme_classic() +
    labs(x = "FCR", y = "", title = taxa_i) +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          axis.text.y = element_blank())
  plot(p)
  file_i <- paste("plot_gamma-regression_missing-dat-predictions_taxa-", taxa_i, ".png", sep = "")
  ggsave(filename = file.path(outdir, file_i), width = 11.5, height = 8)
}


