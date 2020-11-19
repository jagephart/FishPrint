# Create visualizations of brms outputs
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
library(brmstools)

# After running brms regressions (e.g., bayes_gamma_regression.R)

# identify model output to be used for visualizations
name_of_fit <- "fit_with_na" 
brms_output <- get(name_of_fit)

summary(brms_output)

pp_check(brms_output, resp = "y", nsamples = 50) + # density plots
ggtitle("Posterior predictive check")
ggsave(filename = file.path(outdir, "plot_gamma-regression_post-pred-checks.png"), width  = 11.5, height = 8)

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
  
# missing responses are predicted by the model and listed as part of the outputed variables
na_predictions <- brms_output %>%
  spread_draws(Ymi_y[y_na]) %>%
  median_qi() # median values of all missing predicted responses

#Get all the non-NAs in brms_data
brms_non_na <- brms_data %>%
  mutate(y_na = row_number()) %>%
  filter(is.na(y)==FALSE) %>%
  select(Ymi_y = y, y_na)

# combine data with predictions
dat_combine <- na_predictions %>%
  bind_rows(brms_non_na) %>%
  arrange(y_na) %>%
  mutate(data_type = if_else(is.na(.point), true = "data", false = "prediction"))

# combine with original lca data to get meta data
full_dat <- dat_combine %>%
  left_join(lca_with_na %>% mutate(y_na = row_number()), by = "y_na")

# plot predictions of missing responses
p <- ggplot() +
  geom_pointinterval(aes(y = y_na, x = Ymi_y, xmin = .lower, xmax = .upper, color = !!sym(study_dat_p)), size = 0.5, data = full_dat) +
  geom_point(aes(y = y_na, x = Ymi_y, color = data), data = full_dat) +
  coord_cartesian(xlim = c(0, 10)) +
  theme_classic() +
  labs(x = "FCR", y = "LCA study") +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 20))
plot(p)
ggsave(filename = file.path(outdir, "plot_gamma-regression_missing-dat-predictions.png"), width = 8, height = 11.5)

# Reorder by response value
p <- ggplot(data = full_dat %>% 
              mutate(y_na = as.factor(y_na)) %>%
              mutate(y_na = fct_reorder(y_na, Ymi_y)) %>%
              arrange(Ymi_y)) +
  geom_pointinterval(aes(y = y_na, x = Ymi_y, xmin = .lower, xmax = .upper), size = 0.5) +
  geom_point(aes(y = y_na, x = Ymi_y, color = data_type)) +
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
                mutate(y_na = as.factor(y_na)) %>%
                mutate(y_na = fct_reorder(y_na, Ymi_y)) %>%
                arrange(Ymi_y)) +
    geom_pointinterval(aes(y = y_na, x = Ymi_y, xmin = .lower, xmax = .upper), size = 0.5) +
    geom_point(aes(y = y_na, x = Ymi_y, color = !!sym(study_dat_i))) +
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
    mutate(y_na = row_number()) %>%
    replace_na(replace = list(taxa = "unknown", intensity = "unknown", system = "unknown"))
  p <- ggplot() +
    geom_pointinterval(aes(y = y_na, x = Ymi_y, xmin = .lower, xmax = .upper, interval_color = system), size = 2, data = dat_taxa_i) +
    geom_point(aes(y = y_na, x = Ymi_y, color = intensity), data = dat_taxa_i, size = 3) +
    scale_interval_color_discrete(drop = FALSE) +
    scale_color_discrete(drop = FALSE) +
    coord_cartesian(xlim = c(0, 10)) +
    theme_classic() +
    labs(x = "FCR", y = "", title = taxa_i) +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 20))
  plot(p)
  file_i <- paste("plot_gamma-regression_missing-dat-predictions_", taxa_i, ".png", sep = "")
  ggsave(filename = file.path(outdir, file_i), width = 11.5, height = 8)
}


