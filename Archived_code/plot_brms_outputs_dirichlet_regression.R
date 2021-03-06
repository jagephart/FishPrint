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

# identify model and data output and data variable to be used for visualizations
name_of_fit <- "fit_feed_no_na" 
name_of_data <- "full_feed_dat"
name_of_var <- "feed_proportion"

compiled_dat_clean <- read.csv(file.path(datadir, "lca_clean_with_groups.csv"))
brms_output <- get(name_of_fit)
full_dat <- get(name_of_data)

# FIX IT - dropping .lower and upper bounds for now for easier merge
full_dat_for_merge <- full_dat %>%
  select(-c(.row, rowname, .lower, .upper, .width, .point, .interval)) %>% 
  pivot_wider(names_from = .category, values_from = feed_proportion) %>%
  group_by(study_id, clean_sci_name, taxa, intensity, system, data_type) %>%
  # SUM to collapse into a single row
  mutate(feed_soy = sum(feed_soy, na.rm = TRUE),
         feed_crops = sum(feed_crops, na.rm = TRUE),
         feed_fmfo = sum(feed_fmfo, na.rm = TRUE),
         feed_animal = sum(feed_animal, na.rm = TRUE)) %>%
  distinct()

# Combine modeled data (both data and predictions) with the full clean LCA dataset and output this
dat_for_si <- compiled_dat_clean %>%
  left_join(full_dat_for_merge, by = c("study_id", "clean_sci_name", "taxa", "Intensity" = "intensity", "Production_system_group" = "system"))
write.csv(dat_for_si, file.path(outdir, paste("lca_clean_with_model_predictions-", name_of_var, ".csv", sep = "")), row.names = FALSE)

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


# SET THEME
plot_theme <- theme(title = element_text(size = 20),
                    axis.title.x = element_text(size = 20),
                    axis.text=element_text(size=20, color = "black"))

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
    plot_theme
  plot(p_custom)
  
  ggsave(filename = file.path(outdir, paste("plot_dirichlet-regression_coeffs_", feed_coeffs[i], ".png", sep = "")), height = 8.5, width = 11)
}


# feed proportion all taxa
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
  ggsave(file.path(outdir, paste("boxplot_feed-prop_", feed_vars[i], ".png", sep = "")), height = 8, width = 11.5)
}

# Feed proportions with facet_wrap
p <- ggplot(data = full_feed_dat, aes(y = taxa, x = feed_proportion)) +
  geom_boxplot(outlier.shape = NA) +
  #geom_violin(aes(color = taxa), scale = "width") +
  geom_jitter(aes(color = data_type), size = 3) +
  theme_classic() +
  plot_theme +
  labs(title = paste("Boxplots of all feed proportions", sep = ""),
       x = "",
       y = "")  +
  facet_wrap(~.category, nrow = 1)
print(p)
ggsave(file.path(outdir, "boxplot_feed-prop_all-feeds.png"), height = 8, width = 11.5)


# Make separate plot of data and predictions for each taxa group
# feed_cols <- c("feed_soy", "feed_crops", "feed_fmfo", "feed_animal")
# for (i in 1:length(unique(full_dat$taxa))){
#   taxa_i <- unique(full_dat$taxa)[i]
#   #feed_j <- feed_cols[j]
#   dat_taxa_i <- full_dat %>%
#     filter(taxa == taxa_i) %>%
#     #filter(.category == feed_j) %>%
#     mutate(y = row_number()) %>%
#     replace_na(replace = list(taxa = "unknown", intensity = "unknown", system = "unknown"))
#   p <- ggplot() +
#     geom_pointinterval(aes(y = y, x = feed_proportion, xmin = .lower, xmax = .upper, shape = system, point_color = intensity), size = 2, data = dat_taxa_i) +
#     geom_point(aes(y = y, x = feed_proportion, shape = system, color = intensity), data = dat_taxa_i, size = 3) +
#     #scale_interval_shape(drop = FALSE) +
#     scale_shape_discrete(drop = FALSE) +
#     #coord_cartesian(xlim = c(0, 10)) +
#     theme_classic() +
#     labs(x = "Feed proportion", y = "", title = taxa_i) +
#     theme(axis.text = element_text(size = 16),
#           axis.title = element_text(size = 20),
#           axis.text.y = element_blank()) +
#     facet_wrap(~.category, ncol = 1)
#   plot(p)
#   file_i <- paste("plot_dirichlet-regression_missing-dat-predictions_taxa-", taxa_i, ".png", sep = "")
#   ggsave(filename = file.path(outdir, file_i), width = 8, height = 11.5)
# }
