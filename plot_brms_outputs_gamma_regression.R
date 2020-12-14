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

# identify model and data output and data variable to be used for visualizations

# FCR
name_of_fit <- "fit_fcr_no_na"
name_of_data <- "full_fcr_dat"
name_of_var <- "fcr"

# Yield
# name_of_fit <- "fit_yield_no_na"
# name_of_data <- "full_yield_dat"
# name_of_var <- "yield"

# Electricity
# name_of_fit <- "fit_electric_no_na"
# name_of_data <- "full_electric_dat"
# name_of_var <- "electric"

# Diesel
# name_of_fit <- "fit_diesel_no_na"
# name_of_data <- "full_diesel_dat"
# name_of_var <- "diesel"

# Petrol
# name_of_fit <- "fit_petrol_no_na"
# name_of_data <- "full_petrol_dat"
# name_of_var <- "petrol"

# Natural gas
# name_of_fit <- "fit_natgas_no_na"
# name_of_data <- "full_natgas_dat"
# name_of_var <- "natgas"

compiled_dat_clean <- read.csv(file.path(datadir, "lca_clean_with_groups.csv"))
brms_output <- get(name_of_fit)
full_dat <- get(name_of_data)

# Combine modeled data (both data and predictions) with the full clean LCA dataset and output this
dat_for_si <- compiled_dat_clean %>%
  left_join(full_dat, by = c("study_id", "clean_sci_name", "taxa", "Intensity" = "intensity", "Production_system_group" = "system"))
write.csv(dat_for_si, file.path(outdir, paste("lca_clean_with_model_predictions-", name_of_var, ".csv", sep = "")), row.names = FALSE)

summary(brms_output)

# Default posterior predictive check is a density plot:
# Specify response variable in resp
pp_check(brms_output, resp = "y", nsamples = 50) + 
ggtitle(paste("Posterior predictive check: ", name_of_var, sep = ""))
ggsave(filename = file.path(outdir, paste("plot_gamma-regression_", name_of_var, "_post-pred-checks_density.png", sep = "")), width  = 11.5, height = 8)
pp_check(brms_output, type = "scatter_avg", nsamples = 1000, resp = "y") +
ggtitle(paste("Posterior predictive check: ", name_of_var, sep = ""))
ggsave(filename = file.path(outdir, paste("plot_gamma-regression_", name_of_var, "_post-pred-checks_scatter.png", sep = "")), width  = 11.5, height = 8)

# Other posterior predictive checks
# pp_check(brms_output, type = "error_hist", nsamples = 5, resp = "y")
# pp_check(brms_output, type = "stat_2d", resp = "y")
# pp_check(brms_output, type = "stat", resp = "y")

get_variables(brms_output)

# Plot coefficients
p <- mcmc_plot(brms_output, pars = c("^b_"))

coeff_data <- p$data %>% 
  mutate(effect = case_when(m > 0 ~ "positive",
                            m < 0 ~ "negative",
                            TRUE ~ "none"))

ggplot(data = coeff_data) +
  geom_segment(data = coeff_data, mapping = aes(xend = hh, yend = parameter, x = ll, y = parameter)) +
  geom_segment(data = coeff_data, mapping = aes(xend = h, yend = parameter, x = l, y = parameter), size = 2) +
  geom_point(aes(x = m, y = parameter, color = effect), size = 3) +
  labs(x = "", y = "", title = paste("Coefficients for ", name_of_var, sep = "")) +
  theme_classic() +
  theme(axis.text = element_text(size = 16))

ggsave(filename = file.path(outdir, paste("plot_gamma-regression_", name_of_var, "_coeffs.png", sep = "")), height = 8.5, width = 11)
  

# plot predictions of missing responses
p <- ggplot() +
  geom_pointinterval(aes(y = rowname, x = !!sym(name_of_var), xmin = .lower, xmax = .upper), size = 0.5, data = full_dat) +
  geom_point(aes(y = rowname, x = !!sym(name_of_var), color = data_type), data = full_dat) +
  #coord_cartesian(xlim = c(0, 10)) +
  theme_classic() +
  labs(x = name_of_var, y = "") +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 20))
plot(p)
ggsave(filename = file.path(outdir, paste("plot_gamma-regression_", name_of_var, "_missing-dat-predictions.png", sep = "")), width = 8, height = 11.5)

# Reorder by response value
p <- ggplot(data = full_dat %>% 
              mutate(rowname = as.factor(rowname)) %>%
              mutate(rowname = fct_reorder(rowname, !!sym(name_of_var))) %>%
              arrange(!!sym(name_of_var))) +
  geom_pointinterval(aes(y = rowname, x = !!sym(name_of_var), xmin = .lower, xmax = .upper), size = 0.5) +
  geom_point(aes(y = rowname, x = !!sym(name_of_var), color = data_type)) +
  #coord_cartesian(xlim = c(0, 10)) +
  theme_classic() +
  labs(x = name_of_var, y = "") +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 20),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())
plot(p)
ggsave(filename = file.path(outdir, paste("plot_gamma-regression_", name_of_var, "_missing-dat-predictions-ordered.png", sep = "")), width = 8, height = 11.5)

# plot predictions of missing responses with meta-data
p <- ggplot(data = full_dat %>% 
              mutate(rowname = as.factor(rowname)) %>%
              mutate(rowname = fct_reorder(rowname, !!sym(name_of_var))) %>%
              arrange(!!sym(name_of_var))) +
  geom_pointinterval(aes(y = rowname, x = !!sym(name_of_var), xmin = .lower, xmax = .upper, shape = system, point_color = intensity), size = 0.5) +
  geom_point(aes(y = rowname, x = !!sym(name_of_var), shape = system, color = intensity)) +
  #coord_cartesian(xlim = c(0, 10)) +
  theme_classic() +
  labs(x = name_of_var, y = "") +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 20),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())
plot(p)
file_i <- paste("plot_gamma-regression_", name_of_var, "_missing-dat-predictions-ordered-with-metadata", ".png", sep = "")
ggsave(filename = file.path(outdir, file_i), width = 8, height = 11.5)


# Make separate plot for each taxa group
for (i in 1:length(unique(full_dat$taxa))){
  taxa_i <- unique(full_dat$taxa)[i]
  dat_taxa_i <- full_dat %>%
    filter(taxa == taxa_i) %>%
    mutate(rowname = row_number()) %>%
    replace_na(replace = list(taxa = "unknown", intensity = "unknown", system = "unknown"))
  p <- ggplot() +
    geom_pointinterval(aes(y = rowname, x = !!sym(name_of_var), xmin = .lower, xmax = .upper, shape = system, point_color = intensity), size = 2, data = dat_taxa_i) +
    geom_point(aes(y = rowname, x = !!sym(name_of_var), shape = system, color = intensity), data = dat_taxa_i, size = 3) +
    #scale_interval_shape(drop = FALSE) +
    scale_shape_discrete(drop = FALSE) +
    #coord_cartesian(xlim = c(0, 10)) +
    theme_classic() +
    labs(x = name_of_var, y = "", title = taxa_i) +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          axis.text.y = element_blank())
  plot(p)
  file_i <- paste("plot_gamma-regression_", name_of_var, "_missing-dat-predictions_taxa-", taxa_i, ".png", sep = "")
  ggsave(filename = file.path(outdir, file_i), width = 11.5, height = 8)
}


