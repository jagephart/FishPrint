# Author: Kelvin Gorospe
# PLOT FIGURE 1 (AND OTHER SI FIGURES)
# Load model outputs individually, save to plot_grid, and create multi-panel plot
# Remove y-axis labels for panels B, C, and F
# Set plot margins for all plots to 0
# Set axis label margins for panels C and F to 0

rm(list=ls())
library(tidyverse)
library(tidybayes)
library(ggplot2)
library(cowplot) ## for plot_grid
library(scales) ## for scales_x_continuous(labels = comma) # add comma for thousands separator
datadir <- "/Volumes/jgephart/BFA Environment 2/Data"
outdir <- "/Volumes/jgephart/BFA Environment 2/Outputs"

# Theme
# NEW tx_plot_theme_final for Nature submission - smaller x-axis font + hjust on x-axis
tx_plot_theme_final <- theme(axis.title = element_text(size = 8),
                       axis.text.x = element_text(hjust = 0.6, size = 7, color = "black"),
                       axis.text.y = element_text(size = 8, color = "black"),
                       legend.position = "none",
                       plot.margin = unit(c(0, 3, 0, 0), "mm")) # increase the right margin so that there can be some pillover from x-axis

###########################################################################################
# OPTIONS: Set filenames
# EDIBLE WEIGHT RESULTS

# NOTE: GHG and WATER results are in a different folder dated August 2021 because these models had to be re-run (last-minute change to natural gas constant and multiplying evaporation data by 12 to get annual value)
# All other filepaths for N, P, Land, Wild still point to their original folder from May 2021
# EDIBLE WEIGHT Mass allocation + FCR Priors
# c_results <- "Nature-submitted-2021-08/RERUN-GHG-and-Water-Bayesian-Means-Edible-Weight/PRIORS/2021-08-17_full-model-posterior_Global warming potential_Mass-allocation.RData"
# n_results <- "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Nitrogen/2021-04-28_full-model-posterior_Marine eutrophication_Mass-allocation.RData"
# p_results <- "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Phosphorus/2021-04-28_full-model-posterior_Freshwater eutrophication_Mass-allocation.RData"
# land_results <- "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Land-FCR-Priors-Only/2021-04-28_full-model-posterior_Land Use_Mass-allocation.RData"
# water_results <- "Nature-submitted-2021-08/RERUN-GHG-and-Water-Bayesian-Means-Edible-Weight/PRIORS/2021-08-17_full-model-posterior_Water Consumption_Mass-allocation.RData"
# wild_results <- "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Wild/2021-04-27_full-model-posterior_Wild-capture-ghg.RData"

# EDIBLE WEIGHT Economic allocation (no wild capture results) + FCR Priors
# c_results <- "Nature-submitted-2021-08/RERUN-GHG-and-Water-Bayesian-Means-Edible-Weight/PRIORS/2021-08-18_full-model-posterior_Global warming potential_Economic-allocation.RData"
# n_results <- "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Nitrogen/2021-04-28_full-model-posterior_Marine eutrophication_Economic-allocation.RData"
# p_results <- "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Phosphorus/2021-04-28_full-model-posterior_Freshwater eutrophication_Economic-allocation.RData"
# land_results <- "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Land-FCR-Priors-Only/2021-04-30_full-model-posterior_Land Use_Economic-allocation.RData"
# water_results <- "Nature-submitted-2021-08/RERUN-GHG-and-Water-Bayesian-Means-Edible-Weight/PRIORS/2021-08-17_full-model-posterior_Water Consumption_Economic-allocation.RData"

# EDIBLE WEIGHT Gross energy allocation (no wild capture results) + FCR Priors
# c_results <- "Nature-submitted-2021-08/RERUN-GHG-and-Water-Bayesian-Means-Edible-Weight/PRIORS/2021-08-18_full-model-posterior_Global warming potential_Gross energy content-allocation.RData"
# n_results <- "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Nitrogen/2021-04-28_full-model-posterior_Marine eutrophication_Gross energy content-allocation.RData"
# p_results <- "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Phosphorus/2021-04-28_full-model-posterior_Freshwater eutrophication_Gross energy content-allocation.RData"
# land_results <- "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Land-FCR-Priors-Only/2021-04-30_full-model-posterior_Land Use_Gross energy content-allocation.RData"
# water_results <- "Nature-submitted-2021-08/RERUN-GHG-and-Water-Bayesian-Means-Edible-Weight/PRIORS/2021-08-18_full-model-posterior_Water Consumption_Gross energy content-allocation.RData"

# EDIBLE WEIGHT Mass allocation + NO Priors
# c_results <- "Nature-submitted-2021-08/RERUN-GHG-and-Water-Bayesian-Means-Edible-Weight/NO PRIORS/2021-08-18_full-model-posterior_Global warming potential_Mass-allocation.RData"
# n_results <- "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/NO PRIORS/Nitrogen/2021-04-28_full-model-posterior_Marine eutrophication_Mass-allocation.RData"
# p_results <- "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/NO PRIORS/Phosphorus/2021-04-28_full-model-posterior_Freshwater eutrophication_Mass-allocation.RData"
# land_results <- "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/NO PRIORS/Land/2021-04-29_full-model-posterior_Land Use_Mass-allocation.RData"
# water_results <- "Nature-submitted-2021-08/RERUN-GHG-and-Water-Bayesian-Means-Edible-Weight/NO PRIORS/2021-08-18_full-model-posterior_Water Consumption_Mass-allocation.RData"
# wild_results <- "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/NO PRIORS/Wild/2021-04-27_full-model-posterior_Wild-capture-ghg.RData"

########################
# LIVE WEIGHT RESULTS

# LIVE WEIGHT Mass allocation + FCR Priors
c_results <- "Nature-submitted-2021-08/RERUN-GHG-and-Water-Bayesian-Means-Live-Weight/PRIORS/2021-08-18_full-model-posterior_Global warming potential_Mass-allocation.RData"
n_results <- "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/Nitrogen/2021-01-15_full-model-posterior_Marine eutrophication_Mass-allocation.RData"
p_results <- "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/Phosphorus/2021-01-16_full-model-posterior_Freshwater eutrophication_Mass-allocation.RData"
land_results <- "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/Land-FCR-priors-only/2021-01-20_full-model-posterior_Land Use_Mass-allocation.RData"
water_results <- "Nature-submitted-2021-08/RERUN-GHG-and-Water-Bayesian-Means-Live-Weight/PRIORS/2021-08-18_full-model-posterior_Water Consumption_Mass-allocation.RData"
wild_results <- "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/Wild/2021-01-16_full-model-posterior_Wild-capture-ghg.RData"

# LIVE WEIGHT Economic allocation (no wild capture results) + FCR Priors
# c_results <- "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/GHG/2021-01-15_full-model-posterior_Global warming potential_Economic-allocation.RData"
# n_results <- "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/Nitrogen/2021-01-20_full-model-posterior_Marine eutrophication_Economic-allocation.RData"
# p_results <- "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/Phosphorus/2021-01-16_full-model-posterior_Freshwater eutrophication_Economic-allocation.RData"
# land_results <- "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/Land-FCR-priors-only/2021-01-16_full-model-posterior_Land Use_Economic-allocation.RData"
# water_results <- "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/Water/2021-01-15_full-model-posterior_Water Consumption_Economic-allocation.RData"

# LIVE WEIGHT Gross energy allocation (no wild capture results) + FCR Priors
# c_results <- "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/GHG/2021-01-15_full-model-posterior_Global warming potential_Gross energy content-allocation.RData"
# n_results <- "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/Nitrogen/2021-01-15_full-model-posterior_Marine eutrophication_Gross energy content-allocation.RData"
# p_results <- "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/Phosphorus/2021-01-16_full-model-posterior_Freshwater eutrophication_Gross energy content-allocation.RData"
# land_results <- "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/Land-FCR-priors-only/2021-01-16_full-model-posterior_Land Use_Gross energy content-allocation.RData"
# water_results <- "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/Water/2021-01-15_full-model-posterior_Water Consumption_Gross energy content-allocation.RData"

# LIVE WEIGHT Mass allocation + NO Priors
# c_results <- "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/NO PRIORS/GHG/2021-01-15_full-model-posterior_Global warming potential_Mass-allocation.RData"
# n_results <- "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/NO PRIORS/Nitrogen/2021-01-15_full-model-posterior_Marine eutrophication_Mass-allocation.RData"
# p_results <- "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/NO PRIORS/Phosphorus/2021-01-15_full-model-posterior_Freshwater eutrophication_Mass-allocation.RData"
# land_results <- "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/NO PRIORS/Land/2021-01-20_full-model-posterior_Land Use_Mass-allocation.RData"
# water_results <- "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/NO PRIORS/Water/2021-01-16_full-model-posterior_Water Consumption_Mass-allocation.RData"
# wild_results <- "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/NO PRIORS/Wild/2021-01-16_full-model-posterior_Wild-capture-ghg.RData"

###########################################################################################
# SET CHICKEN NUMBERS (need this to produce both live weight and edible weight options below):
chx_ghg_xmin_live <- 3127
chx_ghg_xmax_live <- 3589
chx_land_xmin_live <- 5794
chx_land_xmax_live <- 6048
chx_n_xmin_live <- 77
chx_n_xmax_live <- 92
chx_p_xmin_live <- 11
chx_p_xmax_live <- 16
chx_water_xmin_live <- 170
chx_water_xmax_live <- 202

# OPTIONS: Choose which environmental stressors for Chicken to display - only applies to Mass allocation
# CF for chicken = 40% i.e., For edible weight chicken numbers, multiply live weight by 1/0.40

# FOR EDIBLE WEIGHT + MASS ALLOCATION
# chx_fill <- "#D6CCB7"
# chx_ghg_xmin <- chx_ghg_xmin_live * (1/0.40)
# chx_ghg_xmax <- chx_ghg_xmax_live * (1/0.40)
# chx_land_xmin <- chx_land_xmin_live * (1/0.40)
# chx_land_xmax <- chx_land_xmax_live * (1/0.40)
# chx_n_xmin <- chx_n_xmin_live * (1/0.40)
# chx_n_xmax <- chx_n_xmax_live * (1/0.40)
# chx_p_xmin <- chx_p_xmin_live * (1/0.40)
# chx_p_xmax <- chx_p_xmax_live * (1/0.40)
# chx_water_xmin <- chx_water_xmin_live * (1/0.40)
# chx_water_xmax <- chx_water_xmax_live * (1/0.40)

# FOR LIVE WEIGHT + MASS ALLOCATION
chx_fill <- "#D6CCB7"
chx_ghg_xmin <- chx_ghg_xmin_live
chx_ghg_xmax <- chx_ghg_xmax_live
chx_land_xmin <- chx_land_xmin_live
chx_land_xmax <- chx_land_xmax_live
chx_n_xmin <- chx_n_xmin_live
chx_n_xmax <- chx_n_xmax_live
chx_p_xmin <- chx_p_xmin_live
chx_p_xmax <- chx_p_xmax_live
chx_water_xmin <- chx_water_xmin_live
chx_water_xmax <- chx_water_xmax_live

# FOR ALL ECONOMIC OR GROSS ENERGY CONTENT ALLOCATIONS (Live or Edible weight): I.E., FULLY TRANSPARENT WHITE
# chx_fill <- "#00FFFFFF"
# chx_ghg_xmin <- 0
# chx_ghg_xmax <- 0
# chx_land_xmin <- 0
# chx_land_xmax <- 0
# chx_n_xmin <- 0
# chx_n_xmax <- 0
# chx_p_xmin <- 0
# chx_p_xmax <- 0
# chx_water_xmin <- 0
# chx_water_xmax <- 0

###########################################################################################
# OPTIONS: LIVE WEIGHT VS EDIBLE WEIGHT MINOR PLOT ADJUSTMENTS
# (1) Set x-axis limits for GHG plots so plots for wild and farmed are aligned
# (2) X-axis label: per tonne 'live weight' vs 'edible weight'

# FOR LIVE WEIGHT
ghg_xlimit <- 12500
weight_type <- 'live weight'

# FOR EDIBLE WEIGHT
# ghg_xlimit <- 26000
# weight_type <- 'edible weight'

# SPECIAL CASES:
# FOR EDIBLE WEIGHT - GHG gross-energy-allocation - misc_diad, upper credible interval gets cutoff
# ghg_xlimit <- 29000 # use this for GHG gross-energy-allocation - misc_diad, upper credible interval gets cutoff
# weight_type <- 'edible weight'

# FOR EDIBLE WEIGHT - Mass-allocation NO PRIORS
# ghg_xlimit <- 32000 # otherwise wild lobster upper credible interval gets cutoff
# weight_type <- 'edible weight'

###########################################################################################
# OPTIONS: Customize what is plotted

# Plot only a subset of the data (e.g., for social media posts)
# filter_taxa <- TRUE
# taxa_list_farmed <- c("salmon", "bivalves", "tilapia", "misc carp", "catfish", "shrimp")
# taxa_list_wild <- c("flounder, etc", "shrimp", "bivalves", "salmon, etc", "cod, etc", "herring, etc")
# # MANUALLY reset ghg_xlimit based on what's being plotted
# ghg_xlimit <- 26000

# Plot all the data (e.g., for publication)
filter_taxa <- FALSE

######################################################################################################
# Carbon

load(file.path(outdir, c_results))
units_for_plot = bquote('kg'~CO[2]*'-eq t'^-1~.(weight_type))
interval_palette <- c("#9EA8B7", "#6A7A90", "#364F6B") # Order: light to dark

# Key for naming taxa levels
# Get full taxa group names back
tx_index_key_final <- lca_model_dat %>%
  group_by(clean_sci_name) %>%
  mutate(n_obs = n()) %>%
  ungroup() %>%
  select(taxa, tx) %>%
  unique() %>%
  arrange(taxa) %>%
  mutate(taxa = as.character(taxa),
         full_taxa_name = case_when(taxa == "hypoph_carp" ~ "silver/bighead",
                                    taxa == "oth_carp" ~ "misc carp",
                                    taxa == "misc_diad" ~ "misc diad",
                                    taxa == "misc_fresh" ~ "misc freshwater",
                                    taxa == "misc_marine" ~ "misc marine",
                                    taxa == "fresh_crust" ~ "freshwater crust",
                                    taxa == "plants" ~ "seaweeds",
                                    TRUE ~ taxa),
         taxa = as.factor(taxa),
         full_taxa_name = as.factor(full_taxa_name)) %>%
  mutate(source = "Farmed") %>%
  {if (filter_taxa == TRUE) filter(., full_taxa_name %in% taxa_list_farmed)
    else .}

# Set taxa order from low to high GHG (keep this order for all plots)
ghg_ordered <- fit_no_na %>%
  spread_draws(tx_total_fp_w[tx]) %>%
  median_qi(.width = c(0.95)) %>%
  right_join(tx_index_key_final, by = "tx") %>% # Join with index key to get sci and taxa names; use RIGHT_JOIN so that only desired taxa list are plotted
  # Filter out tx_total_fp_w for plants and bivalves
  filter(taxa %in% c("bivalves", "plants")==FALSE) %>%
  # Substitute on_farm data as total impact for plants and bivalves
  bind_rows(fit_no_na %>% 
              spread_draws(tx_farm_fp_w[tx]) %>% 
              median_qi(.width = c(0.95)) %>% 
              left_join(tx_index_key_final, by = "tx") %>% 
              filter(taxa %in% c("bivalves", "plants")) %>% 
              rename(tx_total_fp_w = tx_farm_fp_w)) %>%
  arrange(tx_total_fp_w) 
  
full_taxa_name_order_final <- ghg_ordered %>%
  pull(full_taxa_name) %>%
  as.character()

# Carbon:
carbon_data <- fit_no_na %>%
  spread_draws(tx_total_fp_w[tx]) %>%
  median_qi(.width = c(0.95, 0.8, 0.5)) %>%
  right_join(tx_index_key_final, by = "tx") %>% # Join with index key to get sci and taxa names; use RIGHT_JOIN so that only desired taxa list are plotted
  # Filter out tx_total_fp_w for plants and bivalves
  filter(taxa %in% c("bivalves", "plants")==FALSE) %>%
  # Substitute on_farm data as total impact for plants and bivalves
  bind_rows(fit_no_na %>% 
              spread_draws(tx_farm_fp_w[tx]) %>% 
              median_qi(.width = c(0.95, 0.8, 0.5)) %>% 
              left_join(tx_index_key_final, by = "tx") %>% 
              filter(taxa %in% c("bivalves", "plants")) %>% 
              rename(tx_total_fp_w = tx_farm_fp_w)) %>%
  # REORDER taxa axis
  mutate(full_taxa_name = fct_relevel(full_taxa_name, full_taxa_name_order_final))

p_carbon <- ggplot(data = carbon_data, aes(y = full_taxa_name, x = tx_total_fp_w, xmin = .lower, xmax = .upper)) +
  scale_x_continuous(labels = comma) +
  # ADD CHICKEN RESULTS:
  geom_rect(aes(ymin = -Inf, ymax = length(tx_index_key_final$full_taxa_name), xmin = chx_ghg_xmin, xmax = chx_ghg_xmax), fill = chx_fill) +
  geom_interval(aes(xmin = .lower, xmax = .upper), interval_size = 2.9) +
  geom_point(aes(y = full_taxa_name, x = tx_total_fp_w)) +
  geom_hline(yintercept = 1:length(tx_index_key_final$full_taxa_name), linetype = "dotted") +
  scale_color_manual(values = interval_palette) +
  coord_cartesian(xlim = c(0, ghg_xlimit)) +
  theme_classic() + 
  tx_plot_theme_final + 
  labs(x = units_for_plot, y = "", title = "")
#p_carbon
#ggsave(filename = file.path(outdir, "plot_Global warming potential_Mass-allocation_TOTAL-IMPACT-TAXA-LEVEL-WEIGHTED.png"), width = 11, height = 8.5)

# Output graphing data for SI: Only need to do this for Figure 1 plots
# carbon_chx_data <- data.frame(full_taxa_name = "chicken",
#                               .lower = chx_ghg_xmin,
#                               .upper = chx_ghg_xmax,
#                               stressor = "ghg")
# carbon_data %>%
#   rename(median = tx_total_fp_w) %>%
#   mutate(stressor = "ghg") %>%
#   select(-c(tx, .point, .interval, taxa)) %>%
#   arrange(full_taxa_name, .width) %>%
#   bind_rows(carbon_chx_data) %>%
#   write.csv(file.path(outdir, "Fig-1-ghg-data.csv"), row.names=FALSE)

chx_info <- ls()[grep("chx", ls())]
# Clear before next model
rm(list=ls()[!(ls() %in% c("outdir", "tx_index_key_final", "full_taxa_name_order_final", "tx_plot_theme_final",
                           "c_results", "p_results", "n_results", "land_results", "water_results", "wild_results",
                           chx_info, "chx_info", "ghg_xlimit", "weight_type",
                           "filter_taxa", "taxa_list_farmed", "taxa_list_wild",
                           "p_carbon"))])

######################################################################################################
# LAND

# Load model outputs individually, save to plot_grid, and create multi-panel plot
load(file.path(outdir, land_results))
units_for_plot = bquote('m'^2*'a t'^-1~.(weight_type))
interval_palette <- c("#57D182", "#42955E", "#2E5C3C")

# Use for land:
land_data <- fit_no_na %>%
  spread_draws(tx_land_total_w[tx]) %>%
  median_qi(.width = c(0.95, 0.8, 0.5)) %>%
  right_join(tx_index_key_final, by = "tx") %>% # Join with index key to get sci and taxa names
  # Set plant and bivalves distributions to 0 (both are 0 for on and off farm)
  mutate(tx_land_total_w = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = tx_land_total_w),
         .lower = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = .lower),
         .upper = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = .upper)) %>%
  # LAND can't be less than zero: (lower bound for misc diad fishes is negative)
  mutate(.lower = if_else(.lower < 0, true = 0, false = .lower)) %>%
  # REORDER taxa axis
  mutate(full_taxa_name = fct_relevel(full_taxa_name, full_taxa_name_order_final))

p_land <- ggplot(land_data, aes(y = full_taxa_name, x = tx_land_total_w, xmin = .lower, xmax = .upper)) +
  scale_x_continuous(labels = comma) +
  # ADD CHICKEN RESULTS:
  geom_rect(aes(ymin = -Inf, ymax = length(tx_index_key_final$full_taxa_name), xmin = chx_land_xmin, xmax = chx_land_xmax), fill = chx_fill) +
  geom_interval(aes(xmin = .lower, xmax = .upper), interval_size = 2.9) +
  geom_point(aes(y = full_taxa_name, x = tx_land_total_w)) +
  geom_point(x = 0, y = "bivalves") +
  geom_point(x = 0, y = "plants") +
  geom_hline(yintercept = 1:length(tx_index_key_final$full_taxa_name), linetype = "dotted") +
  scale_color_manual(values = interval_palette) +
  theme_classic() + 
  tx_plot_theme_final + 
  #theme(axis.text.y=element_blank()) +
  labs(x = units_for_plot, y = "", title = "")
#p_land
#ggsave(filename = file.path(outdir, "plot_Land Use_Mass-allocation_TOTAL-IMPACT-TAXA-LEVEL-WEIGHTED-formatted.png"), width = 11, height = 8.5)

# Output graphing data for SI: Only need to do this for Figure 1 plots
# land_chx_data <- data.frame(full_taxa_name = "chicken",
#                               .lower = chx_land_xmin,
#                               .upper = chx_land_xmax,
#                               stressor = "land")
# land_data %>%
#   rename(median = tx_land_total_w) %>%
#   mutate(stressor = "land") %>%
#   select(-c(tx, .point, .interval, taxa)) %>%
#   arrange(full_taxa_name, .width) %>%
#   bind_rows(land_chx_data) %>%
#   write.csv(file.path(outdir, "Fig-1-land-data.csv"), row.names=FALSE)

rm(list=ls()[!(ls() %in% c("outdir", "tx_index_key_final", "full_taxa_name_order_final", "tx_plot_theme_final",
                           "c_results", "p_results", "n_results", "land_results", "water_results", "wild_results",
                           chx_info, "chx_info", "ghg_xlimit", "weight_type",
                           "filter_taxa", "taxa_list_farmed", "taxa_list_wild",
                           "p_carbon",
                           "p_land"))])

######################################################################################################
# NITROGEN 

load(file.path(outdir, n_results))
units_for_plot = bquote('kg N-eq t'^-1~.(weight_type))
interval_palette <- c("#FFA647", "#B57736", "#704B25")

nitrogen_data <- fit_no_na %>%
  spread_draws(tx_total_fp_w[tx]) %>%
  median_qi(.width = c(0.95, 0.8, 0.5)) %>%
  right_join(tx_index_key_final, by = "tx") %>% # Join with index key to get sci and taxa names
  # REORDER taxa axis
  mutate(full_taxa_name = fct_relevel(full_taxa_name, full_taxa_name_order_final))

p_nitrogen <- ggplot(nitrogen_data, aes(y = full_taxa_name, x = tx_total_fp_w)) +
  # ADD CHICKEN RESULTS:
  geom_rect(aes(ymin = -Inf, ymax = length(tx_index_key_final$full_taxa_name), xmin = chx_n_xmin, xmax = chx_n_xmax), fill = chx_fill) +
  geom_interval(aes(xmin = .lower, xmax = .upper), interval_size = 2.9) +
  geom_point(aes(y = full_taxa_name, x = tx_total_fp_w)) +
  geom_hline(yintercept = 1:length(tx_index_key_final$full_taxa_name), linetype = "dotted") +
  scale_color_manual(values = interval_palette) +
  theme_classic() + 
  tx_plot_theme_final + 
  labs(x = units_for_plot, y = "", title = "")
#p_nitrogen
#ggsave(filename = file.path(outdir, "plot_Marine eutrophication_Mass-allocation_TOTAL-IMPACT-TAXA-LEVEL-WEIGHTED-formatted.png"), width = 11, height = 8.5)

# Output graphing data for SI: Only need to do this for Figure 1 plots
# nitrogen_chx_data <- data.frame(full_taxa_name = "chicken",
#                               .lower = chx_n_xmin,
#                               .upper = chx_n_xmax,
#                               stressor = "nitrogen")
# nitrogen_data %>%
#   rename(median = tx_total_fp_w) %>%
#   mutate(stressor = "nitrogen") %>%
#   select(-c(tx, .point, .interval, taxa)) %>%
#   arrange(full_taxa_name, .width) %>%
#   bind_rows(nitrogen_chx_data) %>%
#   write.csv(file.path(outdir, "Fig-1-nitrogen-data.csv"), row.names=FALSE)

rm(list=ls()[!(ls() %in% c("outdir", "tx_index_key_final", "full_taxa_name_order_final", "tx_plot_theme_final",
                           "c_results", "p_results", "n_results", "land_results", "water_results", "wild_results",
                           chx_info, "chx_info", "ghg_xlimit", "weight_type",
                           "filter_taxa", "taxa_list_farmed", "taxa_list_wild",
                           "p_carbon",
                           "p_land",
                           "p_nitrogen"))])

######################################################################################################
# PHOSPHORUS

load(file.path(outdir, p_results))
units_for_plot = bquote('kg P-eq t'^-1~.(weight_type))
interval_palette <- c("#FFB4C4", "#FF86A4", "#FC5185") # sample palette base to "white"

phosphorus_data <- fit_no_na %>%
  spread_draws(tx_total_fp_w[tx]) %>%
  median_qi(.width = c(0.95, 0.8, 0.5)) %>%
  right_join(tx_index_key_final, by = "tx") %>% # Join with index key to get sci and taxa names
  # REORDER taxa axis
  mutate(full_taxa_name = fct_relevel(full_taxa_name, full_taxa_name_order_final)) 

p_phosphorus <- ggplot(phosphorus_data, aes(y = full_taxa_name, x = tx_total_fp_w)) +
  # ADD CHICKEN RESULTS:
  geom_rect(aes(ymin = -Inf, ymax = length(tx_index_key_final$full_taxa_name), xmin = chx_p_xmin, xmax = chx_p_xmax), fill = chx_fill) +
  geom_interval(aes(xmin = .lower, xmax = .upper), interval_size = 2.9) +
  geom_point(aes(y = full_taxa_name, x = tx_total_fp_w)) +
  geom_hline(yintercept = 1:length(tx_index_key_final$full_taxa_name), linetype = "dotted") +
  scale_color_manual(values = interval_palette) +
  theme_classic() + 
  tx_plot_theme_final + 
  labs(x = units_for_plot, y = "", title = "")

# Output graphing data for SI: Only need to do this for Figure 1 plots
# phosphorus_chx_data <- data.frame(full_taxa_name = "chicken",
#                               .lower = chx_p_xmin,
#                               .upper = chx_p_xmax,
#                               stressor = "phosphorus")
# phosphorus_data %>%
#   rename(median = tx_total_fp_w) %>%
#   mutate(stressor = "phosphorus") %>%
#   select(-c(tx, .point, .interval, taxa)) %>%
#   arrange(full_taxa_name, .width) %>%
#   bind_rows(phosphorus_chx_data) %>%
#   write.csv(file.path(outdir, "Fig-1-phosphorus-data.csv"), row.names=FALSE)

rm(list=ls()[!(ls() %in% c("outdir", "tx_index_key_final", "full_taxa_name_order_final", "tx_plot_theme_final",
                           "c_results", "p_results", "n_results", "land_results", "water_results", "wild_results",
                           chx_info, "chx_info", "ghg_xlimit", "weight_type",
                           "filter_taxa", "taxa_list_farmed", "taxa_list_wild",
                           "p_carbon",
                           "p_land",
                           "p_nitrogen",
                           "p_phosphorus"))])

######################################################################################################
# WATER

load(file.path(outdir, water_results))
units_for_plot <- bquote('m'^3*' t'^-1~.(weight_type))
interval_palette <- c("#3FC1C9", "#348A8F", "#275659") # sample palette base to "black"

# Get off-farm impact for all non-freshwater taxa
# For non-freshwater taxa, on-farm impacts are ZERO, but because uncertainty around the zero is so large, the total impacts are negatie
# Therefore, use off-farm impacts as a substitute for total impacts
nonfresh_taxa <- c("misc_diad", "salmon", "milkfish", "misc_marine", "shrimp")
nonfresh_data <- fit_no_na %>%
  spread_draws(tx_feed_fp_w[tx]) %>%
  median_qi(.width = c(0.95, 0.8, 0.5)) %>%
  right_join(tx_index_key_final, by = "tx") %>% # Join with index key to get sci and taxa names
  # Only keep non-freshwater taxa - i.e., misc diadromous fishes, salmon, milkfish, misc marine fishes, and shrimp
  filter(taxa %in% nonfresh_taxa) %>%
  # rename off-farm impacts as total impacts so it can merge with the rest of the dataset
  rename(tx_total_fp_w = tx_feed_fp_w)

# Mean total impact taxa-level
water_data <- fit_no_na %>%
  spread_draws(tx_total_fp_w[tx]) %>%
  median_qi(.width = c(0.95, 0.8, 0.5)) %>%
  right_join(tx_index_key_final, by = "tx") %>% # Join with index key to get sci and taxa names
  # Remove non-freshwater taxa from previous section for which we're using off-farm impacts only as the total impact
  filter(taxa %in% nonfresh_taxa == FALSE) %>%
  # Set plant and bivalves distributions to 0 (both are 0 for on and off farm)
  mutate(tx_total_fp_w = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = tx_total_fp_w),
         .lower = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = .lower),
         .upper = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = .upper)) %>%
  # Rowbind with nonfresh_data
  bind_rows(nonfresh_data) %>%
  # WATER can't be less than zero: (lower bound for trout is negative)
  mutate(.lower = if_else(.lower < 0, true = 0, false = .lower)) %>%
  # REORDER taxa axis
  mutate(full_taxa_name = fct_relevel(full_taxa_name, full_taxa_name_order_final)) 

p_water <- ggplot(water_data, aes(y = full_taxa_name, x = tx_total_fp_w)) +
  # ADD CHICKEN RESULTS:
  geom_rect(aes(ymin = -Inf, ymax = length(tx_index_key_final$full_taxa_name), xmin = chx_water_xmin, xmax = chx_water_xmax), fill = chx_fill) +
  geom_interval(aes(xmin = .lower, xmax = .upper), interval_size = 2.9) +
  geom_point(x = 0, y = "bivalves") +
  geom_point(x = 0, y = "plants") +
  geom_point(aes(y = full_taxa_name, x = tx_total_fp_w)) +
  geom_hline(yintercept = 1:length(tx_index_key_final$full_taxa_name), linetype = "dotted") +
  scale_color_manual(values = interval_palette) +
  theme_classic() + 
  tx_plot_theme_final + 
  labs(x = units_for_plot, y = "", title = "")

# Output graphing data for SI: Only need to do this for Figure 1 plots
# water_chx_data <- data.frame(full_taxa_name = "chicken",
#                                   .lower = chx_water_xmin,
#                                   .upper = chx_water_xmax,
#                                   stressor = "water")
# water_data %>%
#   rename(median = tx_total_fp_w) %>%
#   mutate(stressor = "water") %>%
#   select(-c(tx, .point, .interval, taxa)) %>%
#   arrange(full_taxa_name, .width) %>%
#   bind_rows(water_chx_data) %>%
#   write.csv(file.path(outdir, "Fig-1-water-data.csv"), row.names=FALSE)

rm(list=ls()[!(ls() %in% c("outdir", "tx_index_key_final", "full_taxa_name_order_final", "tx_plot_theme_final",
                           "c_results", "p_results", "n_results", "land_results", "water_results", "wild_results",
                           chx_info, "chx_info", "ghg_xlimit", "weight_type",
                           "filter_taxa", "taxa_list_farmed", "taxa_list_wild",
                           "p_carbon",
                           "p_land",
                           "p_nitrogen",
                           "p_phosphorus", 
                           "p_water"))])

######################################################################################################
# WILD CATCH - (plot LAST because it has a different tx_index_key_final)

# Load model outputs individually, save to plot_grid, and create multi-panel plot
load(file.path(outdir, wild_results))
units_for_plot = bquote('kg'~CO[2]*'-eq t'^-1~.(weight_type))
interval_palette <- c("#9EA8B7", "#6A7A90", "#364F6B") # Order: light to dark

wild_index_key_final <- wild_dat_new_weights %>%
  group_by(clean_sci_name) %>%
  mutate(n_obs = n()) %>%
  ungroup() %>%
  select(taxa, tx) %>%
  unique() %>%
  arrange(taxa) %>%
  mutate(taxa = tolower(taxa)) %>%
  mutate(plot_taxa = case_when(taxa == "bivalves" ~ "bivalves",
                               taxa == "cephalopods" ~ "squid, etc",
                               taxa == "flatfishes" ~ "flounder, etc",
                               taxa == "gadiformes" ~ "cod, etc",
                               taxa == "jacks, mullets, sauries" ~ "jack, etc",
                               taxa == "large pelagic fishes" ~ "tuna, etc",
                               taxa == "lobsters" ~ "lobster",
                               taxa == "redfishes, basses, congers" ~ "redfish, etc",
                               taxa == "salmonids" ~ "salmon, etc",
                               taxa == "shrimps" ~ "shrimp",
                               taxa == "small pelagic fishes" ~ "herring, etc"
                               )) %>%
  mutate(source = "Wild")  %>%
  {if (filter_taxa == TRUE) filter(., plot_taxa %in% taxa_list_wild)
    else .}

# With out etc
# mutate(plot_taxa = case_when(taxa == "bivalves" ~ "bivalves",
#                              taxa == "cephalopods" ~ "squid",
#                              taxa == "flatfishes" ~ "flounder",
#                              taxa == "gadiformes" ~ "cod",
#                              taxa == "jacks, mullets, sauries" ~ "jack",
#                              taxa == "large pelagic fishes" ~ "tuna",
#                              taxa == "lobsters" ~ "lobster",
#                              taxa == "redfishes, basses, congers" ~ "redfish",
#                              taxa == "salmonids" ~ "salmon",
#                              taxa == "shrimps" ~ "shrimp",
#                              taxa == "small pelagic fishes" ~ "herringactually, "
# ))

# Wild:
wild_data <- fit_no_na %>%
  spread_draws(tx_ghg_w[tx]) %>%
  median_qi(tx_ghg_w, .width = c(0.95, 0.8, 0.5)) %>%
  right_join(wild_index_key_final, by = "tx") %>% # Join with index key to get sci and taxa names
  # Reorder from highest to lowest
  mutate(plot_taxa = as.factor(plot_taxa)) %>%
  mutate(plot_taxa = fct_reorder(plot_taxa, tx_ghg_w))

p_wild <- ggplot(wild_data, aes(y = plot_taxa, x = tx_ghg_w)) +
  scale_x_continuous(labels = comma) +
  # ADD CHICKEN RESULTS:
  geom_rect(aes(ymin = -Inf, ymax = length(wild_index_key_final$plot_taxa), xmin = chx_ghg_xmin, xmax = chx_ghg_xmax), fill = chx_fill) +
  geom_interval(aes(xmin = .lower, xmax = .upper), interval_size = 3) +
  geom_point(aes(y = plot_taxa, x = tx_ghg_w)) +
  geom_hline(yintercept = 1:length(wild_index_key_final$plot_taxa), linetype = "dotted") +
  scale_color_manual(values = interval_palette) +
  theme_classic() + 
  coord_cartesian(xlim = c(0, ghg_xlimit)) +
  tx_plot_theme_final + 
  labs(x = units_for_plot, y = "", title = "")
# p_wild
# ggsave(filename = file.path(outdir, "plot_WILD-GHG-TAXA-LEVEL-WEIGHTED-formatted.png"), width = 11, height = 8.5)

# Output graphing data for SI: Only need to do this for Figure 1 plots
# wild_data %>%
#   rename(median = tx_ghg_w,
#          full_taxa_name = plot_taxa) %>%
#   mutate(stressor = "ghg") %>%
#   select(-c(tx, .point, .interval, taxa)) %>%
#   arrange(full_taxa_name, .width) %>%
#   write.csv(file.path(outdir, "Fig-1-wild-data.csv"), row.names=FALSE)

rm(list=ls()[!(ls() %in% c("outdir", "tx_index_key_final", "full_taxa_name_order_final", "tx_plot_theme_final", "wild_index_key_final",
                           "c_results", "p_results", "n_results", "land_results", "water_results", "wild_results",
                           chx_info, "chx_info", "ghg_xlimit", "weight_type",
                           "filter_taxa", "taxa_list_farmed", "taxa_list_wild",
                           "p_carbon",
                           "p_land",
                           "p_nitrogen",
                           "p_phosphorus", 
                           "p_water",
                           "p_wild"))])

# COMBINE INTO PANEL WITH FACET LABELS:
# USE: facet_grid(switch = "y") # "switch" moves y-axis label from right to left side
# OPTION: Choose 5 panel or 6 panel plot

# 6 panel plot:
p_carbon_facet <- p_carbon + facet_grid(rows = "source", switch = "y") + theme(strip.placement = "outside")
p_wild_facet <- p_wild + facet_grid(rows = "source", switch = "y") + theme(strip.placement = "outside")
p_left <- plot_grid(p_carbon_facet, p_wild_facet, ncol = 1, nrow = 2, align = "hv", labels = c("a", "d"), axis = "l")
p_water_facet <- p_water + facet_grid(rows = "source", switch = "y") + theme(strip.placement = "outside")
p_right <- plot_grid(p_nitrogen + theme(axis.text.y=element_blank()),
                     p_phosphorus + theme(axis.text.y=element_blank()),
                     p_water_facet,
                     p_land + theme(axis.text.y=element_blank()), ncol = 2, nrow = 2, align = "h", axis = "b", labels = c("b", "c", "e", "f"), rel_widths = c(1.3, 1))

# NOTE: this conforms to Nature figure specs (equivalent to 183mm for two-column width)
plot_grid(p_left, p_right, ncol = 2, nrow = 1, rel_widths = c(0.55, 1))
ggsave(filename = file.path(outdir, "plot_Figure-X.png"), width = 183, height = 90, units = "mm") # manually change filename after save

pdf(file = file.path(outdir, "plot_Figure-X.pdf"), width = 7.2, height = 3.5) # equivalent to 183 x 90 mm (Nature figure specs)
plot_grid(p_left, p_right, ncol = 2, nrow = 1, rel_widths = c(0.55, 1))
dev.off()

# Use NO FACET version below for 5-panel plot
# 5 panel plot
# NOTE: this conforms to Nature figure specs (equivalent to 183mm for two-column width)
# p_left <- plot_grid(p_carbon, p_water, nrow = 2, ncol = 1, labels = c("a", "d"))
# p_right <- plot_grid(p_nitrogen + theme(axis.text.y=element_blank()),
#                      p_phosphorus + theme(axis.text.y=element_blank()),
#                      p_land + theme(axis.text.y=element_blank()), ncol = 2, nrow = 2, align = "h", labels = c("b", "c", "e"))
# plot_grid(p_left, p_right, ncol = 2, nrow = 1, align = "v", rel_widths = c(0.6, 1))
# ggsave(filename = file.path(outdir, "plot_Figure-X.png"), width = 183, height = 90, units = "mm")

# Version for social media:
# PDF
# Equal sized X-axes (can be 6 separate plots or one combined plot)
# Farmed: salmon, bivalves, tilapia, misc carps, catfish, shrimp
# Capture: flounder, shrimp, bivalves, salmon, cod, herring 
# pdf(file = file.path(outdir, "plot_Figure-1_Azote-version-squished.pdf"), width = 7.2, height = 2.4) # equivalent to 183 x 90 mm (Nature figure specs)
# plot_grid(p_carbon, p_nitrogen, p_phosphorus,
#           p_wild, p_water, p_land,
#           nrow = 2, ncol = 3, align = "h", labels = c("a", "b", "c", "d", "e", "f"))
# dev.off()
