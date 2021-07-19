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
datadir <- "/Volumes/jgephart/BFA Environment 2/Data"
outdir <- "/Volumes/jgephart/BFA Environment 2/Outputs"

# Theme
# ORIGINAL tx_plot_theme_final for first Nature submission
# tx_plot_theme_final <- theme(axis.title = element_text(size = 9),
#                        axis.text.x = element_text(size = 9, color = "black"),
#                        axis.text.y = element_text(size = 8, color = "black"),
#                        legend.position = "none",
#                        plot.margin = unit(c(0, 3, 0, 0), "mm")) # increase the right margin so that there can be some pillover from x-axis
# NEW tx_plot_theme_final for first Nature submission - smaller x-axis font + hjust on x-axis
tx_plot_theme_final <- theme(axis.title = element_text(size = 9),
                       axis.text.x = element_text(hjust = 0.6, size = 7, color = "black"),
                       axis.text.y = element_text(size = 8, color = "black"),
                       legend.position = "none",
                       plot.margin = unit(c(0, 3, 0, 0), "mm")) # increase the right margin so that there can be some pillover from x-axis

###########################################################################################
# OPTIONS: Set filenames
# EDIBLE WEIGHT RESULTS

# EDIBLE WEIGHT Mass allocation + FCR Priors
c_results <- "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/GHG/2021-04-27_full-model-posterior_Global warming potential_Mass-allocation.RData"
n_results <- "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Nitrogen/2021-04-28_full-model-posterior_Marine eutrophication_Mass-allocation.RData"
p_results <- "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Phosphorus/2021-04-28_full-model-posterior_Freshwater eutrophication_Mass-allocation.RData"
land_results <- "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Land-FCR-Priors-Only/2021-04-28_full-model-posterior_Land Use_Mass-allocation.RData"
water_results <- "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Water/2021-04-28_full-model-posterior_Water Consumption_Mass-allocation.RData"
wild_results <- "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Wild/2021-04-27_full-model-posterior_Wild-capture-ghg.RData"

# EDIBLE WEIGHT Economic allocation (no wild capture results) + FCR Priors
# c_results <- "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/GHG/2021-04-27_full-model-posterior_Global warming potential_Economic-allocation.RData"
# n_results <- "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Nitrogen/2021-04-28_full-model-posterior_Marine eutrophication_Economic-allocation.RData"
# p_results <- "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Phosphorus/2021-04-28_full-model-posterior_Freshwater eutrophication_Economic-allocation.RData"
# land_results <- "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Land-FCR-Priors-Only/2021-04-30_full-model-posterior_Land Use_Economic-allocation.RData"
# water_results <- "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Water/2021-04-28_full-model-posterior_Water Consumption_Economic-allocation.RData"

# EDIBLE WEIGHT Gross energy allocation (no wild capture results) + FCR Priors
# c_results <- "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/GHG/2021-04-27_full-model-posterior_Global warming potential_Gross energy content-allocation.RData"
# n_results <- "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Nitrogen/2021-04-28_full-model-posterior_Marine eutrophication_Gross energy content-allocation.RData"
# p_results <- "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Phosphorus/2021-04-28_full-model-posterior_Freshwater eutrophication_Gross energy content-allocation.RData"
# land_results <- "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Land-FCR-Priors-Only/2021-04-30_full-model-posterior_Land Use_Gross energy content-allocation.RData"
# water_results <- "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Water/2021-04-28_full-model-posterior_Water Consumption_Gross energy content-allocation.RData"

# EDIBLE WEIGHT Mass allocation + NO Priors
# c_results <- "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/NO PRIORS/GHG/2021-04-27_full-model-posterior_Global warming potential_Mass-allocation.RData"
# n_results <- "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/NO PRIORS/Nitrogen/2021-04-28_full-model-posterior_Marine eutrophication_Mass-allocation.RData"
# p_results <- "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/NO PRIORS/Phosphorus/2021-04-28_full-model-posterior_Freshwater eutrophication_Mass-allocation.RData"
# land_results <- "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/NO PRIORS/Land/2021-04-29_full-model-posterior_Land Use_Mass-allocation.RData"
# water_results <- "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/NO PRIORS/Water/2021-04-28_full-model-posterior_Water Consumption_Mass-allocation.RData"
# wild_results <- "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/NO PRIORS/Wild/2021-04-27_full-model-posterior_Wild-capture-ghg.RData"

########################
# LIVE WEIGHT RESULTS

# LIVE WEIGHT Mass allocation + FCR Priors
# c_results <- "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/GHG/2021-01-15_full-model-posterior_Global warming potential_Mass-allocation.RData"
# n_results <- "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/Nitrogen/2021-01-15_full-model-posterior_Marine eutrophication_Mass-allocation.RData"
# p_results <- "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/Phosphorus/2021-01-16_full-model-posterior_Freshwater eutrophication_Mass-allocation.RData"
# land_results <- "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/Land-FCR-priors-only/2021-01-20_full-model-posterior_Land Use_Mass-allocation.RData"
# water_results <- "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/Water/2021-01-15_full-model-posterior_Water Consumption_Mass-allocation.RData"
# wild_results <- "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/Wild/2021-01-16_full-model-posterior_Wild-capture-ghg.RData"

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
chx_ghg_xmin_live <- 3165
chx_ghg_xmax_live <- 3627
chx_land_xmin_live <- 5794
chx_land_xmax_live <- 6047
chx_n_xmin_live <- 77
chx_n_xmax_live <- 92
chx_p_xmin_live <- 11
chx_p_xmax_live <- 16
chx_water_xmin_live <- 170
chx_water_xmax_live <- 202

# OPTIONS: Choose which environmental stressors for Chicken to display - only applies to Mass allocation
# CF for chicken = 40% i.e., For edible weight chicken numbers, multiply live weight by 1/0.40

# FOR EDIBLE WEIGHT + MASS ALLOCATION
chx_fill <- "#D6CCB7"
chx_ghg_xmin <- chx_ghg_xmin_live * (1/0.40)
chx_ghg_xmax <- chx_ghg_xmax_live * (1/0.40)
chx_land_xmin <- chx_land_xmin_live * (1/0.40)
chx_land_xmax <- chx_land_xmax_live * (1/0.40)
chx_n_xmin <- chx_n_xmin_live * (1/0.40)
chx_n_xmax <- chx_n_xmax_live * (1/0.40)
chx_p_xmin <- chx_p_xmin_live * (1/0.40)
chx_p_xmax <- chx_p_xmax_live * (1/0.40)
chx_water_xmin <- chx_water_xmin_live * (1/0.40)
chx_water_xmax <- chx_water_xmax_live * (1/0.40)

# FOR LIVE WEIGHT + MASS ALLOCATION
# chx_fill <- "#D6CCB7"
# chx_ghg_xmin <- chx_ghg_xmin_live
# chx_ghg_xmax <- chx_ghg_xmax_live
# chx_land_xmin <- chx_land_xmin_live
# chx_land_xmax <- chx_land_xmax_live
# chx_n_xmin <- chx_n_xmin_live
# chx_n_xmax <- chx_n_xmax_live
# chx_p_xmin <- chx_p_xmin_live
# chx_p_xmax <- chx_p_xmax_live
# chx_water_xmin <- chx_water_xmin_live
# chx_water_xmax <- chx_water_xmax_live

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
# ghg_xlimit <- 12500
# weight_type <- 'live weight'

# FOR EDIBLE WEIGHT
ghg_xlimit <- 26000
#ghg_xlimit <- 29000 # use this for GHG gross-energy-allocation - misc_diad, upper credible interval gets cutoff
weight_type <- 'edible weight'

###########################################################################################
# OPTIONS: Customize what is plotted

# Plot only a subset of the data (e.g., for social media posts)
filter_taxa <- TRUE
taxa_list_farmed <- c("salmon", "bivalves", "tilapia", "misc carp", "catfish", "shrimp")
taxa_list_wild <- c("flounder, etc", "shrimp", "bivalves", "salmon, etc", "cod, etc", "herring, etc")
# MANUALLY reset ghg_xlimit based on what's being plotted
ghg_xlimit <- 26000

# Plot all the data (e.g., for publication)
# filter_taxa <- FALSE

######################################################################################################
# Carbon

load(file.path(outdir, c_results))
#units_for_plot = "kg CO2-eq per tonne"
#units_for_plot = bquote('kg'~CO[2]*'-eq per tonne')
units_for_plot = bquote('kg'~CO[2]*' t'^-1~.(weight_type))
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
p_carbon <- fit_no_na %>%
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
  mutate(full_taxa_name = fct_relevel(full_taxa_name, full_taxa_name_order_final)) %>%
  #mutate(full_taxa_name = fct_reorder(full_taxa_name, tx_total_fp_w)) %>%
  ggplot(aes(y = full_taxa_name, x = tx_total_fp_w, xmin = .lower, xmax = .upper)) +
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
#units_for_plot = bquote('m'^2*'a per tonne')
units_for_plot = bquote('m'^2*' t'^-1~.(weight_type))

# x <- seq(0, 1, length.out = 16)
# base_color <- "#57D182"
# library(scales)
# show_col(seq_gradient_pal(base_color, "black")(x)) # Get hexadecimals for other colors

# COMPLEMENTARY GREEN (Azote guidelines) sample base to "black"
interval_palette <- c("#57D182", "#42955E", "#2E5C3C")
# COMPLEMENTARY GREEN (Azote guidelines) sample base to "white"
# interval_palette <- c("#B7EBC4", "#8BDEA3", "#57D182")
# PRIMARY GREEN:
#interval_palette <- c("#D9EAB2", "#C2DD86", "#A9D158") # Order: light to dark

# Use for land:
p_land <- fit_no_na %>%
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
  mutate(full_taxa_name = fct_relevel(full_taxa_name, full_taxa_name_order_final)) %>%
  #mutate(full_taxa_name = fct_reorder(full_taxa_name, tx_feed_fp_w)) %>%
  ggplot(., aes(y = full_taxa_name, x = tx_land_total_w, xmin = .lower, xmax = .upper)) +
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

rm(list=ls()[!(ls() %in% c("outdir", "tx_index_key_final", "full_taxa_name_order_final", "tx_plot_theme_final",
                           "c_results", "p_results", "n_results", "land_results", "water_results", "wild_results",
                           chx_info, "chx_info", "ghg_xlimit", "weight_type",
                           "filter_taxa", "taxa_list_farmed", "taxa_list_wild",
                           "p_carbon",
                           "p_land"))])

######################################################################################################
# NITROGEN 

load(file.path(outdir, n_results))
#units_for_plot = "kg N-eq per tonne"
units_for_plot = bquote('kg N-eq t'^-1~.(weight_type))

# x <- seq(0, 1, length.out = 16)
# base_color <- "#FFA647"
# library(scales)
# show_col(seq_gradient_pal(base_color, "black")(x)) # Get hexadecimals for other colors

# PRIMARY COLOR (Azote guidelines)
# interval_palette <- c("#FFEDAD", "#FFE37D", "#FFD947") # Order: light to dark
# COMPLEMENTARY COLOR sample base to "white"
# interval_palette <- c("#FFD5A9", "#FFBD79", "#FFA647")
# COMPLEMENTARY COLOR sample base to "black"
interval_palette <- c("#FFA647", "#B57736", "#704B25")

p_nitrogen <- fit_no_na %>%
  spread_draws(tx_total_fp_w[tx]) %>%
  median_qi(.width = c(0.95, 0.8, 0.5)) %>%
  right_join(tx_index_key_final, by = "tx") %>% # Join with index key to get sci and taxa names
  # Set .lower limit of plants and bivalves to be 0
  #mutate(.lower = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = .lower)) %>% 
  # REORDER taxa axis
  mutate(full_taxa_name = fct_relevel(full_taxa_name, full_taxa_name_order_final)) %>%
  #mutate(full_taxa_name = fct_reorder(full_taxa_name, tx_feed_fp_w)) %>%
  ggplot(aes(y = full_taxa_name, x = tx_total_fp_w)) +
  # ADD CHICKEN RESULTS:
  geom_rect(aes(ymin = -Inf, ymax = length(tx_index_key_final$full_taxa_name), xmin = chx_n_xmin, xmax = chx_n_xmax), fill = chx_fill) +
  geom_interval(aes(xmin = .lower, xmax = .upper), interval_size = 2.9) +
  geom_point(aes(y = full_taxa_name, x = tx_total_fp_w)) +
  geom_hline(yintercept = 1:length(tx_index_key_final$full_taxa_name), linetype = "dotted") +
  scale_color_manual(values = interval_palette) +
  theme_classic() + 
  tx_plot_theme_final + 
  #theme(axis.text.y=element_blank()) +
  #scale_y_discrete(labels = rep("", times = length(tx_index_key_final$full_taxa_name))) +
  # theme(axis.text.y=element_blank(),
  #       axis.title.y=element_blank(),
  #       panel.grid = element_blank(),
  #       panel.border = element_blank()) +
  labs(x = units_for_plot, y = "", title = "")
#p_nitrogen
#ggsave(filename = file.path(outdir, "plot_Marine eutrophication_Mass-allocation_TOTAL-IMPACT-TAXA-LEVEL-WEIGHTED-formatted.png"), width = 11, height = 8.5)


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
#units_for_plot = "kg P-eq per tonne"
units_for_plot = bquote('kg P-eq t'^-1~.(weight_type))

# Order: light to dark
#interval_palette <- c("#FC5185", "#B33E60", "#6F2B3D") # sample palette base to "black"
interval_palette <- c("#FFB4C4", "#FF86A4", "#FC5185") # sample palette base to "white"

p_phosphorus <- fit_no_na %>%
  spread_draws(tx_total_fp_w[tx]) %>%
  median_qi(.width = c(0.95, 0.8, 0.5)) %>%
  right_join(tx_index_key_final, by = "tx") %>% # Join with index key to get sci and taxa names
  # Set .lower limit of plants and bivalves to be 0
  #mutate(.lower = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = .lower)) %>% 
  # REORDER taxa axis
  mutate(full_taxa_name = fct_relevel(full_taxa_name, full_taxa_name_order_final)) %>%
  #mutate(full_taxa_name = fct_reorder(full_taxa_name, tx_feed_fp_w)) %>%
  ggplot(aes(y = full_taxa_name, x = tx_total_fp_w)) +
  # ADD CHICKEN RESULTS:
  geom_rect(aes(ymin = -Inf, ymax = length(tx_index_key_final$full_taxa_name), xmin = chx_p_xmin, xmax = chx_p_xmax), fill = chx_fill) +
  geom_interval(aes(xmin = .lower, xmax = .upper), interval_size = 2.9) +
  geom_point(aes(y = full_taxa_name, x = tx_total_fp_w)) +
  geom_hline(yintercept = 1:length(tx_index_key_final$full_taxa_name), linetype = "dotted") +
  scale_color_manual(values = interval_palette) +
  theme_classic() + 
  tx_plot_theme_final + 
  #theme(axis.text.y=element_blank()) +
  #scale_y_discrete(labels = rep("", times = length(tx_index_key_final$full_taxa_name))) +
  # theme(axis.text.y=element_blank(),
  #       axis.title.y=element_blank(),
  #       panel.grid = element_blank(),
  #       panel.border = element_blank()) +
  labs(x = units_for_plot, y = "", title = "")

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
#units_for_plot <- bquote('m'^3*'per tonne')
units_for_plot <- bquote('m'^3*' t'^-1~.(weight_type))

# x <- seq(0, 1, length.out = 16)
# base_color <- "#3FC1C9"
# library(scales)
# show_col(seq_gradient_pal(base_color, "black")(x)) # Get hexadecimals for other colors

# Order: light to dark
#interval_palette <- c("#B2E2E6", "#80D2D7", "#3FC1C9") # sample palette base to "white"
interval_palette <- c("#3FC1C9", "#348A8F", "#275659") # sample palette base to "black"

# Mean total impact taxa-level
p_water <- fit_no_na %>%
  spread_draws(tx_total_fp_w[tx]) %>%
  median_qi(.width = c(0.95, 0.8, 0.5)) %>%
  right_join(tx_index_key_final, by = "tx") %>% # Join with index key to get sci and taxa names
  # Set plant and bivalves distributions to 0 (both are 0 for on and off farm)
  mutate(tx_total_fp_w = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = tx_total_fp_w),
         .lower = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = .lower),
         .upper = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = .upper)) %>%
  # REORDER taxa axis
  mutate(full_taxa_name = fct_relevel(full_taxa_name, full_taxa_name_order_final)) %>%
  #mutate(full_taxa_name = fct_reorder(full_taxa_name, tx_feed_fp_w)) %>%
  ggplot(aes(y = full_taxa_name, x = tx_total_fp_w)) +
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
#units_for_plot = bquote('kg'~CO[2]*'-eq per tonne')
units_for_plot = bquote('kg'~CO[2]*' t'^-1~.(weight_type))

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
p_wild <- fit_no_na %>%
  spread_draws(tx_ghg_w[tx]) %>%
  median_qi(tx_ghg_w, .width = c(0.95, 0.8, 0.5)) %>%
  right_join(wild_index_key_final, by = "tx") %>% # Join with index key to get sci and taxa names
  # Reorder from highest to lowest
  mutate(plot_taxa = as.factor(plot_taxa)) %>%
  mutate(plot_taxa = fct_reorder(plot_taxa, tx_ghg_w)) %>%
  ggplot(aes(y = plot_taxa, x = tx_ghg_w)) +
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

rm(list=ls()[!(ls() %in% c("outdir", "tx_index_key_final", "full_taxa_name_order_final", "tx_plot_theme_final",
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
ggsave(filename = file.path(outdir, "plot_Figure-X.tiff"), device = "tiff", width = 183, height = 90, units = "mm")  # manually change filename after save

# For PDF
# pdf(file = file.path(outdir, "plot_Figure-X.pdf"), width = 7.2, height = 3.5) # equivalent to 183 x 90 mm (Nature figure specs)
# plot_grid(p_left, p_right, ncol = 2, nrow = 1, rel_widths = c(0.55, 1))
# dev.off()

# Use NO FACET version below for 5-panel plot
# 5 panel plot
# p_left <- plot_grid(p_carbon, p_water, nrow = 2, ncol = 1, labels = c("a", "d"))
# p_right <- plot_grid(p_nitrogen + theme(axis.text.y=element_blank()), 
#                      p_phosphorus + theme(axis.text.y=element_blank()), 
#                      p_land + theme(axis.text.y=element_blank()), ncol = 2, nrow = 2, align = "h", labels = c("b", "c", "e"))
# plot_grid(p_left, p_right, ncol = 2, nrow = 1, align = "v", rel_widths = c(0.6, 1))
# ggsave(filename = file.path(outdir, "plot_Figure-X.png"), width = 183, height = 90, units = "mm")
# ggsave(filename = file.path(outdir, "plot_Figure-X.tiff"), device = "tiff", width = 183, height = 90, units = "mm")

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
