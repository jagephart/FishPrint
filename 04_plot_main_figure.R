# PLOT MAIN FIGURE
# Load model outputs individually, save to plot_grid, and create multi-panel plot
# Remove y-axis labels for panels B, C, and F
# Set plot margins for all plots to 0
# Set axis label margins for panels C and F to 0


#### FIX IT - need to comment out Chicken info (gemo_rect) for energy and economic allocation (the chicken calculations only apply to mass allocation)

rm(list=ls())
library(tidyverse)
library(tidybayes)
library(ggplot2)
library(cowplot) ## for plot_grid
datadir <- "/Volumes/jgephart/BFA Environment 2/Data"
outdir <- "/Volumes/jgephart/BFA Environment 2/Outputs"

# Set filenames:
# Mass allocation
# c_results <- "PRIORS/GHG/2021-01-15_full-model-posterior_Global warming potential_Mass-allocation.RData"
# n_results <- "PRIORS/Nitrogen/2021-01-15_full-model-posterior_Marine eutrophication_Mass-allocation.RData"
# p_results <- "PRIORS/Phosphorus/2021-01-16_full-model-posterior_Freshwater eutrophication_Mass-allocation.RData"
# land_results <- "PRIORS/Land/2021-01-20_full-model-posterior_Land Use_Mass-allocation.RData"
# water_results <- "PRIORS/Water/2021-01-15_full-model-posterior_Water Consumption_Mass-allocation.RData"
# wild_results <- "PRIORS/Wild/2021-01-16_full-model-posterior_Wild-Capture-ghg.RData"

# Gross energy allocation (no wild capture results)
# c_results <- "PRIORS/GHG/2021-01-15_full-model-posterior_Global warming potential_Gross energy content-allocation.RData"
# n_results <- "PRIORS/Nitrogen/2021-01-15_full-model-posterior_Marine eutrophication_Gross energy content-allocation.RData"
# p_results <- "PRIORS/Phosphorus/2021-01-16_full-model-posterior_Freshwater eutrophication_Gross energy content-allocation.RData"
# land_results <- "PRIORS/Land/2021-01-16_full-model-posterior_Land Use_Gross energy content-allocation.RData"
# water_results <- "PRIORS/Water/2021-01-15_full-model-posterior_Water Consumption_Gross energy content-allocation.RData"

# Economic allocation (no wild capture results)
c_results <- "PRIORS/GHG/2021-01-15_full-model-posterior_Global warming potential_Economic-allocation.RData"
n_results <- "PRIORS/Nitrogen/2021-01-20_full-model-posterior_Marine eutrophication_Economic-allocation.RData"
p_results <- "PRIORS/Phosphorus/2021-01-16_full-model-posterior_Freshwater eutrophication_Economic-allocation.RData"
land_results <- "PRIORS/Land/2021-01-16_full-model-posterior_Land Use_Economic-allocation.RData"
water_results <- "PRIORS/Water/2021-01-15_full-model-posterior_Water Consumption_Economic-allocation.RData"

# Mass allocation NO PRIORS
# c_results <- "NO PRIORS/GHG/2021-01-15_full-model-posterior_Global warming potential_Mass-allocation.RData"
# n_results <- "NO PRIORS/Nitrogen/2021-01-15_full-model-posterior_Marine eutrophication_Mass-allocation.RData"
# p_results <- "NO PRIORS/Phosphorus/2021-01-15_full-model-posterior_Freshwater eutrophication_Mass-allocation.RData"
# land_results <- "NO PRIORS/Land/2021-01-20_full-model-posterior_Land Use_Mass-allocation.RData"
# water_results <- "NO PRIORS/Water/2021-01-16_full-model-posterior_Water Consumption_Mass-allocation.RData"
# wild_results <- "NO PRIORS/Wild/2021-01-16_full-model-posterior_Wild-capture-ghg.RData"

######################################################################################################
# Carbon

# Get color for chicken:
# x <- seq(0, 1, length.out = 16)
# base_color <- "#A69569"
# library(scales)
# show_col(seq_gradient_pal(base_color, "white")(x)) # Get hexadecimals for other colors

load(file.path(outdir, c_results))
#units_for_plot = "kg CO2-eq per tonne"
units_for_plot = bquote('kg'~CO[2]*'-eq per tonne')
interval_palette <- c("#9EA8B7", "#6A7A90", "#364F6B") # Order: light to dark

# Theme
tx_plot_theme <- theme(axis.title = element_text(size = 9),
                       axis.text.x = element_text(size = 9, color = "black"),
                       axis.text.y = element_text(size = 8, color = "black"),
                       legend.position = "none",
                       plot.margin = unit(c(0, 3, 0, 0), "mm")) # increase the right margin so that there can be some pillover from x-axis

# Key for naming taxa levels
# Get full taxa group names back
tx_index_key <- lca_model_dat %>%
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
  mutate(source = "Farmed")


# Set taxa order from low to high GHG (keep this order for all plots)
full_taxa_name_order <- fit_no_na %>%
  spread_draws(tx_total_fp_w[tx]) %>%
  median_qi(.width = c(0.95)) %>%
  left_join(tx_index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  # Filter out tx_total_fp_w for plants and bivalves
  filter(taxa %in% c("bivalves", "plants")==FALSE) %>%
  # Substitute on_farm data as total impact for plants and bivalves
  bind_rows(fit_no_na %>% 
              spread_draws(tx_farm_fp_w[tx]) %>% 
              median_qi(.width = c(0.95)) %>% 
              left_join(tx_index_key, by = "tx") %>% 
              filter(taxa %in% c("bivalves", "plants")) %>% 
              rename(tx_total_fp_w = tx_farm_fp_w)) %>%
  arrange(tx_total_fp_w) %>%
  pull(full_taxa_name) %>%
  as.character()

# Carbon:
p_carbon <- fit_no_na %>%
  spread_draws(tx_total_fp_w[tx]) %>%
  median_qi(.width = c(0.95, 0.8, 0.5)) %>%
  left_join(tx_index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  # Filter out tx_total_fp_w for plants and bivalves
  filter(taxa %in% c("bivalves", "plants")==FALSE) %>%
  # Substitute on_farm data as total impact for plants and bivalves
  bind_rows(fit_no_na %>% 
              spread_draws(tx_farm_fp_w[tx]) %>% 
              median_qi(.width = c(0.95, 0.8, 0.5)) %>% 
              left_join(tx_index_key, by = "tx") %>% 
              filter(taxa %in% c("bivalves", "plants")) %>% 
              rename(tx_total_fp_w = tx_farm_fp_w)) %>%
  # REORDER taxa axis
  mutate(full_taxa_name = fct_relevel(full_taxa_name, full_taxa_name_order)) %>%
  #mutate(full_taxa_name = fct_reorder(full_taxa_name, tx_total_fp_w)) %>%
  ggplot(aes(y = full_taxa_name, x = tx_total_fp_w, xmin = .lower, xmax = .upper)) +
  # ADD CHICKEN RESULTS:
  #geom_rect(aes(ymin = -Inf, ymax = nlevels(full_taxa_name), xmin = 3165, xmax = 3627), fill = "#D6CCB7") +
  geom_interval(aes(xmin = .lower, xmax = .upper), interval_size = 2.9) +
  geom_point(aes(y = full_taxa_name, x = tx_total_fp_w)) +
  geom_hline(yintercept = 1:12, linetype = "dotted") +
  scale_color_manual(values = interval_palette) +
  coord_cartesian(xlim = c(0, 12500)) +
  theme_classic() + 
  tx_plot_theme + 
  labs(x = units_for_plot, y = "", title = "")


# Clear before next model
rm(list=ls()[!(ls() %in% c("outdir", "tx_index_key", "full_taxa_name_order", "tx_plot_theme",
                           "c_results", "p_results", "n_results", "land_results", "water_results", "wild_results",
                           "p_carbon"))])

######################################################################################################
# LAND

# Load model outputs individually, save to plot_grid, and create multi-panel plot
load(file.path(outdir, land_results))
units_for_plot = bquote('m'^2*'a per tonne')
# COMPLEMENTARY GREEN (Azote guidelines)
interval_palette <- c("#B7EBC4", "#8BDEA3", "#57D182")
# PRIMARY GREEN:
#interval_palette <- c("#D9EAB2", "#C2DD86", "#A9D158") # Order: light to dark


# Use for land:
p_land <- fit_no_na %>%
  spread_draws(tx_land_total_w[tx]) %>%
  median_qi(.width = c(0.95, 0.8, 0.5)) %>%
  left_join(tx_index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  # Set plant and bivalves distributions to 0 (both are 0 for on and off farm)
  mutate(tx_land_total_w = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = tx_land_total_w),
         .lower = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = .lower),
         .upper = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = .upper)) %>%
  # LAND can't be less than zero: (lower bound for misc diad fishes is negative)
  mutate(.lower = if_else(.lower < 0, true = 0, false = .lower)) %>%
  # REORDER taxa axis
  mutate(full_taxa_name = fct_relevel(full_taxa_name, full_taxa_name_order)) %>%
  #mutate(full_taxa_name = fct_reorder(full_taxa_name, tx_feed_fp_w)) %>%
  ggplot(., aes(y = full_taxa_name, x = tx_land_total_w, xmin = .lower, xmax = .upper)) +
  # ADD CHICKEN RESULTS:
  #geom_rect(aes(ymin = -Inf, ymax = nlevels(full_taxa_name), xmin = 5794, xmax = 6047), fill = "#D6CCB7") +
  geom_interval(aes(xmin = .lower, xmax = .upper), interval_size = 2.9) +
  geom_point(aes(y = full_taxa_name, x = tx_land_total_w)) +
  geom_point(x = 0, y = "bivalves") +
  geom_point(x = 0, y = "plants") +
  geom_hline(yintercept = 1:12, linetype = "dotted") +
  scale_color_manual(values = interval_palette) +
  theme_classic() + 
  tx_plot_theme + 
  theme(axis.text.y=element_blank()) +
  #scale_y_discrete(labels = rep("", times = length(tx_index_key$full_taxa_name))) +
  labs(x = units_for_plot, y = "", title = "")

rm(list=ls()[!(ls() %in% c("outdir", "tx_index_key", "full_taxa_name_order", "tx_plot_theme",
                           "c_results", "p_results", "n_results", "land_results", "water_results", "wild_results",
                           "p_carbon",
                           "p_land"))])

######################################################################################################
# NITROGEN 

load(file.path(outdir, n_results))
units_for_plot = "kg N-eq per tonne"
# PRIMARY COLOR (Azote guidelines)
# interval_palette <- c("#FFEDAD", "#FFE37D", "#FFD947") # Order: light to dark
# COMPLEMENTARY COLOR
interval_palette <- c("#FFD5A9", "#FFBD79", "#FFA647")

p_nitrogen <- fit_no_na %>%
  spread_draws(tx_total_fp_w[tx]) %>%
  median_qi(.width = c(0.95, 0.8, 0.5)) %>%
  left_join(tx_index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  # Set .lower limit of plants and bivalves to be 0
  #mutate(.lower = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = .lower)) %>% 
  # REORDER taxa axis
  mutate(full_taxa_name = fct_relevel(full_taxa_name, full_taxa_name_order)) %>%
  #mutate(full_taxa_name = fct_reorder(full_taxa_name, tx_feed_fp_w)) %>%
  ggplot(aes(y = full_taxa_name, x = tx_total_fp_w)) +
  # ADD CHICKEN RESULTS:
  #geom_rect(aes(ymin = -Inf, ymax = nlevels(full_taxa_name), xmin = 77, xmax = 92), fill = "#D6CCB7") +
  geom_interval(aes(xmin = .lower, xmax = .upper), interval_size = 2.9) +
  geom_point(aes(y = full_taxa_name, x = tx_total_fp_w)) +
  geom_hline(yintercept = 1:12, linetype = "dotted") +
  scale_color_manual(values = interval_palette) +
  theme_classic() + 
  tx_plot_theme + 
  #theme(axis.text.y=element_blank()) +
  scale_y_discrete(labels = rep("", times = length(tx_index_key$full_taxa_name))) +
  # theme(axis.text.y=element_blank(),
  #       axis.title.y=element_blank(),
  #       panel.grid = element_blank(),
  #       panel.border = element_blank()) +
  labs(x = units_for_plot, y = "", title = "")

rm(list=ls()[!(ls() %in% c("outdir", "tx_index_key", "full_taxa_name_order", "tx_plot_theme",
                           "c_results", "p_results", "n_results", "land_results", "water_results", "wild_results",
                           "p_carbon",
                           "p_land",
                           "p_nitrogen"))])

######################################################################################################
# PHOSPHORUS

load(file.path(outdir, p_results))
units_for_plot = "kg P-eq per tonne"
interval_palette <- c("#FFB4C4", "#FF86A4", "#FC5185") # Order: light to dark

p_phosphorus <- fit_no_na %>%
  spread_draws(tx_total_fp_w[tx]) %>%
  median_qi(.width = c(0.95, 0.8, 0.5)) %>%
  left_join(tx_index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  # Set .lower limit of plants and bivalves to be 0
  #mutate(.lower = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = .lower)) %>% 
  # REORDER taxa axis
  mutate(full_taxa_name = fct_relevel(full_taxa_name, full_taxa_name_order)) %>%
  #mutate(full_taxa_name = fct_reorder(full_taxa_name, tx_feed_fp_w)) %>%
  ggplot(aes(y = full_taxa_name, x = tx_total_fp_w)) +
  # ADD CHICKEN RESULTS:
  #geom_rect(aes(ymin = -Inf, ymax = nlevels(full_taxa_name), xmin = 11, xmax = 16), fill = "#D6CCB7") +
  geom_interval(aes(xmin = .lower, xmax = .upper), interval_size = 2.9) +
  geom_point(aes(y = full_taxa_name, x = tx_total_fp_w)) +
  geom_hline(yintercept = 1:12, linetype = "dotted") +
  scale_color_manual(values = interval_palette) +
  theme_classic() + 
  tx_plot_theme + 
  #theme(axis.text.y=element_blank()) +
  scale_y_discrete(labels = rep("", times = length(tx_index_key$full_taxa_name))) +
  # theme(axis.text.y=element_blank(),
  #       axis.title.y=element_blank(),
  #       panel.grid = element_blank(),
  #       panel.border = element_blank()) +
  labs(x = units_for_plot, y = "", title = "")

rm(list=ls()[!(ls() %in% c("outdir", "tx_index_key", "full_taxa_name_order", "tx_plot_theme",
                           "c_results", "p_results", "n_results", "land_results", "water_results", "wild_results",
                           "p_carbon",
                           "p_land",
                           "p_nitrogen",
                           "p_phosphorus"))])

######################################################################################################
# WATER

load(file.path(outdir, water_results))
units_for_plot <- bquote('m'^3*'per tonne')
interval_palette <- c("#B2E2E6", "#80D2D7", "#3FC1C9") # Order: light to dark

# Mean total impact taxa-level
p_water <- fit_no_na %>%
  spread_draws(tx_total_fp_w[tx]) %>%
  median_qi(.width = c(0.95, 0.8, 0.5)) %>%
  left_join(tx_index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  # Set plant and bivalves distributions to 0 (both are 0 for on and off farm)
  mutate(tx_total_fp_w = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = tx_total_fp_w),
         .lower = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = .lower),
         .upper = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = .upper)) %>%
  # REORDER taxa axis
  mutate(full_taxa_name = fct_relevel(full_taxa_name, full_taxa_name_order)) %>%
  #mutate(full_taxa_name = fct_reorder(full_taxa_name, tx_feed_fp_w)) %>%
  ggplot(aes(y = full_taxa_name, x = tx_total_fp_w)) +
  # ADD CHICKEN RESULTS:
  #geom_rect(aes(ymin = -Inf, ymax = nlevels(full_taxa_name), xmin = 170, xmax = 202), fill = "#D6CCB7") +
  geom_interval(aes(xmin = .lower, xmax = .upper), interval_size = 2.9) +
  geom_point(x = 0, y = "bivalves") +
  geom_point(x = 0, y = "plants") +
  geom_point(aes(y = full_taxa_name, x = tx_total_fp_w)) +
  geom_hline(yintercept = 1:12, linetype = "dotted") +
  scale_color_manual(values = interval_palette) +
  theme_classic() + 
  tx_plot_theme + 
  labs(x = units_for_plot, y = "", title = "")

rm(list=ls()[!(ls() %in% c("outdir", "tx_index_key", "full_taxa_name_order", "tx_plot_theme",
                           "c_results", "p_results", "n_results", "land_results", "water_results", "wild_results",
                           "p_carbon",
                           "p_land",
                           "p_nitrogen",
                           "p_phosphorus", 
                           "p_water"))])

######################################################################################################
# WILD CATCH - (plot LAST because it has a different tx_index_key)

# Load model outputs individually, save to plot_grid, and create multi-panel plot
load(file.path(outdir, wild_results))
units_for_plot = bquote('kg'~CO[2]*'-eq per tonne')
interval_palette <- c("#9EA8B7", "#6A7A90", "#364F6B") # Order: light to dark

wild_index_key <- wild_dat_new_weights %>%
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
  mutate(source = "Wild")

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
  left_join(wild_index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  # Reorder from highest to lowest
  mutate(plot_taxa = as.factor(plot_taxa)) %>%
  mutate(plot_taxa = fct_reorder(plot_taxa, tx_ghg_w)) %>%
  ggplot(aes(y = plot_taxa, x = tx_ghg_w)) +
  # ADD CHICKEN RESULTS:
  geom_rect(aes(ymin = -Inf, ymax = nlevels(plot_taxa), xmin = 3165, xmax = 3627), fill = "#D6CCB7") +
  geom_interval(aes(xmin = .lower, xmax = .upper), interval_size = 3) +
  geom_point(aes(y = plot_taxa, x = tx_ghg_w)) +
  geom_hline(yintercept = 1:11, linetype = "dotted") +
  scale_color_manual(values = interval_palette) +
  theme_classic() + 
  coord_cartesian(xlim = c(0, 12500)) +
  tx_plot_theme + 
  labs(x = units_for_plot, y = "", title = "")

rm(list=ls()[!(ls() %in% c("outdir", "tx_index_key", "full_taxa_name_order", "tx_plot_theme",
                           "c_results", "p_results", "n_results", "land_results", "water_results", "wild_results",
                           "p_carbon",
                           "p_land",
                           "p_nitrogen",
                           "p_phosphorus", 
                           "p_water",
                           "p_wild"))])

# COMBINE INTO PANEL WITH FACET LABELS:
# USE: facet_grid(switch = "y") # "switch" moves y-axis label from right to left side

# 6 panel plot:
# p_carbon_facet <- p_carbon + facet_grid(rows = "source", switch = "y") + theme(strip.placement = "outside")
# p_wild_facet <- p_wild + facet_grid(rows = "source", switch = "y") + theme(strip.placement = "outside")
# p_left <- plot_grid(p_carbon_facet, p_wild_facet, ncol = 1, nrow = 2, align = "hv", labels = c("a", "d"), axis = "l")
# p_water_facet <- p_water + facet_grid(rows = "source", switch = "y") + theme(strip.placement = "outside")
# p_right <- plot_grid(p_nitrogen, p_phosphorus, p_water_facet, p_land, ncol = 2, nrow = 2, align = "h", axis = "b", labels = c("b", "c", "e", "f"), rel_widths = c(1.3, 1))
# plot_grid(p_left, p_right, ncol = 2, nrow = 1, rel_widths = c(0.55, 1))
# ggsave(filename = file.path(outdir, "plot_Figure-X_facet-labels.png"), width = 183, height = 90, units = "mm")
# ggsave(filename = file.path(outdir, "plot_Figure-X_facet-labels.tiff"), device = "tiff", width = 183, height = 90, units = "mm")

# Use NO FACET version below for 5-panel plot
# 5 panel plot
p_left <- plot_grid(p_carbon, p_water, nrow = 2, ncol = 1, labels = c("a", "d"))
p_right <- plot_grid(p_nitrogen, p_phosphorus, p_land, ncol = 2, nrow = 2, align = "h", labels = c("b", "c", "e"))
plot_grid(p_left, p_right, ncol = 2, nrow = 1, align = "v", rel_widths = c(0.6, 1))
ggsave(filename = file.path(outdir, "plot_Figure-X.png"), width = 183, height = 90, units = "mm")
ggsave(filename = file.path(outdir, "plot_Figure-X.tiff"), device = "tiff", width = 183, height = 90, units = "mm")


