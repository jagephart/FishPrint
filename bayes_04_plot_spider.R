# Author: Kelvin Gorospe
# Create spiderplot to compare aquaculture taxa across environmental impacts

rm(list=ls())
library(tidyverse)
library(tidybayes)
library(fmsb) # for spider (aka radial) plots
library(ggplot2)
library(cowplot) ## for plot_grid
library(scales) # show_col()
library(colorspace) # adjust_transparency
datadir <- "/Volumes/jgephart/BFA Environment 2/Data"
outdir <- "/Volumes/jgephart/BFA Environment 2/Outputs"

# Set filenames:
# Mass allocation
c_results <- "Nature-submitted-2021-08/RERUN-GHG-and-Water-Bayesian-Means-Edible-Weight/PRIORS/2021-08-17_full-model-posterior_Global warming potential_Mass-allocation.RData"
n_results <- "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Nitrogen/2021-04-28_full-model-posterior_Marine eutrophication_Mass-allocation.RData"
p_results <- "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Phosphorus/2021-04-28_full-model-posterior_Freshwater eutrophication_Mass-allocation.RData"
land_results <- "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Land-FCR-Priors-Only/2021-04-28_full-model-posterior_Land Use_Mass-allocation.RData"
water_results <- "Nature-submitted-2021-08/RERUN-GHG-and-Water-Bayesian-Means-Edible-Weight/PRIORS/2021-08-17_full-model-posterior_Water Consumption_Mass-allocation.RData"
wild_results <- "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Wild/2021-04-27_full-model-posterior_Wild-capture-ghg.RData"

################################################################################################################
# GHG
################################################################################################################
load(file.path(outdir, c_results))

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
         full_taxa_name = as.factor(full_taxa_name))

# Get median of mean GHG impact
p_carbon <- fit_no_na %>%
  spread_draws(tx_total_fp_w[tx]) %>%
  #median_qi(.width = c(0.95, 0.8, 0.5)) %>%
  median_qi() %>%
  left_join(tx_index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  # Filter out tx_total_fp_w for plants and bivalves
  filter(taxa %in% c("bivalves", "plants")==FALSE) %>%
  # Substitute on_farm data as total impact for plants and bivalves
  bind_rows(fit_no_na %>% 
              spread_draws(tx_farm_fp_w[tx]) %>% 
              #median_qi(.width = c(0.95, 0.8, 0.5)) %>% 
              median_qi() %>%
              left_join(tx_index_key, by = "tx") %>% 
              filter(taxa %in% c("bivalves", "plants")) %>% 
              rename(tx_total_fp_w = tx_farm_fp_w)) %>%
  select(full_taxa_name, tx_total_fp_w) %>%
  rename(carbon = tx_total_fp_w)
  
# Clear before next model
rm(list=ls()[!(ls() %in% c("outdir", "tx_index_key", "full_taxa_name_order", "tx_plot_theme",
                           "c_results", "p_results", "n_results", "land_results", "water_results", "wild_results",
                           "p_carbon"))])

################################################################################################################
# LAND
################################################################################################################
load(file.path(outdir, land_results))

p_land <- fit_no_na %>%
  spread_draws(tx_land_total_w[tx]) %>%
  #median_qi(.width = c(0.95, 0.8, 0.5)) %>%
  median_qi() %>%
  left_join(tx_index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  # Set plant and bivalves distributions to 0 (both are 0 for on and off farm)
  mutate(tx_land_total_w = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = tx_land_total_w),
         .lower = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = .lower),
         .upper = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = .upper)) %>%
  # LAND can't be less than zero: (lower bound for misc diad fishes is negative)
  mutate(.lower = if_else(.lower < 0, true = 0, false = .lower)) %>%
  select(full_taxa_name, tx_land_total_w) %>%
  rename(land = tx_land_total_w)

rm(list=ls()[!(ls() %in% c("outdir", "tx_index_key", "full_taxa_name_order", "tx_plot_theme",
                           "c_results", "p_results", "n_results", "land_results", "water_results", "wild_results",
                           "p_carbon",
                           "p_land"))])

################################################################################################################
# NITROGEN 
################################################################################################################
load(file.path(outdir, n_results))

p_nitrogen <- fit_no_na %>%
  spread_draws(tx_total_fp_w[tx]) %>%
  #median_qi(.width = c(0.95, 0.8, 0.5)) %>%
  median_qi() %>%
  left_join(tx_index_key, by = "tx") %>%
  select(full_taxa_name, tx_total_fp_w) %>%
  rename(nitrogen = tx_total_fp_w)

rm(list=ls()[!(ls() %in% c("outdir", "tx_index_key", "full_taxa_name_order", "tx_plot_theme",
                           "c_results", "p_results", "n_results", "land_results", "water_results", "wild_results",
                           "p_carbon",
                           "p_land",
                           "p_nitrogen"))])

################################################################################################################
# PHOSPHORUS
################################################################################################################
load(file.path(outdir, p_results))

p_phosphorus <- fit_no_na %>%
  spread_draws(tx_total_fp_w[tx]) %>%
  #median_qi(.width = c(0.95, 0.8, 0.5)) %>%
  median_qi() %>%
  left_join(tx_index_key, by = "tx") %>%
  select(full_taxa_name, tx_total_fp_w) %>%
  rename(phosphorus = tx_total_fp_w)

rm(list=ls()[!(ls() %in% c("outdir", "tx_index_key", "full_taxa_name_order", "tx_plot_theme",
                           "c_results", "p_results", "n_results", "land_results", "water_results", "wild_results",
                           "p_carbon",
                           "p_land",
                           "p_nitrogen",
                           "p_phosphorus"))])

################################################################################################################
# WATER
################################################################################################################
load(file.path(outdir, water_results))

# Mean total impact taxa-level
p_water <- fit_no_na %>%
  spread_draws(tx_total_fp_w[tx]) %>%
  #median_qi(.width = c(0.95, 0.8, 0.5)) %>%
  median_qi() %>%
  left_join(tx_index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  # Set plant and bivalves distributions to 0 (both are 0 for on and off farm)
  mutate(tx_total_fp_w = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = tx_total_fp_w),
         .lower = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = .lower),
         .upper = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = .upper)) %>%
  select(full_taxa_name, tx_total_fp_w) %>%
  rename(water = tx_total_fp_w)

rm(list=ls()[!(ls() %in% c("outdir", "tx_index_key", "full_taxa_name_order", "tx_plot_theme",
                           "c_results", "p_results", "n_results", "land_results", "water_results", "wild_results",
                           "p_carbon",
                           "p_land",
                           "p_nitrogen",
                           "p_phosphorus", 
                           "p_water"))])

################################################################################################################
# Format data for spider plot
################################################################################################################

# Raw data:
# raw_impacts <- p_carbon %>%
#   left_join(p_land) %>%
#   left_join(p_nitrogen) %>%
#   left_join(p_phosphorus) %>%
#   left_join(p_water)
# radarchart(raw_impacts)

# Normalize vs Scale data?
# NORMALIZE:
# Note for nitrogen and phosphorus, since there are negative values, do a "min-max" normalization
normalize_carbon <- p_carbon %>% 
  mutate(carbon = carbon/max(carbon)) 
normalize_land <- p_land %>% 
  mutate(land = land/max(land)) 
normalize_nitrogen <- p_nitrogen %>% 
  mutate(nitrogen = (nitrogen - min(nitrogen)) / (max(nitrogen) - min(nitrogen)))
normalize_phosphorus <- p_phosphorus %>% 
  mutate(phosphorus = (phosphorus - min(phosphorus)) / (max(phosphorus) - min(phosphorus)))
normalize_water <- p_water %>% 
  mutate(water = water/max(water))

normalize_impacts <- normalize_carbon %>%
  left_join(normalize_land) %>%
  left_join(normalize_nitrogen) %>%
  left_join(normalize_phosphorus) %>%
  left_join(normalize_water) %>%
  # RENAME columns to final axis labels
  rename_with(~str_to_sentence(.), everything()) %>%
  rename(GHG = Carbon)

# To use the fmsb package, I have to add 2 lines to the dataframe: the max and min of each variable to show on the plot!
normalize_impacts_ordered <- rbind(rep(1,5), rep(0, 5), (normalize_impacts %>%
                                                           mutate(total = rowSums(normalize_impacts[,-1])) %>% # arrange in descending total impacts to maximize visibility of all polygons
                                                           arrange(desc(total)) %>%
                                                           column_to_rownames(names(.)[1]) %>% # Get the first column and move to rownames
                                                           select(-total)))

radarchart(normalize_impacts_ordered)
title("Normalized impacts")

# SCALE:
scale_carbon <- p_carbon %>%
  mutate(carbon = scale(carbon, scale = 2*sd(p_carbon$carbon)))
scale_land <- p_land %>%
  mutate(land = scale(land, scale = 2*sd(p_land$land)))
scale_nitrogen <- p_nitrogen %>%
  mutate(nitrogen = scale(nitrogen, scale = 2*sd(p_nitrogen$nitrogen)))
scale_phosphorus <- p_phosphorus %>%
  mutate(phosphorus = scale(phosphorus, scale = 2*sd(p_phosphorus$phosphorus)))
scale_water <- p_water %>%
  mutate(water = scale(water, scale = 2*sd(p_water$water)))

scale_impacts <- scale_carbon %>%
  left_join(scale_land) %>%
  left_join(scale_nitrogen) %>%
  left_join(scale_phosphorus) %>%
  left_join(scale_water) %>%
  # RENAME columns to final axis labels
  rename_with(~str_to_sentence(.), everything()) %>%
  rename(GHG = Carbon)

# To use the fmsb package, I have to add 2 lines to the dataframe: the max and min of each variable to show on the plot!
scale_impacts_ordered <- rbind(rep(1,5), rep(-1, 5), (scale_impacts %>%
                                                        mutate(total = rowSums(scale_impacts[,-1])) %>% # arrange in descending total impacts to maximize visibility of all polygons
                                                        arrange(desc(total)) %>%
                                                        column_to_rownames(names(.)[1]) %>% # Get the first column and move to rownames
                                                        select(-total)))

radarchart(scale_impacts_ordered)
title("Scaled impacts")

################################################################################################################
# Final plots
################################################################################################################

# Specify colors
#364F6B # dark blue
#3FC1C9 # light blue
#A9D158 # light-green
#FFD947 # yellow
#FC5185 # pink
#A69569 # tan
#FFA647 # orange
#C93F3F # red
#70468C # purple
#B389ED # lavender
#57D182 # seagreen
#808080 # gray

color_border <- c("#364F6B", # dark blue
  "#3FC1C9", # light blue
  "#A9D158", # light-green
  "#FFD947", # yellow
  "#FC5185", # pink
  "#A69569", # tan
  "#FFA647", # orange
  "#C93F3F", # red
  "#70468C", # purple
  "#B389ED", # lavender
  "#57D182", # seagreen
  "#808080" # gray
)


color_fill <- adjust_transparency(color_border, alpha = 0.4)

# Normalized impacts
radarchart(normalize_impacts_ordered,
           # customize polygons
           pcol = color_border,
           #pfcol = color_fill, # Remove fill (plot too busy)
           plwd = 4,
           plty = 1,
           # customize grid
           cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,20,5), cglwd=0.8,
           # customize labels
           vlcex=0.8)
# Add legend
legend(x=0.7, y=1.3, legend = rownames(normalize_impacts_ordered[-c(1,2),]), bty = "n", pch=20 , col=color_border , text.col = "black", cex=0.7, pt.cex=3)
#dev.off()
# FIX IT - ggsave not working (save from Plot console)
#ggsave(filename = file.path(outdir, "plot_radar_normalized.png"), width = 11, height = 8.5)


# Scaled impacts
radarchart(scale_impacts_ordered,
           # customize polygons
           pcol = color_border,
           #pfcol = color_fill,
           plwd = 4,
           plty = 1,
           # customize grid
           cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,20,5), cglwd=0.8,
           # customize labels
           vlcex=0.8)
# Add legend
legend(x=0.7, y=1.3, legend = rownames(scale_impacts_ordered[-c(1,2),]), bty = "n", pch=20 , col=color_border , text.col = "black", cex=0.7, pt.cex=3)
#dev.off()
# FIX IT - ggsave not working (save from Plot console)
#ggsave(filename = file.path(outdir, "plot_radar_scaled.png"), width = 11, height = 8.5)
