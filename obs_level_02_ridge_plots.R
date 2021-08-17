# Non-Bayesian estimate figures
rm(list=ls())
library(tidyverse)
library(ggplot2)
library(ggthemes)
library(ggridges)
library(cowplot)
library(egg)

#_______________________________________________________________________________________________________________________#
# Read in data
#_______________________________________________________________________________________________________________________#
# Set data directories
datadir <- "/Volumes/jgephart/BFA Environment 2/Data"
outdir <- "/Volumes/jgephart/BFA Environment 2/Outputs"

# edible weight
stressor_taxa_summary <- read.csv(file.path(outdir,"non-bayes-stressors_farmed_taxa-level_edible-weight.csv"))
stressor_species_summary <- read.csv(file.path(outdir, "non-bayes-stressors_farmed_observation-level_edible-weight.csv"))
df_capture_ghg <- read.csv(file.path(outdir, "non-bayes-stressors_capture_observation-level_edible-weight.csv")) 

#_______________________________________________________________________________________________________________________#
# Summary plots
#_______________________________________________________________________________________________________________________#

# SET UNITS:
weight_type <- "edible weight"

# Single line options
# units_for_ghg <- bquote('kg'~CO[2]*' t'^-1~.(weight_type)) 
# units_for_land <- bquote('m'^2*' t'^-1~.(weight_type))
# units_for_nitrogen <- bquote('kg N-eq t'^-1~.(weight_type))
# units_for_phosphorus <- bquote('kg P-eq t'^-1~.(weight_type))
# units_for_water <- bquote('m'^3*' t'^-1~.(weight_type))

units_for_ghg <- bquote(atop('kg'~CO[2]*' t'^-1~phantom(), .(weight_type))) # atop + phantom() to create line break
units_for_land <- bquote(atop('m'^2*' t'^-1~phantom(), .(weight_type)))
units_for_nitrogen <- bquote(atop('kg N-eq t'^-1~phantom(), .(weight_type)))
units_for_phosphorus <- bquote(atop('kg P-eq t'^-1~phantom(), .(weight_type)))
units_for_water <- bquote(atop('m'^3*' t'^-1~phantom(), .(weight_type)))

# SET THEME:
ridge_plot_theme <- theme(plot.title = element_text(hjust = 0),
                          axis.title.x = element_text(size = 8, hjust = 0.5),
                          axis.text.x = element_text(hjust = 0.5, size = 8, color = "black"),
                          axis.text.y = element_text(size = 12, color = "black"),
                          legend.position = "none",
                          plot.margin = unit(c(0, 1, 0, 0), "mm")) # top right bottom left

# theme(axis.text.x = element_text(angle = 90, hjust = 1),
#       plot.margin=unit(plot_margins, "cm"))

# Data distributions by taxa
stressor_species_summary_plot <- stressor_species_summary %>% 
  pivot_longer(feed_GHG:onfarm_water, names_sep = "_", 
               names_to = c("source", "stressor")) %>%
  filter(!is.na(value)) 

stressor_species_summary_plot$stressor[stressor_species_summary_plot$stressor == "GHG"] <- "ghg"


# Total stressor ridge plots 
total_stressor_species_summary_plot <- stressor_species_summary_plot %>%
  group_by(study_id, taxa,  intensity, system,  clean_sci_name, stressor) %>%
  summarise(total_value = sum(value))

# Set factor levels

# Order taxa levels to same order as Figure 1
full_taxa_order <- c("seaweeds", 
                          "bivalves", 
                          "silver/bighead", 
                          "salmon", 
                          "trout", 
                          "misc carp", 
                          "catfish", 
                          "milkfish", 
                          "shrimp", 
                          "tilapia", 
                          "misc marine", 
                          "misc diad")
full_intensity_order <- c("extensive",
                          "improved extensive",
                          "semi-intensive",
                          "intensive")
full_system_order <- c("on- & off-bottom", "cages & pens", "ponds", "recirculating & tanks")

total_stressor_species_summary_plot <- total_stressor_species_summary_plot %>%
  # Create taxa names for plot
  mutate(taxa = as.character(taxa),
         full_taxa = case_when(taxa == "hypoph_carp" ~ "silver/bighead",
                               taxa == "oth_carp" ~ "misc carp",
                               taxa == "misc_diad" ~ "misc diad",
                               taxa == "misc_fresh" ~ "misc freshwater",
                               taxa == "misc_marine" ~ "misc marine",
                               taxa == "fresh_crust" ~ "freshwater crust",
                               taxa == "plants" ~ "seaweeds",
                               TRUE ~ taxa),
         taxa = as.factor(taxa),
         full_taxa = as.factor(full_taxa)) %>%
  # Set taxa order for plot
  mutate(full_taxa = fct_relevel(full_taxa, full_taxa_order)) %>%
  # Create intesnity names for plot
  mutate(full_intensity = case_when(intensity == "Imp. extensive" ~ "improved extensive",
                                    TRUE ~ intensity),
         full_intensity = str_to_lower(full_intensity),
         full_intensity = as.factor(full_intensity)) %>%
  # Set intensity order for plot
  mutate(full_intensity = fct_relevel(full_intensity, full_intensity_order)) %>%
  # Create system names for plot
  mutate(full_system = case_when(str_detect(system, "and") ~ str_replace(system, pattern = "and", replacement = "&"),
                                 TRUE ~ system),
         full_system = str_to_lower(full_system),
         full_system = as.factor(full_system)) %>%
  # Set system order for plot
  mutate(full_system = fct_relevel(full_system, full_system_order))

# ghg
ghg_ridge_intensity <- ggplot(total_stressor_species_summary_plot %>% 
                                filter(stressor == "ghg" & !is.na(intensity)),
       aes(x = total_value, y = full_intensity)) +
  #geom_density_ridges(alpha = 0.7, fill = "#364F6B") +
  scale_x_continuous(labels = c("0", "1e4", "2e4", "3e4"), breaks = c(0, 10000, 20000, 30000)) +
  geom_density_ridges(alpha = 0.7, fill = "#364F6B", scale = 0.9) +
  labs(title = "GHG", x = units_for_ghg, y = "") +
  theme_ridges() +
  ridge_plot_theme

ghg_ridge_system <- ggplot(total_stressor_species_summary_plot %>% 
                                filter(stressor == "ghg" & !is.na(system)),
                              aes(x = total_value, y = full_system)) +
  #geom_density_ridges(alpha = 0.7, fill = "#364F6B") +
  geom_density_ridges(alpha = 0.7, fill = "#364F6B", scale = 0.9) +
  scale_x_continuous(labels = c("0", "1e4", "2e4", "3e4"), breaks = c(0, 10000, 20000, 30000)) +
  labs(title = "GHG", x = units_for_ghg, y = "") +
  theme_ridges() +
  ridge_plot_theme

ghg_ridge_taxa <- ggplot(total_stressor_species_summary_plot %>% 
                             filter(stressor == "ghg" & !is.na(system)),
                           aes(x = total_value, y = full_taxa)) +
  #geom_density_ridges(alpha = 0.7, fill = "#364F6B") +
  geom_density_ridges(alpha = 0.7, fill = "#364F6B", scale = 0.9) +
  scale_x_continuous(labels = c("0", "1e4", "2e4", "3e4"), breaks = c(0, 10000, 20000, 30000)) +
  scale_y_discrete(labels = rev(c("misc diad", "misc marine", "shrimp", "milkfish",
                                  "tilapia", "catfish", "misc carp", "trout", "salmon",
                                  "silver/bighead", "seaweeds", "bivalves"))) +
  labs(title = "GHG", x = units_for_ghg, y = "") +
  theme_ridges() +
  ridge_plot_theme

# N
N_ridge_intensity <- ggplot(total_stressor_species_summary_plot %>% 
                                filter(stressor == "N" & !is.na(intensity)),
                              aes(x = total_value, y = full_intensity)) +
  #geom_density_ridges(alpha = 0.7, fill = "#FFA647") +
  geom_density_ridges(alpha = 0.7, fill = "#FFA647", scale = 0.9) +
  labs(title = "Nitrogen", x = units_for_nitrogen, y = "") +
  theme_ridges() +
  ridge_plot_theme +
  theme(axis.text.y=element_blank())

N_ridge_system <- ggplot(total_stressor_species_summary_plot %>% 
                             filter(stressor == "N" & !is.na(system)),
                           aes(x = total_value, y = full_system)) +
  #geom_density_ridges(alpha = 0.7, fill = "#FFA647") +
  geom_density_ridges(alpha = 0.7, fill = "#FFA647", scale = 0.9) +
  labs(title = "Nitrogen", x = units_for_nitrogen, y = "") +
  theme_ridges() +
  ridge_plot_theme +
  theme(axis.text.y=element_blank())

N_ridge_taxa <- ggplot(total_stressor_species_summary_plot %>% 
                           filter(stressor == "N" & !is.na(system)),
                         aes(x = total_value, y = full_taxa)) +
  #geom_density_ridges(alpha = 0.7, fill = "#FFA647") +
  geom_density_ridges(alpha = 0.7, fill = "#FFA647", scale = 0.9) +
  labs(title = "Nitrogen", x = units_for_nitrogen, y = "") +
  theme_ridges() +
  ridge_plot_theme +
  theme(axis.text.y=element_blank())

# P
P_ridge_intensity <- ggplot(total_stressor_species_summary_plot %>% 
                                filter(stressor == "P" & !is.na(intensity)),
                              aes(x = total_value, y = full_intensity)) +
  #geom_density_ridges(alpha = 0.7, fill = "#FC5185") +
  geom_density_ridges(alpha = 0.7, fill = "#FC5185", scale = 0.9) +
  labs(title = "Phosphorus", x = units_for_phosphorus, y = "") +
  theme_ridges() +
  ridge_plot_theme +
  theme(axis.text.y=element_blank())

P_ridge_system <- ggplot(total_stressor_species_summary_plot %>% 
                             filter(stressor == "P" & !is.na(system)),
                           aes(x = total_value, y = full_system)) +
  #geom_density_ridges(alpha = 0.7, fill = "#FC5185") +
  geom_density_ridges(alpha = 0.7, fill = "#FC5185", scale = 0.9) +
  labs(title = "Phosphorus", x = units_for_phosphorus, y = "") +
  theme_ridges() +
  ridge_plot_theme +
  theme(axis.text.y=element_blank())

P_ridge_taxa <- ggplot(total_stressor_species_summary_plot %>% 
                           filter(stressor == "P" & !is.na(system)),
                         aes(x = total_value, y = full_taxa)) +
  #geom_density_ridges(alpha = 0.7, fill = "#FC5185") +
  geom_density_ridges(alpha = 0.7, fill = "#FC5185", scale = 0.9) +
  labs(title = "Phosphorus", x = units_for_phosphorus, y = "") +
  theme_ridges() +
  ridge_plot_theme +
  theme(axis.text.y=element_blank())

# land
land_ridge_intensity <- ggplot(total_stressor_species_summary_plot %>% 
                                filter(stressor == "land" & !is.na(intensity)),
                              aes(x = total_value, y = full_intensity)) +
  #geom_density_ridges(alpha = 0.7, fill = "#57D182") +
  geom_density_ridges(alpha = 0.7, fill = "#57D182", scale = 0.9) +
  scale_x_continuous(labels = c("0", "5e4", "1e5"), breaks = c(0, 50000, 100000)) +
  labs(title = "Land", x = units_for_land, y = "") +
  theme_ridges() +
  ridge_plot_theme +
  theme(axis.text.y=element_blank())

land_ridge_system <- ggplot(total_stressor_species_summary_plot %>% 
                             filter(stressor == "land" & !is.na(system)),
                           aes(x = total_value, y = full_system)) +
  #geom_density_ridges(alpha = 0.7, fill = "#57D182") +
  geom_density_ridges(alpha = 0.7, fill = "#57D182", scale = 0.9) +
  labs(title = "Land", x = units_for_land, y = "") +
  scale_x_continuous(labels = c("0", "5e4", "1e5"), breaks = c(0, 50000, 100000)) +
  theme_ridges() +
  ridge_plot_theme +
  theme(axis.text.y=element_blank())

land_ridge_taxa <- ggplot(total_stressor_species_summary_plot %>% 
                              filter(stressor == "land" & !is.na(system)),
                            aes(x = total_value, y = full_taxa)) +
  #geom_density_ridges(alpha = 0.7, fill = "#57D182") +
  geom_density_ridges(alpha = 0.7, fill = "#57D182", scale = 0.9) +
  scale_x_continuous(labels = c("0", "5e4", "1e5"), breaks = c(0, 50000, 100000)) +
  labs(title = "Land", x = units_for_land, y = "") +
  theme_ridges() +
  ridge_plot_theme +
  theme(axis.text.y=element_blank())

# water
water_ridge_intensity <- ggplot(total_stressor_species_summary_plot %>% 
                                filter(stressor == "water" & !is.na(intensity)),
                              aes(x = total_value, y = full_intensity)) +
  #geom_density_ridges(alpha = 0.7, fill ="#3FC1C9") +
  geom_density_ridges(alpha = 0.7, fill ="#3FC1C9", scale = 0.9) +
  scale_x_continuous(labels = c("0", "1e4", "2e4", "3e4"), breaks = c(0, 10000, 20000, 30000)) +
  labs(title = "Water", x = units_for_water, y = "") +
  theme_ridges() +
  ridge_plot_theme +
  theme(axis.text.y=element_blank(),
        axis.text.x=element_text(hjust = 0.7))

water_ridge_system <- ggplot(total_stressor_species_summary_plot %>% 
                             filter(stressor == "water" & !is.na(system)),
                           aes(x = total_value, y = full_system)) +
  #geom_density_ridges(alpha = 0.7, fill ="#3FC1C9") +
  geom_density_ridges(alpha = 0.7, fill ="#3FC1C9", scale = 0.9) +
  scale_x_continuous(labels = c("0", "1e4", "2e4"), breaks = c(0, 10000, 20000)) +
  labs(title = "Water", x = units_for_water, y = "") +
  theme_ridges() +
  ridge_plot_theme +
  theme(axis.text.y=element_blank())

water_ridge_taxa <- ggplot(total_stressor_species_summary_plot %>% 
                               filter(stressor == "water" & !is.na(system)),
                             aes(x = total_value, y = full_taxa)) +
  #geom_density_ridges(alpha = 0.7, fill ="#3FC1C9") +
  geom_density_ridges(alpha = 0.7, fill ="#3FC1C9", scale = 0.9) +
  scale_x_continuous(labels = c("0", "1e4", "2e4"), breaks = c(0, 10000, 20000)) +
  labs(title = "Water", x = units_for_water, y = "") +
  theme_ridges() +
  ridge_plot_theme +
  theme(axis.text.y=element_blank())

png(file.path(outdir, "scaled-ridge-plots_intensity.png"), width = 8, height = 4, units = "in", res = 300)
#pdf(file = file.path(outdir, "scaled-ridge-plots_intensity.pdf"), width = 8, height = 4) 
ggarrange(ghg_ridge_intensity, N_ridge_intensity, P_ridge_intensity,
          land_ridge_intensity, water_ridge_intensity, nrow = 1)
dev.off()

png(file.path(outdir, "scaled-ridge-plots_system.png"), width = 8, height = 4, units = "in", res = 300)
#pdf(file = file.path(outdir, "scaled-ridge-plots_system.pdf"), width = 8, height = 4) 
ggarrange(ghg_ridge_system, N_ridge_system, P_ridge_system,
          land_ridge_system, water_ridge_system, nrow = 1)
dev.off()

png(file.path(outdir, "scaled-ridge-plots_taxa.png"), width = 8, height = 4, units = "in", res = 300)
#pdf(file = file.path(outdir, "scaled-ridge-plots_taxa.pdf"), width = 8, height = 4) 
ggarrange(ghg_ridge_taxa, N_ridge_taxa, P_ridge_taxa,
          land_ridge_taxa, water_ridge_taxa, nrow = 1)
dev.off()


#_______________________________________________________________________________________________________________________#
# Plot capture GHG
#_______________________________________________________________________________________________________________________#
ggplot(df_capture_ghg, aes(x = ghg_kg_t, y = species_group)) +
  geom_bar(stat="identity") +
  labs(title = "Capture Fishery GHGs", x = "CO2-eq (kg/t)", y = "") +
  theme_clean()

