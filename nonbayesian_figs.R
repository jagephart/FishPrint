# Non-Bayesian estimate figures
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

stressor_taxa_summary <- read.csv(file.path(datadir,"20210107_aquaculture_stressors_nonbayes.csv"))
stressor_species_summary <- read.csv(file.path(datadir, "20210107_stressor_species_summary.csv"))
df_capture_ghg <- read.csv(file.path(datadir, "20210107_capture_stressors_nonbayes.csv"))

#_______________________________________________________________________________________________________________________#
# Summary plots
#_______________________________________________________________________________________________________________________#
# Data distributions by taxa
stressor_species_summary_plot <- stressor_species_summary %>% 
  pivot_longer(feed_GHG:onfarm_water, names_sep = "_", 
               names_to = c("source", "stressor")) %>%
  filter(!is.na(value)) 

stressor_species_summary_plot$stressor[stressor_species_summary_plot$stressor == "GHG"] <- "ghg"

# Plot GHG by taxa
ggplot(stressor_species_summary_plot %>% filter(stressor %in% c("GHG", "ghg")),
       aes(x = value, y = taxa)) +
  geom_density_ridges() +
  theme_ridges() + 
  labs(x = "GHG (kg CO2-eq per t)", y = "") +
  facet_wrap(~source, scales = "free")

# Plot GHG by system
ggplot(stressor_species_summary_plot %>% filter(stressor %in% c("GHG", "ghg"), !is.na(system)),
       aes(x = value, y = system)) +
  geom_density_ridges() +
  theme_ridges() + 
  labs(x = "GHG (kg CO2-eq per t)", y = "") +
  facet_wrap(~source, scales = "free")

# Plot GHG by system
ggplot(stressor_species_summary_plot %>% filter(stressor %in% c("GHG", "ghg"), !is.na(intensity)),
       aes(x = value, y = intensity)) +
  geom_density_ridges() +
  theme_ridges() + 
  labs(x = "GHG (kg CO2-eq per t)", y = "") +
  facet_wrap(~source, scales = "free")

# Plot N by taxa
ggplot(stressor_species_summary_plot %>% filter(stressor %in% c("N")),
       aes(x = value, y = taxa)) +
  geom_density_ridges() +
  theme_ridges() + 
  labs(x = "N (kg N per t)", y = "") +
  facet_wrap(~source, scales = "free")

# Plot N by system
ggplot(stressor_species_summary_plot %>% filter(stressor %in% c("N"), !is.na(system)),
       aes(x = value, y = system)) +
  geom_density_ridges() +
  theme_ridges() + 
  labs(x = "N (kg N per t)", y = "") +
  facet_wrap(~source, scales = "free")

# Plot N by taxa
ggplot(stressor_species_summary_plot %>% filter(stressor %in% c("N"), !is.na(intensity)),
       aes(x = value, y = intensity)) +
  geom_density_ridges() +
  theme_ridges() + 
  labs(x = "N (kg N per t)", y = "") +
  facet_wrap(~source, scales = "free")

# Plot P by taxa
ggplot(stressor_species_summary_plot %>% filter(stressor %in% c("P")),
       aes(x = value, y = taxa)) +
  geom_density_ridges() +
  theme_ridges() + 
  labs(x = "P (kg P per t)", y = "") +
  facet_wrap(~source, scales = "free")

# Plot P by system
ggplot(stressor_species_summary_plot %>% filter(stressor %in% c("P"), !is.na(system)),
       aes(x = value, y = system)) +
  geom_density_ridges() +
  theme_ridges() + 
  labs(x = "P (kg P per t)", y = "") +
  facet_wrap(~source, scales = "free")

# Plot P by intensity
ggplot(stressor_species_summary_plot %>% filter(stressor %in% c("P"), !is.na(intensity)),
       aes(x = value, y = intensity)) +
  geom_density_ridges() +
  theme_ridges() + 
  labs(x = "P (kg P per t)", y = "") +
  facet_wrap(~source, scales = "free")

# Plot land by taxa
ggplot(stressor_species_summary_plot %>% filter(stressor %in% c("land")),
       aes(x = value, y = taxa)) +
  geom_density_ridges() +
  theme_ridges() + 
  labs(x = "Land (m2 per t)", y = "") +
  facet_wrap(~source, scales = "free")

# Plot land by system
ggplot(stressor_species_summary_plot %>% filter(stressor %in% c("land"), !is.na(system)),
       aes(x = value, y = system)) +
  geom_density_ridges() +
  theme_ridges() + 
  labs(x = "Land (m2 per t)", y = "") +
  facet_wrap(~source, scales = "free")

# Plot land by intensity
ggplot(stressor_species_summary_plot %>% filter(stressor %in% c("land"), !is.na(intensity)),
       aes(x = value, y = intensity)) +
  geom_density_ridges() +
  theme_ridges() + 
  labs(x = "Land (m2 per t)", y = "") +
  facet_wrap(~source, scales = "free")

# Plot water by taxa
ggplot(stressor_species_summary_plot %>% filter(stressor %in% c("water")),
       aes(x = value, y = taxa)) +
  geom_density_ridges() +
  labs(x = "Water (m3 per t)", y = "") +
  theme_ridges() + 
  facet_wrap(~source, scales = "free")

# Plot water by system
ggplot(stressor_species_summary_plot %>% filter(stressor %in% c("water"), !is.na(system)),
       aes(x = value, y = system)) +
  geom_density_ridges() +
  labs(x = "Water (m3 per t)", y = "") +
  theme_ridges() + 
  facet_wrap(~source, scales = "free")

# Plot water by intensity
ggplot(stressor_species_summary_plot %>% filter(stressor %in% c("water"), !is.na(intensity)),
       aes(x = value, y = intensity)) +
  geom_density_ridges() +
  labs(x = "Water (m3 per t)", y = "") +
  theme_ridges() + 
  facet_wrap(~source, scales = "free")

# Mean stressor by taxa
stressor_taxa_summary_plot <- stressor_taxa_summary %>%
  select(-c("total_ghg", "total_N", "total_P", "total_land", "total_water", 
            "onfarm_GHG_unweighted_stressor", "onfarm_N_unweighted_stressor",
            "onfarm_P_unweighted_stressor", "onfarm_land_unweighted_stressor", "onfarm_water_unweighted_stressor")) %>%
  pivot_longer(feed_P:onfarm_water_weighted_stressor, names_sep = "_", 
               names_to = c("source", "stressor", "drop_1", "drop_2")) %>%
  select(-c("drop_1", "drop_2"))
  
  
ggplot(stressor_taxa_summary_plot, aes(x = value, y = taxa, fill = source)) + 
  geom_bar(position="stack", stat="identity") +
  facet_wrap(~stressor, scales = "free") +
  theme_clean()

# Total stressor ridge plots 
total_stressor_species_summary_plot <- stressor_species_summary_plot %>%
  group_by(study_id, taxa,  intensity, system,  clean_sci_name, stressor) %>%
  summarise(total_value = sum(value))

# Set factor levels
stressor_species_summary_plot$intensity <- factor(stressor_species_summary_plot$intensity, 
                                                  levels = c("Extensive", "Imp. extensive", 
                                                             "Semi-intensive", "Intensive", "NA"))
stressor_species_summary_plot$system <- factor(stressor_species_summary_plot$system,
                                               levels = c("On- and off-bottom", "Cages & pens", 
                                                          "Ponds", "Recirculating and tanks", "NA"))

# ghg
ghg_ridge_intensity <- ggplot(total_stressor_species_summary_plot %>% 
                                filter(stressor == "ghg" & !is.na(intensity)),
       aes(x = total_value, y = intensity)) +
  geom_density_ridges() +
  labs(title = "GHG", x = "kg CO2-eq/t", y = "") +
  theme_ridges() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ghg_ridge_system <- ggplot(total_stressor_species_summary_plot %>% 
                                filter(stressor == "ghg" & !is.na(system)),
                              aes(x = total_value, y = system)) +
  geom_density_ridges() +
  labs(title = "GHG", x = "kg CO2-eq/t", y = "") +
  theme_ridges() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ghg_ridge_taxa <- ggplot(total_stressor_species_summary_plot %>% 
                             filter(stressor == "ghg" & !is.na(system)),
                           aes(x = total_value, y = taxa)) +
  geom_density_ridges() +
  labs(title = "GHG", x = "kg CO2-eq/t", y = "") +
  theme_ridges() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# N
N_ridge_intensity <- ggplot(total_stressor_species_summary_plot %>% 
                                filter(stressor == "N" & !is.na(intensity)),
                              aes(x = total_value, y = intensity)) +
  geom_density_ridges() +
  labs(title = "Nitrogen", x = "kg N/t", y = "") +
  theme_ridges() +
  theme(axis.text.y=element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))

N_ridge_system <- ggplot(total_stressor_species_summary_plot %>% 
                             filter(stressor == "N" & !is.na(system)),
                           aes(x = total_value, y = system)) +
  geom_density_ridges() +
  labs(title = "Nitrogen", x = "kg N/t", y = "") +
  theme_ridges() +
  theme(axis.text.y=element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))

N_ridge_taxa <- ggplot(total_stressor_species_summary_plot %>% 
                           filter(stressor == "N" & !is.na(system)),
                         aes(x = total_value, y = taxa)) +
  geom_density_ridges() +
  labs(title = "Nitrogen", x = "kg N/t", y = "") +
  theme_ridges() +
  theme(axis.text.y=element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))

# P
P_ridge_intensity <- ggplot(total_stressor_species_summary_plot %>% 
                                filter(stressor == "P" & !is.na(intensity)),
                              aes(x = total_value, y = intensity)) +
  geom_density_ridges() +
  labs(title = "Phosphorus", x = "kg P/t", y = "") +
  theme_ridges() +
  theme(axis.text.y=element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))

P_ridge_system <- ggplot(total_stressor_species_summary_plot %>% 
                             filter(stressor == "P" & !is.na(system)),
                           aes(x = total_value, y = system)) +
  geom_density_ridges() +
  labs(title = "Phosphorus", x = "kg P/t", y = "") +
  theme_ridges() +
  theme(axis.text.y=element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))

P_ridge_taxa <- ggplot(total_stressor_species_summary_plot %>% 
                           filter(stressor == "P" & !is.na(system)),
                         aes(x = total_value, y = taxa)) +
  geom_density_ridges() +
  labs(title = "Phosphorus", x = "kg P/t", y = "") +
  theme_ridges() +
  theme(axis.text.y=element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))

# land
land_ridge_intensity <- ggplot(total_stressor_species_summary_plot %>% 
                                filter(stressor == "land" & !is.na(intensity)),
                              aes(x = total_value, y = intensity)) +
  geom_density_ridges() +
  labs(title = "Land", x = "m2/t", y = "") +
  theme_ridges() +
  theme(axis.text.y=element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))

land_ridge_system <- ggplot(total_stressor_species_summary_plot %>% 
                             filter(stressor == "land" & !is.na(system)),
                           aes(x = total_value, y = system)) +
  geom_density_ridges() +
  labs(title = "Land", x = "m2/t", y = "") +
  theme_ridges() +
  theme(axis.text.y=element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))

land_ridge_taxa <- ggplot(total_stressor_species_summary_plot %>% 
                              filter(stressor == "land" & !is.na(system)),
                            aes(x = total_value, y = taxa)) +
  geom_density_ridges() +
  labs(title = "Land", x = "m2/t", y = "") +
  theme_ridges() +
  theme(axis.text.y=element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))

# water
water_ridge_intensity <- ggplot(total_stressor_species_summary_plot %>% 
                                filter(stressor == "water" & !is.na(intensity)),
                              aes(x = total_value, y = intensity)) +
  geom_density_ridges() +
  labs(title = "Water", x = "m3/t", y = "") +
  theme_ridges() +
  theme(axis.text.y=element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))

water_ridge_system <- ggplot(total_stressor_species_summary_plot %>% 
                             filter(stressor == "water" & !is.na(system)),
                           aes(x = total_value, y = system)) +
  geom_density_ridges() +
  labs(title = "Water", x = "m3/t", y = "") +
  theme_ridges() +
  theme(axis.text.y=element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))

water_ridge_taxa <- ggplot(total_stressor_species_summary_plot %>% 
                               filter(stressor == "water" & !is.na(system)),
                             aes(x = total_value, y = taxa)) +
  geom_density_ridges() +
  labs(title = "Water", x = "m3/t", y = "") +
  theme_ridges() +
  theme(axis.text.y=element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1))

png("intensity_ridge_plots.png", width = 8, height = 4, units = "in", res = 300)
ggarrange(ghg_ridge_intensity, N_ridge_intensity, P_ridge_intensity,
          land_ridge_intensity, water_ridge_intensity, nrow = 1)
dev.off()

png("system_ridge_plots.png", width = 8, height = 4, units = "in", res = 300)
ggarrange(ghg_ridge_system, N_ridge_system, P_ridge_system,
          land_ridge_system, water_ridge_system, nrow = 1)
dev.off()

png("taxa_ridge_plots.png", width = 8, height = 4, units = "in", res = 300)
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

