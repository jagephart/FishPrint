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


# Total stressor ridge plots 
total_stressor_species_summary_plot <- stressor_species_summary_plot %>%
  group_by(study_id, taxa,  intensity, system,  clean_sci_name, stressor) %>%
  summarise(total_value = sum(value))

# Set factor levels
stressor_species_summary_plot$taxa <- factor(stressor_species_summary_plot$taxa, levels = rev(c("misc_diad", "misc_marine", "shrimp", "milkfish",
                                                                                                "tilapia", "catfish", "oth_carp", "trout", "salmon",
                                                                                                "hypoph_carp", "plants", "bivalves")))
stressor_species_summary_plot$intensity <- factor(stressor_species_summary_plot$intensity, 
                                                  levels = c("Extensive", "Imp. extensive", 
                                                             "Semi-intensive", "Intensive", "NA"))
stressor_species_summary_plot$system <- factor(stressor_species_summary_plot$system,
                                               levels = c("On- and off-bottom", "Cages & pens", 
                                                          "Ponds", "Recirculating and tanks", "NA"))
# Set plot margins
plot_margins <- c(0.1,0,0,0) # top, right, bottom, left

# ghg
ghg_ridge_intensity <- ggplot(total_stressor_species_summary_plot %>% 
                                filter(stressor == "ghg" & !is.na(intensity)),
       aes(x = total_value, y = intensity)) +
  geom_density_ridges(alpha = 0.7, fill = "#364F6B") +
  labs(title = "GHG", x = expression("kg CO"[2]~"t"^"-1"), y = "") +
  theme_ridges() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        plot.margin=unit(plot_margins, "cm"))

ghg_ridge_system <- ggplot(total_stressor_species_summary_plot %>% 
                                filter(stressor == "ghg" & !is.na(system)),
                              aes(x = total_value, y = system)) +
  geom_density_ridges(alpha = 0.7, fill = "#364F6B") +
  labs(title = "GHG", x = expression("kg CO"[2]~"t"^"-1"), y = "") +
  theme_ridges() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        plot.margin=unit(plot_margins, "cm"))

ghg_ridge_taxa <- ggplot(total_stressor_species_summary_plot %>% 
                             filter(stressor == "ghg" & !is.na(system)),
                           aes(x = total_value, y = taxa)) +
  geom_density_ridges(alpha = 0.7, fill = "#364F6B") +
  scale_y_discrete(labels = rev(c("misc diad", "misc marine", "shrimp", "milkfish",
                                  "tilapia", "catfish", "misc carp", "trout", "salmon",
                                  "silver/bighead", "seaweeds", "bivalves"))) +
  labs(title = "GHG", x = expression("kg CO"[2]~"t"^"-1"), y = "") +
  theme_ridges() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        plot.margin=unit(plot_margins, "cm"))

# N
N_ridge_intensity <- ggplot(total_stressor_species_summary_plot %>% 
                                filter(stressor == "N" & !is.na(intensity)),
                              aes(x = total_value, y = intensity)) +
  geom_density_ridges(alpha = 0.7, fill = "#FFA647") +
  labs(title = "Nitrogen", x = expression("kg N-eq t"^"-1"), y = "") +
  theme_ridges() +
  theme(axis.text.y=element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1),
        plot.margin=unit(plot_margins, "cm"))

N_ridge_system <- ggplot(total_stressor_species_summary_plot %>% 
                             filter(stressor == "N" & !is.na(system)),
                           aes(x = total_value, y = system)) +
  geom_density_ridges(alpha = 0.7, fill = "#FFA647") +
  labs(title = "Nitrogen", x = expression("kg N-eq t"^"-1"), y = "") +
  theme_ridges() +
  theme(axis.text.y=element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1),
        plot.margin=unit(plot_margins, "cm"))

N_ridge_taxa <- ggplot(total_stressor_species_summary_plot %>% 
                           filter(stressor == "N" & !is.na(system)),
                         aes(x = total_value, y = taxa)) +
  geom_density_ridges(alpha = 0.7, fill = "#FFA647") +
  labs(title = "Nitrogen", x = expression("kg N-eq t"^"-1"), y = "") +
  theme_ridges() +
  theme(axis.text.y=element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1),
        plot.margin=unit(plot_margins, "cm"))

# P
P_ridge_intensity <- ggplot(total_stressor_species_summary_plot %>% 
                                filter(stressor == "P" & !is.na(intensity)),
                              aes(x = total_value, y = intensity)) +
  geom_density_ridges(alpha = 0.7, fill = "#FC5185") +
  labs(title = "Phosphorus", x = expression("kg P-eq t"^"-1"), y = "") +
  theme_ridges() +
  theme(axis.text.y=element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1),
        plot.margin=unit(plot_margins, "cm"))

P_ridge_system <- ggplot(total_stressor_species_summary_plot %>% 
                             filter(stressor == "P" & !is.na(system)),
                           aes(x = total_value, y = system)) +
  geom_density_ridges(alpha = 0.7, fill = "#FC5185") +
  labs(title = "Phosphorus", x = expression("kg P-eq t"^"-1"), y = "") +
  theme_ridges() +
  theme(axis.text.y=element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1),
        plot.margin=unit(plot_margins, "cm"))

P_ridge_taxa <- ggplot(total_stressor_species_summary_plot %>% 
                           filter(stressor == "P" & !is.na(system)),
                         aes(x = total_value, y = taxa)) +
  geom_density_ridges(alpha = 0.7, fill = "#FC5185") +
  labs(title = "Phosphorus", x = expression("kg P-eq t"^"-1"), y = "") +
  theme_ridges() +
  theme(axis.text.y=element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1),
        plot.margin=unit(plot_margins, "cm"))

# land
land_ridge_intensity <- ggplot(total_stressor_species_summary_plot %>% 
                                filter(stressor == "land" & !is.na(intensity)),
                              aes(x = total_value, y = intensity)) +
  geom_density_ridges(alpha = 0.7, fill = "#57D182") +
  labs(title = "Land", x = expression("m"^2~"t"^"-1"), y = "") +
  theme_ridges() +
  theme(axis.text.y=element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1),
        plot.margin=unit(plot_margins, "cm"))

land_ridge_system <- ggplot(total_stressor_species_summary_plot %>% 
                             filter(stressor == "land" & !is.na(system)),
                           aes(x = total_value, y = system)) +
  geom_density_ridges(alpha = 0.7, fill = "#57D182") +
  labs(title = "Land", x = expression("m"^2~"t"^"-1"), y = "") +
  theme_ridges() +
  theme(axis.text.y=element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1),
        plot.margin=unit(plot_margins, "cm"))

land_ridge_taxa <- ggplot(total_stressor_species_summary_plot %>% 
                              filter(stressor == "land" & !is.na(system)),
                            aes(x = total_value, y = taxa)) +
  geom_density_ridges(alpha = 0.7, fill = "#57D182") +
  labs(title = "Land", x = expression("m"^2~"t"^"-1"), y = "") +
  theme_ridges() +
  theme(axis.text.y=element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1),
        plot.margin=unit(plot_margins, "cm"))

# water
water_ridge_intensity <- ggplot(total_stressor_species_summary_plot %>% 
                                filter(stressor == "water" & !is.na(intensity)),
                              aes(x = total_value, y = intensity)) +
  geom_density_ridges(alpha = 0.7, fill ="#3FC1C9") +
  labs(title = "Water", x = expression("m"^"3"~"t"^"-1"), y = "") +
  theme_ridges() +
  theme(axis.text.y=element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1),
        plot.margin=unit(plot_margins + c(0, 1, 0, 0), "cm"))

water_ridge_system <- ggplot(total_stressor_species_summary_plot %>% 
                             filter(stressor == "water" & !is.na(system)),
                           aes(x = total_value, y = system)) +
  geom_density_ridges(alpha = 0.7, fill ="#3FC1C9") +
  labs(title = "Water", x = expression("m"^"3"~"t"^"-1"), y = "") +
  theme_ridges() +
  theme(axis.text.y=element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1),
        plot.margin=unit(plot_margins + c(0, 1, 0, 0), "cm"))

water_ridge_taxa <- ggplot(total_stressor_species_summary_plot %>% 
                               filter(stressor == "water" & !is.na(system)),
                             aes(x = total_value, y = taxa)) +
  geom_density_ridges(alpha = 0.7, fill ="#3FC1C9") +
  labs(title = "Water", x = expression("m"^"3"~"t"^"-1"), y = "") +
  theme_ridges() +
  theme(axis.text.y=element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1),
        plot.margin=unit(plot_margins + c(0, 1, 0, 0), "cm"))

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

