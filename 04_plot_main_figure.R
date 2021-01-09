# PLOT MAIN FIGURE
# Load model outputs individually, save to plot_grid, and create multi-panel plot
# Remove y-axis labels for panels B, C, and F
# Set plot margins for all plots to 0
# Set axis label margins for panels C and F to 0

# ATTEMPTS TO REMOVE AXIS LABEL MARGINS:
#in tx_plot_theme: theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

#in creation of p object
#theme(axis.text.y=element_blank()) +
#scale_y_discrete(labels = rep("", times = length(tx_index_key$full_taxa_name))) +
# theme(axis.text.y=element_blank(),
#       panel.grid = element_blank(),
#       panel.border = element_blank()) +

rm(list=ls())
library(tidyverse)
library(tidybayes)
library(ggplot2)
library(cowplot) ## for plot_grid
datadir <- "/Volumes/jgephart/BFA Environment 2/Data"
outdir <- "/Volumes/jgephart/BFA Environment 2/Outputs"


######################################################################################################
# Carbon
load(file.path(outdir, "Carbon/2021-01-07_full-model-posterior_Global warming potential_Mass-allocation.RData"))
#units_for_plot = "kg CO2-eq per tonne"
units_for_plot = bquote('kg'~CO[2]*'-eq per tonne')
interval_palette <- c("#9EA8B7", "#6A7A90", "#364F6B") # Order: light to dark

# Theme
tx_plot_theme <- theme(axis.title = element_text(size = 9),
                       axis.text.x = element_text(size = 9, color = "black"),
                       axis.text.y = element_text(size = 9, color = "black"),
                       legend.position = "none",
                       plot.margin = unit(c(0, 3, 0, 0), "mm")) # increase the right margin so that there can be some pillover from x-axis

# Key for naming taxa levels

# FIX IT - later remove the section for creating lca_model_dat; after re-doing model with priors, lca_model_dat will be part of the saved .RData loaded back into this working space)
# Load Data
library(countrycode)
datadir <- "/Volumes/jgephart/BFA Environment 2/Data"
lca_full_dat <- read.csv(file.path(datadir, "2021-01-06_lca-dat-imputed-vars_rep-sqrt-n-farms.csv"), fileEncoding="UTF-8-BOM")
# Get farm-associated carbon footprints data
# Add iso3c
electric_fp_dat <- read.csv(file.path(datadir, "electricity_GWP.csv")) %>%
  mutate(iso3c = countrycode(Country, origin = "country.name", destination = "iso3c"))
other_energy_fp_dat <- read.csv(file.path(datadir, "energy_carriers_impact_factors.csv"))
# units for diesel and petrol constants are in kg CO2-eq / L
diesel_fp <- other_energy_fp_dat %>% filter(Impact.category == "Global warming potential" & Input == "Diesel") %>% pull(Value)
petrol_fp <- other_energy_fp_dat %>% filter(Impact.category == "Global warming potential" & Input == "Petrol") %>% pull(Value)
# NOTE: natural gas is in kg CO2-eq / m3, need to convert to kg CO2-eq / L - multiply by 0.001 m3 / L
natgas_fp <- other_energy_fp_dat %>% filter(Impact.category == "Global warming potential" & Input == "Natural gas") %>% pull(Value) * 0.001
# Format data for model:
lca_model_dat <- lca_full_dat %>%
  select(study_id, iso3c, clean_sci_name, taxa, intensity, system, 
         feed_soy, feed_crops, feed_fmfo, feed_animal, 
         fcr, 
         electric = Electricity_kwh,
         diesel = Diesel_L,
         petrol = Petrol_L,
         natgas = NaturalGas_L) %>%
  drop_na() %>% # just in case
  # OPTION 1: TEST MODEL ON SPECIES THAT ARE FED 
  #filter(fcr != 0)  %>% 
  # OPTION 2: INCLUDE FED AND NON-FED SPECIES BUT mutate feed proportions to be the average within it's clean_sci_name; otherwise, give it an arbitrary simplex (0.25 per component) to avoid STAN error for simplexes that don't sum to 1
  # FIX IT - in terms of on-farm footprint, this is OK because feed proportions are multiplied by FCR == 0 so on-farm footprint for these studies will be 0
  # BUT this will affect the pooled sci and taxa level feed proportions since they enter as 0.25
  group_by(clean_sci_name) %>%
  mutate(feed_soy = if_else(fcr==0, true = mean(feed_soy), false = feed_soy),
         feed_crops = if_else(fcr==0, true = mean(feed_crops), false = feed_crops),
         feed_fmfo = if_else(fcr==0, true = mean(feed_fmfo), false = feed_fmfo),
         feed_animal = if_else(fcr==0, true = mean(feed_animal), false = feed_animal)) %>%
  ungroup() %>%
  # Need to re-normalize everything (using the mean value does not guarantee that proportion sums to 1)
  mutate(feed_sum = feed_soy + feed_crops + feed_fmfo + feed_animal) %>%
  mutate(feed_soy = if_else(feed_sum != 0, true = feed_soy / feed_sum, false = feed_soy), # need conditional feed_sum != 0 - otherwise species with FCR == 0 and no sci-level average will renormalize using denominator of 0
         feed_crops = if_else(feed_sum != 0, true = feed_crops / feed_sum, false = feed_crops),
         feed_fmfo = if_else(feed_sum != 0, true = feed_fmfo / feed_sum, false = feed_fmfo),
         feed_animal = if_else(feed_sum != 0, true = feed_animal / feed_sum, false = feed_animal)) %>%
  # Give the remaining things 0.25
  mutate(feed_soy = if_else(fcr == 0 & feed_soy == 0 & feed_crops == 0 & feed_fmfo == 0 & feed_animal == 0, true = 0.25, false = feed_soy),
         feed_crops = if_else(fcr == 0 & feed_crops == 0 & feed_fmfo == 0 & feed_animal == 0, true = 0.25, false = feed_crops),
         feed_fmfo = if_else(fcr == 0 & feed_fmfo == 0 & feed_animal == 0, true = 0.25, false = feed_fmfo),
         feed_animal = if_else(fcr == 0 & feed_animal == 0, true = 0.25, false = feed_animal)) %>%
  # Multiply energy input by their Carbon (i.e., GHG) footprint, then sum across energy inputs
  # (Electricity use * country-specific GHG of electricity) + (Diesel * GHG of diesel) + (Petrol * GHG of petrol) + (Natural gas * GHG of natural gas)
  # Calculate electriciy GHG footprint
  left_join(electric_fp_dat %>% select(-Country), by = "iso3c") %>% 
  mutate(GWP_perkWh_kgCO2eq = if_else(is.na(GWP_perkWh_kgCO2eq), true = mean(GWP_perkWh_kgCO2eq, na.rm = TRUE), false = GWP_perkWh_kgCO2eq)) %>% # for studies with no country data, just use the average across countries
  mutate(electric_ghg = electric * GWP_perkWh_kgCO2eq) %>%
  # Calculate diesel GHG footprint
  mutate(diesel_ghg = diesel * diesel_fp) %>%
  # Calculate petrol GHG footprint
  mutate(petrol_ghg = petrol * petrol_fp) %>%
  # Calculate natural gas GHG footprint
  mutate(natgas_ghg = natgas * natgas_fp) %>%
  # Calculate sum total of GHG footprint
  mutate(total_ghg = electric_ghg + diesel_ghg + petrol_ghg + natgas_ghg) %>%
  # LAST FORMATING STEP - always arrange by clean_sci_name
  arrange(clean_sci_name) %>%
  mutate(clean_sci_name = as.factor(clean_sci_name),
         sci = as.numeric(clean_sci_name),
         taxa = as.factor(taxa),
         tx = as.numeric(taxa)) 


# Get full taxa group names back
tx_index_key <- lca_model_dat %>%
  group_by(clean_sci_name) %>%
  mutate(n_obs = n()) %>%
  ungroup() %>%
  select(taxa, tx) %>%
  unique() %>%
  arrange(taxa) %>%
  mutate(taxa = as.character(taxa),
         full_taxa_name = case_when(taxa == "hypoph_carp" ~ "big/silverhead carp",
                                    taxa == "misc_marine" ~ "misc marine fishes",
                                    taxa == "misc_fresh" ~ "misc freshwater fishes",
                                    taxa == "misc_diad" ~ "misc diad fishes",
                                    taxa == "oth_carp" ~ "misc carps",
                                    taxa == "fresh_crust" ~ "freshwater crustaceans",
                                    TRUE ~ taxa),
         taxa = as.factor(taxa),
         full_taxa_name = as.factor(full_taxa_name))


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
                           "p_carbon"))])

######################################################################################################
# LAND

# Load model outputs individually, save to plot_grid, and create multi-panel plot
load(file.path(outdir, "Land/2021-01-07_full-model-posterior_Land use_Mass-allocation.RData"))
units_for_plot = bquote('m'^2*'a per tonne')
interval_palette <- c("#D9EAB2", "#C2DD86", "#A9D158") # Order: light to dark

# Use for land:
p_land <- fit_no_na %>%
  spread_draws(tx_land_total_w[tx]) %>%
  median_qi(.width = c(0.95, 0.8, 0.5)) %>%
  left_join(tx_index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  # Set plant and bivalves distributions to 0 (both are 0 for on and off farm)
  mutate(tx_land_total_w = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = tx_land_total_w),
         .lower = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = .lower),
         .upper = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = .upper)) %>%
  # REORDER taxa axis
  mutate(full_taxa_name = fct_relevel(full_taxa_name, full_taxa_name_order)) %>%
  #mutate(full_taxa_name = fct_reorder(full_taxa_name, tx_feed_fp_w)) %>%
  ggplot(., aes(y = full_taxa_name, x = tx_land_total_w, xmin = .lower, xmax = .upper)) +
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
                           "p_carbon",
                           "p_land"))])

######################################################################################################
# NITROGEN 
load(file.path(outdir, "Nitrogen/2021-01-07_full-model-posterior_Marine eutrophication_Mass-allocation.RData"))
units_for_plot = "kg N-eq per tonne"
interval_palette <- c("#FFEDAD", "#FFE37D", "#FFD947") # Order: light to dark

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
                           "p_carbon",
                           "p_land",
                           "p_nitrogen"))])

######################################################################################################
# PHOSPHORUS
load(file.path(outdir, "Phosphorus/2021-01-07_full-model-posterior_Freshwater eutrophication_Mass-allocation.RData"))
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
                           "p_carbon",
                           "p_land",
                           "p_nitrogen",
                           "p_phosphorus"))])

######################################################################################################
# WATER
load(file.path(outdir, "Water/2021-01-08_full-model-posterior_Water consumption_Mass-allocation.RData"))
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
                           "p_carbon",
                           "p_land",
                           "p_nitrogen",
                           "p_phosphorus", 
                           "p_water"))])

######################################################################################################
# WILD CATCH - (plot LAST because it has a different tx_index_key)

# Load model outputs individually, save to plot_grid, and create multi-panel plot
load(file.path(outdir, "Wild/2021-01-07_full-model-posterior_Wild-capture-ghg.RData"))
units_for_plot = bquote('kg'~CO[2]*'-eq per tonne')
interval_palette <- c("#9EA8B7", "#6A7A90", "#364F6B") # Order: light to dark


# FIX IT - can remove this after re-doing model with priors (model data will be part of the saved .RData loaded back into this working space)
datadir <- "/Volumes/jgephart/BFA Environment 2/Data"
wild_dat <- read.csv(file.path(datadir, "fisheries_fuel_use.csv")) %>% as_tibble()
wild_dat_new_weights <- wild_dat %>%
  filter(species_group != "Finfish") %>%
  # Remove mixed gear and nei observations
  filter(!str_detect(pattern = " nei", species)) %>%
  filter(gear != "Other, Mixed, or Unknown") %>%
  # Remove observations with 0 gear, species, or consumption weighting
  filter(gear_weighting > 0 & species_weighting > 0 & consumption_weighting > 0) %>%
  # Re-weight gear within species
  group_by(species) %>%
  mutate(gear_weights_new = gear_weighting/sum(gear_weighting)) %>%
  ungroup() %>%
  # Re-weight species and consumption within species_groups
  group_by(species_group) %>%
  mutate(species_weights_new = species_weighting/sum(species_weighting),
         consumption_weights_new = consumption_weighting/sum(consumption_weighting)) %>%
  ungroup() %>%
  # Calculate overall weights and re-weight
  mutate(prod_of_weights = gear_weights_new * species_weights_new * consumption_weights_new) %>%
  group_by(species_group) %>%
  mutate(overall_weights = prod_of_weights/sum(prod_of_weights)) %>%
  ungroup() %>%
  # LAST FORMATING STEP - always arrange by clean_sci_name
  arrange(species) %>%
  mutate(clean_sci_name = as.factor(species),
         sci = as.numeric(clean_sci_name),
         taxa = as.factor(species_group),
         tx = as.numeric(taxa)) 

wild_index_key <- wild_dat_new_weights %>%
  group_by(clean_sci_name) %>%
  mutate(n_obs = n()) %>%
  ungroup() %>%
  select(taxa, tx) %>%
  unique() %>%
  arrange(taxa) 

# Wild:
p_wild <- fit_no_na %>%
  spread_draws(tx_ghg_w[tx]) %>%
  median_qi(tx_ghg_w, .width = c(0.95, 0.8, 0.5)) %>%
  left_join(wild_index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  # Reorder from highest to lowest
  mutate(taxa = tolower(taxa)) %>%
  mutate(taxa = fct_reorder(taxa, tx_ghg_w)) %>%
  ggplot(aes(y = taxa, x = tx_ghg_w)) +
  geom_interval(aes(xmin = .lower, xmax = .upper), interval_size = 3) +
  geom_point(aes(y = taxa, x = tx_ghg_w)) +
  geom_hline(yintercept = 1:11, linetype = "dotted") +
  scale_color_manual(values = interval_palette) +
  theme_classic() + 
  coord_cartesian(xlim = c(0, 12500)) +
  tx_plot_theme + 
  labs(x = units_for_plot, y = "", title = "")

rm(list=ls()[!(ls() %in% c("outdir", "tx_index_key", "full_taxa_name_order", "tx_plot_theme",
                           "p_carbon",
                           "p_land",
                           "p_nitrogen",
                           "p_phosphorus", 
                           "p_water",
                           "p_wild"))])


# COMBINE INTO PANEL:
p_left <- plot_grid(p_carbon, p_wild, ncol = 1, nrow = 2, align = "hv")
p_right <- plot_grid(p_land, p_nitrogen, p_water, p_phosphorus, ncol = 2, nrow = 2, align = "h")
plot_grid(p_left, p_right, ncol = 2, nrow = 1, rel_widths = c(0.7, 1))
ggsave(filename = file.path(outdir, "plot_Figure-X.png"), width = 183, height = 90, units = "mm")
ggsave(filename = file.path(outdir, "plot_Figure-X.tiff"), device = "tiff", width = 183, height = 90, units = "mm")


