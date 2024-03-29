# Analysis of levers using alternate feed constants
# Do this for both Mass allocation and Economic allocation versions

#_______________________________________________________________________________________________________________________#
# Load packages and source functions
#_______________________________________________________________________________________________________________________#
rm(list=ls())
library(tidyverse)
library(countrycode)
library(ggplot2)
library(hrbrthemes)
library(viridis)
library(RColorBrewer)
library(ggthemes)
library(ggpubr)
library(cowplot) ## for plot_grid

source("Functions.R")

# Set data directories
datadir <- "/Volumes/jgephart/BFA Environment 2/Data"
outdir <- "/Volumes/jgephart/BFA Environment 2/Outputs"

##########################################################################################################################
# OPTION: Choose allocation method (used for all analyses and plotting below)
#allocation_method <- "Mass"
allocation_method <- "Economic"
##########################################################################################################################

#_______________________________________________________________________________________________________________________#
# Load data
#_______________________________________________________________________________________________________________________#
# Load full lca data with predicted parameters
df <- read.csv(file.path(outdir, "lca-dat-imputed-vars_rep-sqrt-n-farms_edible-weight.csv"))

# Load and join weightings
prod_weightings <- read.csv(file.path(outdir, "aqua_prod_weightings.csv"))

df <- df %>% 
  left_join(prod_weightings, by = c("clean_sci_name", "taxa")) 

# Load and join electricity GHG constants
electricity_gwp <- read.csv(file.path(datadir, "electricity_GWP.csv"))
electricity_gwp$iso3c <- countrycode(electricity_gwp$Country, origin = "country.name", destination = "iso3c")

df <- df %>% 
  left_join(electricity_gwp %>% select(-Country), by = "iso3c")

# Load diesel, petrol, natural gas constants
energy_gwp <- read.csv(file.path(datadir, "energy_carriers_impact_factors.csv"))
energy_gwp <- energy_gwp %>%
  filter(Impact.category == "Global warming potential") %>%
  select(Input, Value) %>%
  mutate(Input = case_when(
    (Input == "Diesel") ~ "Diesel_L",
    (Input == "Petrol") ~ "Petrol_L", 
    (Input == "Natural gas") ~"NaturalGas_L"
  ))

# Load and join evaporative loss constants
evap <- read.csv(file.path(outdir, "clim_summarise_by_country.csv"))
evap$iso3c <- countrycode(evap$admin, origin = "country.name", destination = "iso3c")
evap <- evap %>%
  mutate(evap_rate_m3_per_m2 = mean_evap_mm/1000) %>%
  filter(!is.na(iso3c)) %>%
  select(iso3c, evap_rate_m3_per_m2)

df <- df %>%
  left_join(evap, by = "iso3c")

# Load feed data
feed_fp <- read.csv(file.path(outdir, "weighted_feed_fp.csv"))
# Change names to match df
feed_fp <- feed_fp %>%
  mutate(feed_type = case_when(
    (Input.type == "Soy") ~ "soy",
    (Input.type == "Crop") ~ "crops",
    (Input.type == "Fishery") ~ "fmfo",
    (Input.type == "Livestock") ~ "animal"
  )) %>%
  select(-Input.type) %>%
  mutate(stressor = case_when(
    Impact.category == "Global warming potential" ~ "ghg",
    Impact.category == "Marine eutrophication" ~ "N",
    Impact.category == "Freshwater eutrophication" ~ "P",
    Impact.category == "Land use" ~ "land",
    Impact.category == "Water consumption" ~ "water"
  )) %>%
  filter(Allocation == allocation_method) %>%
  select(feed_type, stressor, ave_stressor)

# Load feed N/P data
feed_NP <- read.csv(file.path(outdir, "feed_NP_clean.csv"))

# Load and join fish N/P data 
fish_NP <- read.csv(file.path(outdir, "fish_NP_clean.csv"))
fish_NP <- fish_NP %>% 
  select(clean_sci_name, N_t_liveweight_t, P_t_liveweight_t)

df <- df %>%
  left_join(fish_NP, by = c("clean_sci_name", "N_t_liveweight_t", "P_t_liveweight_t"))

#_______________________________________________________________________________________________________________________#
# Format data into taxa means and SDs for each model parameter
#_______________________________________________________________________________________________________________________#
df_taxa <- df %>%
  # Only count land for ponds and recirculating systems - add zeros elsewhere
  mutate(Yield_m2_per_t = ifelse(system %in% c("Ponds", "Recirculating and tanks"), Yield_m2_per_t, 0)) %>%
  # Create species means and SDs
  group_by(taxa, clean_sci_name, prod_weighting) %>% # Note to self: do we need to retain intensity or system info?
  summarise(feed_soy_mean = mean(feed_soy, na.rm = TRUE), feed_soy_sd = sd(feed_soy, na.rm = TRUE),
            feed_crops_mean = mean(feed_crops, na.rm = TRUE), feed_crops_sd = sd(feed_crops, na.rm = TRUE),
            feed_fmfo_mean = mean(feed_fmfo, na.rm = TRUE), feed_fmfo_sd = sd(feed_fmfo, na.rm = TRUE),
            feed_animal_mean = mean(feed_animal, na.rm = TRUE), feed_animal_sd = sd(feed_animal, na.rm = TRUE),
            fcr_mean = mean(fcr, na.rm = TRUE), fcr_sd = sd(fcr, na.rm = TRUE),
            electricity_mean = mean(Electricity_kwh, na.rm = TRUE), electricity_sd = sd(Electricity_kwh, na.rm = TRUE),
            diesel_mean = mean(Diesel_L, na.rm = TRUE), diesel_sd = sd(Diesel_L, na.rm = TRUE),
            petrol_mean = mean(Petrol_L, na.rm = TRUE), petrol_sd = sd(Petrol_L, na.rm = TRUE),
            naturalgas_mean = mean(NaturalGas_L, na.rm = TRUE), naturalgas_sd = sd(NaturalGas_L, na.rm = TRUE),
            yield_mean = mean(Yield_m2_per_t, na.rm = TRUE), yield_sd = sd(Yield_m2_per_t, na.rm = TRUE), 
            electricity_ghg_mean = mean(GWP_perkWh_kgCO2eq, na.rm = TRUE), electricity_ghg_sd = sd(GWP_perkWh_kgCO2eq, na.rm = TRUE),
            evap_mean = mean(evap_rate_m3_per_m2, na.rm = TRUE), evap_sd = sd(evap_rate_m3_per_m2, na.rm = TRUE),
            fish_N_mean = mean(N_t_liveweight_t, na.rm = TRUE), fish_N_sd = sd(N_t_liveweight_t, na.rm = TRUE), 
            fish_P_mean = mean(P_t_liveweight_t, na.rm = TRUE), fish_P_sd = sd(P_t_liveweight_t, na.rm = TRUE)
            ) %>% 
    # Create weighted taxa means and SDs
    group_by(taxa) %>%
    summarise(feed_soy_weighted_mean = weighted.mean(feed_soy_mean, prod_weighting, na.rm = TRUE), feed_soy_weighted_sd = sum(prod_weighting * (feed_soy_mean - feed_soy_weighted_mean)^2),
              feed_crops_weighted_mean = weighted.mean(feed_crops_mean, prod_weighting, na.rm = TRUE), feed_crops_weighted_sd = sum(prod_weighting * (feed_crops_mean - feed_crops_weighted_mean)^2),
              feed_fmfo_weighted_mean = weighted.mean(feed_fmfo_mean, prod_weighting, na.rm = TRUE), feed_fmfo_weighted_sd = sum(prod_weighting * (feed_fmfo_mean - feed_fmfo_weighted_mean)^2),
              feed_animal_weighted_mean = weighted.mean(feed_animal_mean, prod_weighting, na.rm = TRUE), feed_animal_weighted_sd = sum(prod_weighting * (feed_animal_mean - feed_animal_weighted_mean)^2),
              fcr_weighted_mean = weighted.mean(fcr_mean, prod_weighting, na.rm = TRUE), fcr_weighted_sd = sum(prod_weighting * (fcr_mean - fcr_weighted_mean)^2),
              electricity_weighted_mean = weighted.mean(electricity_mean, prod_weighting, na.rm = TRUE), electricity_weighted_sd = sum(prod_weighting * (electricity_mean - electricity_weighted_mean)^2),
              diesel_weighted_mean = weighted.mean(diesel_mean, prod_weighting, na.rm = TRUE), diesel_weighted_sd = sum(prod_weighting * (diesel_mean - diesel_weighted_mean)^2),
              petrol_weighted_mean = weighted.mean(petrol_mean, prod_weighting, na.rm = TRUE), petrol_weighted_sd = sum(prod_weighting * (petrol_mean - petrol_weighted_mean)^2),
              naturalgas_weighted_mean = weighted.mean(naturalgas_mean, prod_weighting, na.rm = TRUE), naturalgas_weighted_sd = sum(prod_weighting * (naturalgas_mean - naturalgas_weighted_mean)^2),
              yield_weighted_mean = weighted.mean(yield_mean, prod_weighting, na.rm = TRUE), yield_weighted_sd = sum(prod_weighting * (yield_mean - yield_weighted_mean)^2),
              electricity_ghg_weighted_mean = weighted.mean(electricity_ghg_mean, prod_weighting, na.rm = TRUE), electricity_ghg_weighted_sd = sum(prod_weighting * (electricity_ghg_mean - electricity_ghg_weighted_mean)^2),
              evap_weighted_mean = weighted.mean(evap_mean, prod_weighting, na.rm = TRUE), evap_weighted_sd = sum(prod_weighting * (evap_mean - evap_weighted_mean)^2),
              fish_N_weighted_mean = weighted.mean(fish_N_mean, prod_weighting, na.rm = TRUE), fish_N_weighted_sd = sum(prod_weighting * (fish_N_mean - fish_N_weighted_mean)^2),
              fish_P_weighted_mean = weighted.mean(fish_P_mean, prod_weighting, na.rm = TRUE), fish_P_weighted_sd = sum(prod_weighting * (fish_P_mean - fish_P_weighted_mean)^2)
              
    ) %>% 
  filter(taxa != "fresh_crust") # Freshwater crustaceans will drop out due to lack of data


#_______________________________________________________________________________________________________________________#
# Test sensitivity of each parameter
#_______________________________________________________________________________________________________________________#
# Estimates without any perturbation
stressor_0 <- stressor_sensitivity(data_lca = df_taxa, data_feed_stressors = feed_fp, data_feed_NP = feed_NP, data_energy = energy_gwp,
                             delta_FCR = 0, delta_soy = 0, delta_crops = 0, delta_animal = 0, delta_fmfo = 0,
                             delta_electricity = 0, delta_diesel = 0, delta_petrol = 0, delta_natgas = 0, delta_yield = 0)
# Test function
stressor_1 <- stressor_sensitivity(data_lca = df_taxa, data_feed_stressors = feed_fp, data_feed_NP = feed_NP, data_energy = energy_gwp,
                                   delta_FCR = -0.1*df_taxa$fcr_weighted_mean, delta_soy = 0, delta_crops = 0, delta_animal = 0, delta_fmfo = 0,
                                   delta_electricity = 0, delta_diesel = 0, delta_petrol = 0, delta_natgas = 0, delta_yield = 0)

# Run n perturbation function
perturbation_mean <- compare_perturbations_mean(n = -0.10)

# Reformat to plot as heatmap
plot_perturbation <- perturbation_mean %>% 
  pivot_longer(cols = fcr_ghg_percent_change:yield_water_percent_change, names_sep = "_", names_to = c("Parameter", "Stressor", "drop1", "drop2")) %>%
  select(-contains("drop")) %>%
  mutate(taxa = case_when(taxa == "hypoph_carp" ~ "silver/bighead",
                          taxa == "oth_carp" ~ "misc carp",
                          taxa == "misc_diad" ~ "misc diad",
                          taxa == "misc_fresh" ~ "misc freshwater",
                          taxa == "misc_marine" ~ "misc marine",
                          taxa == "fresh_crust" ~ "freshwater crust",
                          taxa == "plants" ~ "seaweeds",
                          TRUE ~ taxa
  ))

plot_perturbation$Parameter <- factor(plot_perturbation$Parameter, levels = c("fcr", "soy", "crops", "animal", "fmfo", "energy", "yield"))
# Set taxa order to match Figure 1
taxa_order <- c("seaweeds",
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
plot_perturbation$taxa <- factor(plot_perturbation$taxa, levels = taxa_order)

base_size <- 10
base_family <- "sans"
color_lim <- max(plot_perturbation$value, na.rm = TRUE)*c(-1,1)

plot_perturbation$Stressor <- factor(plot_perturbation$Stressor, levels = c("ghg", "land", "water", "N", "P"))
label_names <- c("GHG", "Land", "Water", "N", "P")
names(label_names) <- c("ghg", "land", "water", "N", "P")


low_color <- "#364F6B" # blue
#low_color <- "#799442" # green i.e., For green light vs red light (reduction vs increase in stressors)
high_color <- "#C93F3F" # red
mid_color <- "white"

# For mass allocation insert NA for > 20 and color cell as black
plot_perturbation$value[is.na(plot_perturbation$value)] <- 0
plot_perturbation$value[plot_perturbation$value > 20] <- NA 
fig_4a <- ggplot(plot_perturbation, aes(x = Parameter, y = taxa, fill = value)) +
  geom_tile() +
  labs(x = "", y = "") +
  facet_wrap(~Stressor, nrow = 1, labeller = labeller(Stressor = label_names)) +
  scale_fill_gradient2(low = low_color,
                       mid = mid_color,
                       high = high_color,
                       midpoint = 0,
                       na.value = "black") +
  guides(fill = guide_colourbar(barwidth = 10, barheight = 0.5)) +
  theme(axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"), 
        axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"), 
        axis.text = element_text(size = ceiling(base_size*0.7), colour = "black"),
        axis.title = element_text(size = ceiling(base_size*0.8)), 
        axis.text.x = element_text(angle = 90, vjust = 0.2, hjust = 1),
        panel.grid.minor = element_blank(), 
        panel.grid.major.y = element_line(colour = "gray", linetype = "dotted"), 
        panel.grid.major.x = element_blank(), 
        panel.background = element_blank(), panel.border = element_blank(),
        strip.background = element_rect(linetype = 1, fill = "white"), strip.text = element_text(), 
        strip.text.x = element_text(vjust = 0, hjust = 0), 
        strip.text.y = element_text(angle = -90), 
        legend.text = element_text(size = ceiling(base_size*0.7), family = "sans"), 
        legend.title = element_blank(), 
        legend.key = element_rect(fill = "white", colour = NA), 
        legend.position="bottom",
        legend.margin=margin(t=-0.5, r=0, b=0, l=0, unit="cm"),
        plot.title = element_text(size = ceiling(base_size*1.1), face = "bold"), 
        plot.subtitle = element_text(size = ceiling(base_size*1.05)))

# Output graphing data for SI: 
write.csv(plot_perturbation, file.path(outdir, "data-to-plot-Fig-4a.csv"), row.names = FALSE)

#_______________________________________________________________________________________________________________________#
# Scenario analysis
#_______________________________________________________________________________________________________________________#
# Baseline
stressor_baseline <- stressor_sensitivity(data_lca = df_taxa, data_feed_stressors = feed_fp, data_feed_NP = feed_NP, data_energy = energy_gwp,
                                          delta_FCR = 0, delta_soy = 0, delta_crops = 0, delta_animal = 0, delta_fmfo = 0,
                                          delta_electricity = 0, delta_diesel = 0, delta_petrol = 0, delta_natgas = 0, delta_yield = 0)
stressor_baseline <- stressor_baseline %>% 
  select(taxa, contains("total")) %>%
  pivot_longer(cols = total_ghg:total_water, names_sep = "_", names_to = c("drop", "Stressor")) %>%
  select(-contains("drop"))
stressor_baseline$scenario <- "Baseline"

# Scenario 1: Take lowest 20% point of FCR from data distribution by taxa
fcr_20 <- df %>%
  group_by(taxa) %>%
  summarise(percentile_20 = quantile(fcr, probs = .2)[[1]])

## Calculate how much to subtract off of the mean to get to the 20th percentile
fcr_20_delta <- df_taxa$fcr_weighted_mean - fcr_20$percentile_20

stressor_s1 <- stressor_sensitivity(data_lca = df_taxa, data_feed_stressors = feed_fp, data_feed_NP = feed_NP, data_energy = energy_gwp,
                                   delta_FCR = -1*fcr_20_delta, delta_soy = 0, delta_crops = 0, delta_animal = 0, delta_fmfo = 0,
                                   delta_electricity = 0, delta_diesel = 0, delta_petrol = 0, delta_natgas = 0, delta_yield = 0)
stressor_s1 <- stressor_s1 %>% 
  select(taxa, contains("total")) %>%
  pivot_longer(cols = total_ghg:total_water, names_sep = "_", names_to = c("drop", "Stressor")) %>%
  select(-contains("drop"))
stressor_s1$scenario <- "FCR lower 20th"

# Scenario 2a: Replace FMFO with all soy
stressor_s2a <- stressor_sensitivity(data_lca = df_taxa, data_feed_stressors = feed_fp, data_feed_NP = feed_NP, data_energy = energy_gwp,
                                    delta_FCR = 0, delta_soy = df_taxa$feed_fmfo_weighted_mean, delta_crops = 0, 
                                    delta_animal = 0, delta_fmfo = -1*df_taxa$feed_fmfo_weighted_mean,
                                    delta_electricity = 0, delta_diesel = 0, delta_petrol = 0, delta_natgas = 0, delta_yield = 0)
stressor_s2a <- stressor_s2a %>% 
  select(taxa, contains("total")) %>%
  pivot_longer(cols = total_ghg:total_water, names_sep = "_", names_to = c("drop", "Stressor")) %>%
  select(-contains("drop"))
stressor_s2a$scenario <- "Replace FMFO with soy"

# Scenario 2b: Replace FMFO with land change-free soy (currently only replacement soy is land change-free)
## Land use change free soy
feed_fp_s7_file <- paste("feed_fp_scenario_7_", str_to_lower(allocation_method), ".csv", sep = "")
feed_fp_s7 <- read.csv(file.path(outdir, feed_fp_s7_file))

stressor_s2b <- stressor_sensitivity(data_lca = df_taxa, data_feed_stressors = feed_fp, data_feed_NP = feed_NP, data_energy = energy_gwp,
                                     delta_FCR = 0, delta_soy = 0, delta_crops = 0, delta_animal = 0, delta_fmfo = 0,
                                     delta_electricity = 0, delta_diesel = 0, delta_petrol = 0, delta_natgas = 0, delta_yield = 0,
                                     fmfo_with_land_change_free_crops = feed_fp_s7)
stressor_s2b <- stressor_s2b %>% 
  select(taxa, contains("total")) %>%
  pivot_longer(cols = total_ghg:total_water, names_sep = "_", names_to = c("drop", "Stressor")) %>%
  select(-contains("drop"))
stressor_s2b$scenario <- "Replace FMFO with deforestation-free soy"

# Scenario 2c: Replace FMFO with fishery by-products
feed_fp_s2c_file <- paste("feed_fp_scenario_2c_", str_to_lower(allocation_method), ".csv", sep = "")
feed_fp_s2c <- read.csv(file.path(outdir, feed_fp_s2c_file))

stressor_s2c <- stressor_sensitivity(data_lca = df_taxa, data_feed_stressors = feed_fp, data_feed_NP = feed_NP, data_energy = energy_gwp,
                                     delta_FCR = 0, delta_soy = 0, delta_crops = 0, delta_animal = 0, delta_fmfo = 0,
                                     delta_electricity = 0, delta_diesel = 0, delta_petrol = 0, delta_natgas = 0, 
                                     delta_yield = 0, fishery_byproducts = feed_fp_s2c)
stressor_s2c <- stressor_s2c %>% 
  select(taxa, contains("total")) %>%
  pivot_longer(cols = total_ghg:total_water, names_sep = "_", names_to = c("drop", "Stressor")) %>%
  select(-contains("drop"))
stressor_s2c$scenario <- "Replace FMFO w/ byproducts"

# Scenario 2d: Replace FMFO with low impact fishery by-products
feed_fp_s2d_file <- paste("feed_fp_scenario_2d_", str_to_lower(allocation_method), ".csv", sep = "")
feed_fp_s2d <- read.csv(file.path(outdir, feed_fp_s2d_file))

stressor_s2d <- stressor_sensitivity(data_lca = df_taxa, data_feed_stressors = feed_fp, data_feed_NP = feed_NP, data_energy = energy_gwp,
                                     delta_FCR = 0, delta_soy = 0, delta_crops = 0, 
                                     delta_animal = 0, delta_fmfo = 0,
                                     delta_electricity = 0, delta_diesel = 0, delta_petrol = 0, delta_natgas = 0, delta_yield = 0,
                                    low_impact_fishery_byproducts = feed_fp_s2d)
stressor_s2d <- stressor_s2d %>% 
  select(taxa, contains("total")) %>%
  pivot_longer(cols = total_ghg:total_water, names_sep = "_", names_to = c("drop", "Stressor")) %>%
  select(-contains("drop"))
stressor_s2d$scenario <- "Replace FMFO with low impact fishery by-products"

# Scenario 3: Look at yield-FCR relationship and simultaneously change both in model; Take highest 20% yield and calculate FCR
yield_20 <- df %>%
  mutate(Yield_m2_per_t = ifelse(system %in% c("Ponds", "Recirculating and tanks"), Yield_m2_per_t, 0)) %>%
  group_by(taxa) %>%
  summarise(percentile_20 = quantile(Yield_m2_per_t, probs = .2, na.rm = TRUE)[[1]])

## Calculate how much to subtract off of the mean to get to the 20th percentile
yield_20_delta <- df_taxa$yield_weighted_mean - yield_20$percentile_20
# One ends up negative (likely due to weighting), so replace with zero
yield_20_delta[yield_20_delta<0] <- 0
  
df_taxa$yield_weighted_mean - yield_20_delta

## Note: These FCR's end up heavily influenced by the intercept, so many take on the average FCR
new_fcr <- lm(fcr ~ Yield_m2_per_t, data = df)$coefficients[1] + 
  lm(fcr ~ Yield_m2_per_t, data = df)$coefficients[2] * yield_20$percentile_20

## For now, just take the new yield
stressor_s3 <- stressor_sensitivity(data_lca = df_taxa, data_feed_stressors = feed_fp, data_feed_NP = feed_NP, data_energy = energy_gwp,
                                     delta_FCR = 0, delta_soy = 0, delta_crops = 0, delta_animal = 0, delta_fmfo = 0,
                                     delta_electricity = 0, delta_diesel = 0, delta_petrol = 0, delta_natgas = 0, 
                                    delta_yield = -1*yield_20_delta)
stressor_s3 <- stressor_s3 %>% 
  select(taxa, contains("total")) %>%
  pivot_longer(cols = total_ghg:total_water, names_sep = "_", names_to = c("drop", "Stressor")) %>%
  select(-contains("drop"))
stressor_s3$scenario <- "Yield upper 20th"

# Scenario 4: Capture fisheries - catch 1.2 x as much fish with 60% less effort
df_capture <- read.csv(file.path(outdir, "non-bayes-stressors_capture_observation-level_edible-weight.csv"))

delta_ghg <- 0.56/1.13 # Confirm with Rob

stressor_s4 <- df_capture %>%
  mutate(ghg_kg_t = delta_ghg*ghg_kg_t)
stressor_s4$scenario <- "13% more catch with 56% of the effort"

df_capture$scenario <- "capture_baseline"

# Scenario 5: Capture fisheries - Take minimum gear GHG per species and apply to the group
df_capture_raw <- read.csv(file.path(datadir, "fisheries_fuel_use.csv"))

stressor_s5 <- df_capture_raw %>%
  # Remove mixed gear and nei observations
  filter(!str_detect(pattern = " nei", species)) %>%
  filter(gear != "Other, Mixed, or Unknown") %>%
  # Remove observations with 0 gear, species, or consumption weighting
  filter(gear_weighting > 0 & species_weighting > 0 & consumption_weighting > 0) %>%
  # Represent each species by the min ghg gear type
  group_by(species_group, species) %>%
  # Create species gear-weighted means
  summarise(species_ghg_kg_t = min(ghg), 
            species_weighting = mean(species_weighting), 
            consumption_weighting = mean(consumption_weighting)) %>%
  # Re-weight species and consumption within taxa group
  ungroup() %>%
  group_by(species_group) %>%
  mutate(species_consumption_weighting = (species_weighting*consumption_weighting)/sum(species_weighting*consumption_weighting)) %>%
  summarise(ghg_kg_t = sum(species_ghg_kg_t*species_consumption_weighting))
stressor_s5$scenario <- "Min GHG gear-type"

# Bind capture scenarios
capture_scenarios <- df_capture %>%
  bind_rows(stressor_s4) %>%
  bind_rows(stressor_s5)

# Scenario 6: All by-products sourced from Alaska Pollock
feed_fp_s6_file <- paste("feed_fp_scenario_6_", str_to_lower(allocation_method), ".csv", sep = "")
feed_fp_s6 <- read.csv(file.path(outdir, feed_fp_s6_file))

stressor_s6 <- stressor_sensitivity(data_lca = df_taxa, data_feed_stressors = feed_fp, data_feed_NP = feed_NP, data_energy = energy_gwp,
                                     delta_FCR = 0, delta_soy = 0, delta_crops = 0, delta_animal = 0, delta_fmfo = 0,
                                     delta_electricity = 0, delta_diesel = 0, delta_petrol = 0, delta_natgas = 0, delta_yield = 0,
                                     fmfo_low_impact_fishery_byproducts = feed_fp_s6)
stressor_s6 <- stressor_s6 %>% 
  select(taxa, contains("total")) %>%
  pivot_longer(cols = total_ghg:total_water, names_sep = "_", names_to = c("drop", "Stressor")) %>%
  select(-contains("drop"))
stressor_s6$scenario <- "All by-products sourced from low impact fisheries"

# Scenario 7: All soy and crops from non-rainforest depleting sources
feed_fp_s7_file <- paste("feed_fp_scenario_7_", str_to_lower(allocation_method), ".csv", sep = "")
feed_fp_s7 <- read.csv(file.path(outdir, feed_fp_s7_file))

stressor_s7 <- stressor_sensitivity(data_lca = df_taxa, data_feed_stressors = feed_fp, data_feed_NP = feed_NP, data_energy = energy_gwp,
                                    delta_FCR = 0, delta_soy = 0, delta_crops = 0, delta_animal = 0, delta_fmfo = 0,
                                    delta_electricity = 0, delta_diesel = 0, delta_petrol = 0, delta_natgas = 0, delta_yield = 0,
                                    land_change_free_soy = feed_fp_s7, land_change_free_crops = feed_fp_s7)
stressor_s7 <- stressor_s7 %>% 
  select(taxa, contains("total")) %>%
  pivot_longer(cols = total_ghg:total_water, names_sep = "_", names_to = c("drop", "Stressor")) %>%
  select(-contains("drop"))
stressor_s7$scenario <- "Deforestation-free soy & crops"

# Scenario 8: Use 0 for electricity GHG constant
stressor_s8 <- stressor_sensitivity(data_lca = df_taxa, data_feed_stressors = feed_fp, data_feed_NP = feed_NP, data_energy = energy_gwp,
                                    delta_FCR = 0, delta_soy = 0, delta_crops = 0, delta_animal = 0, delta_fmfo = 0,
                                    delta_electricity = 0, delta_diesel = 0, delta_petrol = 0, delta_natgas = 0, delta_yield = 0,
                                    zero_emissions_electric = TRUE)
stressor_s8 <- stressor_s8 %>% 
  select(taxa, contains("total")) %>%
  pivot_longer(cols = total_ghg:total_water, names_sep = "_", names_to = c("drop", "Stressor")) %>%
  select(-contains("drop"))
stressor_s8$scenario <- "Zero emission electricity"

# Bind all scenarios
scenarios <- stressor_baseline %>%
  bind_rows(stressor_s1) %>%
  bind_rows(stressor_s2a) %>%
  bind_rows(stressor_s2b) %>%
  bind_rows(stressor_s2c) %>%
  bind_rows(stressor_s2d) %>%
  bind_rows(stressor_s3) %>%
  bind_rows(stressor_s6) %>%
  bind_rows(stressor_s7) %>%
  bind_rows(stressor_s8) %>%
  filter(!(taxa %in% c("bivalves", "plants"))) %>%
  mutate(taxa = case_when(taxa == "hypoph_carp" ~ "silver/bighead",
                          taxa == "oth_carp" ~ "misc carp",
                          taxa == "misc_diad" ~ "misc diad",
                          taxa == "misc_fresh" ~ "misc freshwater",
                          taxa == "misc_marine" ~ "misc marine",
                          taxa == "fresh_crust" ~ "freshwater crust",
                          taxa == "plants" ~ "seaweeds",
                          TRUE ~ taxa
  ))

scenarios$Stressor <- factor(scenarios$Stressor, levels = c("ghg", "land", "water", "N", "P"))
scenarios$taxa <- factor(scenarios$taxa, levels = taxa_order)

scenarios_diff <- scenarios %>%
  filter(scenario != "Baseline") %>%
  left_join(scenarios %>% filter(scenario == "Baseline") %>% select(taxa, Stressor, "baseline_value" = "value"), 
            by = c("taxa", "Stressor")) %>%
  mutate(value_change = value - baseline_value, value_percent_change = round(100*(value - baseline_value)/baseline_value, digits = 0))

lollipop_theme <- theme(axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"), 
                        axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"), 
                        axis.text = element_text(size = ceiling(base_size*0.7), colour = "black"),
                        axis.title = element_text(size = ceiling(base_size*0.8)), 
                        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                        panel.grid.minor = element_blank(), 
                        panel.grid.major.y = element_line(colour = "gray", linetype = "dotted"), 
                        panel.grid.major.x = element_blank(), 
                        panel.background = element_blank(), panel.border = element_blank(),
                        strip.background = element_rect(linetype = 1, fill = "white"), 
                        strip.text = element_text(), 
                        strip.text.x = element_text(hjust = 0, vjust = 0), 
                        strip.text.y = element_text(angle = -90, size = base_size*0.7), 
                        strip.placement.y = "outside",
                        legend.text = element_blank(), 
                        legend.title = element_blank(), 
                        legend.key = element_blank(), 
                        legend.position="none",
                        plot.title = element_text(size = ceiling(base_size*1.1), face = "bold"), 
                        plot.subtitle = element_text(size = ceiling(base_size*1.05)))

# Use a 50% cutoff for lollipop plot (use arrowheads if beyond cutoff point)
lollipop_dat <- scenarios_diff %>%
  mutate(plot_shape = ifelse(value_percent_change > -50, "under", "over")) %>%
  mutate(value_percent_change = ifelse(value_percent_change < -50, -50, value_percent_change))

fig_4b_dat <- lollipop_dat %>%
  filter(scenario %in% c("FCR lower 20th", "Replace FMFO w/ byproducts",
                         "Yield upper 20th", "Deforestation-free soy & crops"))

# Output graphing data for SI: 
write.csv(fig_4b_dat %>% select(-plot_shape), file.path(outdir, "data-to-plot-Fig-4b.csv"), row.names = FALSE)

fig_4b <- ggplot(fig_4b_dat, aes(x = value_percent_change, y = taxa, shape = plot_shape)) +
  geom_point(size = 2) +
  scale_shape_manual(values=c(60, 20)) +
  geom_segment(aes(x = 0, xend=value_percent_change, yend = taxa)) +
  geom_vline(xintercept = 0) +
  labs(x = "% Change", y = "") +
  scale_x_continuous(limits = c(-50, 25), labels = c(-50, -25, 0), breaks = c(-50, -25, 0)) +
  facet_grid(rows = vars(scenario), cols = vars(Stressor), switch = "y",
             labeller = labeller(Stressor = label_names, scenario =label_wrap_gen(15))) +
  lollipop_theme

margin_theme <- theme(plot.margin = unit(c(0, 0, 0, -3), "mm"),
                                 panel.spacing = unit(1, "mm"),
                                 axis.title.y = element_text(margin = unit(c(0, 0, 0, -1), "mm")),
                                 axis.text.y = element_text(margin = unit(c(0, 0, 0, -3), "mm")))

fig_4_file <- paste("plot_Figure-4_", str_to_lower(allocation_method), sep = "")
# NOTE: this conforms to Nature figure specs (89 mm for one-column width)
#png(file = file.path(outdir, paste(fig_4_file, ".png", sep = "")), width = 89, height = 189, units = "mm", res = 300)
pdf(file = file.path(outdir, paste(fig_4_file, ".pdf", sep = "")), width = 3.5, height = 7.44) # convert 89 x 189 mm to inches
# Adjust spacing between and around plots
plot_grid(fig_4a + margin_theme,
          fig_4b + margin_theme,
          nrow = 2, align = "v", labels = c("a", "b"), axis = "l", rel_heights = c(0.4, 1))
dev.off()


#_______________________________________________________________________________________________________________________#
# Additional aquaculture intervention figures for SI
#_______________________________________________________________________________________________________________________#
                   
if (allocation_method == "Mass") {
  SI_fig_4b_dat <- lollipop_dat %>%
    filter(scenario %in% c("Replace FMFO with deforestation-free soy",
                           "Replace FMFO with low impact fishery by-products",
                           #"Replace FMFO with soy & crops", 
                           #"All by-products sourced from low impact fisheries",
                           "Zero emission electricity" ))
} else if (allocation_method == "Economic") {
  SI_fig_4b_dat <- lollipop_dat %>%
    filter(scenario %in% c("All by-products sourced from low impact fisheries",
                           "Replace FMFO with deforestation-free soy",
                           "Replace FMFO with low impact fishery by-products",
                           #"Replace FMFO with soy & crops", 
                           "Zero emission electricity" ))
}

 
SI_fig_4b_other_scenarios <- ggplot(SI_fig_4b_dat, aes(x = value_percent_change, y = taxa, shape = plot_shape)) +
  scale_shape_manual(values=c(60, 20)) +
  geom_point() +
  geom_segment(aes(x = 0, xend=value_percent_change, yend = taxa)) +
  geom_vline(xintercept = 0) +
  labs(x = "% Change", y = "") +
  scale_x_continuous(limits = c(-50, 25), labels = c(-50, -25, 0), breaks = c(-50, -25, 0)) +
  facet_grid(rows = vars(scenario), cols = vars(Stressor), switch = "y",
             labeller = labeller(Stressor = label_names, scenario = label_wrap_gen(20))) +
  lollipop_theme

SI_file <- paste("plot_Figure-SI-X_other-farmed-lever-scenarios_", str_to_lower(allocation_method), sep = "")
png(file.path(outdir, paste(SI_file, ".png", sep = "")), width = 89, height = 189*(1/1.4)*(3/4), units = "mm", res = 300) # Match relative height of main figure
#pdf(file = file.path(outdir, paste(SI_file, ".pdf", sep = "")), width = 3.5, height = 7.44*(1/1.4)*(3/4))
SI_fig_4b_other_scenarios + margin_theme
dev.off()

#_______________________________________________________________________________________________________________________#
# Capture fishery interventions figures for SI
# Note: allocation_method option does not apply (only Mass allocation available for capture)
#_______________________________________________________________________________________________________________________#

# Capture fishery intervention figures
capture_scenarios_diff <- capture_scenarios %>%
  filter(scenario != "capture_baseline") %>%
  rename("value" = "ghg_kg_t") %>%
  left_join(capture_scenarios %>% filter(scenario == "capture_baseline") %>% select(species_group, "baseline_value" = "ghg_kg_t"), 
            by = c("species_group")) %>%
  mutate(value_change = value - baseline_value, value_percent_change = 100*(value - baseline_value)/baseline_value) %>%
  # Format names to match Figure 1
  mutate(species_group = tolower(species_group)) %>%
  mutate(plot_taxa = case_when(species_group == "bivalves" ~ "bivalves",
                               species_group == "cephalopods" ~ "squid, etc",
                               species_group == "flatfishes" ~ "flounder, etc",
                               species_group == "gadiformes" ~ "cod, etc",
                               species_group == "jacks, mullets, sauries" ~ "jack, etc",
                               species_group == "large pelagic fishes" ~ "tuna, etc",
                               species_group == "lobsters" ~ "lobster",
                               species_group == "redfishes, basses, congers" ~ "redfish, etc",
                               species_group == "salmonids" ~ "salmon, etc",
                               species_group == "shrimps" ~ "shrimp",
                               species_group == "small pelagic fishes" ~ "herring, etc")) 

wild_taxa_order <- c("herring, etc",
                     "cod, etc",
                     "salmon, etc",
                     "tuna, etc",
                     "squid, etc",
                     "jack, etc",
                     "redfish, etc",
                     "bivalves",
                     "shrimp",
                     "lobster",
                     "flounder, etc")
capture_scenarios_diff$plot_taxa <- factor(capture_scenarios_diff$plot_taxa, levels = wild_taxa_order)

# DON'T ADD 50% CUTOFF FOR WILD (results in all taxa having more than 50% for scenario: 13% catch with 56% effort )

# SI Fig for the capture scenarios
SI_fig_4b_capture_scenarios <- ggplot(capture_scenarios_diff, aes(x = value_percent_change, y = plot_taxa)) +
  geom_point() + 
  geom_segment(aes(x = 0, xend=value_percent_change, yend = plot_taxa)) +
  geom_vline(xintercept = 0) +
  labs(x = "% Change", y = "") +
  facet_grid(rows = vars(scenario), switch = "y", labeller = label_wrap_gen(20)) +
  lollipop_theme

png(file.path(outdir, "plot_Figure-SI-X_capture-lever-scenarios.png"), width = 89, height = 189*(1/1.4)*(1/2), units = "mm", res = 300)
SI_fig_4b_capture_scenarios # Margin theme not needed
dev.off()

#_______________________________________________________________________________________________________________________#
# Identify common features of low stressor systems for 
#_______________________________________________________________________________________________________________________#
df_species <- read.csv(file.path(outdir, "non-bayes-stressors_farmed_observation-level_edible-weight.csv"))
df_species <- df_species %>%
  filter(taxa %in% c("salmon", "hypoph_carp", "oth_carp", "catfish", "tilapia", "shrimp")) %>%
  mutate(total_ghg = feed_GHG + onfarm_ghg, total_N = feed_N + onfarm_N, total_P = feed_P + onfarm_P,
         total_water = feed_water + onfarm_water, total_land = feed_land + onfarm_land) 

df_species_lower20 <- df_species %>%
  group_by(taxa) %>%
  summarise(ghg_lower_20 = quantile(total_ghg, probs = .2, na.rm = TRUE)[[1]], 
            N_lower_20 = quantile(total_N, probs = .2, na.rm = TRUE)[[1]], 
            P_lower_20 = quantile(total_P, probs = .2, na.rm = TRUE)[[1]], 
            water_lower_20 = quantile(total_water, probs = .2, na.rm = TRUE)[[1]], 
            land_lower_20 = quantile(total_land, probs = .2, na.rm = TRUE)[[1]])

df_species <- df_species %>%
  left_join(df_species_lower20, by = "taxa") %>%
  mutate(ghg_lower_20_logic = ifelse(total_ghg <= ghg_lower_20, 1, 0),
         N_lower_20_logic = ifelse(total_N <= N_lower_20, 1, 0),
         P_lower_20_logic = ifelse(total_P <= P_lower_20, 1, 0),
         water_lower_20_logic = ifelse(total_water <= water_lower_20, 1, 0),
         land_lower_20_logic = ifelse(total_land <= land_lower_20, 1, 0)) %>%
  mutate(lower_total = ghg_lower_20_logic + N_lower_20_logic + P_lower_20_logic +
                           water_lower_20_logic + land_lower_20_logic)




