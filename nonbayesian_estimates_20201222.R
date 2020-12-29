# Non-Bayesian estimates 2

#_______________________________________________________________________________________________________________________#
# Load packages and source functions
#_______________________________________________________________________________________________________________________#
library(tidyverse)
library(countrycode)

source("Functions.R")

# Set data directories
datadir <- "/Volumes/jgephart/BFA Environment 2/Data"
outdir <- "/Volumes/jgephart/BFA Environment 2/Outputs"

#_______________________________________________________________________________________________________________________#
# Load lca data
#_______________________________________________________________________________________________________________________#
# Load full data with predicted parameters
#df <- read.csv(file.path(datadir, "lca-dat-imputed-vars_rep-n-farms.csv"))
df <- read.csv(file.path(datadir, "lca-dat-imputed-vars_rep-sqrt-n-farms.csv"))

prod_weightings <- read.csv(file.path(datadir, "aqua_prod_weightings.csv"))

df <- df %>% 
  left_join(prod_weightings, by = c("clean_sci_name", "taxa")) 

#_______________________________________________________________________________________________________________________#
# Calculate feed-associated stressors 
#_______________________________________________________________________________________________________________________#
feed_fp <- read.csv(file.path(datadir, "20201217_weighted_feed_fp.csv"))
# Change names to match df
feed_fp <- feed_fp %>%
  mutate(feed_type = case_when(
    (Input.type == "Soy") ~ "soy",
    (Input.type == "Crop") ~ "crops",
    (Input.type == "Fishery") ~ "fmfo",
    (Input.type == "Livestock") ~ "animal"
  )) %>%
  select(-Input.type)

df_feed <- df %>% 
  pivot_longer(feed_soy:feed_animal, names_to = c("drop", "feed_type"), names_sep = "_", values_to = "feed_proportion") %>%
  select(-drop) %>%
  left_join(feed_fp, by = "feed_type") %>%
  mutate(stressor = fcr*feed_proportion*ave_stressor) %>% 
  group_by(study_id, clean_sci_name, taxa, intensity, system, Impact.category, Allocation, Units, prod_weighting) %>%
  summarise(stressor = sum(stressor, na.rm = TRUE))

# Check that all taxa weightings sum to 1
df_feed %>% 
  ungroup() %>%
  select(taxa, clean_sci_name, prod_weighting) %>% 
  distinct() %>%
  group_by(taxa) %>%
  summarise(total_weighting = sum(prod_weighting))

df_feed_taxa <- df_feed %>%
  ungroup() %>%
  group_by(taxa, clean_sci_name, Impact.category, Allocation, Units, prod_weighting) %>%
  summarise(stressor = mean(stressor, na.rm = TRUE)) %>%
  # Convert units to per tonne
  mutate(stressor = 1000*stressor) %>% 
  # Calculate taxa mean
  group_by(taxa, Impact.category, Allocation, Units) %>%
  summarise(weighted_stressor = sum(stressor*prod_weighting), 
            unweighted_stressor = mean(stressor, na.rm = TRUE))

#_______________________________________________________________________________________________________________________#
# Calculate on-farm ghg
#_______________________________________________________________________________________________________________________#
# Load on-farm energy data
electricity_gwp <- read.csv(file.path(datadir, "electricity_GWP.csv"))
electricity_gwp$iso3c <- countrycode(electricity_gwp$Country, origin = "country.name", destination = "iso3c")
energy_gwp <- read.csv(file.path(datadir, "energy_carriers_impact_factors.csv"))
energy_gwp <- energy_gwp %>%
  filter(Impact.category == "Global warming potential") %>%
  select(Input, Value) %>%
  mutate(Input = case_when(
    (Input == "Diesel") ~ "Diesel_L",
    (Input == "Petrol") ~ "Petrol_L", 
    (Input == "Natural gas") ~"NaturalGas_L"
  ))

df_onfarm_ghg <- df %>%
  # Multiply diesel, petrol and natgas by associated GHG emissions
  pivot_longer(c("Diesel_L", "Petrol_L", "NaturalGas_L"), names_to = "Input", values_to = "fuel_L") %>%
  left_join(energy_gwp, by = "Input") %>%
  mutate(fuel_ghg_kgCO2 = fuel_L*Value) %>%
  select(-c("fuel_L", "Value")) %>%
  pivot_wider(names_from = Input, values_from = fuel_ghg_kgCO2) %>%
  # Multiply electricity by country electricity GHG
  left_join(electricity_gwp, by = "iso3c") %>%
  mutate(electricity_ghg_kgCO2 = Electricity_kwh*GWP_perkWh_kgCO2eq) %>% 
  mutate(onfarm_ghg_kgCO2 = Diesel_L + Petrol_L + NaturalGas_L + electricity_ghg_kgCO2) %>% 
  select("study_id", "clean_sci_name", "taxa", "intensity", "system", 
         "Country" = "Country.x", "iso3c", "onfarm_ghg_kgCO2", "prod_weighting")

df_onfarm_ghg_taxa <- df_onfarm_ghg %>% 
  # Calculate species means
  group_by(taxa, clean_sci_name, prod_weighting) %>% 
  summarise(mean_onfarm_ghg = mean(onfarm_ghg_kgCO2, na.rm = TRUE)) %>%
  # Calculate taxa means
  group_by(taxa) %>%
  summarise(onfarmGHG_weighted_stressor = sum(mean_onfarm_ghg*prod_weighting), 
            onfarmGHG_unweighted_stressor = mean(mean_onfarm_ghg, na.rm = TRUE))

#_______________________________________________________________________________________________________________________#
# Calculate on-farm N and P
#_______________________________________________________________________________________________________________________#
feed_NP <-clean_feedNutrition(feedNutrition_data = 
                                read.csv(file.path(datadir, "United-States-Canadian-Tables-of-Feed-1982-pages-68-921-with_CrudeProtein.csv"),
                                         stringsAsFactors = FALSE))
feed_NP <- feed_NP %>%
  mutate(feed_type = case_when(
    (ingredient == "Soy") ~ "soy",
    (ingredient == "Crop") ~ "crops",
    (ingredient == "Fishery") ~ "fmfo",
    (ingredient == "Animal by-products") ~ "animal"
  )) %>%
  select(-c("ingredient", "sd")) %>%
  pivot_wider(names_from = "element", values_from = "value")

fish_NP <- read.csv(file.path(datadir, "fish_NP_clean.csv"))

df_onfarm_NP <- df %>%
  # Calculate N and P from each feed component
  pivot_longer(feed_soy:feed_animal, names_to = c("drop", "feed_type"), names_sep = "_", values_to = "feed_proportion") %>%
  select(-drop) %>%
  left_join(feed_NP, by = "feed_type") %>%
  mutate(feed_N = feed_proportion*(N/100), 
         feed_P = feed_proportion*(P/100)) %>%
  select(-c("feed_proportion", "N", "P")) %>%
  group_by(study_id,clean_sci_name, taxa, intensity, system, fcr, prod_weighting) %>% # removed Country, iso3c since they are temporarily missing from df
  summarise(feed_N = fcr*sum(feed_N, na.rm = TRUE), 
            feed_P = fcr*sum(feed_P, na.rm = TRUE)) %>%
  # Subtract product N and P
  left_join(fish_NP, by = c("clean_sci_name")) %>%
  mutate(N_emissions_kg_per_t = 1000*(feed_N - N_t_liveweight_t), 
         P_emissions_kg_per_t = 1000*(feed_P - P_t_liveweight_t))

df_onfarm_NP_taxa <- df_onfarm_NP %>% 
  # Calculate species means
  group_by(taxa, clean_sci_name, prod_weighting) %>% 
  summarise(mean_onfarm_N = mean(N_emissions_kg_per_t, na.rm = TRUE), 
            mean_onfarm_P = mean(P_emissions_kg_per_t, na.rm = TRUE)) %>% 
  # Calculate taxa means
  group_by(taxa) %>%
  summarise(onfarmN_weighted_stressor = sum(mean_onfarm_N*prod_weighting), 
            onfarmN_unweighted_stressor = mean(mean_onfarm_N, na.rm = TRUE),
            onfarmP_weighted_stressor = sum(mean_onfarm_P*prod_weighting), 
            onfarmP_unweighted_stressor = mean(mean_onfarm_P, na.rm = TRUE))

#_______________________________________________________________________________________________________________________#
# Calculate on-farm land
#_______________________________________________________________________________________________________________________#
df_onfarm_land_taxa <- df %>% 
  # Only count land for ponds and recirculating systems
  mutate(Yield_m2_per_t = ifelse(system %in% c("Ponds", "Recirculating and tanks"), Yield_m2_per_t, 0)) %>%
  # Calculate species means
  group_by(taxa, clean_sci_name, prod_weighting) %>% 
  summarise(mean_onfarm_land = mean(Yield_m2_per_t, na.rm = TRUE)) %>% 
  # Calculate taxa means
  group_by(taxa) %>%
  summarise(onfarmland_weighted_stressor = sum(mean_onfarm_land*prod_weighting), 
            onfarmland_unweighted_stressor = mean(mean_onfarm_land, na.rm = TRUE))

#_______________________________________________________________________________________________________________________#
# Calculate on-farm water
#_______________________________________________________________________________________________________________________#
evap <- read.csv(file.path(datadir, "20201222_clim_summarise_by_country.csv"))
evap$iso3c <- countrycode(evap$admin, origin = "country.name", destination = "iso3c")
evap <- evap %>%
  mutate(evap_rate_m3_per_m2 = mean_evap_mm/1000) %>%
  filter(!is.na(iso3c))

fw_taxa <- c("oth_carp", "catfish", "hypoph_carp", "tilapia", "trouts", "fresh_crust")

df_onfarm_water <- df %>%
  left_join(evap, by = "iso3c") %>%
  # Apply evap only to freshwater ponds
  mutate(on_farm_water = ifelse(taxa %in% fw_taxa & system %in% c("Ponds", "Recirculating and tanks"), 
                                Yield_m2_per_t*evap_rate_m3_per_m2, 0))

df_onfarm_water_taxa <- df_onfarm_water %>% 
  # Calculate species means
  group_by(taxa, clean_sci_name, prod_weighting) %>% 
  summarise(mean_onfarm_water = mean(on_farm_water, na.rm = TRUE)) %>%
  # Calculate taxa means
  group_by(taxa) %>%
  summarise(onfarmwater_weighted_stressor = sum(mean_onfarm_water*prod_weighting), 
            onfarmwater_unweighted_stressor = mean(mean_onfarm_water, na.rm = TRUE))

df_onfarm_water_taxa$onfarmwater_weighted_stressor[is.na(df_onfarm_water_taxa$onfarmwater_weighted_stressor)] <- 0
df_onfarm_water_taxa$onfarmwater_unweighted_stressor[is.na(df_onfarm_water_taxa$onfarmwater_unweighted_stressor)] <- 0

#_______________________________________________________________________________________________________________________#
# Calculate capture GHGs
#_______________________________________________________________________________________________________________________#
df_capture <- read.csv(file.path(datadir, "fisheries_fuel_use.csv"))

df_capture_ghg <- df_capture %>%
  # Remove mixed gear and nei observations
  filter(!str_detect(pattern = " nei", species)) %>%
  filter(gear != "Other, Mixed, or Unknown") %>%
  # Remove observations with 0 gear, species, or consumption weighting
  filter(gear_weighting > 0 & species_weighting > 0 & consumption_weighting > 0) %>%
  # Re-weight gear within each species
  group_by(species_group, species) %>%
  mutate(gear_weighting_new = gear_weighting/sum(gear_weighting)) %>%
  # Create species gear-weighted means
  summarise(species_ghg_kg_t = sum(ghg*gear_weighting_new), 
            species_weighting = mean(species_weighting), 
            consumption_weighting = mean(consumption_weighting)) %>%
  # Re-weight species and consumption within taxa group
  ungroup() %>%
  group_by(species_group) %>%
  mutate(species_consumption_weighting = (species_weighting*consumption_weighting)/sum(species_weighting*consumption_weighting)) %>%
  summarise(ghg_kg_t = sum(species_ghg_kg_t*species_consumption_weighting))

#_______________________________________________________________________________________________________________________#
# Summarize all by taxa
#_______________________________________________________________________________________________________________________#
df_feed_taxa <- df_feed_taxa %>%
  ungroup() %>%
  filter(Allocation == "Mass") %>%
  mutate(spread_col = case_when(
    Impact.category == "Global warming potential" ~ "feed_GHG",
    Impact.category == "Freshwater eutrophication" ~ "feed_P",
    Impact.category == "Marine eutrophication" ~ "feed_N",
    Impact.category == "Land use" ~ "feed_land",
    Impact.category == "Water consumption" ~ "feed_water"
  )) %>%
  select(taxa, spread_col, weighted_stressor) %>%
  pivot_wider(names_from = spread_col, values_from = weighted_stressor) 

stressor_taxa_summary <- df_feed_taxa %>%
  left_join(df_onfarm_ghg_taxa, by = "taxa") %>%
  left_join(df_onfarm_NP_taxa, by = "taxa") %>%
  left_join(df_onfarm_land_taxa, by = "taxa") %>%
  left_join(df_onfarm_water_taxa, by = "taxa") %>%
  mutate(total_ghg = feed_GHG + onfarmGHG_weighted_stressor, 
         total_N = feed_N + onfarmN_weighted_stressor, 
         total_P = feed_P + onfarmP_weighted_stressor, 
         total_land = feed_land + onfarmland_weighted_stressor, 
         total_water = feed_water + onfarmwater_weighted_stressor)

