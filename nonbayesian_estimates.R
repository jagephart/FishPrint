# Non-Bayesian comparison estimates

#_______________________________________________________________________________________________________________________#
# Load packages and source functions
#_______________________________________________________________________________________________________________________#
library(tidyverse)
library(ggplot2)
library(ggthemes)
library(countrycode)

source("Functions.R")

# Set data directories
datadir <- "/Volumes/jgephart/BFA Environment 2/Data"
outdir <- "/Volumes/jgephart/BFA Environment 2/Outputs"

#_______________________________________________________________________________________________________________________#
# Load lca data
#_______________________________________________________________________________________________________________________#
# Load full data
df <- read.csv(file.path(outdir, "dat_for_bayesian-carbon.csv"))

taxa_weights <- read.csv(file.path(datadir, "aqua_prod_weightings.csv"))

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
  group_by(study_id, clean_sci_name, taxa, intensity, system, Impact.category, Allocation, Units) %>%
  summarise(stressor = sum(stressor))


df_feed_taxa_unweighted <- df_feed %>%
  group_by(taxa, Impact.category, Allocation, Units) %>%
  summarise(weighted_stressor = mean(stressor)) %>%
  # Convert units to per tonne
  mutate(weighted_stressor = 1000*weighted_stressor)

# Calculate production-weighted average
# Not currently working since weightings no longer sum to 1 within a taxa group after the merge
# df_feed_taxa <- df_feed %>%
#   ungroup() %>%
#   group_by(clean_sci_name, taxa, Impact.category, Allocation, Units) %>%
#   # calculate mean species stressor values
#   summarise(stressor = mean(stressor), .groups = 'drop') %>%
#   left_join(taxa_weights, by = c("clean_sci_name", "taxa")) %>%
#   mutate(weighted_stressor = stressor * weighting) %>%
#   group_by(taxa, Impact.category, Allocation, Units) %>%
#   summarise(taxa_weighted_stressor = sum(weighted_stressor)) %>%
#   # Convert units to per tonne
#   mutate(taxa_weighted_stressor = 1000*taxa_weighted_stressor)

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
    (Input == "Diesel") ~ "diesel",
    (Input == "Petrol") ~ "petrol", 
    (Input == "Natural gas") ~"natgas"
  ))

df_onfarm_ghg <- df %>%
  # Multiply diesel, petrol and natgas by associated GHG emissions
  pivot_longer(c("diesel", "petrol", "natgas"), names_to = "Input", values_to = "fuel_L") %>%
  left_join(energy_gwp, by = "Input") %>%
  mutate(fuel_ghg_kgCO2 = fuel_L*Value) %>%
  select(-c("fuel_L", "Value")) %>%
  pivot_wider(names_from = Input, values_from = fuel_ghg_kgCO2) %>%
  # Multiply electricity by country electricity GHG
  left_join(electricity_gwp, by = "iso3c") %>%
  mutate(electricity_ghg_kgCO2 = electric*GWP_perkWh_kgCO2eq) %>% 
  mutate(onfarm_ghg_kgCO2 = diesel + petrol + natgas + electricity_ghg_kgCO2) %>% 
  select("study_id", "clean_sci_name", "taxa", "intensity", "system", "Country" = "Country.x", "iso3c", "onfarm_ghg_kgCO2")

df_onfarm_ghg_taxa <- df_onfarm_ghg %>% 
  group_by(taxa) %>% 
  summarise(mean_onfarm_ghg = mean(onfarm_ghg_kgCO2))

#_______________________________________________________________________________________________________________________#
# Calculate on-farm N and P
#_______________________________________________________________________________________________________________________#
feed_NP <-clean_feedNutrition(feedNutrition_data = 
                                read.csv(file.path(datadir, "United-States-Canadian-Tables-of-Feed-1982-pages-68-921-with_CrudeProtein.csv"),
                                         stringsAsFactors = FALSE))
# FIX IT: Check with Alon on N and P units - currently assuming % N and P for dry matter

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
  group_by(study_id,Country, iso3c, clean_sci_name, taxa, intensity, system, fcr) %>%
  summarise(feed_N = fcr*sum(feed_N), 
            feed_P = fcr*sum(feed_P)) %>%
  # Subtract product N and P
  left_join(fish_NP, by = c("clean_sci_name" = "Scientific.Name")) %>%
  mutate(N_emissions_kg_per_t = 1000*(feed_N - N_t_liveweight_t), 
         P_emissions_kg_per_t = 1000*(feed_P - P_t_liveweight_t))

# FIX IT: These values seem quite high. Talk to Alon about going through feed N and P code and for someone to double check all units

df_onfarm_NP_taxa <- df_onfarm_NP %>% 
  group_by(taxa) %>% 
  summarise(mean_onfarm_N = mean(N_emissions_kg_per_t), 
            mean_onfarm_P = mean(P_emissions_kg_per_t))

#_______________________________________________________________________________________________________________________#
# Calculate on-farm land
#_______________________________________________________________________________________________________________________#
df_onfarm_land <- read.csv(file.path(outdir, "dat_for_bayesian-land.csv"))

df_onfarm_land_taxa <- df_onfarm_land %>% 
  group_by(taxa) %>% 
  summarise(mean_onfarm_land = mean(yield))

#_______________________________________________________________________________________________________________________#
# Calculate on-farm water
#_______________________________________________________________________________________________________________________#
evap <- read.csv(file.path(datadir, "20201222_clim_summarise_by_country.csv"))
evap$iso3c <- countrycode(evap$admin, origin = "country.name", destination = "iso3c")
evap <- evap %>%
  mutate(evap_rate_m3_per_m2 = mean_evap_mm/1000) %>%
  filter(!is.na(iso3c))

# Bring back country codes
cc <- df %>%
  select(study_id, iso3c)

fw_taxa <- c("oth_carp", "catfish", "hypoph_carp", "tilapia", "trouts", "fresh_crust")
  
# FIX IT: Units currently unknown
df_onfarm_water <- df_onfarm_land %>%
  left_join(cc, by = "study_id") %>%
  left_join(evap, by = "iso3c") %>%
  # Apply evap only to freshwater ponds
  mutate(on_farm_water = ifelse(taxa %in% fw_taxa & system %in% c("Ponds", "Recirculating and tanks"), 
                                yield*evap_rate_m3_per_m2, 0))

df_onfarm_water_taxa <- df_onfarm_water %>% 
  group_by(taxa) %>% 
  summarise(mean_onfarm_water = mean(on_farm_water, na.rm = TRUE))

df_onfarm_water_taxa$mean_onfarm_water[is.na(df_onfarm_water_taxa$mean_onfarm_water)] <- 0
  
#_______________________________________________________________________________________________________________________#
# Summarize all by taxa
#_______________________________________________________________________________________________________________________#
df_feed_taxa <- df_feed_taxa_unweighted %>%
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
  mutate(total_ghg = feed_GHG + mean_onfarm_ghg, 
            total_N = feed_N + mean_onfarm_N, 
            total_P = feed_P + mean_onfarm_P, 
            total_land = feed_land + mean_onfarm_land, 
            total_water = feed_water + mean_onfarm_water)
 