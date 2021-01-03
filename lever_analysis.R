# Analysis of levers

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
# Load data
#_______________________________________________________________________________________________________________________#
# Load full lca data with predicted parameters
df <- read.csv(file.path(datadir, "lca-dat-imputed-vars_rep-sqrt-n-farms.csv"))

# Load and join weightings
prod_weightings <- read.csv(file.path(datadir, "aqua_prod_weightings.csv"))

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
evap <- read.csv(file.path(datadir, "20201222_clim_summarise_by_country.csv"))
evap$iso3c <- countrycode(evap$admin, origin = "country.name", destination = "iso3c")
evap <- evap %>%
  mutate(evap_rate_m3_per_m2 = mean_evap_mm/1000) %>%
  filter(!is.na(iso3c)) %>%
  select(iso3c, evap_rate_m3_per_m2)

df <- df %>%
  left_join(evap, by = "iso3c")

# Load feed data
feed_fp <- read.csv(file.path(datadir, "20201217_weighted_feed_fp.csv"))
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
  filter(Allocation == "Mass") %>%
  select(feed_type, stressor, ave_stressor)

# Load feed N/P data
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
  pivot_wider(names_from = "element", values_from = "value") %>%
  mutate(N = N/100, P = P/100) # NOTE: Check why for units

# Load and join fish N/P data 
fish_NP <- read.csv(file.path(datadir, "fish_NP_clean.csv"))
fish_NP <- fish_NP %>% 
  select(clean_sci_name, N_t_liveweight_t, P_t_liveweight_t)

df <- df %>%
  left_join(fish_NP, by = "clean_sci_name")

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
              fcr_weighted_mean = weighted.mean(fcr_mean, prod_weighting, na.rm = TRUE), fcr_sd = sum(prod_weighting * (fcr_mean - fcr_weighted_mean)^2),
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
# Change in stressors with change in parameters function
#_______________________________________________________________________________________________________________________#

stressor_sensitivity <- function(data_lca, data_feed_stressors, data_feed_NP, data_energy,
                             delta_FCR = 0, delta_soy = 0, delta_crops = 0, delta_animal = 0, delta_fmfo = 0,
                             delta_electricity = 0, delta_diesel = 0, delta_petrol = 0, delta_natgas = 0, delta_yield = 0){
  # Set feed component stressor constants
  soy_stressor_ghg <- data_feed_stressors %>% filter(stressor == "ghg" & feed_type == "soy") %>% pull(ave_stressor)
  soy_stressor_N <- data_feed_stressors %>% filter(stressor == "N" & feed_type == "soy") %>% pull(ave_stressor)
  soy_stressor_P <- data_feed_stressors %>% filter(stressor == "P" & feed_type == "soy") %>% pull(ave_stressor)
  soy_stressor_land <- data_feed_stressors %>% filter(stressor == "land" & feed_type == "soy") %>% pull(ave_stressor)
  soy_stressor_water <- data_feed_stressors %>% filter(stressor == "water" & feed_type == "soy") %>% pull(ave_stressor)
  
  crops_stressor_ghg <- data_feed_stressors %>% filter(stressor == "ghg" & feed_type == "crops") %>% pull(ave_stressor)
  crops_stressor_N <- data_feed_stressors %>% filter(stressor == "N" & feed_type == "crops") %>% pull(ave_stressor)
  crops_stressor_P <- data_feed_stressors %>% filter(stressor == "P" & feed_type == "crops") %>% pull(ave_stressor)
  crops_stressor_land <- data_feed_stressors %>% filter(stressor == "land" & feed_type == "crops") %>% pull(ave_stressor)
  crops_stressor_water <- data_feed_stressors %>% filter(stressor == "water" & feed_type == "crops") %>% pull(ave_stressor)
  
  animal_stressor_ghg <- data_feed_stressors %>% filter(stressor == "ghg" & feed_type == "animal") %>% pull(ave_stressor)
  animal_stressor_N <- data_feed_stressors %>% filter(stressor == "N" & feed_type == "animal") %>% pull(ave_stressor)
  animal_stressor_P <- data_feed_stressors %>% filter(stressor == "P" & feed_type == "animal") %>% pull(ave_stressor)
  animal_stressor_land <- data_feed_stressors %>% filter(stressor == "land" & feed_type == "animal") %>% pull(ave_stressor)
  animal_stressor_water <- data_feed_stressors %>% filter(stressor == "water" & feed_type == "animal") %>% pull(ave_stressor)
  
  fmfo_stressor_ghg <- data_feed_stressors %>% filter(stressor == "ghg" & feed_type == "fmfo") %>% pull(ave_stressor)
  fmfo_stressor_N <- data_feed_stressors %>% filter(stressor == "N" & feed_type == "fmfo") %>% pull(ave_stressor)
  fmfo_stressor_P <- data_feed_stressors %>% filter(stressor == "P" & feed_type == "fmfo") %>% pull(ave_stressor)
  fmfo_stressor_land <- data_feed_stressors %>% filter(stressor == "land" & feed_type == "fmfo") %>% pull(ave_stressor)
  fmfo_stressor_water <- data_feed_stressors %>% filter(stressor == "water" & feed_type == "fmfo") %>% pull(ave_stressor)
  
  # Set energy constants
  diesel_ghg <- data_energy %>% filter(Input == "Diesel_L") %>% pull(Value)
  petrol_ghg <- data_energy %>% filter(Input == "Petrol_L") %>% pull(Value)
  natgas_ghg <- data_energy %>% filter(Input == "NaturalGas_L") %>% pull(Value)
  
  # Set feed N and P constants
  N_content_soy <- data_feed_NP %>% filter(feed_type == "soy") %>% pull(N)
  P_content_soy <- data_feed_NP %>% filter(feed_type == "soy") %>% pull(P)
  
  N_content_crops <- data_feed_NP %>% filter(feed_type == "crops") %>% pull(N)
  P_content_crops <- data_feed_NP %>% filter(feed_type == "crops") %>% pull(P)
  
  N_content_animal <- data_feed_NP %>% filter(feed_type == "animal") %>% pull(N)
  P_content_animal <- data_feed_NP %>% filter(feed_type == "animal") %>% pull(P)
  
  N_content_fmfo <- data_feed_NP %>% filter(feed_type == "fmfo") %>% pull(N)
  P_content_fmfo <- data_feed_NP %>% filter(feed_type == "fmfo") %>% pull(P)

  data_lca <- data_lca %>%
    # Add deltas to calculate perturbation
    mutate(fcr_weighted_mean = fcr_weighted_mean + delta_FCR,
           feed_soy_weighted_mean = feed_soy_weighted_mean + delta_soy, 
           feed_crops_weighted_mean = feed_crops_weighted_mean + delta_crops,
           feed_animal_weighted_mean = feed_animal_weighted_mean + delta_animal,
           feed_fmfo_weighted_mean = feed_fmfo_weighted_mean + delta_fmfo,
           electricity_weighted_mean = electricity_weighted_mean + delta_electricity,
           diesel_weighted_mean = diesel_weighted_mean + delta_diesel, 
           petrol_weighted_mean = petrol_weighted_mean + delta_petrol, 
           naturalgas_weighted_mean = naturalgas_weighted_mean + delta_natgas, 
           yield_weighted_mean = yield_weighted_mean + delta_yield) %>% 
    # Reweight feed components to ensure they still sum to 1
    mutate(feed_soy_weighted_mean_new = feed_soy_weighted_mean/(feed_soy_weighted_mean + feed_crops_weighted_mean +
                                                               feed_animal_weighted_mean + feed_fmfo_weighted_mean),
           feed_crops_weighted_mean_new = feed_crops_weighted_mean/(feed_soy_weighted_mean + feed_crops_weighted_mean +
                                                                      feed_animal_weighted_mean + feed_fmfo_weighted_mean),
           feed_animal_weighted_mean_new = feed_animal_weighted_mean/(feed_soy_weighted_mean + feed_crops_weighted_mean +
                                                                        feed_animal_weighted_mean + feed_fmfo_weighted_mean),
           feed_fmfo_weighted_mean_new = feed_fmfo_weighted_mean/(feed_soy_weighted_mean + feed_crops_weighted_mean +
                                                                    feed_animal_weighted_mean + feed_fmfo_weighted_mean)) %>%
    mutate(feed_soy_weighted_mean_new = ifelse(is.na(feed_soy_weighted_mean_new), 0, feed_soy_weighted_mean_new),
           feed_crops_weighted_mean_new = ifelse(is.na(feed_crops_weighted_mean_new), 0, feed_crops_weighted_mean_new),
           feed_animal_weighted_mean_new = ifelse(is.na(feed_animal_weighted_mean_new), 0, feed_animal_weighted_mean_new),
           feed_fmfo_weighted_mean_new = ifelse(is.na(feed_fmfo_weighted_mean_new), 0, feed_fmfo_weighted_mean_new)) %>%
    # Calculate all stressors
    mutate(
      # Calculate feed-associate components
      feed_ghg = 1000*fcr_weighted_mean*(feed_soy_weighted_mean_new*soy_stressor_ghg + 
                             feed_crops_weighted_mean_new*crops_stressor_ghg + 
                             feed_animal_weighted_mean_new*animal_stressor_ghg + 
                             feed_fmfo_weighted_mean_new*fmfo_stressor_ghg),
      
      feed_N = 1000*fcr_weighted_mean*(feed_soy_weighted_mean_new*soy_stressor_N + 
                             feed_crops_weighted_mean_new*crops_stressor_N + 
                             feed_animal_weighted_mean_new*animal_stressor_N + 
                             feed_fmfo_weighted_mean_new*fmfo_stressor_N),
      
      feed_P = 1000*fcr_weighted_mean*(feed_soy_weighted_mean_new*soy_stressor_P + 
                             feed_crops_weighted_mean_new*crops_stressor_P + 
                             feed_animal_weighted_mean_new*animal_stressor_P + 
                             feed_fmfo_weighted_mean_new*fmfo_stressor_P),
      
      feed_land = 1000*fcr_weighted_mean*(feed_soy_weighted_mean_new*soy_stressor_land + 
                             feed_crops_weighted_mean_new*crops_stressor_land + 
                             feed_animal_weighted_mean_new*animal_stressor_land + 
                             feed_fmfo_weighted_mean_new*fmfo_stressor_land),
      
      feed_water = 1000*fcr_weighted_mean*(feed_soy_weighted_mean_new*soy_stressor_water + 
                             feed_crops_weighted_mean_new*crops_stressor_water + 
                             feed_animal_weighted_mean_new*animal_stressor_water + 
                             feed_fmfo_weighted_mean_new*fmfo_stressor_water), 
      
      # Calculate on farm GHG
      onfarm_ghg = electricity_weighted_mean*electricity_ghg_weighted_mean +
        diesel_weighted_mean*diesel_ghg + petrol_weighted_mean*petrol_ghg + naturalgas_weighted_mean*natgas_ghg,
      
      # Calculate on farm N and P
      onfarm_N = 1000*((fcr_weighted_mean*(feed_soy_weighted_mean_new*N_content_soy + 
                        feed_crops_weighted_mean_new*N_content_crops + 
                        feed_animal_weighted_mean_new*N_content_animal + 
                        feed_fmfo_weighted_mean_new*N_content_fmfo)) - fish_N_weighted_mean),
      
      onfarm_P = 1000*((fcr_weighted_mean*(feed_soy_weighted_mean_new*P_content_soy + 
                        feed_crops_weighted_mean_new*P_content_crops + 
                        feed_animal_weighted_mean_new*P_content_animal + 
                        feed_fmfo_weighted_mean_new*P_content_fmfo)) - fish_P_weighted_mean),
      
      # Calculate on farm land
      onfarm_land = yield_weighted_mean,
      
      # Calculate on farm water (for freshwater species only)
      onfarm_water =  ifelse(taxa %in% c("oth_carp", "catfish", "hypoph_carp", "tilapia", "trouts", "fresh_crust"), 
                             onfarm_land*evap_weighted_mean, 0),
      
      # Calculate totals
      total_ghg = feed_ghg + onfarm_ghg,
      total_N = feed_N + onfarm_N,
      total_P = feed_P + onfarm_P,
      total_land = feed_land + onfarm_land,
      total_water = feed_water + onfarm_water
    ) %>%
    select(taxa, feed_ghg:total_water)
  
  return(data_lca)
}

#_______________________________________________________________________________________________________________________#
# Loop through and perturb each parameter
#_______________________________________________________________________________________________________________________#
# Estimates without any perturbation
stressor_0 <- stressor_sensitivity(data_lca = df_taxa, data_feed_stressors = feed_fp, data_feed_NP = feed_NP, data_energy = energy_gwp,
                             delta_FCR = 0, delta_soy = 0, delta_crops = 0, delta_animal = 0, delta_fmfo = 0,
                             delta_electricity = 0, delta_diesel = 0, delta_petrol = 0, delta_natgas = 0, delta_yield = 0)
 
stressor_1 <- stressor_sensitivity(data_lca = df_taxa, data_feed_stressors = feed_fp, data_feed_NP = feed_NP, data_energy = energy_gwp,
                                   delta_FCR = -0.1*df_taxa$fcr_weighted_mean, delta_soy = 0, delta_crops = 0, delta_animal = 0, delta_fmfo = 0,
                                   delta_electricity = 0, delta_diesel = 0, delta_petrol = 0, delta_natgas = 0, delta_yield = 0)



  
  
  

