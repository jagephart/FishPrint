# Calculate chicken stressors

#_______________________________________________________________________________________________________________________#
# Load packages and source functions
#_______________________________________________________________________________________________________________________#
library(tidyverse)
library(countrycode)
library(ggpubr)

source("Functions.R")

# Set data directories
datadir <- "/Volumes/jgephart/BFA Environment 2/Data"
outdir <- "/Volumes/jgephart/BFA Environment 2/Outputs"

#_______________________________________________________________________________________________________________________#
# Load lca data
#_______________________________________________________________________________________________________________________#
# Load full data with predicted parameters
df <- read.csv(file.path(datadir, "poultry_parameters.csv"), header = TRUE)
df <- df %>% 
  filter(Country != "NA") %>%
  mutate(NaturalGas_L = NaturalGas_MJ/0.039) %>%
  mutate(iso3c = countrycode(Country, origin = "country.name", destination = "iso3c")) %>% 
  # Rename columes to us same function
  select(iso3c, taxa, "fcr_weighted_mean" = "FCR", "feed_soy_weighted_mean" = "soy", "feed_crops_weighted_mean" = "crops",
         "feed_animal_weighted_mean" = "animal", "feed_fmfo_weighted_mean" = "fmfo", 
         "electricity_weighted_mean" = "Electricity_kWh", "naturalgas_weighted_mean" = "NaturalGas_L",
         diesel_weighted_mean, petrol_weighted_mean, yield_weighted_mean, fish_N_weighted_mean, fish_P_weighted_mean)

# Load and join electricity GHG constants
electricity_gwp <- read.csv(file.path(datadir, "electricity_GWP.csv"))
electricity_gwp$iso3c <- countrycode(electricity_gwp$Country, origin = "country.name", destination = "iso3c")

df <- df %>% 
  left_join(electricity_gwp %>% select(-Country), by = "iso3c") %>% 
  rename("electricity_ghg_weighted_mean" = "GWP_perkWh_kgCO2eq")

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
  filter(Allocation == "Mass") %>%
  select(feed_type, stressor, ave_stressor)

# Load feed N/P data
feed_NP <- read.csv(file.path(outdir, "feed_NP_clean.csv"))

#_______________________________________________________________________________________________________________________#
# Stressor function
#_______________________________________________________________________________________________________________________#

stressor_sensitivity <- function(data_lca, data_feed_stressors, data_feed_NP, data_energy,
                                 delta_FCR = 0, delta_soy = 0, delta_crops = 0, delta_animal = 0, delta_fmfo = 0,
                                 delta_electricity = 0, delta_diesel = 0, delta_petrol = 0, delta_natgas = 0, delta_yield = 0,
                                 # Optional changes to constants 
                                 land_change_free_soy = NULL, land_change_free_crops = NULL, 
                                 fmfo_with_land_change_free_crops = NULL,
                                 fishery_byproducts = NULL, low_impact_fishery_byproducts = NULL,
                                 zero_emissions_electric = FALSE){
  # Set feed component stressor constants
  soy_stressor_ghg <- data_feed_stressors %>% filter(stressor == "ghg" & feed_type == "soy") %>% pull(ave_stressor)
  soy_stressor_N <- data_feed_stressors %>% filter(stressor == "N" & feed_type == "soy") %>% pull(ave_stressor)
  soy_stressor_P <- data_feed_stressors %>% filter(stressor == "P" & feed_type == "soy") %>% pull(ave_stressor)
  soy_stressor_land <- data_feed_stressors %>% filter(stressor == "land" & feed_type == "soy") %>% pull(ave_stressor)
  soy_stressor_water <- data_feed_stressors %>% filter(stressor == "water" & feed_type == "soy") %>% pull(ave_stressor)
  
  if(is.null(land_change_free_soy) == FALSE){
    soy_stressor_ghg <- land_change_free_soy %>% filter(stressor == "ghg" & feed_type == "soy") %>% pull(ave_stressor)
    soy_stressor_N <- land_change_free_soy %>% filter(stressor == "N" & feed_type == "soy") %>% pull(ave_stressor)
    soy_stressor_P <- land_change_free_soy %>% filter(stressor == "P" & feed_type == "soy") %>% pull(ave_stressor)
    soy_stressor_land <- land_change_free_soy %>% filter(stressor == "land" & feed_type == "soy") %>% pull(ave_stressor)
    soy_stressor_water <- land_change_free_soy %>% filter(stressor == "water" & feed_type == "soy") %>% pull(ave_stressor)
  }else{}
  
  crops_stressor_ghg <- data_feed_stressors %>% filter(stressor == "ghg" & feed_type == "crops") %>% pull(ave_stressor)
  crops_stressor_N <- data_feed_stressors %>% filter(stressor == "N" & feed_type == "crops") %>% pull(ave_stressor)
  crops_stressor_P <- data_feed_stressors %>% filter(stressor == "P" & feed_type == "crops") %>% pull(ave_stressor)
  crops_stressor_land <- data_feed_stressors %>% filter(stressor == "land" & feed_type == "crops") %>% pull(ave_stressor)
  crops_stressor_water <- data_feed_stressors %>% filter(stressor == "water" & feed_type == "crops") %>% pull(ave_stressor)
  
  if(is.null(land_change_free_crops) == FALSE){
    crops_stressor_ghg <- land_change_free_crops %>% filter(stressor == "ghg" & feed_type == "crops") %>% pull(ave_stressor)
    crops_stressor_N <- land_change_free_crops %>% filter(stressor == "N" & feed_type == "crops") %>% pull(ave_stressor)
    crops_stressor_P <- land_change_free_crops %>% filter(stressor == "P" & feed_type == "crops") %>% pull(ave_stressor)
    crops_stressor_land <- land_change_free_crops %>% filter(stressor == "land" & feed_type == "crops") %>% pull(ave_stressor)
    crops_stressor_water <- land_change_free_crops %>% filter(stressor == "water" & feed_type == "crops") %>% pull(ave_stressor)
  }else{}
  
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
  
  if(is.null(fishery_byproducts) == FALSE){
    fmfo_stressor_ghg <- fishery_byproducts %>% filter(stressor == "ghg" & feed_type == "fmfo") %>% pull(ave_stressor)
    fmfo_stressor_N <- fishery_byproducts %>% filter(stressor == "N" & feed_type == "fmfo") %>% pull(ave_stressor)
    fmfo_stressor_P <- fishery_byproducts %>% filter(stressor == "P" & feed_type == "fmfo") %>% pull(ave_stressor)
    fmfo_stressor_land <- fishery_byproducts %>% filter(stressor == "land" & feed_type == "fmfo") %>% pull(ave_stressor)
    fmfo_stressor_water <- fishery_byproducts %>% filter(stressor == "water" & feed_type == "fmfo") %>% pull(ave_stressor)
  }else{}
  
  if(is.null(low_impact_fishery_byproducts) == FALSE){
    fmfo_stressor_ghg <- low_impact_fishery_byproducts %>% filter(stressor == "ghg" & feed_type == "fmfo") %>% pull(ave_stressor)
    fmfo_stressor_N <- low_impact_fishery_byproducts %>% filter(stressor == "N" & feed_type == "fmfo") %>% pull(ave_stressor)
    fmfo_stressor_P <- low_impact_fishery_byproducts %>% filter(stressor == "P" & feed_type == "fmfo") %>% pull(ave_stressor)
    fmfo_stressor_land <- low_impact_fishery_byproducts %>% filter(stressor == "land" & feed_type == "fmfo") %>% pull(ave_stressor)
    fmfo_stressor_water <- low_impact_fishery_byproducts %>% filter(stressor == "water" & feed_type == "fmfo") %>% pull(ave_stressor)
  }else{}
  
  if(is.null(fmfo_with_land_change_free_crops) == FALSE){
    fmfo_stressor_ghg <- fmfo_with_land_change_free_crops %>% filter(stressor == "ghg" & feed_type == "fmfo") %>% pull(ave_stressor)
    fmfo_stressor_N <- fmfo_with_land_change_free_crops %>% filter(stressor == "N" & feed_type == "fmfo") %>% pull(ave_stressor)
    fmfo_stressor_P <- fmfo_with_land_change_free_crops %>% filter(stressor == "P" & feed_type == "fmfo") %>% pull(ave_stressor)
    fmfo_stressor_land <- fmfo_with_land_change_free_crops %>% filter(stressor == "land" & feed_type == "fmfo") %>% pull(ave_stressor)
    fmfo_stressor_water <- fmfo_with_land_change_free_crops %>% filter(stressor == "water" & feed_type == "fmfo") %>% pull(ave_stressor)
  }else{}
  
  # Set energy constants
  diesel_ghg <- data_energy %>% filter(Input == "Diesel_L") %>% pull(Value)
  petrol_ghg <- data_energy %>% filter(Input == "Petrol_L") %>% pull(Value)
  natgas_ghg <- data_energy %>% filter(Input == "NaturalGas_L") %>% pull(Value)
  
  if(zero_emissions_electric == TRUE){
    data_lca$electricity_ghg_weighted_mean <- 0
  }else{}
  
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
    # Add grow out period constants
    mutate(grow_out_yr_prop = case_when(
      taxa == "oth_carp" ~ 300/365,
      taxa == "hypoph_carp" ~ 300/365,
      taxa == "catfish" ~ 210/365,
      taxa == "tilapia" ~ 200/365,
      taxa == "trouts" ~ 365/365,
      taxa == "fresh_crust" ~ 240/365
    )) %>%
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
                             onfarm_land*evap_weighted_mean*grow_out_yr_prop, 0),
      
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
# Calculate poultry stressors
#_______________________________________________________________________________________________________________________#
stressor_summary <- stressor_sensitivity(data_lca = df, data_feed_stressors = feed_fp, data_feed_NP = feed_NP, data_energy = energy_gwp,
                                   delta_FCR = 0, delta_soy = 0, delta_crops = 0, delta_animal = 0, delta_fmfo = 0,
                                   delta_electricity = 0, delta_diesel = 0, delta_petrol = 0, delta_natgas = 0, delta_yield = 0)


df_plot <- stressor_summary %>%
  select(-c("total_ghg", "total_N", "total_P", "total_land", "total_water")) %>%
  group_by(taxa) %>%
  summarise_all(mean) %>%
  pivot_longer(feed_ghg:onfarm_water, names_sep = "_", 
               names_to = c("source", "stressor")) %>% 
  mutate(median = value*(1/0.40)) %>% # Convert to edible portion
  rename("plot_taxa_name" = "taxa")

base_size <- 12
base_family <- "sans"

# GHG plot
ghg_plot <- ggplot(df_plot %>% filter(stressor == "ghg"), 
                   aes(x = median, y = plot_taxa_name, fill = source)) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(values = c("#70468C", "#57D182")) +
  labs(title = "GHG", x = expression("kg CO"[2]~"t"^"-1"), y = "") +
  theme(axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"), 
        axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"), 
        axis.text = element_text(size = ceiling(base_size*0.7), colour = "black"),
        axis.title = element_text(size = ceiling(base_size*0.8)), 
        axis.text.x = element_text(angle = 90),
        panel.grid.minor = element_blank(), 
        panel.grid.major.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        panel.background = element_blank(), panel.border = element_blank(),
        strip.background = element_rect(linetype = 1, fill = "white"), 
        strip.text = element_text(), 
        strip.text.x = element_text(vjust = 0.5), 
        strip.text.y = element_text(angle = -90), 
        strip.placement.y = "outside",
        legend.position="bottom",
        plot.title = element_text(size = ceiling(base_size*1.1), face = "bold"), 
        plot.subtitle = element_text(size = ceiling(base_size*1.05)))

# N plot
N_plot <- ggplot(df_plot %>% filter(stressor == "N"), 
                 aes(x = median, y = plot_taxa_name, fill = source)) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(values = c("#70468C", "#57D182")) +
  labs(title = "Nitrogen", x = expression("kg N-eq t"^"-1"), y = "") +
  theme(axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"), 
        axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"), 
        axis.text = element_blank(),
        axis.title = element_text(size = ceiling(base_size*0.8)), 
        axis.text.x = element_text(angle = 90),
        panel.grid.minor = element_blank(), 
        panel.grid.major.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        panel.background = element_blank(), panel.border = element_blank(),
        strip.background = element_rect(linetype = 1, fill = "white"), 
        strip.text = element_text(), 
        strip.text.x = element_text(vjust = 0.5), 
        strip.text.y = element_text(angle = -90), 
        strip.placement.y = "outside",
        legend.position="bottom",
        plot.title = element_text(size = ceiling(base_size*1.1), face = "bold"), 
        plot.subtitle = element_text(size = ceiling(base_size*1.05)))

# P plot
P_plot <- ggplot(df_plot %>% filter(stressor == "P"), 
                 aes(x = median, y = plot_taxa_name, fill = source)) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(values = c("#70468C", "#57D182")) +
  labs(title = "Phosphorus", x = expression("kg P-eq t"^"-1"), y = "") +
  theme(axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"), 
        axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"), 
        axis.text = element_blank(),
        axis.title = element_text(size = ceiling(base_size*0.8)), 
        axis.text.x = element_text(angle = 90),
        panel.grid.minor = element_blank(), 
        panel.grid.major.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        panel.background = element_blank(), panel.border = element_blank(),
        strip.background = element_rect(linetype = 1, fill = "white"), 
        strip.text = element_text(), 
        strip.text.x = element_text(vjust = 0.5), 
        strip.text.y = element_text(angle = -90), 
        strip.placement.y = "outside",
        legend.position="bottom",
        plot.title = element_text(size = ceiling(base_size*1.1), face = "bold"), 
        plot.subtitle = element_text(size = ceiling(base_size*1.05)))

# Land plot
land_plot <- ggplot(df_plot %>% filter(stressor == "land"), 
                    aes(x = median, y = plot_taxa_name, fill = source)) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(values = c("#70468C", "#57D182")) +
  labs(title = "Land", x = expression("m"^2~"t"^"-1"), y = "") +
  theme(axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"), 
        axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"), 
        axis.text = element_blank(),
        axis.title = element_text(size = ceiling(base_size*0.8)), 
        axis.text.x = element_text(angle = 90),
        panel.grid.minor = element_blank(), 
        panel.grid.major.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        panel.background = element_blank(), panel.border = element_blank(),
        strip.background = element_rect(linetype = 1, fill = "white"), 
        strip.text = element_text(), 
        strip.text.x = element_text(vjust = 0.5), 
        strip.text.y = element_text(angle = -90), 
        strip.placement.y = "outside",
        legend.position="bottom",
        plot.title = element_text(size = ceiling(base_size*1.1), face = "bold"), 
        plot.subtitle = element_text(size = ceiling(base_size*1.05)))

# Water plot
water_plot <- ggplot(df_plot %>% filter(stressor == "water"), 
                     aes(x = median, y = plot_taxa_name, fill = source)) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(values = c("#70468C", "#57D182")) +
  labs(title = "Water", x = expression("m"^"3"~"t"^"-1"), y = "") +
  theme(axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"), 
        axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"), 
        axis.text = element_blank(),
        axis.title = element_text(size = ceiling(base_size*0.8)), 
        axis.text.x = element_text(angle = 90),
        panel.grid.minor = element_blank(), 
        panel.grid.major.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        panel.background = element_blank(), panel.border = element_blank(),
        strip.background = element_rect(linetype = 1, fill = "white"), 
        strip.text = element_text(), 
        strip.text.x = element_text(vjust = 0.5), 
        strip.text.y = element_text(angle = -90), 
        strip.placement.y = "outside",
        legend.position="bottom",
        plot.title = element_text(size = ceiling(base_size*1.1), face = "bold"), 
        plot.subtitle = element_text(size = ceiling(base_size*1.05)))

png(file.path(outdir, "plot_poultry_stressor_by_source.png"), width = 8.5, height = 2, units = "in", res = 300)
ggarrange(ghg_plot, N_plot, P_plot, water_plot, land_plot, nrow = 1,
          common.legend = TRUE, legend = "bottom", align = "h", widths = c(1.1, 0.7, 0.7, 0.7, 0.7))
dev.off()


# Stats for comparison
stressor_summary_table <- stressor_summary %>%
  select(c("taxa", "total_ghg", "total_N", "total_P", "total_land", "total_water")) %>%
  mutate(total_ghg = total_ghg*(1/0.4), total_N = total_N*(1/0.4),
         total_P = total_P*(1/0.4), total_land = total_land*(1/0.4),
         total_water = total_water*(1/0.4)) %>%
  group_by(taxa) %>%
  summarise_all(list("mean", "min", "max"))

# Summary for Table S6
summary(stressor_summary %>% 
          select(c("total_ghg", "total_N", "total_P", "total_land", "total_water")) %>% 
          mutate(total_ghg = total_ghg*(1/0.4), total_N = total_N*(1/0.4),
                 total_P = total_P*(1/0.4), total_land = total_land*(1/0.4),
                 total_water = total_water*(1/0.4)))
