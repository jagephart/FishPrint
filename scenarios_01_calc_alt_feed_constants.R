# Generate alternate feed constants for lever scenarios
# Do this for both Mass allocation and Economic allocation versions

# Libraries for processing and analyses
rm(list=ls())
library(tidyverse)
library(countrycode) # part of clean.lca
source("Functions.R")

#_________________________________________________________________________________________________________________________________#
# Load and clean data
#_________________________________________________________________________________________________________________________________#
# Mac
datadir <- "/Volumes/jgephart/BFA Environment 2/Data"
outdir <- "/Volumes/jgephart/BFA Environment 2/Outputs"

feed_fp <- read.csv(file.path(datadir, "Feed_impact_factors_20201203.csv"))
feed_fp$iso3c <- countrycode(feed_fp$Country.Region, origin = "country.name", destination = "iso3c")
feed_fp$iso3c[is.na(feed_fp$iso3c)] <- "Other"

feed_fp <- feed_fp %>% 
  mutate(iso3c = ifelse(Country.Region == "Europe" & Input %in% c("Chicken by-product meal", "Chicken by-product oil"),
                        "EUR", iso3c)) 

# DATA DOCUMENTATION: FAO export data should be downloaded from FAOSTAT
# http://www.fao.org/faostat/en/#data/TP
# Query terms for downloading: SELECT All countries, Export quantity, All items
# File can also be provided by corresponding author upon request
faostat <- read.csv(file.path(datadir, "FAOSTAT_data_12-9-2020.csv"))
faostat$iso3c <- countrycode(faostat$Area, origin = "country.name", destination = "iso3c")

##########################################################################################################################
# OPTION: Choose allocation method (used for all analyses below and filenames)
# allocation_method <- "Mass"
allocation_method <- "Economic"
##########################################################################################################################

#_________________________________________________________________________________________________________________________________#
# Load baseline feed data
#_________________________________________________________________________________________________________________________________#
feed_fp_baseline <- read.csv(file.path(datadir, "weighted_feed_fp.csv"))

#_________________________________________________________________________________________________________________________________#
# No land use soy and crops - create mass and economic allocation versions of this
#_________________________________________________________________________________________________________________________________#

weighted_soy <- calc_soy_weights(faostat, feed_fp, deforestation_free = TRUE)
weighted_crop <- calc_crop_weights(faostat, feed_fp, deforestation_free = TRUE)

# Replace baseline soy and crops with new soy and crops
feed_fp_noland <- feed_fp_baseline %>%
  filter(Input.type != "Soy") %>%
  filter(Input.type != "Crop") %>%
  bind_rows(weighted_soy) %>%
  bind_rows(weighted_crop)

# Format for analysis
feed_fp_noland <- feed_fp_noland %>%
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

feed_fp_noland_filename <- paste("feed_fp_scenario_7_", str_to_lower(allocation_method), ".csv", sep = "")
write.csv(feed_fp_noland, file.path(datadir, feed_fp_noland_filename), row.names = FALSE)

#_________________________________________________________________________________________________________________________________#
# Replace FMFO with fishery by-products - create mass and economic allocation versions of this
#_________________________________________________________________________________________________________________________________#

# Calculate fish byproduct weights
fmfo_prod <- read.csv(file.path(datadir, "fish_weightings.csv"))
weighted_fishbyproduct <- calc_byproduct_weights(fmfo_prod, feed_fp)
weighted_fishbyproduct$Input.type <- "Fishery"

feed_fp_fish_bp <- feed_fp_baseline %>%
  filter(Input.type != "Fishery") %>%
  bind_rows(weighted_fishbyproduct)

# Format for analysis
feed_fp_fish_bp <- feed_fp_fish_bp %>%
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

feed_fp_fish_bp_filename <- paste("feed_fp_scenario_2c_", str_to_lower(allocation_method), ".csv", sep = "")
write.csv(feed_fp_fish_bp, file.path(datadir, feed_fp_fish_bp_filename), row.names = FALSE)

#_________________________________________________________________________________________________________________________________#
# Replace FMFO with low impact (Alaska Pollock) fishery by-products - create mass and economic allocation versions of this
#_________________________________________________________________________________________________________________________________#

# Modified function calc_byproduct_weights for just Alaska Pollock
weighted_fishbyproduct_pollock <- feed_fp %>%
  filter(Input.type == c("Fishery by-product")) %>% 
  filter(str_detect(Input, pattern = "pollock")) %>%
  left_join(fmfo_prod, by = c("Input" = "Name")) %>%
  group_by(Input.type, Impact.category, Allocation, Units) %>%
  # Normalize weightings to sum to 1
  mutate(reweighting = Weighting/sum(Weighting, na.rm = TRUE)) %>%
  summarise(Value = sum(Value * reweighting)) %>%
  # If weighting ingredient types, do here along with country weightings
  group_by(Impact.category, Allocation, Units) %>%
  summarise(ave_stressor = mean(Value, na.rm = TRUE))

weighted_fishbyproduct_pollock$Input.type <- "Fishery"

feed_fp_fish_low_impact_bp <- feed_fp_baseline %>%
  filter(Input.type != "Fishery") %>%
  bind_rows(weighted_fishbyproduct_pollock)

# Format for analysis
# Change names to match df
feed_fp_fish_low_impact_bp <- feed_fp_fish_low_impact_bp %>%
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

feed_fp_fish_low_impact_bp_filename <- paste("feed_fp_scenario_2d_", str_to_lower(allocation_method), ".csv", sep = "")
write.csv(feed_fp_fish_low_impact_bp, file.path(datadir, feed_fp_fish_low_impact_bp_filename), row.names = FALSE)

#_________________________________________________________________________________________________________________________________#
# Low impact (Alaska Pollock) fishery by-products in FMFO - create mass and economic allocation versions of this
#_________________________________________________________________________________________________________________________________#

weighted_fishery <- calc_fishery_weights(fmfo_prod, feed_fp)
weighted_fish <- combine_fish_weights(weighted_fishery, weighted_fishbyproduct_pollock)
weighted_fish$Input.type <- "Fishery"

feed_fp_fmfo_low_impact_bp <- feed_fp_baseline %>%
  filter(Input.type != "Fishery") %>%
  bind_rows(weighted_fish)

# Format for analysis
# Change names to match df
feed_fp_fmfo_low_impact_bp <- feed_fp_fmfo_low_impact_bp %>%
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

feed_fp_fmfo_low_impact_bp_filename <- paste("feed_fp_scenario_6_", str_to_lower(allocation_method), ".csv", sep = "")
write.csv(feed_fp_fmfo_low_impact_bp, file.path(datadir, feed_fp_fmfo_low_impact_bp_filename), row.names = FALSE)
