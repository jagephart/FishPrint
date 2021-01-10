# Generate alternate feed constants for scenarios

# Libraries for processing and analyses
library(tidyverse)
library(countrycode) # part of clean.lca

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

faostat <- read.csv(file.path(datadir, "FAOSTAT_data_12-9-2020.csv"))
faostat$iso3c <- countrycode(faostat$Area, origin = "country.name", destination = "iso3c")

# Set allocation method
allocation_method <- "Mass"

#_________________________________________________________________________________________________________________________________#
# Load baseline feed data
#_________________________________________________________________________________________________________________________________#
feed_fp_baseline <- read.csv(file.path(datadir, "20201217_weighted_feed_fp.csv"))

#_________________________________________________________________________________________________________________________________#
# No land use soy and crops
#_________________________________________________________________________________________________________________________________#
# Calculate soy weightings
# Since the specific soy products don't map onto the FAO items, use all soy exports for all soy products
weightings <-  faostat %>% 
  filter(Unit == "tonnes") %>% 
  filter(Item %in% c("Soybeans", "Cake, soybeans", "Oil, soybean")) %>%
  filter(iso3c != "BRA") %>%
  group_by(iso3c) %>%
  summarise(Exports = sum(Value, na.rm = TRUE)) %>%
  filter(Exports > 0) %>% 
  left_join(feed_fp %>% filter(Input.type == "Soy"), by = c("iso3c")) %>%
  filter(is.na(Input.type) == FALSE) %>%
  group_by(iso3c) %>%
  summarise(Exports = sum(Exports, na.rm = TRUE)) %>%
  mutate(weighting = Exports/sum(Exports)) %>%
  select(iso3c, weighting)

weighted_soy <- feed_fp %>% 
  filter(Input.type == "Soy") %>%
  left_join(weightings, by = "iso3c") %>%
  filter(is.na(weighting) == FALSE) %>%
  group_by(Input.type, Input, Impact.category, Allocation, Units) %>%
  # Normalize weightings to sum to 1
  mutate(reweighting = weighting/sum(weighting, na.rm = TRUE)) %>%
  summarise(Value = sum(Value * reweighting)) %>%
  # If weighting ingredient types, do here along with country weightings
  group_by(Input.type, Impact.category, Allocation, Units) %>%
  summarise(ave_stressor = mean(Value, na.rm = TRUE))

# Calculate crop weightings
weightings <- faostat %>% 
  filter(Unit == "tonnes") %>% 
  filter(iso3c != "ARG") %>%
  group_by(iso3c, Item) %>%
  summarise(Exports = sum(Value, na.rm = TRUE)) %>%
  filter(Exports > 0) %>% 
  mutate(Input = case_when(
    (Item %in% c("Cassava Equivalent")) ~ "Cassava",
    (Item %in% c("Maize")) ~ "Maize",
    (Item %in% c("Cake, maize")) ~ "Corn gluten meal",
    (Item %in% c("Cake, groundnuts")) ~ "Peanut meal",
    (Item %in% c("Rape and Mustard Oils")) ~ "Rapeseed oil",
    (Item %in% c("Cake, rapeseed")) ~ "Rapeseed meal",
    (Item %in% c("Wheat")) ~ "Wheat",
    (Item %in% c("Bran, wheat")) ~ "Wheat bran",
    (Item %in% c("Rice")) ~ "Rice bran",
    (Item %in% c("Cake, sunflower")) ~ "Sunflower meal"
  )) %>%
  filter(is.na(Input) == FALSE) %>% 
  left_join(feed_fp %>% filter(Input.type == "Crop"), by = c("iso3c", "Input")) %>%
  filter(is.na(Input.type) == FALSE) %>%
  ungroup() %>%
  group_by(Input, iso3c) %>%
  summarise(Exports = sum(Exports, na.rm = TRUE)) %>%
  mutate(weighting = Exports/sum(Exports)) %>%
  select(iso3c, Input, weighting)

weighted_crop <- feed_fp %>% 
  filter(Input.type == "Crop") %>%
  left_join(weightings, by = c("iso3c", "Input")) %>%
  filter(is.na(weighting) == FALSE) %>% 
  # If weighting soy ingredient types, do here along with country weightings
  group_by(Input.type, Input, Impact.category, Allocation, Units) %>%
  # Normalize weightings to sum to 1
  mutate(reweighting = weighting/sum(weighting, na.rm = TRUE)) %>%
  summarise(Value = sum(Value * reweighting)) %>%
  # If weighting ingredient types, do here along with country weightings
  group_by(Input.type, Impact.category, Allocation, Units) %>%
  summarise(ave_stressor = mean(Value, na.rm = TRUE))

# Replace baseline soy and crops with new soy and crops
feed_fp_noland <- feed_fp_baseline %>%
  filter(Input.type != "Soy") %>%
  filter(Input.type != "Crop") %>%
  bind_rows(weighted_soy) %>%
  bind_rows(weighted_crop)

# Format for analysis
# Change names to match df
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

write.csv(feed_fp_noland, file.path(datadir, "feed_fp_scenario_7_mass.csv"), row.names = FALSE)

#_________________________________________________________________________________________________________________________________#
# Replace FMFO with fishery by-products
#_________________________________________________________________________________________________________________________________#
fmfo_trade <- read.csv(file.path(datadir, "FMFO_trade.csv"))
weightings <- fmfo_trade %>% 
  filter(!(Importer %in% c("Other Asia, nes", "Areas, nes", "Other Europe, nes", "Free Zones"))) %>%
  filter(!(Exporter %in% c("Other Asia, nes", "Areas, nes", "Other Europe, nes", "Free Zones")))%>% 
  group_by(Exporter.ISO) %>% 
  summarise(Exports = sum(Max.Weight.Live, na.rm = TRUE)) %>% 
  mutate(weighting = Exports/sum(Exports)) %>%
  select("iso3c" = "Exporter.ISO", weighting)

weighted_fishbyproduct <- feed_fp %>%
  filter(Input.type == c("Fishery by-product")) %>% 
  left_join(weightings, by = "iso3c") %>%
  filter(is.na(weighting) == FALSE) %>% 
  group_by(Input.type, Input, Impact.category, Allocation, Units) %>%
  # Normalize weightings to sum to 1
  mutate(reweighting = weighting/sum(weighting, na.rm = TRUE)) %>%
  summarise(Value = sum(Value * reweighting)) %>%
  # If weighting ingredient types, do here along with country weightings
  group_by(Impact.category, Allocation, Units) %>%
  summarise(ave_stressor = mean(Value, na.rm = TRUE))

weighted_fishbyproduct$Input.type <- "Fishery"

feed_fp_fish_bp <- feed_fp_baseline %>%
  filter(Input.type != "Fishery") %>%
  bind_rows(weighted_fishbyproduct)

# Format for analysis
# Change names to match df
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

write.csv(feed_fp_fish_bp, file.path(datadir, "feed_fp_scenario_2c_mass.csv"), row.names = FALSE)

#_________________________________________________________________________________________________________________________________#
# Replace FMFO with low impact (Alaska Pollock) fishery by-products
#_________________________________________________________________________________________________________________________________#
weightings <- fmfo_trade %>% 
  filter(!(Importer %in% c("Other Asia, nes", "Areas, nes", "Other Europe, nes", "Free Zones"))) %>%
  filter(!(Exporter %in% c("Other Asia, nes", "Areas, nes", "Other Europe, nes", "Free Zones")))%>% 
  group_by(Exporter.ISO) %>% 
  summarise(Exports = sum(Max.Weight.Live, na.rm = TRUE)) %>% 
  mutate(weighting = Exports/sum(Exports)) %>%
  select("iso3c" = "Exporter.ISO", weighting)


weighted_fishbyproduct <- feed_fp %>%
  filter(Input.type == c("Fishery by-product")) %>% 
  filter(str_detect(Input, pattern = "pollock")) %>%
  left_join(weightings, by = "iso3c") %>%
  filter(is.na(weighting) == FALSE) %>% 
  group_by(Input.type, Input, Impact.category, Allocation, Units) %>%
  # Normalize weightings to sum to 1
  mutate(reweighting = weighting/sum(weighting, na.rm = TRUE)) %>%
  summarise(Value = sum(Value * reweighting)) %>%
  # If weighting ingredient types, do here along with country weightings
  group_by(Impact.category, Allocation, Units) %>%
  summarise(ave_stressor = mean(Value, na.rm = TRUE))
weighted_fishbyproduct$Input.type <- "Fishery"

feed_fp_fish_low_impact_bp <- feed_fp_baseline %>%
  filter(Input.type != "Fishery") %>%
  bind_rows(weighted_fishbyproduct)

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

write.csv(feed_fp_fish_low_impact_bp, file.path(datadir, "feed_fp_scenario_2d_mass.csv"), row.names = FALSE)

#_________________________________________________________________________________________________________________________________#
# Low impact (Alaska Pollock) fishery by-products in FMFO
#_________________________________________________________________________________________________________________________________#
weightings <- fmfo_trade %>% 
  filter(!(Importer %in% c("Other Asia, nes", "Areas, nes", "Other Europe, nes", "Free Zones"))) %>%
  filter(!(Exporter %in% c("Other Asia, nes", "Areas, nes", "Other Europe, nes", "Free Zones")))%>% 
  group_by(Exporter.ISO) %>% 
  summarise(Exports = sum(Max.Weight.Live, na.rm = TRUE)) %>% 
  mutate(weighting = Exports/sum(Exports)) %>%
  select("iso3c" = "Exporter.ISO", weighting)

weighted_fishery <- feed_fp %>%
  filter(Input.type == c("Fishery")) %>% 
  left_join(weightings, by = "iso3c") %>%
  filter(is.na(weighting) == FALSE) %>% 
  group_by(Input.type, Input, Impact.category, Allocation, Units) %>%
  # Normalize weightings to sum to 1
  mutate(reweighting = weighting/sum(weighting, na.rm = TRUE)) %>%
  summarise(Value = sum(Value * reweighting)) %>%
  # If weighting ingredient types, do here along with country weightings
  group_by(Impact.category, Allocation, Units) %>%
  summarise(ave_stressor = mean(Value, na.rm = TRUE))

weighted_fishbyproduct <- feed_fp %>%
  filter(Input.type == c("Fishery by-product")) %>% 
  filter(str_detect(Input, pattern = "pollock")) %>%
  left_join(weightings, by = "iso3c") %>%
  filter(is.na(weighting) == FALSE) %>% 
  group_by(Input.type, Input, Impact.category, Allocation, Units) %>%
  # Normalize weightings to sum to 1
  mutate(reweighting = weighting/sum(weighting, na.rm = TRUE)) %>%
  summarise(Value = sum(Value * reweighting)) %>%
  # If weighting ingredient types, do here along with country weightings
  group_by(Impact.category, Allocation, Units) %>%
  summarise(ave_stressor = mean(Value, na.rm = TRUE))
weighted_fishbyproduct$Input.type <- "Fishery"

weighted_fish <- weighted_fishery %>% 
  left_join(weighted_fishbyproduct, by = c("Impact.category", "Allocation", "Units")) %>%
  mutate(ave_stressor = 0.675*ave_stressor.x + 0.325*ave_stressor.y) %>%
  select(-ave_stressor.x, -ave_stressor.y)
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

write.csv(feed_fp_fmfo_low_impact_bp, file.path(datadir, "feed_fp_scenario_6_mass.csv"), row.names = FALSE)
