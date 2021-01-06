# Clean compiled LCA data needed for all models

rm(list = ls())

# Libraries for processing and analyses
library(tidyverse)
library(rstan)
library(data.table)
library(countrycode) # part of clean.lca
library(bayesplot) # for mcmc_areas_ridges
library(shinystan)
library(brms)
library(tidybayes)

# Mac
datadir <- "/Volumes/jgephart/BFA Environment 2/Data"
outdir <- "/Volumes/jgephart/BFA Environment 2/Outputs"
# Windows
# datadir <- "K:/BFA Environment 2/Data"
# outdir <- "K:BFA Environment 2/Outputs"

lca_dat <- read.csv(file.path(datadir, "LCA_compiled_20201222.csv"), fileEncoding="UTF-8-BOM") #fileEncoding needed when reading in file from windows computer (suppresses BOM hidden characters)
source("Functions.R")

# Clean LCA data
lca_dat_clean <- clean.lca(LCA_data = lca_dat)

# Rebuild FAO fish production from zip file
fishstat_dat <- rebuild_fish("/Volumes/jgephart/FishStatR/Data/Production-Global/ZippedFiles/GlobalProduction_2019.1.0.zip")

# Classify species into taxa groupings
# Use ISSCAAP grouping (in FAO data) to help with classification
lca_dat_clean_groups <- add_taxa_group(lca_dat_clean, fishstat_dat) %>%
  group_by(clean_sci_name, taxa_group_name, taxa) %>%
  mutate(n_in_sci = n()) %>%
  ungroup()

# Check after add_taxa_group, there should be no taxa == "unassigned"
sort(unique(lca_dat_clean_groups$taxa))

# Output clean data with groups (later will read this back in to join with model predictions)
write.csv(lca_dat_clean_groups, file.path(datadir, "lca_clean_with_groups.csv"), row.names = FALSE)

# Output taxa groupings and sample sizes:
n_and_composition <- lca_dat_clean_groups %>% select(taxa_group_name, Source, clean_sci_name) %>% distinct() %>% group_by(taxa_group_name, clean_sci_name) %>% tally() # number of studies/sources per taxa
write.csv(n_and_composition, file.path(outdir, "taxa_group_n_and_composition.csv"))

#________________________________________________________________________________________________________________________________________________________________#
# Calculate production weightings for each taxa group
#________________________________________________________________________________________________________________________________________________________________#
prod_weightings <- fishstat_dat %>% 
  filter(year > 2012) %>%
  filter(unit == "t") %>%
  filter(source_name_en %in% c("Aquaculture production (marine)", "Aquaculture production (brackishwater)", "Aquaculture production (freshwater)")) %>%
  group_by(isscaap_group, species_scientific_name) %>%
  summarise(species_prod = sum(quantity, na.rm = TRUE))

# Check for species names in lca data that does not match fao data
unique(lca_dat_clean_groups$clean_sci_name[!(lca_dat_clean_groups$clean_sci_name %in% prod_weightings$species_scientific_name)])

# Correct species names to merge with fao production data

lca_dat_clean_groups <- lca_dat_clean_groups %>%
  mutate(fao_species = case_when(
    (clean_sci_name == "Litopenaeus vannamei") ~ "Penaeus vannamei",
    (clean_sci_name == "Scophthalmidae") ~ "",
    (clean_sci_name == "Cynoscion spp") ~ "Sciaenidae", # Replace with the family, but this may underestimate production
    (clean_sci_name == "Laminaria digitata") ~ "Macrocystis pyrifera", # Or use Saccharina latissima to weight evenly with other brown seaweeds
    (clean_sci_name == "Gracilaria chilensis") ~ "Gracilaria spp",  # Use as representative for red algae in general
    (clean_sci_name == "Mixed Hypophthalmichthys molitrix and H. nobilis") ~ "Hypophthalmichthys molitrix",
    (clean_sci_name == "Ctenopharyngodon idella") ~ "Ctenopharyngodon idellus",
    (clean_sci_name == "Anoplopoma fimbria") ~ "",
    (clean_sci_name == "Macrobrachium amazonicum") ~ "Macrobrachium spp"
  ))
lca_dat_clean_groups$fao_species[is.na(lca_dat_clean_groups$fao_species)] <- lca_dat_clean_groups$clean_sci_name[is.na(lca_dat_clean_groups$fao_species)]

prod_weightings <-  prod_weightings %>% 
  right_join(lca_dat_clean_groups, by = c("species_scientific_name" = "fao_species", "isscaap_group")) %>%
  select(isscaap_group, taxa_group_name, taxa, clean_sci_name, species_prod) %>%
  distinct()

# Add small amount to na's 
prod_weightings$species_prod[is.na(prod_weightings$species_prod)] <- 10

# Calculate proportions
prod_weightings <-  prod_weightings %>% 
  group_by(taxa_group_name, taxa) %>%
  # For taxa levels higher than species, use production volumes equal to the average in the taxa group
  mutate(species_prod_reweight = ifelse(str_detect(clean_sci_name, pattern = "spp")|length(str_split(string = clean_sci_name, pattern = "\\s")[[1]]) == 1, 
                                        mean(species_prod, na.rm = TRUE), species_prod),
    prod_weighting = species_prod_reweight/sum(species_prod_reweight, na.rm = TRUE))

# Check that species sum to 1
prod_weightings %>% group_by(taxa_group_name, taxa) %>% summarise(total_weight = sum(prod_weighting))

write.csv(prod_weightings, file.path(datadir, "20210106_aqua_prod_weightings.csv"), row.names = FALSE)
#________________________________________________________________________________________________________________________________________________________________#
# Calculate the weighted averages for the feed components # FIX IT: Need to add in groupings (not finished yet)
#________________________________________________________________________________________________________________________________________________________________#
feed_fp <- read.csv(file.path(datadir, "Feed_impact_factors_20201203.csv"))
feed_fp$iso3c <- countrycode(feed_fp$Country.Region, origin = "country.name", destination = "iso3c")
feed_fp$iso3c[is.na(feed_fp$iso3c)] <- "Other"

faostat <- read.csv(file.path(datadir, "FAOSTAT_data_12-9-2020.csv"))
faostat$iso3c <- countrycode(faostat$Area, origin = "country.name", destination = "iso3c")

# Calculate soy weightings
# Since the specific soy products don't map onto the FAO items, use all soy exports for all soy products
weightings <-  faostat %>% 
  filter(Unit == "tonnes") %>% 
  filter(Item %in% c("Soybeans", "Cake, soybeans", "Oil, soybean")) %>%
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
  group_by(Input.type, Impact.category, Allocation, Units) %>%
  # If weighting soy ingredient types, do here along with country weightings
  summarise(ave_stressor = sum(Value * weighting))
  
# Calculate crop weightings
weightings <- faostat %>% 
  filter(Unit == "tonnes") %>% 
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
  group_by(Input.type, Impact.category, Allocation, Units) %>%
  # Currently filtering out NAs weightings (no trade matched for corn gluten meal)
  filter(is.na(weighting) == FALSE) %>% 
  # If weighting soy ingredient types, do here along with country weightings
  summarise(ave_stressor = sum(Value * weighting))

# Weightings for animal by products
# Using weightings for all pigmeat and chicken exports
weightings <-  faostat %>% 
  filter(Unit == "tonnes") %>% 
  filter(Item %in% c("Pigmeat")) %>%
  group_by(iso3c) %>%
  summarise(Exports = sum(Value, na.rm = TRUE)) %>%
  filter(Exports > 0) %>% 
  left_join(feed_fp %>% filter(Input == "Pork blood meal"), by = c("iso3c")) %>%
  filter(is.na(Input.type) == FALSE) %>%
  group_by(iso3c) %>%
  summarise(Exports = sum(Exports, na.rm = TRUE)) %>%
  mutate(weighting = Exports/sum(Exports)) %>%
  select(iso3c, weighting)

weighted_pig <- feed_fp %>% 
  filter(Input == "Pork blood meal") %>%
  left_join(weightings, by = "iso3c") %>%
  group_by(Input.type, Impact.category, Allocation, Units) %>%
  # If weighting soy ingredient types, do here along with country weightings
  summarise(ave_stressor = sum(Value * weighting))

weightings <-  faostat %>% 
  filter(Unit == "tonnes") %>% 
  filter(Item %in% c("Poultry Meat")) %>%
  group_by(iso3c) %>%
  summarise(Exports = sum(Value, na.rm = TRUE)) %>%
  filter(Exports > 0) %>% 
  left_join(feed_fp %>% 
              filter(Input %in% c("Chicken by-product meal", "Chicken by-product oil")), by = c("iso3c")) %>%
  filter(is.na(Input.type) == FALSE) %>%
  group_by(iso3c) %>%
  summarise(Exports = sum(Exports, na.rm = TRUE)) %>%
  mutate(weighting = Exports/sum(Exports)) %>%
  select(iso3c, weighting)

weighted_chicken <- feed_fp %>% 
  filter(Input %in% c("Chicken by-product meal", "Chicken by-product oil")) %>%
  left_join(weightings, by = "iso3c") %>%
  group_by(Input.type, Impact.category, Allocation, Units) %>%
  # If weighting soy ingredient types, do here along with country weightings
  summarise(ave_stressor = sum(Value * weighting, na.rm = TRUE))

weighted_livestock <- weighted_pig %>%
  bind_rows(weighted_chicken) %>%
  group_by(Input.type, Impact.category, Allocation, Units) %>%
  summarise(ave_stressor = mean(ave_stressor))

# Fishery products (currently unweighted since fish products are not in FAOSTAT trade)
weighted_fish <- feed_fp %>%
  filter(Input.type == "Fishery") %>% 
  group_by(Input.type, Impact.category, Allocation, Units) %>%
  # If weighting soy ingredient types, do here along with country weightings
  summarise(ave_stressor = mean(Value))

# Combine data frames
weighted_fp <- rbind(weighted_soy, weighted_crop)
weighted_fp <- rbind(weighted_fp, weighted_livestock)
weighted_fp <- rbind(weighted_fp, weighted_fish)

write.csv(weighted_fp, file.path(datadir, "20210106_weighted_feed_fp.csv"), row.names = FALSE)

