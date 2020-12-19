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

lca_dat <- read.csv(file.path(datadir, "LCA_compiled_20201214.csv"), fileEncoding="UTF-8-BOM") #fileEncoding needed when reading in file from windows computer (suppresses BOM hidden characters)
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
  ungroup() #%>%
  #filter(n_in_sci > 2)
# doign this removes whole taxa groups: crabs, aquatic plants get filtered out

# Check after add_taxa_group, there should be no taxa == "unassigned"
sort(unique(lca_dat_clean_groups$taxa))

# Output clean data with groups (later will read this back in to join with model predictions)
write.csv(lca_dat_clean_groups, file.path(datadir, "lca_clean_with_groups.csv"), row.names = FALSE)

# Output taxa groupings and sample sizes:
# data.frame(table(lca_dat_clean_groups$taxa_group_name))
# data.frame(table(lca_dat_clean_groups$taxa)) # abbreviated version of taxa_group_name for writing models
# lca_dat_clean_groups %>% select(taxa_group_name, Source) %>% distinct() %>% group_by(taxa_group_name) %>% tally() # number of studies/sources per taxa
# lca_dat_clean_groups %>% select(taxa_group_name, clean_sci_name) %>% group_by(taxa_group_name, clean_sci_name) %>% mutate(n_obs = n()) %>% unique() %>% arrange(taxa_group_name) %>% print(n=50)
# write.csv(data.frame(table(lca_dat_clean_groups$taxa_group_name)), file.path(outdir, "taxa_group_sample_size.csv"))
# write.csv(lca_dat_clean_groups %>% select(taxa_group_name, clean_sci_name) %>% unique() %>% arrange(taxa_group_name), file.path(outdir, "taxa_group_composition.csv"))

# Add fish N and P content
#calculate C,N,P content of fish (based on energy and fat content ) and compute discharge based on Czamanski et al 2011, Marine Biology,
#Carbon, nitrogen and phosphorus elemental stoichiometry in aquacultured and wild-caught fish and consequences
#for pelagic nutrient dynamics
#And also a simplified NPZ method approach

#N and P content of fish in % of DM
#read Zach's data and make calculations using built functions
fishNutrition <- read.csv(file.path(datadir, "AFCD_live.csv"), stringsAsFactors = FALSE)
fishNutrition <-fishNutrition %>%filter(Processing=='r') %>%  #filter out all non raw observations (e.g. dried, cooked)
  slice(c(1:1928))    #cut rows without names of species 
fishNutrition1<-fishNutrition %>% mutate(N_C=fishN_via_C(Water,Energy.total.metabolizable.calculated.from.the.energy.producing.food.components.original.as.from.source.kcal)
                                         ,P_C=fishP_via_C(Water,Energy.total.metabolizable.calculated.from.the.energy.producing.food.components.original.as.from.source.kcal)
                                         ,N_fat=fishN_viaFat(Water,Fat.total)
                                         ,P_fat=fishP_viaFat(Water,Fat.total)
                                         ,N_built_in=Nitrogen.total*100/(100-Water),P_built_in=Phosphorus/1000*100/(100-Water)) #conversion to units of percentages in DM. Nitrogen.total is in units of g/100 edible gram; Phosphorous in the dataset is in mg/100 edible gram

fishNutrition2<-fishNutrition1 %>% group_by(species) %>% summarise(P_byC=mean(P_C,na.rm = TRUE),N_byC=mean(N_C,na.rm = TRUE),
                                                                   P_fat_1=mean(P_fat,na.rm = TRUE),N_fat_1=mean(N_fat,na.rm = TRUE),
                                                                   P_built_in1=mean(P_built_in,na.rm = TRUE),N_built_in1=mean(N_built_in,na.rm = TRUE)) %>%                                 
  rowwise() %>%                                                    
  mutate(N_avg=mean(c(N_byC,N_fat_1,N_built_in1),na.rm = TRUE))%>%
  mutate(P_avg=mean(c(P_byC,P_fat_1,P_built_in1),na.rm = TRUE))%>%
  select(c("species","N_avg","P_avg"))

# FIX IT: Currently there are 25 unmatched species. Emailed Alon Dec 17 to see if he can fix
lca_dat_clean_groups_tmp <- lca_dat_clean_groups %>% 
  left_join(fishNutrition2, by = c("Scientific.Name" = "species"))

#________________________________________________________________________________________________________________________________________________________________#
# Calculate production weightings for each taxa group
#________________________________________________________________________________________________________________________________________________________________#
prod_weightings <- fishstat_dat %>% 
  filter(year > 2012) %>%
  filter(unit == "t") %>%
  filter(source_name_en %in% c("Aquaculture production (marine)", "Aquaculture production (brackishwater)", "Aquaculture production (freshwater)")) %>%
  group_by(isscaap_group, species_scientific_name) %>%
  summarise(total = sum(quantity, na.rm = TRUE))

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
  select(isscaap_group, taxa_group_name, taxa, clean_sci_name, total) %>%
  distinct()

# Add small amount to na's 
prod_weightings$total[is.na(prod_weightings$total)] <- 10

prod_weightings <-  prod_weightings %>% 
  group_by(taxa_group_name, taxa) %>%
  mutate(weighting = total/sum(total, na.rm = TRUE))

write.csv(prod_weightings, file.path(datadir, "aqua_prod_weightings.csv"), row.names = FALSE)
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

write.csv(weighted_fp, file.path(datadir, "20201217_weighted_feed_fp.csv"), row.names = FALSE)

