# Clean compiled LCA data needed for all models

rm(list = ls())
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

lca_dat <- read.csv(file.path(datadir, "LCA_compiled_20201109.csv"), fileEncoding="UTF-8-BOM") #fileEncoding needed when reading in file from windows computer (suppresses BOM hidden characters)
source("Functions.R")

# Clean LCA data
lca_dat_clean <- clean.lca(LCA_data = lca_dat)

# Rebuild FAO fish production from zip file
fishstat_dat <- rebuild_fish("/Volumes/jgephart/FishStatR/Data/Production-Global/ZippedFiles/GlobalProduction_2019.1.0.zip")

# Classify species into taxa groupings
# Use ISSCAAP grouping (in FAO data) to help with classification
lca_dat_clean_groups <- add_taxa_group(lca_dat_clean, fishstat_dat)

# Output clean data with groups (later will read this back in to join with model predictions)
write.csv(lca_dat_clean_groups, file.path(datadir, "lca_clean_with_groups.csv"), row.names = FALSE)

# Output taxa groupings and sample sizes:
#data.frame(table(lca_dat_clean_groups$taxa_group_name))
#data.frame(table(lca_dat_clean_groups$taxa)) # abbreviated version of taxa_group_name for writing models
#lca_dat_clean_groups %>% select(taxa_group_name, clean_sci_name) %>% unique() %>% arrange(taxa_group_name)
#write.csv(data.frame(table(lca_dat_clean_groups$taxa_group_name)), file.path(outdir, "taxa_group_sample_size.csv"))
#write.csv(lca_dat_clean_groups %>% select(taxa_group_name, clean_sci_name) %>% unique() %>% arrange(taxa_group_name), file.path(outdir, "taxa_group_composition.csv"))

# Calculate production weightings for each taxa group
prod_weightings <- fishstat_dat %>% 
  filter(year > 2012) %>%
  filter(unit == "t") %>%
  filter(source_name_en %in% c("Aquaculture production (marine)", "Aquaculture production (brackishwater)", "Aquaculture production (freshwater)")) %>%
  group_by(isscaap_group, species_scientific_name) %>%
  summarise(total = sum(quantity, na.rm = TRUE))

# FIX IT: the species merge isn't working. Need to do some species name cleaning 
prod_weightings <-  prod_weightings %>% 
  right_join(lca_dat_clean_groups, by = c("species_scientific_name" = "clean_sci_name", "isscaap_group")) 

# Calculate the weighted averages for the feed components # FIX IT: Need to add in groupings (not finished yet)
feed_fp <- read.csv(file.path(datadir, "Feed_impact_factors_20201203.csv"))

faostat <- read.csv(file.path(datadir, "FAOSTAT_data_12-8-2020.csv"))
faostat_summary <- faostat %>% 
  filter(Unit == "tonnes") %>% 
  group_by(Area, Item) %>%
  summarise(Value = sum(Value, na.rm = TRUE)) %>%
  filter(Value > 0) %>% 
  mutate(Item_group = case_when(
    (Item %in% c("Soybeans", "Cake, soybeans", )) ~ "Soy",
    (Item %in% c("Flour, cassava", "Starch, cassava")) ~ "Cassava",
    (Item %in% c()) ~ "Corn gluten meal",
    (Item %in% c()) ~ "Maize",
    (Item %in% c()) ~ "Peanut meal",
    (Item %in% c()) ~ "Rapeseed meal",
    (Item %in% c("Oil, rapeseed")) ~ "Rapeseed oil",
    (Item %in% c()) ~ "Wheat",
    (Item %in% c()) ~ "Rice bran",
    (Item %in% c()) ~ "Sunflower meal",
    (Item %in% c()) ~ "Wheat bran",
    (Item %in% c()) ~ "Soy protein concentrate",
    (Item %in% c()) ~ "Soy protein isolate",
    (Item %in% c()) ~ "Soybean meal"
  ))


 "Chicken by-product meal"              "Chicken by-product oil"              
[25] "Pork blood meal"                                                     
[28]    


