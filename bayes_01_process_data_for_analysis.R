# Author: Kelvin Gorospe
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
source("Functions.R")

# Full raw dataset used in published analysis
# lca_dat <- read.csv(file.path(datadir, "LCI_compiled_FINAL.csv")) 
# Clean and replicate LCA data 
# lca_dat_clean <- clean.lca(LCA_data = lca_dat, replicate_dat = TRUE)

# Read in cleaned, aggregated, and replicated dataset:
# DATA DOCUMENTATION: 80% of the data points from the original analysis are identical in the aggregated file, the remaining 20% were aggregated as they come from unpublished studies
# ORIGINAL RAW FILE CAN BE REQUESTED: CONTACT THE CORRESPONDING AUTHOR FOR INFORMATION
lca_dat_clean <- read.csv(file.path(datadir, "LCI_compiled_for_SI.csv")) 

# Rebuild FAO fish production from zip file
#fishstat_dat <- rebuild_fish("/Volumes/jgephart/FishStatR/Data/Production-Global/ZippedFiles/GlobalProduction_2019.1.0.zip")
fishstat_dat <- rebuild_fish("/Volumes/jgephart/FishStatR/Data/Production-Global/ZippedFiles/GlobalProduction_2020.1.0.zip")
# DATA DOCUMENTATION: Zip file can be downloaded from FAO FishStat: http://www.fao.org/fishery/static/Data/GlobalProduction_2020.1.0.zip
# File can also be provided by corresponding author upon request

# Classify species into taxa groupings
# Use ISSCAAP grouping (in FAO data) to help with classification
lca_dat_clean_groups <- add_taxa_group(lca_dat_clean, fishstat_dat) %>%
  group_by(clean_sci_name, taxa_group_name, taxa) %>%
  mutate(n_in_sci = n()) %>%
  ungroup()

# Check after add_taxa_group, there should be no taxa == "unassigned"
sort(unique(lca_dat_clean_groups$taxa))

# Output clean data with groups (later will read this back in to join with model predictions)
write.csv(lca_dat_clean_groups, file.path(outdir, "lca_clean_with_groups.csv"), row.names = FALSE)

# Output taxa groupings and sample sizes:
n_and_composition <- lca_dat_clean_groups %>% 
  select(taxa_group_name, Sample_size_n_farms, Country, iso3c, Source, clean_sci_name) %>% 
  distinct() %>% 
  group_by(taxa_group_name, clean_sci_name, Country, iso3c) %>% 
  summarize(n_farms = sum(Sample_size_n_farms),
            n_studies = n()) %>% # each row is a distinct study
  ungroup() # number of studies/sources and farms per taxa
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

write.csv(prod_weightings, file.path(outdir, "aqua_prod_weightings.csv"), row.names = FALSE)

#________________________________________________________________________________________________________________________________________________________________#
# Calculate the weighted averages for the feed components
#________________________________________________________________________________________________________________________________________________________________#
feed_fp <- read.csv(file.path(datadir, "feed_impact_factors.csv"))
feed_fp$iso3c <- countrycode(feed_fp$Country.Region, origin = "country.name", destination = "iso3c")
feed_fp$iso3c[is.na(feed_fp$iso3c)] <- "Other"

feed_fp <- feed_fp %>% 
  mutate(iso3c = ifelse(Country.Region == "Europe" & Input %in% c("Chicken by-product meal", "Chicken by-product oil"),
                        "EUR", iso3c)) 

# DATA DOCUMENTATION: These data come from FAO and can be downloaded from FAOSTAT
# http://www.fao.org/faostat/en/#data/TP
# Query terms for downloading: SELECT All countries, Export quantity, All items
# File can also be provided by corresponding author upon request
faostat <- read.csv(file.path(datadir, "FAOSTAT_data_12-9-2020.csv"))
faostat$iso3c <- countrycode(faostat$Area, origin = "country.name", destination = "iso3c")

# Calculate soy weights
weighted_soy <- calc_soy_weights(faostat, feed_fp)

# Calculate crop weights
weighted_crop <- calc_crop_weights(faostat, feed_fp)
  
# Weights for animal by products - only using chicken
weighted_livestock <- calc_chicken_weights(faostat, feed_fp)

# Fishery products
# # UN Comtrade data for 2015
# fmfo_trade <- read.csv(file.path(datadir, "FMFO_trade.csv"))
# Use landings weightings instead of trade weightings
fmfo_prod <- read.csv(file.path(datadir, "fish_weightings.csv"))

# Calculate fishery weights
weighted_fishery <- calc_fishery_weights(fmfo_prod, feed_fp)

# Calculate fishery byproduct weights 
weighted_fishbyproduct <- calc_byproduct_weights(fmfo_prod, feed_fp)

# Combine fishery and byproduct weights
weighted_fish <- combine_fish_weights(weighted_fishery, weighted_fishbyproduct)
weighted_fish$Input.type <- "Fishery"

# Combine data frames
weighted_fp <- rbind(weighted_soy, weighted_crop)
weighted_fp <- rbind(weighted_fp, weighted_livestock)
weighted_fp <- rbind(weighted_fp, weighted_fish)

write.csv(weighted_fp, file.path(outdir, "weighted_feed_fp.csv"), row.names = FALSE)

