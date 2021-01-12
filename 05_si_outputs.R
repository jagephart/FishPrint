# Get additional plots/summaries for SI

rm(list=ls())
library(tidyverse)
source("Functions.R") # for rebuild_fish

datadir <- "/Volumes/jgephart/BFA Environment 2/Data"
outdir <- "/Volumes/jgephart/BFA Environment 2/Outputs"




# Rebuild FAO fish production from zip file
fishstat_dat <- rebuild_fish("/Volumes/jgephart/FishStatR/Data/Production-Global/ZippedFiles/GlobalProduction_2019.1.0.zip")

# NOTES on matching Taxa group to ISSCAAP group(s)
# Silver/bighead carp - remove from Carps, barbels and other cyprinids
# Bivalves - add Mussels, Oysters, Scallops, pectens, and Clams, cockles, arkshells
# Catfish
# Milkfish
# Misc carps - all of Carps barbels and other cyprinids except big/silverhead carp
# Misc diad fishes
# Misc marine fishes
# Plants
# Salmon
# Shrimp
# Tilapia
# Trout

# 2012 - 2019
prod_by_taxa <- fishstat_dat %>% 
  filter(year > 2012) %>%
  filter(unit == "t") %>%
  filter(source_name_en %in% c("Aquaculture production (marine)", "Aquaculture production (brackishwater)", "Aquaculture production (freshwater)")) %>%
  group_by(isscaap_group, order, family, species_major_group, species_scientific_name) %>%
  summarise(species_prod = sum(quantity, na.rm = TRUE)) %>%
  ungroup() %>%
  # Identify rows that are part of this study
  mutate(taxa = case_when(isscaap_group %in% c("Oysters", "Mussels") ~ "bivalves",
                          order == "SILURIFORMES" ~ "catfish", 
                          isscaap_group == "Carps, barbels and other cyprinids" & str_detect(species_scientific_name, "Hypophthalmichthys") ~ "hypoph_carp",
                          species_scientific_name == "Chanos chanos" ~ "milkfish",
                          isscaap_group == "Carps, barbels and other cyprinids" & str_detect(species_scientific_name, "Hypophthalmichthys")==FALSE ~ "misc_carp",
                          isscaap_group == "Miscellaneous diadromous fishes" ~ "misc_diad",
                          isscaap_group %in% c("Miscellaneous coastal fishes", "Miscellaneous demersal fishes", "Miscellaneous pelagic fishes") ~ "misc_marine",
                          species_major_group == "PLANTAE AQUATICAE" ~ "plants",
                          family == "SALMONIDAE" ~ "salmon",
                          isscaap_group == "Shrimps, prawns" ~ "shrimp",
                          isscaap_group == "Tilapias and other cichlids" ~ "tilapia",
                          str_detect(species_scientific_name, "Oncorhynchus|Salmo|Salvelinus") ~ "trout",
                          TRUE ~ "other_taxa")) %>%
  group_by(taxa) %>%
  summarise(taxa_prod = sum(species_prod, na.rm = TRUE)) %>%
  ungroup()

# Check total sums are equal
sum(prod_by_taxa$taxa_prod)
fishstat_dat %>% 
  filter(year > 2012) %>%
  filter(unit == "t") %>%
  filter(source_name_en %in% c("Aquaculture production (marine)", "Aquaculture production (brackishwater)", "Aquaculture production (freshwater)")) %>%
  pull(quantity) %>%
  sum()
  

# How much of total production do our taxa groupings represent?
prod_by_taxa %>%
  filter(taxa != "other_taxa") %>%
  pull(taxa_prod) %>% sum() / sum(prod_by_taxa$taxa_prod)

# Proportion of global aquaculture production represented by analysis: 0.8245676
  
  
  
  