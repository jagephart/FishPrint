# TITLE: explore_edible_portion.R
# AUTHOR: Jessica
# DATE: 12-Apr-21

# Load packages
library(tidyverse)
library(ggthemes)

# Load data
# Peter's compiled data formatted for R
df_all <- read.csv("Data/edible_muscle_raw.csv")

# Matched data
df_matched <- read.csv("Data/species_edible_portion.csv")

# Exploratory plots
## Edible fraction versus fillet yield
ggplot(df_all, aes(x = max_edible_fraction_from_live_weight_percent,
                   y = mean_fillet_yield_from_live_weight_percent, 
                   color = Group_coarse,
                   alpha = 0.7)) +
  labs(x = "Max edible fraction from live weight (%)",
       y = "Mean fillet yield from live weight (%)") +
  geom_point() + 
  theme_minimal()

ggplot(df_all %>% filter(Group_coarse == "finfish"), 
       aes(x = max_edible_fraction_from_live_weight_percent,
           y = mean_fillet_yield_from_live_weight_percent, 
           alpha = 0.7)) +
  geom_point() +
  lims(x = c(30, 70)) + 
  labs(x = "Max edible fraction from live weight (%)",
       y = "Mean fillet yield from live weight (%)") +
  theme_minimal() +
  theme(legend.position = "none") 

df_finfish <- df_all %>%
  filter(Group_coarse == "finfish")

summary(lm(df_finfish$max_edible_fraction_from_live_weight_percent ~
             df_finfish$mean_fillet_yield_from_live_weight_percent))

## Edible fraction versus protein yield
ggplot(df_all, aes(x = max_edible_fraction_from_live_weight_percent,
                   y = mean_protein_fraction_of_edible_portion_percent, 
                   color = Group_coarse,
                   alpha = 0.7)) +
  geom_point() + 
  theme_minimal()

ggplot(df_all %>% filter(Group_coarse == "finfish"), 
       aes(x = max_edible_fraction_from_live_weight_percent,
           y = mean_protein_fraction_of_edible_portion_percent, 
           alpha = 0.7)) +
  geom_point() +
  lims(x = c(30, 70)) + 
  labs(x = "Max edible fraction from live weight (%)",
       y = "Mean protein fraction from live weight (%)") +
  theme_minimal() +
  theme(legend.position = "none") 

summary(lm(df_finfish$max_edible_fraction_from_live_weight_percent ~
             df_finfish$mean_protein_fraction_of_edible_portion_percent))

# Summary of group values from matched data
## Number of matches per species
species_summary <- df_matched %>% 
  select(fishprint_species, TaxonName) %>%
  filter(!is.na(TaxonName)) %>%
  filter(TaxonName != "") %>%
  group_by(fishprint_species) %>%
  tally()


## Number of matches per taxa group
taxa_summary <- df_matched %>% 
  select(fishprint_taxa, TaxonName) %>%
  filter(!is.na(TaxonName)) %>%
  filter(TaxonName != "") %>%
  group_by(fishprint_taxa) %>%
  tally()

## Mean edible by group
means_summary <- df_matched %>%
  group_by(fishprint_taxa) %>%
  summarise(mean_maxedible = mean(max_edible_fraction_from_live_weight_percent, na.rm = TRUE),
            mean_fillet = mean(mean_fillet_yield_from_live_weight_percent, na.rm = TRUE),
            mean_protein = mean(mean_protein_fraction_of_edible_portion_percent, na.rm = TRUE))

