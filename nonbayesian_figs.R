# Non-Bayesian estimate figures
library(tidyverse)
library(ggplot2)
library(ggthemes)

source("nonbayesian_estimates.R")

#_______________________________________________________________________________________________________________________#
# Plot feed-associated stressors
#_______________________________________________________________________________________________________________________#

ggplot(df_feed_taxa_unweighted %>% filter(Impact.category == "Global warming potential"), 
       aes(x = weighted_stressor, y = taxa, fill = Allocation)) +
  geom_bar(position="dodge", stat="identity") + 
  labs(x = "GHG kg CO2-eq/t", y = "", title = "GHG") +
  theme_clean()

ggplot(df_feed_taxa_unweighted %>% filter(Impact.category == "Marine eutrophication"), 
       aes(x = weighted_stressor, y = taxa, fill = Allocation)) +
  geom_bar(position="dodge", stat="identity") + 
  labs(x = "kg N-eq/t", y = "", title = "Nitrogen") +
  theme_clean()

ggplot(df_feed_taxa_unweighted %>% filter(Impact.category == "Freshwater eutrophication"), 
       aes(x = weighted_stressor, y = taxa, fill = Allocation)) +
  geom_bar(position="dodge", stat="identity") + 
  labs(x = "kg P-eq/t", y = "", title = "Phosphorus") +
  theme_clean()

ggplot(df_feed_taxa_unweighted %>% filter(Impact.category == "Land use"), 
       aes(x = weighted_stressor, y = taxa, fill = Allocation)) +
  geom_bar(position="dodge", stat="identity") + 
  labs(x = "m2a/t", y = "", title = "Land") +
  theme_clean()

ggplot(df_feed_taxa_unweighted %>% filter(Impact.category == "Water consumption"), 
       aes(x = weighted_stressor, y = taxa, fill = Allocation)) +
  geom_bar(position="dodge", stat="identity") + 
  labs(x = "m3/t", y = "", title = "Water") +
  theme_clean()

#_______________________________________________________________________________________________________________________#
# Plot on farm GHG
#_______________________________________________________________________________________________________________________#
# Plot unweighted taxa means
ggplot(df_onfarm_ghg %>% group_by(taxa) %>% summarise(mean_onfarm_ghg = mean(onfarm_ghg_kgCO2)), 
       aes(y = taxa, x = mean_onfarm_ghg)) +
  geom_bar(stat = "identity") + 
  labs(x = "On farm GHG (kg CO2-eq)", y = "") +
  theme_clean()


#_______________________________________________________________________________________________________________________#
# Plot on farm N and P
#_______________________________________________________________________________________________________________________#
# Plot unweighted taxa means
ggplot(df_onfarm_NP %>% group_by(taxa) %>% summarise(mean_onfarm_N = mean(N_emissions_kg_per_t)), 
       aes(y = taxa, x = mean_onfarm_N)) +
  geom_bar(stat = "identity") + 
  labs(x = "On farm N (kg N per t)", y = "") +
  theme_clean()

ggplot(df_onfarm_NP %>% group_by(taxa) %>% summarise(mean_onfarm_P = mean(P_emissions_kg_per_t)), 
       aes(y = taxa, x = mean_onfarm_P)) +
  geom_bar(stat = "identity") + 
  labs(x = "On farm P (kg P per t)", y = "") +
  theme_clean()

#_______________________________________________________________________________________________________________________#
# Plot on farm land
#_______________________________________________________________________________________________________________________#
# Plot unweighted taxa means
ggplot(df_onfarm_land %>% group_by(taxa) %>% summarise(mean_onfarm_land = mean(yield)), 
       aes(y = taxa, x = mean_onfarm_land)) +
  geom_bar(stat = "identity") + 
  labs(x = "On farm land (m2 per t)", y = "") +
  theme_clean()
