# Box 1 figure

library(tidyverse)
library(ggplot2)
library(ggthemes)
library(plotrix)


# Load data
mm_risk <- read.csv("Data/marine_mammal_risk.csv")
fui <- read.csv("Data/marine_mammal_risk_fui.csv")

# Tidy data and join FUI data
mm_risk <- mm_risk %>%
  pivot_longer(high:low, names_to = "risk", values_to = "n.species") %>%
  filter(!(is.na(n.species))) %>%
  select("mm_species" = species, "species" = fui, gear, Region, risk, n.species) %>%
  left_join(fui, by = c("species", "gear")) %>%
  filter(!is.na(ghg))

mm_risk$risk <- ordered(mm_risk$risk, levels = c("low", "medium", "high"))

# Gillnets and Entangling Nets
ggplot(mm_risk, aes(x = n.species, y = ghg, colour = factor(risk))) +
  geom_boxplot() +
  facet_grid(rows = vars(gear), cols = vars(mm_species))

# Plot average ghg with se bars
mm_risk_aveghg <- mm_risk %>%
  group_by(mm_species, gear, risk, n.species) %>%
  summarise(ghg.ave = mean(ghg), ghg.se = std.error(ghg, na.rm = TRUE)) %>%
  mutate(species_short = case_when(
    (mm_species == "shrimps and prawns in northern Australia") ~ "N. Aus. Shrimps & Prawns", 
    (mm_species == "Demersal species in NE Atlantic (Irish EEZ)") ~ "NE Atl. Demersals", 
    (mm_species == "Pelagic species in NE Atlantic (Irish EEZ)") ~ "NE Atl. Pelagics", 
    (mm_species == "Crustaceans in NE Atlantic (Irish EEZ)") ~ "NE Atl. Crustaceans", 
    (mm_species == "small-scale fishing of lobsters") ~ "SSF Lobsters", 
    (mm_species == "Tunas in the Pacific") ~ "Pac. Tunas"
  )) %>%
  mutate(fishery = paste(species_short, ", ", gear, sep = ""))

mm_risk_aveghg$fishery <- ordered(mm_risk_aveghg$fishery, levels = 
                                    c("NE Atl. Crustaceans, Traps and Lift Nets" , 
                                      "SSF Lobsters, Traps and Lift Nets", 
                                      "NE Atl. Demersals, Bottom Trawls",
                                      "N. Aus. Shrimps & Prawns, Bottom Trawls",
                                      "NE Atl. Pelagics, Midwater Trawls",
                                      "NE Atl. Demersals, Gillnets and Entangling Nets",
                                      "Pac. Tunas, Surrounding Nets"))

ggplot(mm_risk_aveghg, aes(x = ghg.ave, y = fishery, colour = factor(risk))) +
  geom_point(aes(size = n.species)) +
  geom_errorbar(aes(xmin=ghg.ave-(1.96*ghg.se), xmax=ghg.ave+(1.96*ghg.se)), width=.1) 

ggplot(mm_risk_aveghg, aes(x = ghg.ave, y = fishery, colour = factor(risk))) +
  geom_point(aes(size = n.species)) +
  geom_errorbar(aes(xmin=ghg.ave-(1.96*ghg.se), xmax=ghg.ave+(1.96*ghg.se)), width=.1) +
  facet_wrap(~risk) +
  labs(y = "", x = "kg CO2-eq per tonne", size = "N. Species", colour = "Risk") +
  theme_clean()

# Plots with combined biodiversity risk index
mm_riskindex_aveghg <- mm_risk_aveghg %>%
  mutate(risk.index = ifelse(risk == "low", 1*n.species,
                             ifelse(risk == "medium", 2*n.species,
                                    ifelse(risk == "high", 3*n.species, NA)))) %>%
  group_by(fishery, mm_species, gear, ghg.ave, ghg.se) %>%
  summarise(risk.index = sum(risk.index))

ggplot(mm_riskindex_aveghg, aes(x = ghg.ave, y = risk.index)) +
  geom_point() +
  geom_errorbar(aes(xmin=ghg.ave-(1.96*ghg.se), xmax=ghg.ave+(1.96*ghg.se)), width=.1) +
  labs(y = "Risk index", x = "kg CO2-eq per tonne", size = "N. Species", colour = "Risk") +
  theme_clean() +
  facet_wrap(~gear, nrow = length(unique(mm_riskindex_aveghg$gear)))

ggplot(mm_riskindex_aveghg, aes(x = ghg.ave, y = risk.index, colour = factor(gear))) +
  geom_point() +
  geom_errorbar(aes(xmin=ghg.ave-(1.96*ghg.se), xmax=ghg.ave+(1.96*ghg.se)), width=.1) +
  labs(y = "Risk index", x = "kg CO2-eq per tonne", size = "N. Species", colour = "Risk") +
  theme_clean() 

