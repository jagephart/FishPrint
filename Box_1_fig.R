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
  filter(str_detect(Reference, pattern = "Brown")| str_detect(Reference, pattern = "Micheli")) %>%
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

base_size <- 6
base_family <- "sans"

cols <- c("#57D182", "#FFD947", "#FFA647", "#70468C")


png("Fig3.png", width = 89, height = 89, units = "mm", res = 300)
ggplot(mm_riskindex_aveghg, aes(x = ghg.ave, y = risk.index, colour = factor(gear))) +
  geom_point() +
  geom_errorbar(aes(xmin=ghg.ave-(1.96*ghg.se), xmax=ghg.ave+(1.96*ghg.se)), width=.1) +
  scale_color_manual(values = cols, labels = function(x) str_wrap(x, width = 12)) +
  geom_text(aes(x = ghg.ave, y = risk.index, label = mm_species), hjust = 0.2, vjust = -0.7, 
            size = 2, colour = "black") +
  labs(y = "Risk index", x = "kg CO2-eq per tonne", size = "N. Species", colour = "Risk") +
  guides(colour=guide_legend(nrow=1, byrow=TRUE)) + 
  theme(axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"), 
        axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"), 
        axis.text = element_text(size = ceiling(base_size), colour = "black"),
        axis.title = element_text(size = ceiling(base_size)), 
        panel.grid.minor = element_blank(), 
        panel.grid.major.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        panel.background = element_blank(), panel.border = element_blank(), 
        strip.background = element_rect(linetype = 0), strip.text = element_text(), 
        strip.text.x = element_text(vjust = 0.5), strip.text.y = element_text(angle = -90), 
        legend.text = element_text(size = ceiling(0.9*base_size), family = "sans"), 
        legend.title = element_blank(), 
        legend.key = element_rect(fill = "white", colour = NA), 
        legend.position="bottom",
        legend.margin=margin(c(1,1,1,1)),
        plot.title = element_text(size = ceiling(base_size*1), face = "bold"), 
        plot.subtitle = element_text(size = ceiling(base_size*1)))
dev.off()
