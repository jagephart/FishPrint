# TITLE: clean_edible_portion.R
# AUTHOR: Jessica
# DATE: 12-Apr-21

# Load packages
library(tidyverse)
library(ggthemes)
library(rfishbase)
library(ggpubr)

# Set data directories
datadir <- "/Volumes/jgephart/BFA Environment 2/Data"
outdir <- "/Volumes/jgephart/BFA Environment 2/Outputs"

# Load data FIXIT: Copy data to drive and read from there
df_all <- read.csv("Data/edible_muscle_raw.csv")
df_aqua <- read.csv("Data/aquaculture_edible_portion.csv")
taxa_isscaap <- read.csv("Data/taxacode_ISSCAAP_matches.csv")

df_stressors <- read.csv(file.path(outdir, "SI_stressor_results.csv"))

# Select edible portion data and filter out NAs
df_all <- df_all %>%
  select(TaxonName, TaxonKey, TaxLevel, CommonName, Group_coarse, 
         max_edible_fraction_from_live_weight_percent, ref_meat_yield) %>%
  filter(!is.na(max_edible_fraction_from_live_weight_percent))

# Standardize scientific names
taxa_synonyms <- synonyms(taxa_isscaap$ASFIS_sci_name)
taxa_synonyms <- taxa_synonyms %>%
  select(synonym, SpecCode, Species)
taxa_isscaap <- taxa_isscaap %>%
  left_join(taxa_synonyms, by = c("ASFIS_sci_name" = "synonym"))

taxa_synonyms <- synonyms(df_all$TaxonName)
taxa_synonyms <- taxa_synonyms %>%
  select(synonym, SpecCode, Species)
df_all <- df_all %>%
  left_join(taxa_synonyms, by = c("TaxonName" = "synonym"))

# Add ISSCAAP groups to data
taxa_isscaap <- taxa_isscaap %>%
  mutate(ASFIS_sci_name = tolower(ASFIS_sci_name), 
         Species = tolower(Species)) %>%
  select("TaxonName" = "ASFIS_sci_name", ISSCAAP_group, "Species_alt" = "Species") 

df_all <- df_all %>%
  mutate(TaxonName = tolower(TaxonName), 
         Species = tolower(Species)) %>%
  left_join(taxa_isscaap, by = "TaxonName")

df_all_unmatched <- df_all %>%
  filter(is.na(ISSCAAP_group)) %>%
  select(-c(Species_alt, ISSCAAP_group)) %>%
  filter(!is.na(Species)) %>%
  left_join(taxa_isscaap, by = c("Species" = "Species_alt")) %>%
  select(-TaxonName.y) %>%
  rename("TaxonName" = "TaxonName.x")

df_all <- df_all %>%
  bind_rows(df_all_unmatched) %>%
  filter(!is.na(ISSCAAP_group)) %>%
  select(-c("SpecCode", "Species", "Species_alt")) %>%
  distinct()

# Summarize mean by ISSCAAP group for capture
capture_summary <- df_all %>% 
  group_by(ISSCAAP_group) %>%
  summarise(edible_mean = mean(max_edible_fraction_from_live_weight_percent),
            edible_median = median(max_edible_fraction_from_live_weight_percent),
            edible_min = min(max_edible_fraction_from_live_weight_percent), 
            edible_max = max(max_edible_fraction_from_live_weight_percent))

df_all %>% 
  group_by(Group_coarse) %>%
  summarise(edible_mean = mean(max_edible_fraction_from_live_weight_percent),
            edible_median = median(max_edible_fraction_from_live_weight_percent),
            edible_min = min(max_edible_fraction_from_live_weight_percent), 
            edible_max = max(max_edible_fraction_from_live_weight_percent))

# Summarize mean by taxa group for aquaculture
aquaculture_summary <- df_aqua %>%
  group_by(fishprint_taxa) %>% 
  summarise(edible_mean = mean(max_edible_fraction_from_live_weight_percent, na.rm = TRUE),
            edible_median = median(max_edible_fraction_from_live_weight_percent, na.rm = TRUE),
            edible_min = min(max_edible_fraction_from_live_weight_percent, na.rm = TRUE), 
            edible_max = max(max_edible_fraction_from_live_weight_percent, na.rm = TRUE))

# Catfish average dress-out percentage from Argue et al. 2003
# (https://www.sciencedirect.com/science/article/pii/S004484860300245X?via%3Dihub#BIB6)
#catfish_ave <- mean(c(61.1, 58.9, 58.3, 58.2, 57.5, 57.3, 57.3, 56.9))

# Replace with dress-out average
#aquaculture_summary$edible_mean[aquaculture_summary$fishprint_taxa == "catfish"] <- catfish_ave

# Or replace with finfish average
aquaculture_summary[aquaculture_summary$fishprint_taxa == "catfish", 2:5] <- df_all %>% 
  group_by(Group_coarse) %>%
  summarise(edible_mean = mean(max_edible_fraction_from_live_weight_percent),
            edible_median = median(max_edible_fraction_from_live_weight_percent),
            edible_min = min(max_edible_fraction_from_live_weight_percent), 
            edible_max = max(max_edible_fraction_from_live_weight_percent)) %>%
  filter(Group_coarse == "finfish") %>%
  select(edible_mean:edible_max)

# Replace misc_diad with finfish average
aquaculture_summary[aquaculture_summary$fishprint_taxa == "misc_diad", 2:5] <- df_all %>% 
  group_by(Group_coarse) %>%
  summarise(edible_mean = mean(max_edible_fraction_from_live_weight_percent),
            edible_median = median(max_edible_fraction_from_live_weight_percent),
            edible_min = min(max_edible_fraction_from_live_weight_percent), 
            edible_max = max(max_edible_fraction_from_live_weight_percent)) %>%
  filter(Group_coarse == "finfish") %>%
  select(edible_mean:edible_max)

# Compare the results fo live weight versus edible portion
# Aquaculture plots
aqua_stressors <- df_stressors %>%
  filter(Production == "Aquaculture") %>%
  right_join(aquaculture_summary, by = c("taxa" = "fishprint_taxa")) %>%
  mutate(median_edible = median/(edible_mean/100)) %>% 
  mutate(median_edible = case_when(
    taxa == "plants" ~ median,
    taxa != "plants" ~ median_edible
  )) %>%
  select(taxa, full_taxa_name, stressor, source, "live" = "median", "edible" = "median_edible") %>%
  pivot_longer(cols = live:edible, names_to = "functional_unit", values_to = "stressor_est")

ggplot(aqua_stressors, aes(x = functional_unit, y = stressor_est, fill = source)) +
  geom_bar(position = "stack", stat = "identity") +
  facet_grid(stressor~full_taxa_name, scales = "free_y") +
  theme(legend.position = "none")

aqua_stressors_total <- aqua_stressors %>%
  group_by(taxa, full_taxa_name, stressor, functional_unit) %>%
  summarise(stressor_est = sum(stressor_est))

ghg_live <- ggplot(aqua_stressors_total %>% filter(functional_unit == "live", stressor == "GHG"), 
             aes(x = reorder(full_taxa_name, stressor_est), y = stressor_est)) +
  geom_bar(stat = "identity") + 
  coord_flip() +
  labs(x = "", y = "GHG") 

ghg_edible <- ggplot(aqua_stressors_total %>% filter(functional_unit == "edible", stressor == "GHG"), 
                aes(x = reorder(full_taxa_name, stressor_est), y = stressor_est)) +
  geom_bar(stat = "identity") + 
  coord_flip() +
  labs(x = "", y = "GHG") 

ggarrange(ghg_live, ghg_edible)

land_live <- ggplot(aqua_stressors_total %>% filter(functional_unit == "live", stressor == "Land"), 
                   aes(x = reorder(full_taxa_name, stressor_est), y = stressor_est)) +
  geom_bar(stat = "identity") + 
  coord_flip() +
  labs(x = "", y = "Land") 

land_edible <- ggplot(aqua_stressors_total %>% filter(functional_unit == "edible", stressor == "Land"), 
                     aes(x = reorder(full_taxa_name, stressor_est), y = stressor_est)) +
  geom_bar(stat = "identity") + 
  coord_flip() +
  labs(x = "", y = "Land") 

ggarrange(land_live, land_edible)

water_live <- ggplot(aqua_stressors_total %>% filter(functional_unit == "live", stressor == "Water"), 
                   aes(x = reorder(full_taxa_name, stressor_est), y = stressor_est)) +
  geom_bar(stat = "identity") + 
  coord_flip() +
  labs(x = "", y = "Water") 

water_edible <- ggplot(aqua_stressors_total %>% filter(functional_unit == "edible", stressor == "Water"), 
                     aes(x = reorder(full_taxa_name, stressor_est), y = stressor_est)) +
  geom_bar(stat = "identity") + 
  coord_flip() +
  labs(x = "", y = "Water") 

ggarrange(water_live, water_edible)

N_live <- ggplot(aqua_stressors_total %>% filter(functional_unit == "live", stressor == "N"), 
                    aes(x = reorder(full_taxa_name, stressor_est), y = stressor_est)) +
  geom_bar(stat = "identity") + 
  coord_flip() +
  labs(x = "", y = "N") 

N_edible <- ggplot(aqua_stressors_total %>% filter(functional_unit == "edible", stressor == "N"), 
                      aes(x = reorder(full_taxa_name, stressor_est), y = stressor_est)) +
  geom_bar(stat = "identity") + 
  coord_flip() +
  labs(x = "", y = "N") 

ggarrange(N_live, N_edible)

P_live <- ggplot(aqua_stressors_total %>% filter(functional_unit == "live", stressor == "P"), 
                 aes(x = reorder(full_taxa_name, stressor_est), y = stressor_est)) +
  geom_bar(stat = "identity") + 
  coord_flip() +
  labs(x = "", y = "P") 

P_edible <- ggplot(aqua_stressors_total %>% filter(functional_unit == "edible", stressor == "P"), 
                   aes(x = reorder(full_taxa_name, stressor_est), y = stressor_est)) +
  geom_bar(stat = "identity") + 
  coord_flip() +
  labs(x = "", y = "P") 

ggarrange(P_live, P_edible)

# Capture plots
capture_summary <- capture_summary %>%
  mutate(full_taxa_name = case_when(
    ISSCAAP_group %in% c("Oysters", "Mussels", "Scallops, pectens", "Clams, cockles, arkshells") ~ "Bivalves",
    ISSCAAP_group %in% c("Squids, cuttlefishes, octopuses") ~ "Cephalopods",
    ISSCAAP_group %in% c("Flounders, halibuts, soles") ~ "Flatfishes",
    ISSCAAP_group %in% c("Cods, hakes, haddocks") ~ "Gadiformes",
    ISSCAAP_group %in% c("Miscellaneous pelagic fishes") ~ "Jacks, mullets, sauries",
    ISSCAAP_group %in% c("Tunas, bonitos, billfishes") ~ "Large pelagic fishes",
    ISSCAAP_group %in% c("Lobsters, spiny-rock lobsters") ~ "Lobsters",
    ISSCAAP_group %in% c("Miscellaneous coastal fishes", "Miscellaneous demersal fishes") ~ "Redfishes, basses, congers",
    ISSCAAP_group %in% c("Salmons, trouts, smelts") ~ "Salmonids",
    ISSCAAP_group %in% c("Shrimps, prawns") ~ "Shrimps",
    ISSCAAP_group %in% c("Herrings, sardines, anchovies") ~ "Small pelagic fishes"
  )) %>% 
  group_by(full_taxa_name) %>%
  summarise(edible_mean = mean(edible_mean))

capture_stressors <- df_stressors %>%
  filter(Production == "Capture") %>%
  right_join(capture_summary, by = c("taxa" = "full_taxa_name")) %>%
  mutate(median_edible = median/(edible_mean/100)) %>% 
  select(taxa, full_taxa_name, stressor, source, "live" = "median", "edible" = "median_edible") %>%
  pivot_longer(cols = live:edible, names_to = "functional_unit", values_to = "stressor_est") %>%
  filter(!is.na(taxa))

ggplot(capture_stressors, aes(x = functional_unit, y = stressor_est)) +
  geom_bar(stat = "identity") +
  facet_wrap(~full_taxa_name, nrow = 1) 


capture_live <- ggplot(capture_stressors %>% filter(functional_unit == "live"), 
                 aes(x = reorder(full_taxa_name, stressor_est), y = stressor_est)) +
  geom_bar(stat = "identity") + 
  coord_flip() +
  labs(x = "", y = "GHG") 

capture_edible <- ggplot(capture_stressors %>% filter(functional_unit == "edible"), 
                   aes(x = reorder(full_taxa_name, stressor_est), y = stressor_est)) +
  geom_bar(stat = "identity") + 
  coord_flip() +
  labs(x = "", y = "GHG") 

ggarrange(capture_live, capture_edible)

# Write out data files
write.csv(capture_summary, "Data/capture_edible_CFs.csv", row.names = FALSE)
write.csv(aquaculture_summary %>% select(fishprint_taxa, edible_mean), "Data/aquaculture_edible_CFs.csv", row.names = FALSE)
  