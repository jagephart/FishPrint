# Get additional plots/summaries for SI

rm(list=ls())
library(tidyverse)
library(countrycode)
source("Functions.R") # for rebuild_fish

datadir <- "/Volumes/jgephart/BFA Environment 2/Data"
outdir <- "/Volumes/jgephart/BFA Environment 2/Outputs"


# FIX IT - move all of this to 01_process data for analysis since it uses the same datasets (fishstat_dat, lca_dat_clean_groups, etc)

# Rebuild FAO fish production from zip file
fishstat_dat <- rebuild_fish("/Volumes/jgephart/FishStatR/Data/Production-Global/ZippedFiles/GlobalProduction_2019.1.0.zip")

# Match taxa group to ISSCAAP group(s)
# 2012 - 2019
prod_clean <- fishstat_dat %>% 
  filter(year > 2012) %>%
  filter(unit == "t") %>%
  filter(source_name_en %in% c("Aquaculture production (marine)", "Aquaculture production (brackishwater)", "Aquaculture production (freshwater)")) %>%
  filter(is.na(quantity)==FALSE) %>%
  # Identify rows that are part of this study
  mutate(taxa_group_name = case_when(isscaap_group %in% c("Brown seaweeds", "Red seaweeds")  ~ "Aquatic plants",
                                     isscaap_group %in% c("Oysters", "Mussels") ~ "Bivalves",
                                     order == "SILURIFORMES" ~ "Catfish", 
                                     species_scientific_name == "Chanos chanos" ~ "Milkfish",
                                     isscaap_group == "Miscellaneous diadromous fishes" | species_name_en %in% c("Arctic char", "Chars nei") ~ "Miscellaneous diadromous fishes",
                                     isscaap_group %in% c("Miscellaneous coastal fishes", "Miscellaneous demersal fishes", "Miscellaneous pelagic fishes") ~ "Miscellaneous marine fishes",
                                     isscaap_group == "Carps, barbels and other cyprinids" & str_detect(species_scientific_name, "Hypophthalmichthys")==FALSE ~ "Other carps, barbels and cyprinids",                      
                                     str_detect(species_name_en, "salmon|Salmon") ~ "Salmon", # Match by common name
                                     isscaap_group == "Shrimps, prawns" ~ "Shrimps, prawns",
                                     isscaap_group == "Carps, barbels and other cyprinids" & str_detect(species_scientific_name, "Hypophthalmichthys") ~ "Silver and bighead carp",
                                     isscaap_group == "Tilapias and other cichlids" ~ "Tilapias and other cichlids",
                                     str_detect(species_name_en, "trout|Trout") ~ "Trout", # Match by common name
                                     TRUE ~ "other_taxa")) %>%
  mutate(plot_name = tolower(taxa_group_name)) %>%
  mutate(plot_name = case_when(plot_name == "aquatic plants" ~ "plants",
                               plot_name == "miscellaneous diadromous fishes" ~ "misc diad",
                               plot_name == "miscellaneous marine fishes" ~ "misc marine",
                               plot_name == "other carps, barbels and cyprinids" ~ "misc carps",
                               plot_name == "shrimps, prawns" ~ "shrimp",
                               plot_name == "silver and bighead carp" ~ "silver/bighead",
                               plot_name == "tilapias and other cichlids" ~ "tilapia",
                               TRUE ~ plot_name)) %>%
  mutate(plot_name = as.factor(plot_name))
  
prod_by_taxa <- prod_clean %>%
  group_by(plot_name, taxa_group_name) %>%
  summarise(taxa_prod = sum(quantity, na.rm = TRUE)) %>%
  ungroup() %>%
  # Adjust taxa production for algae (only 6.8% used for human consumption according to 2018 FAO food balance sheet)
  mutate(taxa_prod = if_else(plot_name == "plants", true = taxa_prod * 0.068, false = taxa_prod))

# How much of total production do our taxa groupings represent?
prod_by_taxa %>%
  filter(plot_name != "other_taxa") %>%
  pull(taxa_prod) %>% sum() / sum(prod_by_taxa$taxa_prod)
# Proportion of global aquaculture production represented by analysis: 0.747532
  
#########################################
# Plot n studies (and n farms) per taxa group vs production (shape = taxa group)
# CHOOSE n_type: sum of n_farms or n_studies
n_type <- "n_farms"
taxa_dat <- read.csv(file.path(outdir, "taxa_group_n_and_composition.csv")) %>%
  group_by(taxa_group_name) %>%
  summarise(n = sum(!!sym(n_type))) %>%
  ungroup()

plot_dat <- prod_by_taxa %>% 
  left_join(taxa_dat, by = "taxa_group_name") %>% 
  filter(taxa_group_name != "other_taxa")
         
summary(lm(taxa_prod ~ n, data = plot_dat))
        
ggplot(data = plot_dat, aes(x = n, y = taxa_prod, shape = plot_name)) +
  geom_point() + 
  scale_shape_manual(values = c(0:11)) +
  geom_smooth(method = "lm", color = "red") +
  theme_classic() +
  labs(x = paste(unlist(str_split(n_type, "_")), collapse = " "), y = "Total Aquaculture Production 2012-2017", shape = "taxa group")
ggsave(filename = file.path(outdir, paste("plot_production_vs_", n_type, ".png", sep = "")), width = 11, height = 8.5)


#########################################
# Plot n farms (and n studies) per country vs national production per taxa group (shape = country)

# Clean LCA data
country_dat <- read.csv(file.path(outdir, "taxa_group_n_and_composition.csv")) %>%
  group_by(taxa_group_name, iso3c) %>%
  summarise(n = sum(!!sym(n_type))) %>%
  ungroup() %>%
  filter(is.na(iso3c)==FALSE) %>%
  mutate(region = countrycode(iso3c, origin = "iso3c", destination = "region")) %>%
  arrange(desc(n))

prod_by_taxa_iso <- prod_clean %>%
  group_by(taxa_group_name, country_iso3_code) %>%
  summarise(taxa_iso_prod = sum(quantity)) %>%
  ungroup()

plot_country_dat <- country_dat %>%
  left_join(prod_by_taxa_iso, by = c("taxa_group_name", "iso3c" = "country_iso3_code")) %>%
  drop_na()

summary(lm(taxa_iso_prod ~ n, data = plot_country_dat))

ggplot(data = plot_country_dat, aes(x = n, y = taxa_iso_prod)) +
  geom_point() + 
  scale_shape_manual(values = c(0:11)) +
  geom_smooth(method = "lm", color = "red") +
  theme_classic() +
  labs(x = paste(unlist(str_split(n_type, "_")), collapse = " "), y = "Total Aquaculture Production 2012-2017")
ggsave(filename = file.path(outdir, paste("plot_national_production_vs_", n_type, ".png", sep = "")), width = 11, height = 8.5)

# Which are the three outlier studies
plot_country_dat %>%
  arrange(desc(taxa_iso_prod)) %>%
  select(taxa_group_name, iso3c)

# same for n_type = "n_studies" or "n_farms"
# A tibble: 51 x 2
# taxa_group_name                    iso3c
# <chr>                              <chr>
# 1 Other carps, barbels and cyprinids CHN  
# 2 Silver and bighead carp            CHN  
# 3 Other carps, barbels and cyprinids IND  