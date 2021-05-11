# Get additional plots/summaries for SI

rm(list=ls())
library(tidyverse)
library(countrycode)
source("Functions.R") # for rebuild_fish

datadir <- "/Volumes/jgephart/BFA Environment 2/Data"
outdir <- "/Volumes/jgephart/BFA Environment 2/Outputs"

#################################################################
# FORMAT DATA for journal SI:

lca_dat <- read.csv(file.path(datadir, "LCI_compiled_FINAL.csv"), fileEncoding="UTF-8-BOM") #fileEncoding needed when reading in file from windows computer (suppresses BOM hidden characters)
source("Functions.R")

# Create clean and aggregated (Patrick's and Wenbo's studies) LCI data for publication
# Data should not be replicated yet (comment out )
lca_dat_for_si <- clean.lca(LCA_data = lca_dat, replicate_dat = TRUE) # for SI, we will publish the cleaned, aggregated, and replicated dataset

# Patrick's unpublished India data:
henriksson_unpub <- lca_dat_for_si %>%
  filter(Source == "Henriksson et al. (unpubl. Data)")

henriksson_unpub_agg <- lca_dat_for_si %>%
  filter(Source == "Henriksson et al. (unpubl. Data)") %>%
  group_by(Source, Country, iso3c, clean_sci_name, Production_system_group, Intensity) %>% # NOT grouping by Common.Name, Scientific.Name, Product creates NA's when merging back with full dataset, but don't need these columns for analysis
  summarise(across(c(Yield_m2_per_t, Grow_out_period_days, FCR, feed_soy_new, feed_crops_new, feed_fmfo_new, feed_animal_new, Electricity_kwh, Diesel_L, Petrol_L, NaturalGas_L), mean, na.rm = TRUE),
            Sample_size_n_farms = sum(Sample_size_n_farms),
            study_id = min(study_id),
            n_study = n()) %>% # just use the study_id for the first row
  ungroup() %>%
  arrange(study_id)

# Wenbo's Chinese carp  data:
chn_carp <- lca_dat_for_si %>%
  filter(Source == "Zhang & Newton (unpubl. Data)")

chn_carp_agg <- lca_dat_for_si %>%
  filter(Source == "Zhang & Newton (unpubl. Data)") %>%
  group_by(Source, Country, iso3c, clean_sci_name, Production_system_group, Intensity) %>%
  summarise(across(c(Yield_m2_per_t, Grow_out_period_days, FCR, feed_soy_new, feed_crops_new, feed_fmfo_new, feed_animal_new, Electricity_kwh, Diesel_L, Petrol_L, NaturalGas_L), mean, na.rm = TRUE),
            Sample_size_n_farms = sum(Sample_size_n_farms),
            study_id = min(study_id),
            n_study = n()) %>%
  ungroup() %>%
  arrange(study_id)

aggregated_source <- c("Henriksson et al. (unpubl. Data)", "Zhang & Newton (unpubl. Data)")

lca_dat_for_si_clean <- lca_dat_for_si %>%
  # Filter out the raw data that was aggregated
  filter(Source != "Henriksson et al. (unpubl. Data)") %>%
  filter(Source != "Zhang & Newton (unpubl. Data)") %>%
  # Add aggregated data back in
  bind_rows(henriksson_unpub_agg %>% select(-n_study)) %>%
  bind_rows(chn_carp_agg %>% select(-n_study)) %>%
  mutate(data_type = if_else(Source %in% aggregated_source, true = "aggregated", false = "raw"))

write.csv(lca_dat_for_si_clean, file = file.path(datadir, "LCI_compiled_for_SI.csv"), row.names = FALSE)

#################################################################
# PRODUCTION CALCULATIONS:

# Rebuild FAO fish production from zip file
#fishstat_dat <- rebuild_fish("/Volumes/jgephart/FishStatR/Data/Production-Global/ZippedFiles/GlobalProduction_2019.1.0.zip")
fishstat_dat <- rebuild_fish("/Volumes/jgephart/FishStatR/Data/Production-Global/ZippedFiles/GlobalProduction_2020.1.0.zip")
# DATA DOCUMENTATION: Zip file can be downloaded from FAO FishStat: http://www.fao.org/fishery/static/Data/GlobalProduction_2020.1.0.zip

# Match taxa group to ISSCAAP group(s)
# 2012 - 2019
prod_clean <- fishstat_dat %>% 
  filter(year > 2011) %>%
  #filter(year > 2013) %>%
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
  mutate(plot_name = as.factor(plot_name)) %>%
  mutate(plot_name_2 = paste(plot_name, " (", country_iso3_code, ")", sep = ""))
  
prod_by_taxa <- prod_clean %>%
  group_by(plot_name, taxa_group_name, year) %>%
  summarise(taxa_prod = sum(quantity, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(plot_name, taxa_group_name) %>%
  summarise(taxa_prod = mean(taxa_prod, na.rm = TRUE)) %>% # For MEAN production
  ungroup() %>%
  # Adjust taxa production for algae (only 6.8% used for human consumption according to 2018 FAO food balance sheet)
  mutate(taxa_prod = if_else(plot_name == "plants", true = taxa_prod * 0.068, false = taxa_prod))

# How much of total production do our taxa groupings represent?
prod_by_taxa %>%
  filter(plot_name != "other_taxa") %>%
  pull(taxa_prod) %>% sum() / sum(prod_by_taxa$taxa_prod)
# Proportion of global aquaculture production represented in our study: 0.7586919

# How much of total production is carps?
prod_by_taxa %>%
  filter(plot_name %in% c("misc carps", "silver/bighead")) %>%
  pull(taxa_prod) %>% sum() / sum(prod_by_taxa$taxa_prod)

# How much of total WILD CAPTURE production do our wild taxa groupings represent?
# First, REMOVE the following isscaap groups (mostly non-human purposes)
non_human_isscaap <- c("Pearls, mother-of-pearl, shells", "Corals", "Sponges", "Brown seaweeds", "Red seaweeds", "Green seaweeds", "Miscellaneous aquatic plants") 

# Match taxa group to ISSCAAP group(s)
# 2012 - 2018
# Aligning Rob's groupings with ISSCAAP: according to http://www.fao.org/3/Y5852E10.htm
# Jacks, mullets, sauries are Misc demersal fishes
# Redfishes, basses, and congers are Misc coastal fishes
# Large pelagic fishes are Tunas, bonitos, billfishes + Misc pelagic fishes
# Small pelagic fishes are Herrings, sardines, anchovies <- all Clupeidae found here, find number to adjust this
wild_prod_clean <- fishstat_dat %>% 
  filter(year > 2011) %>%
  #filter(year > 2013) %>%
  filter(unit == "t") %>%
  filter(source_name_en == "Capture production") %>%
  filter(isscaap_group %in% non_human_isscaap == FALSE) %>%
  filter(is.na(quantity)==FALSE) %>%
  filter(species_scientific_name != "Osteichthyes") %>%
  # Identify rows that are part of this study
  mutate(taxa_group_name = case_when(isscaap_group %in% c("Mussels", "Oysters", "Scallops, pectens", "Clams, cockles, arkshells")  ~ "bivalves",
                                     isscaap_group == "Squids, cuttlefishes, octopuses"  ~ "cephalopods",
                                     isscaap_group == "Flouders, halibuts, soles"  ~ "flatfishes",
                                     isscaap_group == "Cods, hakes, haddocks"  ~ "gadiformes",
                                     isscaap_group == "Miscellaneous diadromous fishes"  ~ "jacks, mullets, sauries",
                                     isscaap_group == "Tunas, bonitos, billfishes" ~ "large pelagic fishes",
                                     isscaap_group == "Lobsters, spiny-rock lobsters" ~ "lobsters",
                                     isscaap_group == "Miscellaneous coastal fishes" ~ "redfishes, basses, congers",
                                     isscaap_group == "Salmons, trouts, smelts" ~ "salmonids",
                                     isscaap_group == "Shrimps, prawns" ~ "shrimps", 
                                     isscaap_group == "Herrings, sardines, anchovies" ~ "small pelagic fishes",
                                     TRUE ~ "other_taxa"))

# From Tacon A, Metian M. Fishing for feed or fishing for food: increasing global competition for small pelagic forage fish. Ambio J Hum Environ (2009) 38(6): 294-302.
# 36.2% of the total catch of small pelagics went to non human use in 2006 (63.8% went to human use)
wild_prod_by_taxa <- wild_prod_clean %>%
  group_by(taxa_group_name) %>%
  summarise(taxa_prod = sum(quantity, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(taxa_prod = if_else(taxa_group_name == "small pelagic fishes", true = taxa_prod * .638, false = taxa_prod))

# How much of total production do our wild capture taxa groupings represent?
wild_prod_by_taxa %>%
  filter(taxa_group_name != "other_taxa") %>%
  pull(taxa_prod) %>% sum() / sum(wild_prod_by_taxa$taxa_prod)
# Removing non-human isscaap groups + filter Osteichthyes + adjust small pelagics):
# 0.6459381 

# CHECK CODE:
# wild_study <- wild_prod_by_taxa %>%
#   filter(taxa_group_name != "other_taxa") %>%
#   pull(taxa_prod) %>% sum()
# 
# wild_global <- wild_prod_by_taxa %>%
#   pull(taxa_prod) %>% sum()
# 
# wild_study / wild_global

# How much of total capture AND aquaculture production do our taxa groupings represent?
global_total <- prod_by_taxa %>%
  bind_rows(wild_prod_by_taxa) %>%
  select(-plot_name) %>%
  pull(taxa_prod) %>% sum()

study_total <- prod_by_taxa %>%
  bind_rows(wild_prod_by_taxa) %>%
  filter(taxa_group_name != "other_taxa") %>%
  select(-plot_name) %>%
  pull(taxa_prod) %>% sum()

study_total / global_total
# 0.6335239

#########################################
# Plot n studies (and n farms) per taxa group vs production (shape = taxa group)
# CHOOSE n_type: sum of "n_farms" or "n_studies"

plot_theme <- theme(axis.title.x = element_text(size = 20),
                    axis.text = element_text(size = 18, color = "black", margin = margin(t = 20)),
                    axis.title.y = element_text(size = 20, margin = margin(t = 0, r = 20, b = 0, l = 0)),
                    legend.title = element_text(size = 16), 
                    legend.text = element_text(size = 14),
                    plot.margin = unit(c(1,1,1,1), "cm"))

#n_type <- "n_farms"
taxa_dat <- read.csv(file.path(outdir, "taxa_group_n_and_composition.csv")) %>%
  group_by(taxa_group_name) %>%
  #summarise(n = sum(!!sym(n_type))) %>%
  summarise(n_farms = sum(n_farms),
            n_studies = sum(n_studies)) %>%
  ungroup()

plot_dat <- prod_by_taxa %>% 
  left_join(taxa_dat, by = "taxa_group_name") %>% 
  filter(taxa_group_name != "other_taxa")
         
#summary(lm(taxa_prod ~ n, data = plot_dat))
      
library(ggrepel)
ggplot(data = plot_dat, aes(x = n_farms, y = taxa_prod, size = n_studies)) +
  geom_point() + 
  geom_text_repel(aes(label=plot_name), box.padding = unit(0.5, "lines"), size = 4) +
  #scale_size_discrete(name = "Number of Studies") +
  #scale_shape_manual(values = c(0:11)) +
  #geom_smooth(method = "lm", color = "red") +
  theme_classic() +
  #scale_y_log10() +
  scale_y_log10(labels = c(1, 10, 100), breaks = c(1e5, 1e6, 1e7)) + 
  #labs(x = "No. of Farms", y = paste("Mean Aquaculture Production \n2014-2018 ", bquote('(10'^5*' tonnes)'), sep = ""), size = "No. of Studies") +
  #labs(x = "No. of Farms", y = expression(Mean~Aquaculture~Production~2014-2018~(10^5~tonnes)), size = "No. of Studies") +
  labs(x = "No. of Farms", y = expression(atop("Mean Aquaculture Production", paste((2014-2018~"in"~10^5~tonnes)))), size = "No. of Studies") +
  theme(legend.position=c(0.90, 0.25),
        legend.background = element_rect(linetype = 1, size = 0.5, colour = 1)) +
  plot_theme 
ggsave(filename = file.path(outdir, paste("plot_production_vs_n_farms.png", sep = "")), width = 11, height = 8.5)
ggsave(filename = file.path(outdir, paste("plot_production_vs_n_farms.tiff", sep = "")), width = 11, height = 8.5)


#########################################
# Plot n farms (and n studies) per country vs national production per taxa group (shape = country)
# CHOOSE n_type: sum of "n_farms" or "n_studies"
#n_type <- "n_studies"

# Clean LCA data
country_dat <- read.csv(file.path(outdir, "taxa_group_n_and_composition.csv")) %>%
  group_by(taxa_group_name, iso3c) %>%
  #summarise(n = sum(!!sym(n_type))) %>%
  summarise(n_farms = sum(n_farms),
            n_studies = sum(n_studies)) %>%
  ungroup() %>%
  filter(is.na(iso3c)==FALSE) %>%
  mutate(region = countrycode(iso3c, origin = "iso3c", destination = "region")) 

prod_by_taxa_iso <- prod_clean %>%
  group_by(taxa_group_name, plot_name_2, country_iso3_code, year) %>%
  summarise(taxa_iso_prod = sum(quantity, na.rm = TRUE)) %>% # total per year
  ungroup() %>%
  group_by(taxa_group_name, plot_name_2, country_iso3_code) %>%
  summarise(taxa_iso_prod = mean(taxa_iso_prod, na.rm = TRUE)) %>% # mean across years
  ungroup()

plot_country_dat <- country_dat %>%
  left_join(prod_by_taxa_iso, by = c("taxa_group_name", "iso3c" = "country_iso3_code")) %>%
  drop_na()

# If interested in non-represented production
plot_country_dat_all <- country_dat %>%
  full_join(prod_by_taxa_iso, by = c("taxa_group_name", "iso3c" = "country_iso3_code")) %>%
  filter(taxa_group_name != "other_taxa") %>%
  filter(is.na(n_studies)) %>%
  arrange(desc(taxa_iso_prod))

ggplot(plot_country_dat_all %>% filter(taxa_group_name != "Aquatic plants") %>% filter(taxa_iso_prod > 50000)) +
  #geom_tile(aes(x = iso3c, y = taxa_group_name, fill = log10(taxa_iso_prod))) +
  geom_tile(aes(x = iso3c, y = taxa_group_name, fill = taxa_iso_prod)) +
  #scale_fill_gradient(name = "mean production", trans = "log") +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "", y = "", fill = "mean production") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, color = "black", vjust = 0.5),
        axis.text.y = element_text(color = "black"))
ggsave(filename = file.path(outdir, paste("plot_prod_iso_missing_LCA_data.png", sep = "")), width = 11, height = 3.5)
ggsave(filename = file.path(outdir, paste("plot_prod_iso_missing_LCA_data.tiff", sep = "")), width = 11, height = 3.5)

# LOG-transformed
ggplot(plot_country_dat_all %>% filter(taxa_group_name != "Aquatic plants") %>% filter(taxa_iso_prod > 50000)) +
  geom_tile(aes(x = iso3c, y = taxa_group_name, fill = log10(taxa_iso_prod))) +
  #geom_tile(aes(x = iso3c, y = taxa_group_name, fill = taxa_iso_prod)) +
  #scale_fill_gradient(name = "mean production", trans = "log") +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "", y = "", fill = "log(mean production)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, color = "black", vjust = 0.5),
        axis.text.y = element_text(color = "black"))
ggsave(filename = file.path(outdir, paste("plot_prod_iso_missing_LCA_data-log.png", sep = "")), width = 11, height = 3.5)
ggsave(filename = file.path(outdir, paste("plot_prod_iso_missing_LCA_data-log.tiff", sep = "")), width = 11, height = 3.5)

# Identify which points to label
# Get top 10 by number of farms
top_n_farms <- plot_country_dat %>%
  arrange(desc(n_farms)) %>%
  slice_head(n=7) %>%
  pull(plot_name_2)

# Get top production examples
top_prod <- plot_country_dat %>%
  arrange(desc(taxa_iso_prod)) %>%
  slice_head(n=10) %>%
  pull(plot_name_2)

# Get cases where n_studies & n_farms is low
# low_n <- plot_country_dat %>%
#   filter(n_studies == 1 & n_farms == 1) %>%
#   pull(plot_name_2)

# Get cases where no. of farms AND no. of studies is low but produciton is high
low_n_farms <- plot_country_dat %>%
  filter(n_studies == 1 & n_farms == 1) %>%
  arrange(desc(taxa_iso_prod)) %>%
  slice_head(n = 5) %>%
  pull(plot_name_2)

# Get cases were production is low but no. of farms is high
over_represented <- plot_country_dat %>%
  arrange(taxa_iso_prod) %>%
  filter(n_farms > 1) %>%
  slice_head(n = 3) %>%
  pull(plot_name_2)

interesting_points <- unique(c(top_n_farms, top_prod, low_n_farms, over_represented))

ggplot(data = plot_country_dat, aes(x = n_farms, y = taxa_iso_prod, size = n_studies)) +
  geom_point() +
  #geom_jitter() + 
  geom_text_repel(data = . %>% mutate(label = if_else(plot_name_2 %in% interesting_points, true = plot_name_2, false = "")), 
                  aes(label=label), box.padding = unit(0.5, "lines"), size = 4) +
  #scale_shape_manual(values = c(0:11)) +
  #geom_smooth(method = "lm", color = "red") +
  theme_classic() +
  #scale_y_log10() +
  scale_y_log10(labels = c(1, 100, 10000), breaks = c(1e3, 1e5, 1e7)) + 
  scale_x_log10() +
  #coord_cartesian(xlim = c(-100, 300)) +
  labs(x = "No. of Farms", y = expression(atop("Mean Aquaculture Production", paste((2014-2018~"in"~10^3~tonnes)))), size = "No. of Studies") +
  theme(legend.position=c(0.90, 0.25),
        legend.background = element_rect(linetype = 1, size = 0.5, colour = 1)) +
  plot_theme 
# PRO TIP: adjust image in Plot Window until happy with geom_text_repel, then save manually with "Export"
ggsave(filename = file.path(outdir, paste("plot_national_production_vs_n_farms.png", sep = "")), width = 11, height = 8.5)
ggsave(filename = file.path(outdir, paste("plot_national_production_vs_n_farms.tiff", sep = "")), width = 11, height = 8.5)

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