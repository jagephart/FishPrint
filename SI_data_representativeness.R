# Get additional plots/summaries for SI

rm(list=ls())
library(tidyverse)
library(countrycode)
source("Functions.R") # for rebuild_fish

datadir <- "/Volumes/jgephart/BFA Environment 2/Data"
outdir <- "/Volumes/jgephart/BFA Environment 2/Outputs"

#################################################################
# FORMAT DATA for journal SI:

# Create new version with cleaned references
# lca_dat <- read.csv(file.path(datadir, "LCA_compiled_20201222.csv"), fileEncoding="UTF-8-BOM") %>% #fileEncoding needed when reading in file from windows computer (suppresses BOM hidden characters)
#   # Clean up references (Source column)
#   mutate(Source = str_replace_all(Source, pattern = "-", replacement = " "))
# write.csv(lca_dat, file.path(datadir, "LCA_compiled_20210405.csv"), row.names = FALSE)
# Then manually clean references (Remove "a" which all seem unnecessary except for Iribarren 2010)

lca_dat <- read.csv(file.path(datadir, "LCA_compiled_20210405.csv"), fileEncoding="UTF-8-BOM") #fileEncoding needed when reading in file from windows computer (suppresses BOM hidden characters)
source("Functions.R")

# Create clean and aggregated (Patrick's and Wenbo's studies) LCA data for Zenodo data repository
# lca_dat_for_si should be able to work with full code starting at bayes_01_process_data_for_analysis starting with the function (add_taxa_group)
lca_dat_for_si <- clean.lca(LCA_data = lca_dat)

# Aggregate data:
# Patrick's Indonesia 
henriksson_indo <- lca_dat_for_si %>%
  filter(Source == "Henriksson et al. 2019" & Country == "Indonesia") 

henriksson_indo_agg <- lca_dat_for_si %>%
  filter(Source == "Henriksson et al. 2019" & Country == "Indonesia") %>%
  group_by(Source, Country, iso3c, clean_sci_name, Production_system_group, Intensity) %>% # NOT grouping by Common.Name, Scientific.Name, Product creates NA's when merging back with full dataset, but don't need these columns for analysis
  summarise(across(c(Yield_m2_per_t, Grow_out_period_days, FCR, feed_soy_new, feed_crops_new, feed_fmfo_new, feed_animal_new, Electricity_kwh, Diesel_L, Petrol_L, NaturalGas_L), mean, na.rm = TRUE),
            Sample_size_n_farms = sum(Sample_size_n_farms),
            study_id = min(study_id)) %>% # just use the study_id for the first row
  ungroup() %>%
  arrange(study_id)

# FIX IT - this doesn't really collapse much; does Patrick want studies collapsed by just the clean_sci_name? but then can't impute data because we won't have Intensity + System info?

# Wenbo's Chinese carp  data:
chn_carp <- lca_dat_for_si %>%
  filter(Source == "Zhang & Newton (unpubl. Data)")

chn_carp_agg <- lca_dat_for_si %>%
  filter(Source == "Zhang & Newton (unpubl. Data)") %>%
  group_by(Source, Country, iso3c, clean_sci_name, Production_system_group, Intensity) %>%
  summarise(across(c(Yield_m2_per_t, Grow_out_period_days, FCR, feed_soy_new, feed_crops_new, feed_fmfo_new, feed_animal_new, Electricity_kwh, Diesel_L, Petrol_L, NaturalGas_L), mean, na.rm = TRUE),
            Sample_size_n_farms = sum(Sample_size_n_farms),
            study_id = min(study_id)) %>%
  ungroup() %>%
  arrange(study_id)

aggregated_source <- c("Henriksson et al. 2019", "Zhang & Newton \\(unpubl. Data\\)")

lca_dat_for_si_clean <- lca_dat_for_si %>%
  # Filter out the raw data that needs to be aggregated
  filter((Source == "Henriksson et al. 2019" & Country == "Indonesia")==FALSE) %>%
  filter(Source != "Zhang & Newton (unpubl. Data)") %>%
  # Add aggregated data back in
  bind_rows(henriksson_indo_agg) %>%
  bind_rows(chn_carp_agg) %>%
  mutate(data_type = if_else(Source %in% aggregated_source, true = "aggregated", false = "raw"))

write.csv(lca_dat_for_si_clean, file = file.path(datadir, "LCA_compiled_for_SI.csv"), row.names = FALSE)


# OLD CODE:
# lca_dat_clean_for_si <- lca_dat_clean %>%
#   filter(str_detect(Source, pattern = "unpubl")==FALSE) # Remove unpublished data
# 
# write.csv(lca_dat_clean_for_si, file.path(outdir, "data_for_si.csv"))

# OLD CODE:
# Load Imputed Data
# lca_full_dat <- read.csv(file.path(datadir, "2021-01-06_lca-dat-imputed-vars_rep-sqrt-n-farms.csv"), fileEncoding="UTF-8-BOM")
# 
# # Load Cleaned Data (before imputation) just to get "Source" and "study_id" columns:
# lca_dat <- read.csv(file.path(datadir, "lca_clean_with_groups.csv"), fileEncoding="UTF-8-BOM") %>% #fileEncoding needed when reading in file from windows computer (suppresses BOM hidden characters)
#   select(study_id, Source) %>%
#   mutate(Source = str_trim(Source))
# 
# # Join back together to get "Source" column back
# lca_dat_for_si <- lca_full_dat %>%
#   #left_join(lca_dat, by = intersect(names(.), names(lca_dat))) %>%
#   left_join(lca_dat, by = "study_id") %>%
#   select(c( "Source", !!names(lca_full_dat))) %>%
#   # remove data replication:
#   select(-study_id) %>%
#   unique() 
#   
# # data that should be aggregated prior to publishing:
# # Patrick's Indonesia 
# henriksson_indo <- lca_dat_for_si %>%
#   filter(Source == "Henriksson et al. 2019" & Country == "Indonesia") %>%
#   group_by(Source, clean_sci_name, taxa, intensity, system, Country, iso3c) %>%
#   summarise(across(is.numeric, mean, na.rm = TRUE)) %>%
#   ungroup()
# 
# # Chinese carp data:
# chn_carp_dat <- lca_dat_for_si %>%
#   filter(Source == "Zhang & Newton (unpubl. Data)") %>%
#   group_by(Source, clean_sci_name, taxa, intensity, system, Country, iso3c) %>%
#   summarise(across(is.numeric, mean, na.rm = TRUE)) %>%
#   ungroup()
# 
# aggregated_source <- c("Henriksson et al. 2019", "Zhang & Newton \\(unpubl. Data\\)")
# 
# lca_dat_for_si_clean <- lca_dat_for_si %>%
#   # Filter out the raw data that needs to be aggregated
#   filter((Source == "Henriksson et al. 2019" & Country == "Indonesia")==FALSE) %>%
#   filter(Source != "Zhang & Newton (unpubl. Data)") %>%
#   # Add aggregated data back in
#   bind_rows(henriksson_indo) %>%
#   bind_rows(chn_carp_dat) %>%
#   mutate(feed_data_type = if_else(Source %in% aggregated_source, true = "aggregated", false = feed_data_type),
#          fcr_data_type = if_else(Source %in% aggregated_source, true = "aggregated", false = fcr_data_type),
#          electric_data_type = if_else(Source %in% aggregated_source, true = "aggregated", false = electric_data_type),
#          diesel_data_type = if_else(Source %in% aggregated_source, true = "aggregated", false = diesel_data_type),
#          petrol_data_type = if_else(Source %in% aggregated_source, true = "aggregated", false = petrol_data_type),
#          natgas_data_type = if_else(Source %in% aggregated_source, true = "aggregated", false = natgas_data_type),
#          yield_data_type = if_else(Source %in% aggregated_source, true = "aggregated", false = yield_data_type))
# 
# write.csv(lca_dat_for_si_clean, file = file.path(outdir, "lca_data_for_si_clean.csv"), row.names = FALSE)

#################################################################
# PRODUCTION CALCULATIONS:

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
  mutate(plot_name = as.factor(plot_name)) %>%
  mutate(plot_name_2 = paste(plot_name, " (", country_iso3_code, ")", sep = ""))
  
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
# Proportion of global aquaculture production represented in our study: 0.747532

# How much of total production is carps?
prod_by_taxa %>%
  filter(plot_name %in% c("misc carps", "silver/bighead")) %>%
  pull(taxa_prod) %>% sum() / sum(prod_by_taxa$taxa_prod)

# How much of total WILD CAPTURE production do our wild taxa groupings represent?
# First, REMOVE the following isscaap groups (mostly non-human purposes)
non_human_isscaap <- c("Pearls, mother-of-pearl, shells", "Corals", "Sponges", "Brown seaweeds", "Red seaweeds", "Green seaweeds", "Miscellaneous aquatic plants") 

# Match taxa group to ISSCAAP group(s)
# 2012 - 2019
# Aligning Rob's groupings with ISSCAAP: according to http://www.fao.org/3/Y5852E10.htm
# Jacks, mullets, sauries are Misc demersal fishes
# Redfishes, basses, and congers are Misc coastal fishes
# Large pelagic fishes are Tunas, bonitos, billfishes + Misc pelagic fishes
# Small pelagic fishes are Herrings, sardines, anchovies <- all Clupeidae found here, find number to adjust this
wild_prod_clean <- fishstat_dat %>% 
  filter(year > 2012) %>%
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
# Including non-human isscaap groups:
# Proportion of global wild capture production represented in our study: 
# 0.6339886 (filter Osteichthyes, adjust small pelagics) 
# 0.6622046 (filter Osteichthyes, no adjustment to small pelagics)
# 0.5234017 (keep Osteichthyes, adjust small pelagics)
# 0.5541475 (keep Osteichthyes, no adjustment to small pelagics)

# Removing non-human isscaap groups:
# 0.6445885 (filter Osteichthyes, adjust small pelagics)
# 0.6724095 (filter Osteichthyes, no adjustment to small pelagics)
# 0.5306052 (keep Osteichthyes, adjust small pelagics)
# 0.5612758 (keep Osteichthyes, no adjustment to small pelagics)

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
  scale_y_log10(labels = c(1, 10, 100), breaks = c(1e6, 1e7, 1e8)) +
  labs(x = "No. of Farms", y = "Total Aquaculture Production \n2012-2017 (million tonnes)", size = "No. of Studies") +
  theme(legend.position=c(0.8, 0.25),
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
  group_by(taxa_group_name, plot_name_2, country_iso3_code) %>%
  summarise(taxa_iso_prod = sum(quantity)) %>%
  ungroup()

plot_country_dat <- country_dat %>%
  left_join(prod_by_taxa_iso, by = c("taxa_group_name", "iso3c" = "country_iso3_code")) %>%
  drop_na()

# Identify which points to label
# Get top 10 by number of farms
top_n_farms <- plot_country_dat %>%
  arrange(desc(n_farms)) %>%
  slice_head(n=3) %>%
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

interesting_points <- unique(c(top_n_farms, low_n_farms, over_represented))
  



#summary(lm(taxa_iso_prod ~ n, data = plot_country_dat))

ggplot(data = plot_country_dat, aes(x = n_farms, y = taxa_iso_prod, size = n_studies)) +
  geom_point() +
  #geom_jitter() + 
  geom_text_repel(data = . %>% mutate(label = if_else(plot_name_2 %in% interesting_points, true = plot_name_2, false = "")), 
                  aes(label=label), box.padding = unit(0.5, "lines"), size = 4) +
  #scale_shape_manual(values = c(0:11)) +
  #geom_smooth(method = "lm", color = "red") +
  theme_classic() +
  scale_y_log10(labels = c(1, 100, 10000), breaks = c(1e3, 1e5, 1e7)) +
  scale_x_log10() +
  #coord_cartesian(xlim = c(-100, 300)) +
  labs(x = "No. of Farms", y = "Total Aquaculture Production \n2012-2017 (thousand tonnes)", size = "No. of Studies") +
  theme(legend.position=c(0.8, 0.25),
        legend.background = element_rect(linetype = 1, size = 0.5, colour = 1)) +
  plot_theme 
# PRO TIP: adjust image in Plot Window until happy with geom_text_repel, then save manually with "Export"
#ggsave(filename = file.path(outdir, paste("plot_national_production_vs_n_farms.png", sep = "")), width = 11, height = 8.5)

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