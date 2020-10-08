rm(list=ls())
library(tidyverse)
library(rstan)
library(taxize)
library(data.table)

datadir <- "/Volumes/jgephart/BFA Environment 2/Data"
outdir <- "/Volumes/jgephart/BFA Environment 2/Outputs"
lca_dat <- read.csv(file.path(datadir, "LCA_compiled.csv"))

# FIX IT - blanks for scientific names; add a taxanomic name that can be recognized by taxize
# FIX IT - remove bivalves from FCR analysis
# Clean lca_dat:
# Clean species names
# Remove unnecessary columns
lca_dat_clean <- lca_dat %>% 
  mutate(clean_sci_name = case_when(str_detect(Species.scientific.name, "spp") ~ str_replace(Species.scientific.name, pattern = " spp\\.| spp", replacement = ""),
                                    Species.scientific.name == "Morone chrysops x M. saxatilis" ~ "Morone",
                                    TRUE ~ Species.scientific.name)) %>%
  filter(Species.scientific.name != "") %>%
  mutate(FCR = case_when(str_detect(Species.scientific.name, "Thunnus") ~ FCR/5,
                         TRUE ~ FCR)) %>%
  select(-c(Date_entered, Description, Note_on_system, Notes, Person_entering, Product, Production_system, Sample_size, SeaWEED.ID, Source, Specific_location, Strain))


################################################################################################################
# First use package taxize to get higher classification levels for each species

# Use NCBI database? more resolution, but have to remove all the different "clade" ranks - can't merge on this name
# classify_ncbi <- classification("Penaeus monodon", db = "ncbi")
# classify_ncbi[[1]]
# classify_worms <- classification("Penaeus monodon", db = "worms")
# classify_worms[[1]]
# 
# classify_ncbi <- classification("Thunnus albacares", db = "ncbi")
# classify_ncbi[[1]]
# classify_worms <- classification("Thunnus albacares", db = "worms")
# classify_worms[[1]]

# Get full species list
species_list <- sort(unique(lca_dat_clean$clean_sci_name))
# TEST: species_list <- c("Penaeus monodon", "Oreochromis niloticus")
# TEST just genus: species_list <- "Macrobrachium"

# Use classification function to get higher ranks for all species
# IMPORTANT: RUN THIS AS A SINGLE LINE, otherwise when function asks for further input Re: osteichthyes, it will accept next line of code as the input
# NOTE: when prompted choose "bony fishes" for osteichthyes
taxa_ranks <- lapply(species_list, classification, db = "ncbi")

# Which taxa did not find a match in the classification database:
classify_test <- data.frame(no_results = unlist(lapply(taxa_ranks, is.na))) %>%
  filter(no_results == TRUE)

# Function to reformat ranks into columns
format_ranks <- function(taxa_ranks) {
  higher_ranks <- taxa_ranks[[1]] %>%
    filter(rank %in% c("clade", "no rank") == FALSE) %>%
    # filter(rank %in% c("Class", "Subclass", "Superorder", "Order", "Suborder", "Superfamily", "Family", "Subfamily", "Genus", "Species")) %>%
    select(name, rank) %>%
    pivot_wider(names_from = rank, values_from = name) %>%
    mutate(name = names(taxa_ranks))
}

# Reformat all classification function outputs into columns and bind into single data table
rank_cols <- lapply(taxa_ranks, format_ranks)
rank_cols_dt <- data.table::rbindlist(rank_cols, use.names = TRUE, fill = TRUE) %>%
  # Order taxonomic ranks in descending order
  select(superkingdom, kingdom, phylum, subphylum, superclass, class, subclass, cohort, superorder, order, suborder, infraorder, superfamily, family, subfamily, tribe, genus, species, name) #%>%

### LEFT OFF HERE: able to get column name that matches the taxa name, but need to vectorize this for the entire data frame
# test:
rank_cols_row <- rank_cols_dt[1,]
colnames(rank_cols_row)[match(rank_cols_row$name, rank_cols_row)]

  
# Join with lca_dat_clean
lca_with_ranks <- lca_dat_clean %>%
  left_join(rank_cols_dt, by = c("clean_sci_name" = "species"))




lca_with_ranks %>%
  select(clean_sci_name, genus) %>%
  unique()

# Create species identifiers but keep FCR = NA when it is the only study entry for that species (remove all other NAs)
# Create species level and other higher group levels based on examining NAs
lca_dat_groups_full <- lca_with_ranks %>%
  group_by(clean_sci_name) %>%
  mutate(n_studies = n()) %>%
  ungroup() %>%
  filter((is.na(FCR) & n_studies != 1)==FALSE) %>%
  mutate(Species.scientific.name = as.factor(Species.scientific.name),
         sp = as.numeric(Species.scientific.name)) #%>%
mutate()
#select(Species.scientific.name, FCR, sp)
