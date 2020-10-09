rm(list=ls())
library(tidyverse)
library(rstan)
library(taxize)
library(data.table)

# Mac
#datadir <- "/Volumes/jgephart/BFA Environment 2/Data"
#outdir <- "/Volumes/jgephart/BFA Environment 2/Outputs"
# Windows
datadir <- "K:/BFA Environment 2/Data"
outdir <- "K:BFA Environment 2/Outputs"

lca_dat <- read.csv(file.path(datadir, "LCA_compiled.csv"), fileEncoding="UTF-8-BOM") #fileEncoding needed when reading in file from windows computer (suppresses BOM hidden characters)

# Clean lca_dat:
# Clean species names
# Manually fill in blank scientific names
# Remove unnecessary columns
lca_fix_names <- lca_dat %>%
  mutate(Species.scientific.name = case_when(Species.common.name == "Freshwater prawn" ~ "Dendrobranchiata",
                                             Species.common.name == "Indo-Pacific swamp crab; Swimming crabs, etc. nei" ~ "Brachyura",
                                             Species.common.name == "Red crayfish" ~ "Astacidea", # crayfish are split into two superfamilies, so go to the next higher-classification, infraorder = Astacidea
                                             Species.common.name == "Salmonids nei" ~ "Salmonidae",
                                             Species.common.name == "Yellowtail_Seriola_Almaco jack" ~ "Seriola rivoliana",
                                             TRUE ~ Species.scientific.name)) %>%
  mutate(clean_sci_name = case_when(str_detect(Species.scientific.name, "spp") ~ str_replace(Species.scientific.name, pattern = " spp\\.| spp", replacement = ""),
                                    Species.scientific.name == "Morone chrysops x M. saxatilis" ~ "Morone",
                                    TRUE ~ Species.scientific.name)) %>%
  #filter(Species.scientific.name != "") %>% # blank scientific names manually filled in above
  mutate(FCR = case_when(str_detect(Species.scientific.name, "Thunnus") ~ FCR/5,
                         TRUE ~ FCR)) %>%
  select(-c(Date_entered, Description, Note_on_system, Notes, Person_entering, Product, Production_system, Sample_size, SeaWEED.ID, Source, Specific_location, Strain))


################################################################################################################
# use package taxize to get higher classification levels for each species

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
species_list <- sort(unique(lca_fix_names$clean_sci_name))
# TEST: species_list <- c("Penaeus monodon", "Oreochromis niloticus")
# TEST just genus: species_list <- "Macrobrachium"

# Use classification function to get higher ranks for all species
# IMPORTANT: RUN THIS AS A SINGLE LINE, otherwise when function asks for further input Re: osteichthyes, it will accept next line of code as the input
# NOTE: when prompted choose "bony fishes" for osteichthyes
taxa_ranks <- lapply(species_list, classification, db = "ncbi")

# Which taxa did not find a match in the classification database:
classify_test <- data.frame(no_results = unlist(lapply(taxa_ranks, is.na))) %>%
  filter(no_results == TRUE)
# classify_test is empty - all were matched in classification database

# Function to reformat ranks into columns
format_ranks <- function(taxa_ranks) {
  
  # pivot taxa_ranks wider after filtering multiple un-named ranks called "clades"
  higher_ranks <- taxa_ranks[[1]] %>%
    filter(rank %in% c("clade", "no rank") == FALSE) %>%
    # filter(rank %in% c("Class", "Subclass", "Superorder", "Order", "Suborder", "Superfamily", "Family", "Subfamily", "Genus", "Species")) %>%
    select(name, rank) %>%
    pivot_wider(names_from = rank, values_from = name)
  
  # get the last column name
  last_col_name <- names(higher_ranks)[ncol(higher_ranks)]
  
  # record the rank of the taxa in a new column sci_name_rank and rename the taxa column as clean_sci_name (this is what will be used to join bank with lca_dat)
  higher_ranks <- higher_ranks %>%
    mutate(sci_name_rank = colnames(.)[ncol(.)]) %>%
    rename(clean_sci_name := !!last_col_name)
  
  return(higher_ranks)
}

# Reformat all classification function outputs into columns and bind into single data table
rank_cols <- lapply(taxa_ranks, format_ranks)
rank_cols_dt <- data.table::rbindlist(rank_cols, use.names = TRUE, fill = TRUE) %>%
  # note: need to recreate species column; lowest ranking for each sci_name was renamed to "clean_sci_name" to allow for left_join with lca_fix_names
  mutate(species = case_when(sci_name_rank == "species" ~ clean_sci_name,
                             TRUE ~ "")) %>%
  mutate(species = na_if(species, "")) %>%
  # Order taxonomic ranks in descending order; 
  select(superkingdom, kingdom, phylum, subphylum, superclass, class, subclass, infraclass, cohort, subcohort, superorder, order, suborder, infraorder, superfamily, family, subfamily, tribe, genus, species, clean_sci_name, sci_name_rank)

# Join with lca_fix_names
lca_dat_clean <- lca_fix_names %>%
  left_join(rank_cols_dt, by = "clean_sci_name")

## OTHER CLEANING? ADD HERE...


# Output lca_with_ranks to use for Bayesian analyses
write.csv(lca_dat_clean, file.path(datadir, "lca_clean_with_ranks.csv"))
  
