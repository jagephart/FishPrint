# Summarize results from posteriors

# Set data directories
datadir <- "/Volumes/jgephart/BFA Environment 2/Data"
outdir <- "/Volumes/jgephart/BFA Environment 2/Outputs"

#_______________________________________________________________________________________________________________________#
# Load packages and source functions
#_______________________________________________________________________________________________________________________#
library(tidyverse)
library(ggplot2)

#_______________________________________________________________________________________________________________________#
# Load results files and merge
#_______________________________________________________________________________________________________________________#

df <- read.csv(file.path(outdir, "PRIORS/GHG/summary_Global warming potential_Economic-allocation_ON-FARM-TAXA-LEVEL-WEIGHTED.csv"))
df <- df %>%
  mutate(source = "on-farm", stressor = "GHG") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "on_farm", "lower_95" = ".lower", "upper_95" = ".upper")

tmp <- read.csv(file.path(outdir, "PRIORS/GHG/summary_Global warming potential_Economic-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>% 
  mutate(source = "off-farm", stressor = "GHG") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "off_farm", "lower_95" = ".lower", "upper_95" = ".upper")
df <- df %>%
  bind_rows(tmp)

tmp <- read.csv(file.path(outdir, "PRIORS/Nitrogen/summary_Marine eutrophication_Economic-allocation_ON-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "on-farm", stressor = "N") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "on_farm", "lower_95" = ".lower", "upper_95" = ".upper")
df <- df %>%
  bind_rows(tmp) 

tmp <- read.csv(file.path(outdir, "PRIORS/Nitrogen/summary_Marine eutrophication_Economic-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "off-farm", stressor = "N") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "off_farm", "lower_95" = ".lower", "upper_95" = ".upper")
df <- df %>%
  bind_rows(tmp) 

tmp <- read.csv(file.path(outdir, "PRIORS/Phosphorus/summary_Freshwater eutrophication_Economic-allocation_ON-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "on-farm", stressor = "P") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "on_farm", "lower_95" = ".lower", "upper_95" = ".upper")
df <- df %>%
  bind_rows(tmp) 

tmp <- read.csv(file.path(outdir, "PRIORS/Phosphorus/summary_Freshwater eutrophication_Economic-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "off-farm", stressor = "P") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "off_farm", "lower_95" = ".lower", "upper_95" = ".upper")
df <- df %>%
  bind_rows(tmp) 

tmp <- read.csv(file.path(outdir, "PRIORS/Water/summary_Water consumption_Economic-allocation_ON-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "on-farm", stressor = "Water") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "on_farm", "lower_95" = ".lower", "upper_95" = ".upper")
df <- df %>%
  bind_rows(tmp) 

tmp <- read.csv(file.path(outdir, "PRIORS/Water/summary_Water consumption_Economic-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "off-farm", stressor = "Water") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "off_farm", "lower_95" = ".lower", "upper_95" = ".upper")
df <- df %>%
  bind_rows(tmp) 

tmp <- read.csv(file.path(outdir, "PRIORS/Land/summary_Land Use_Economic-allocation_ON-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "on-farm", stressor = "Land") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "on_farm", "lower_95" = ".lower", "upper_95" = ".upper")
df <- df %>%
  bind_rows(tmp) 

tmp <- read.csv(file.path(outdir, "PRIORS/Land/summary_Land Use_Economic-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "off-farm", stressor = "Land") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "off_farm", "lower_95" = ".lower", "upper_95" = ".upper")
df <- df %>%
  bind_rows(tmp) 

df_total <- df %>% 
  group_by(taxa, full_taxa_name, stressor) %>% 
  summarise(total = sum(median))

#_______________________________________________________________________________________________________________________#
# Write out results tables
#_______________________________________________________________________________________________________________________#


write.csv(df, file.path(outdir, "SI_stressor_results_econ.csv"), row.names = FALSE)
