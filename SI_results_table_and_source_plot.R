# Bind results from posteriors for SI

# Set data directories
rm(list=ls())
datadir <- "/Volumes/jgephart/BFA Environment 2/Data"
outdir <- "/Volumes/jgephart/BFA Environment 2/Outputs"

#_______________________________________________________________________________________________________________________#
# Load packages and source functions
#_______________________________________________________________________________________________________________________#
library(tidyverse)
library(ggplot2)
library(ggpubr)

#_______________________________________________________________________________________________________________________#
# Combine tables for on and off farm stressors
#_______________________________________________________________________________________________________________________#
#_______________________________________________________________________________________________________________________#
# Load main results files and merge (Priors + Mass allocation + Edible Weight)
#_______________________________________________________________________________________________________________________#

# ON-FARM GHG
df <- read.csv(file.path(outdir, "Nature-submitted-2021-08/RERUN-GHG-and-Water-Bayesian-Means-Edible-Weight/PRIORS/summary_Global warming potential_Mass-allocation_ON-FARM-TAXA-LEVEL-WEIGHTED.csv"))
df <- df %>%
  mutate(source = "on-farm", stressor = "GHG", weight = "edible", allocation = "mass") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "on_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)

# OFF-FARM GHG
tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-08/RERUN-GHG-and-Water-Bayesian-Means-Edible-Weight/PRIORS/summary_Global warming potential_Mass-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>% 
  mutate(source = "off-farm", stressor = "GHG", weight = "edible", allocation = "mass") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "off_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df <- df %>%
  bind_rows(tmp)

# ON-FARM NITROGEN
tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Nitrogen/summary_Marine eutrophication_Mass-allocation_ON-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "on-farm", stressor = "N", weight = "edible", allocation = "mass") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "on_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df <- df %>%
  bind_rows(tmp) 

# OFF-FARM NITROGEN
tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Nitrogen/summary_Marine eutrophication_Mass-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "off-farm", stressor = "N", weight = "edible", allocation = "mass") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "off_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df <- df %>%
  bind_rows(tmp) 

# ON-FARM PHOSPHORUS
tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Phosphorus/summary_Freshwater eutrophication_Mass-allocation_ON-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "on-farm", stressor = "P", weight = "edible", allocation = "mass") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "on_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df <- df %>%
  bind_rows(tmp) 

# OFF-FARM PHOSPHORUS
tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Phosphorus/summary_Freshwater eutrophication_Mass-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "off-farm", stressor = "P", weight = "edible", allocation = "mass") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "off_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df <- df %>%
  bind_rows(tmp) 

# ON-FARM WATER
tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-08/RERUN-GHG-and-Water-Bayesian-Means-Edible-Weight/PRIORS/summary_Water Consumption_Mass-allocation_ON-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "on-farm", stressor = "Water", weight = "edible", allocation = "mass") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "on_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df <- df %>%
  bind_rows(tmp) 

# OFF-FARM WATER
tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-08/RERUN-GHG-and-Water-Bayesian-Means-Edible-Weight/PRIORS/summary_Water Consumption_Mass-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "off-farm", stressor = "Water", weight = "edible", allocation = "mass") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "off_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df <- df %>%
  bind_rows(tmp) 

# ON-FARM LAND
tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Land-FCR-Priors-Only-Rerun/summary_Land Use_Mass-allocation_ON-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "on-farm", stressor = "Land", weight = "edible", allocation = "mass") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "on_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df <- df %>%
  bind_rows(tmp) 

# OFF-FARM LAND
tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Land-FCR-Priors-Only-Rerun/summary_Land Use_Mass-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "off-farm", stressor = "Land", weight = "edible", allocation = "mass") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "off_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df <- df %>%
  bind_rows(tmp) 

# CAPTURE
# Capture fishery results (edible weight)
df_capture_edible <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Wild/summary_WILD-GHG-TAXA-LEVEL-WEIGHTED.csv"))

df_capture_edible <- df_capture_edible %>%
  filter(.width == 0.95)

df_capture_edible <- df_capture_edible %>% 
  mutate(stressor = "GHG", full_taxa_name = taxa) %>%
  mutate(source = "all", production = "capture", weight = "edible", allocation = "mass") %>%
  select(production, taxa, full_taxa_name, stressor, source,  "median" = "total_stressor", 
         "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)

#_______________________________________________________________________________________________________________________#
# Load SI results files and merge 
# Priors + Mass allocation + Live Weight
#_______________________________________________________________________________________________________________________#

# ON-FARM GHG
df_mass_live <- read.csv(file.path(outdir, "Nature-submitted-2021-08/RERUN-GHG-and-Water-Bayesian-Means-Live-Weight/PRIORS/summary_Global warming potential_Mass-allocation_ON-FARM-TAXA-LEVEL-WEIGHTED.csv"))
df_mass_live <- df_mass_live %>%
  mutate(source = "on-farm", stressor = "GHG", weight = "live", allocation = "mass") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "on_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)

# OFF-FARM GHG
tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-08/RERUN-GHG-and-Water-Bayesian-Means-Live-Weight/PRIORS/summary_Global warming potential_Mass-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>% 
  mutate(source = "off-farm", stressor = "GHG", weight = "live", allocation = "mass") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "off_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_mass_live <- df_mass_live %>%
  bind_rows(tmp)

# ON-FARM NITROGEN
tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/Nitrogen/summary_Marine eutrophication_Mass-allocation_ON-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "on-farm", stressor = "N", weight = "live", allocation = "mass") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "on_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_mass_live <- df_mass_live %>%
  bind_rows(tmp) 

# OFF-FARM NITROGEN
tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/Nitrogen/summary_Marine eutrophication_Mass-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "off-farm", stressor = "N", weight = "live", allocation = "mass") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "off_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_mass_live <- df_mass_live %>%
  bind_rows(tmp) 

# ON-FARM PHOSPHORUS
tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/Phosphorus/summary_Freshwater eutrophication_Mass-allocation_ON-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "on-farm", stressor = "P", weight = "live", allocation = "mass") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "on_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_mass_live <- df_mass_live %>%
  bind_rows(tmp) 

# OFF-FARM PHOSPHORUS
tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/Phosphorus/summary_Freshwater eutrophication_Mass-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "off-farm", stressor = "P", weight = "live", allocation = "mass") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "off_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_mass_live <- df_mass_live %>%
  bind_rows(tmp) 

# ON-FARM WATER
tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-08/RERUN-GHG-and-Water-Bayesian-Means-Live-Weight/PRIORS/summary_Water Consumption_Mass-allocation_ON-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "on-farm", stressor = "Water", weight = "live", allocation = "mass") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "on_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_mass_live <- df_mass_live %>%
  bind_rows(tmp) 

# OFF-FARM WATER
tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-08/RERUN-GHG-and-Water-Bayesian-Means-Live-Weight/PRIORS/summary_Water Consumption_Mass-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "off-farm", stressor = "Water", weight = "live", allocation = "mass") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "off_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_mass_live <- df_mass_live %>%
  bind_rows(tmp) 

# ON-FARM LAND
tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/Land-FCR-Priors-Only/summary_Land Use_Mass-allocation_ON-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "on-farm", stressor = "Land", weight = "live", allocation = "mass") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "on_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_mass_live <- df_mass_live %>%
  bind_rows(tmp) 

# OFF-FARM LAND
tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/Land-FCR-Priors-Only/summary_Land Use_Mass-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "off-farm", stressor = "Land", weight = "live", allocation = "mass") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "off_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_mass_live <- df_mass_live %>%
  bind_rows(tmp) 

# Capture fishery results (live weight)
df_capture_live <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/Wild/summary_WILD-GHG-TAXA-LEVEL-WEIGHTED.csv"))

df_capture_live <- df_capture_live %>%
  filter(.width == 0.95)

df_capture_live <- df_capture_live %>% 
  mutate(stressor = "GHG", full_taxa_name = taxa) %>%
  mutate(source = "all", production = "capture", weight = "live", allocation = "mass") %>%
  select(production, taxa, full_taxa_name, stressor, source,  "median" = "total_stressor", 
         "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)

#_______________________________________________________________________________________________________________________#
# Load SI results files and merge 
# Priors + Economic allocation + Edible Weight
#_______________________________________________________________________________________________________________________#

# ON-FARM GHG
df_economic_edible <- read.csv(file.path(outdir, "Nature-submitted-2021-08/RERUN-GHG-and-Water-Bayesian-Means-Edible-Weight/PRIORS/summary_Global warming potential_Economic-allocation_ON-FARM-TAXA-LEVEL-WEIGHTED.csv"))
df_economic_edible <- df_economic_edible %>%
  mutate(source = "on-farm", stressor = "GHG", weight = "edible", allocation = "economic") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "on_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)

# OFF-FARM GHG
tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-08/RERUN-GHG-and-Water-Bayesian-Means-Edible-Weight/PRIORS/summary_Global warming potential_Economic-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>% 
  mutate(source = "off-farm", stressor = "GHG", weight = "edible", allocation = "economic") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "off_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_economic_edible <- df_economic_edible %>%
  bind_rows(tmp)

# ON-FARM NITROGEN
tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Nitrogen/summary_Marine eutrophication_Economic-allocation_ON-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "on-farm", stressor = "N", weight = "edible", allocation = "economic") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "on_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_economic_edible <- df_economic_edible %>%
  bind_rows(tmp) 

# OFF-FARM NITROGEN
tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Nitrogen/summary_Marine eutrophication_Economic-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "off-farm", stressor = "N", weight = "edible", allocation = "economic") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "off_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_economic_edible <- df_economic_edible %>%
  bind_rows(tmp) 

# ON-FARM PHOSPHORUS
tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Phosphorus/summary_Freshwater eutrophication_Economic-allocation_ON-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "on-farm", stressor = "P", weight = "edible", allocation = "economic") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "on_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_economic_edible <- df_economic_edible %>%
  bind_rows(tmp) 

# OFF-FARM PHOSPHORUS
tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Phosphorus/summary_Freshwater eutrophication_Economic-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "off-farm", stressor = "P", weight = "edible", allocation = "economic") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "off_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_economic_edible <- df_economic_edible %>%
  bind_rows(tmp) 

# ON-FARM WATER
tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-08/RERUN-GHG-and-Water-Bayesian-Means-Edible-Weight/PRIORS/summary_Water Consumption_Economic-allocation_ON-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "on-farm", stressor = "Water", weight = "edible", allocation = "economic") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "on_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_economic_edible <- df_economic_edible %>%
  bind_rows(tmp) 

# OFF-FARM WATER
tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-08/RERUN-GHG-and-Water-Bayesian-Means-Edible-Weight/PRIORS/summary_Water Consumption_Economic-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "off-farm", stressor = "Water", weight = "edible", allocation = "economic") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "off_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_economic_edible <- df_economic_edible %>%
  bind_rows(tmp) 

# ON-FARM LAND
tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Land-FCR-Priors-Only-Rerun/summary_Land Use_Economic-allocation_ON-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "on-farm", stressor = "Land", weight = "edible", allocation = "economic") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "on_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_economic_edible <- df_economic_edible %>%
  bind_rows(tmp) 

# OFF-FARM LAND
tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Land-FCR-Priors-Only-Rerun/summary_Land Use_Economic-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "off-farm", stressor = "Land", weight = "edible", allocation = "economic") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "off_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_economic_edible <- df_economic_edible %>%
  bind_rows(tmp) 


#_______________________________________________________________________________________________________________________#
# Load SI results files and merge 
# Priors + Economic allocation + Live Weight
#_______________________________________________________________________________________________________________________#

# GHG ON-FARM
df_economic_live <- read.csv(file.path(outdir, "Nature-submitted-2021-08/RERUN-GHG-and-Water-Bayesian-Means-Live-Weight/PRIORS/summary_Global warming potential_Economic-allocation_ON-FARM-TAXA-LEVEL-WEIGHTED.csv"))
df_economic_live <- df_economic_live %>%
  mutate(source = "on-farm", stressor = "GHG", weight = "live", allocation = "economic") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "on_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)

# GHG OFF-FARM
tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-08/RERUN-GHG-and-Water-Bayesian-Means-Live-Weight/PRIORS/summary_Global warming potential_Economic-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>% 
  mutate(source = "off-farm", stressor = "GHG", weight = "live", allocation = "economic") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "off_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_economic_live <- df_economic_live %>%
  bind_rows(tmp)

# NITROGEN ON-FARM
tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/Nitrogen/summary_Marine eutrophication_Economic-allocation_ON-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "on-farm", stressor = "N", weight = "live", allocation = "economic") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "on_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_economic_live <- df_economic_live %>%
  bind_rows(tmp) 

# NITROGEN OFF-FARM
tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/Nitrogen/summary_Marine eutrophication_Economic-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "off-farm", stressor = "N", weight = "live", allocation = "economic") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "off_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_economic_live <- df_economic_live %>%
  bind_rows(tmp) 

# PHOSPHORUS ON-FARM
tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/Phosphorus/summary_Freshwater eutrophication_Economic-allocation_ON-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "on-farm", stressor = "P", weight = "live", allocation = "economic") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "on_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_economic_live <- df_economic_live %>%
  bind_rows(tmp) 

# PHOSPHORUS OFF-FARM
tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/Phosphorus/summary_Freshwater eutrophication_Economic-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "off-farm", stressor = "P", weight = "live", allocation = "economic") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "off_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_economic_live <- df_economic_live %>%
  bind_rows(tmp) 

# WATER ON-FARM
tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-08/RERUN-GHG-and-Water-Bayesian-Means-Live-Weight/PRIORS/summary_Water Consumption_Economic-allocation_ON-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "on-farm", stressor = "Water", weight = "live", allocation = "economic") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "on_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_economic_live <- df_economic_live %>%
  bind_rows(tmp) 

# WATER OFF-FARM
tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-08/RERUN-GHG-and-Water-Bayesian-Means-Live-Weight/PRIORS/summary_Water Consumption_Economic-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "off-farm", stressor = "Water", weight = "live", allocation = "economic") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "off_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_economic_live <- df_economic_live %>%
  bind_rows(tmp) 

# LAND ON-FARM
tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/Land-FCR-Priors-Only/summary_Land Use_Economic-allocation_ON-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "on-farm", stressor = "Land", weight = "live", allocation = "economic") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "on_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_economic_live <- df_economic_live %>%
  bind_rows(tmp) 

# LAND OFF-FARM
tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/Land-FCR-Priors-Only/summary_Land Use_Economic-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "off-farm", stressor = "Land", weight = "live", allocation = "economic") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "off_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_economic_live <- df_economic_live %>%
  bind_rows(tmp) 

#_______________________________________________________________________________________________________________________#
# Write out results tables
#_______________________________________________________________________________________________________________________#
  
df_all <- df %>%
  bind_rows(df_mass_live, df_economic_edible, df_economic_live) %>%
  mutate(production = "aquaculture") %>%
  mutate(full_taxa_name = if_else(taxa == "hypoph_carp", true = "silver/bighead", false = full_taxa_name)) %>% # Fix full taxa name for hypoph/carp
  bind_rows(df_capture_edible) %>%
  bind_rows(df_capture_live)

write.csv(df_all, file.path(outdir, "SI_stressor_results_on_off_farm.csv"), row.names = FALSE)

#_______________________________________________________________________________________________________________________#
# Combine tables for total stressors (i.e., on + off farm)
#_______________________________________________________________________________________________________________________#
#_______________________________________________________________________________________________________________________#
# Load main results files and merge (Priors + Mass allocation + Edible Weight)
#_______________________________________________________________________________________________________________________#

# TOTAL GHG
df_total <- read.csv(file.path(outdir, "Nature-submitted-2021-08/RERUN-GHG-and-Water-Bayesian-Means-Edible-Weight/PRIORS/summary_Global warming potential_Mass-allocation_TOTAL-IMPACT-TAXA-LEVEL-WEIGHTED.csv"))
df_total <- df_total %>%
  mutate(source = "total on- and off-farm", stressor = "GHG", weight = "edible", allocation = "mass") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "total_stressor", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)

# TOTAL NITROGEN
tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Nitrogen/summary_Marine eutrophication_Mass-allocation_TOTAL-IMPACT-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "total on- and off-farm", stressor = "N", weight = "edible", allocation = "mass") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "total_stressor", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_total <- df_total %>%
  bind_rows(tmp) 

# TOTAL PHOSPHORUS
tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Phosphorus/summary_Freshwater eutrophication_Mass-allocation_TOTAL-IMPACT-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "total on- and off-farm", stressor = "P", weight = "edible", allocation = "mass") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "total_stressor", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_total <- df_total %>%
  bind_rows(tmp) 

# TOTAL WATER
tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-08/RERUN-GHG-and-Water-Bayesian-Means-Edible-Weight/PRIORS/summary_Water Consumption_Mass-allocation_TOTAL-IMPACT-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "total on- and off-farm", stressor = "Water", weight = "edible", allocation = "mass") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "total_stressor", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_total <- df_total %>%
  bind_rows(tmp) 

# TOTAL LAND
tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Land-FCR-Priors-Only-Rerun/summary_Land Use_Mass-allocation_TOTAL-IMPACT-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "total on- and off-farm", stressor = "Land", weight = "edible", allocation = "mass") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "total_stressor", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_total <- df_total %>%
  bind_rows(tmp) 

#_______________________________________________________________________________________________________________________#
# Load SI results files and merge 
# Priors + Mass allocation + Live Weight
#_______________________________________________________________________________________________________________________#

# TOTAL GHG
df_total_mass_live <- read.csv(file.path(outdir, "Nature-submitted-2021-08/RERUN-GHG-and-Water-Bayesian-Means-Live-Weight/PRIORS/summary_Global warming potential_Mass-allocation_TOTAL-IMPACT-TAXA-LEVEL-WEIGHTED.csv"))
df_total_mass_live <- df_total_mass_live %>%
  mutate(source = "total on- and off-farm", stressor = "GHG", weight = "live", allocation = "mass") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "total_stressor", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)

# TOTAL NITROGEN
tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/Nitrogen/summary_Marine eutrophication_Mass-allocation_TOTAL-IMPACT-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "total on- and off-farm", stressor = "N", weight = "live", allocation = "mass") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "total_stressor", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_total_mass_live <- df_total_mass_live %>%
  bind_rows(tmp) 

# TOTAL PHOSPHORUS
tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/Phosphorus/summary_Freshwater eutrophication_Mass-allocation_TOTAL-IMPACT-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "total on- and off-farm", stressor = "P", weight = "live", allocation = "mass") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "total_stressor", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_total_mass_live <- df_total_mass_live %>%
  bind_rows(tmp) 

# TOTAL WATER
tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-08/RERUN-GHG-and-Water-Bayesian-Means-Live-Weight/PRIORS/summary_Water Consumption_Mass-allocation_TOTAL-IMPACT-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "total on- and off-farm", stressor = "Water", weight = "live", allocation = "mass") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "total_stressor", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_total_mass_live <- df_total_mass_live %>%
  bind_rows(tmp) 

# TOTAL LAND
tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/Land-FCR-Priors-Only/summary_Land Use_Mass-allocation_TOTAL-IMPACT-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "total on- and off-farm", stressor = "Land", weight = "live", allocation = "mass") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "total_stressor", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_total_mass_live <- df_total_mass_live %>%
  bind_rows(tmp) 

#_______________________________________________________________________________________________________________________#
# Load SI results files and merge 
# Priors + Economic allocation + Edible Weight
#_______________________________________________________________________________________________________________________#

# TOTAL GHG
df_total_economic_edible <- read.csv(file.path(outdir, "Nature-submitted-2021-08/RERUN-GHG-and-Water-Bayesian-Means-Edible-Weight/PRIORS/summary_Global warming potential_Economic-allocation_TOTAL-IMPACT-TAXA-LEVEL-WEIGHTED.csv"))
df_total_economic_edible <- df_total_economic_edible %>%
  mutate(source = "total on- and off-farm", stressor = "GHG", weight = "edible", allocation = "economic") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "total_stressor", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)

# TOTAL NITROGEN
tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Nitrogen/summary_Marine eutrophication_Economic-allocation_TOTAL-IMPACT-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "total on- and off-farm", stressor = "N", weight = "edible", allocation = "economic") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "total_stressor", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_total_economic_edible <- df_total_economic_edible %>%
  bind_rows(tmp) 

# TOTAL PHOSPHORUS
tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Phosphorus/summary_Freshwater eutrophication_Economic-allocation_TOTAL-IMPACT-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "total on- and off-farm", stressor = "P", weight = "edible", allocation = "economic") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "total_stressor", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_total_economic_edible <- df_total_economic_edible %>%
  bind_rows(tmp) 

# TOTAL WATER
tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-08/RERUN-GHG-and-Water-Bayesian-Means-Edible-Weight/PRIORS/summary_Water Consumption_Economic-allocation_TOTAL-IMPACT-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "total on- and off-farm", stressor = "Water", weight = "edible", allocation = "economic") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "total_stressor", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_total_economic_edible <- df_total_economic_edible %>%
  bind_rows(tmp) 

# TOTAL LAND
tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Land-FCR-Priors-Only-Rerun/summary_Land Use_Economic-allocation_TOTAL-IMPACT-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "total on- and off-farm", stressor = "Land", weight = "edible", allocation = "economic") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "total_stressor", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_total_economic_edible <- df_total_economic_edible %>%
  bind_rows(tmp) 

#_______________________________________________________________________________________________________________________#
# Load SI results files and merge 
# Priors + Economic allocation + Live Weight
#_______________________________________________________________________________________________________________________#

# TOTAL GHG
df_total_economic_live <- read.csv(file.path(outdir, "Nature-submitted-2021-08/RERUN-GHG-and-Water-Bayesian-Means-Live-Weight/PRIORS/summary_Global warming potential_Economic-allocation_TOTAL-IMPACT-TAXA-LEVEL-WEIGHTED.csv"))
df_total_economic_live <- df_total_economic_live %>%
  mutate(source = "total on- and off-farm", stressor = "GHG", weight = "live", allocation = "economic") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "total_stressor", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)

# TOTAL NITROGEN
tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/Nitrogen/summary_Marine eutrophication_Economic-allocation_TOTAL-IMPACT-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "total on- and off-farm", stressor = "N", weight = "live", allocation = "economic") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "total_stressor", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_total_economic_live <- df_total_economic_live %>%
  bind_rows(tmp) 

# TOTAL PHOSPHORUS
tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/Phosphorus/summary_Freshwater eutrophication_Economic-allocation_TOTAL-IMPACT-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "total on- and off-farm", stressor = "P", weight = "live", allocation = "economic") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "total_stressor", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_total_economic_live <- df_total_economic_live %>%
  bind_rows(tmp) 

# TOTAL WATER
tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-08/RERUN-GHG-and-Water-Bayesian-Means-Live-Weight/PRIORS/summary_Water Consumption_Economic-allocation_TOTAL-IMPACT-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "total on- and off-farm", stressor = "Water", weight = "live", allocation = "economic") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "total_stressor", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_total_economic_live <- df_total_economic_live %>%
  bind_rows(tmp) 

# TOTAL LAND
tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/Land-FCR-Priors-Only/summary_Land Use_Economic-allocation_TOTAL-IMPACT-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "total on- and off-farm", stressor = "Land", weight = "live", allocation = "economic") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "total_stressor", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_total_economic_live <- df_total_economic_live %>%
  bind_rows(tmp) 

#_______________________________________________________________________________________________________________________#
# Write out results for total stressors
#_______________________________________________________________________________________________________________________#

df_total_all <- df_total %>%
  bind_rows(df_total_mass_live, df_total_economic_edible, df_total_economic_live) %>%
  mutate(production = "aquaculture") %>%
  mutate(full_taxa_name = if_else(taxa == "hypoph_carp", true = "silver/bighead", false = full_taxa_name)) %>% # Fix full taxa name for hypoph/carp
  bind_rows(df_capture_live) %>% # capture stressor source is designated as "all"
  bind_rows(df_capture_edible)

write.csv(df_total_all, file.path(outdir, "SI_stressor_results_total.csv"), row.names = FALSE)

#_______________________________________________________________________________________________________________________#
# Summarize patterns for paper text
#_______________________________________________________________________________________________________________________#

# Load totals from above summary files
df_total <- read.csv("Data/Table S11.csv")
df_total_mass_lw <- df_total %>%
  filter(weight == "edible", allocation == "mass")

df_source <- read.csv("Data/SI_stressor_results_on_off_farm.csv")
df_source_mass_lw <- df_source %>%
  filter(weight == "edible", allocation == "mass")  

# Lowest taxa across stressors
df_total_mass_lw %>% 
  group_by(stressor) %>%
  slice_min(median, n = 4)

df_total_mass_lw %>% 
  group_by(stressor) %>%
  slice_max(median, n = 3)

# Stressor correlations 
tmp <- df_total_mass_lw %>% 
  filter(production == "aquaculture") %>%
  select(taxa, stressor, median) %>%
  pivot_wider(names_from = stressor, values_from = median)
cor(tmp[,2:6])

# Percent of on- versus off-farm 
source_percent <- df_source_mass_lw %>%
  select(taxa, full_taxa_name, stressor, source, median) %>%
  pivot_wider(names_from = source, values_from = median) %>%
  mutate(total = `on-farm` + `off-farm`) %>%
  mutate(percent_onfarm = 100*(`on-farm`/total),
         percent_offfarm = 100*(`off-farm`/total))

# GHG results
source_percent %>% 
  filter(stressor == "GHG") %>%
  arrange(percent_onfarm)

source_percent %>% 
  filter(stressor == "GHG") %>%
  arrange(`on-farm`)

source_percent %>% 
  filter(stressor == "GHG") %>%
  arrange(total)

df_total_mass_lw %>%
  filter(taxa %in% c("bivalves", "Bivalves"), stressor == "GHG")

df_total_mass_lw %>%
  filter(taxa %in% c("shrimp", "Shrimps"), stressor == "GHG")

df_total_mass_lw %>%
  filter(taxa %in% c("salmon", "trout", "Salmonids"), stressor == "GHG")


# Land results
source_percent %>% 
  filter(stressor == "Land") %>%
  arrange(`on-farm`)

source_percent %>% 
  filter(stressor == "Land") %>%
  arrange(`off-farm`)

source_percent %>% 
  filter(stressor == "Land") %>%
  arrange(percent_offfarm)

source_percent %>% 
  filter(stressor == "Land") %>%
  arrange(desc(total))

# Water results
source_percent %>% 
  filter(stressor == "Water") %>%
  arrange(percent_offfarm)

source_percent %>% 
  filter(stressor == "Water") %>%
  arrange(`off-farm`)

source_percent %>% 
  filter(stressor == "Water") %>%
  arrange(total)

# N and P results
source_percent %>% 
  filter(stressor == "N") %>%
  arrange(percent_offfarm)

source_percent %>% 
  filter(stressor == "P") %>%
  arrange(percent_offfarm)

source_percent %>% 
  filter(stressor == "N") %>%
  arrange(desc(total))

source_percent %>% 
  filter(stressor == "P") %>%
  arrange(desc(total))

#_______________________________________________________________________________________________________________________#
# Stacked bar of impact by source
#_______________________________________________________________________________________________________________________#

df_plot <- df %>%
  mutate(plot_taxa_name = case_when(taxa == "hypoph_carp" ~ "silver/bighead",
                                    taxa == "misc_diad" ~ "misc diad",
                                    taxa == "misc_marine" ~ "misc marine",
                                    taxa == "oth_carp" ~ "misc carp",
                                    taxa == "plants" ~ "seaweeds",
                                    TRUE ~ taxa),
         plot_taxa_name = as.factor(plot_taxa_name)) # Create new column for plotting taxa names

plot_taxa_name_order <- c("seaweeds",
                          "bivalves",
                          "silver/bighead",
                          "salmon",
                          "trout",
                          "misc carp",
                          "catfish",
                          "milkfish",
                          "shrimp",
                          "tilapia",
                          "misc marine",
                          "misc diad")

weight_type <- "edible weight"
units_for_ghg <- bquote(atop('kg'~CO[2]*' t'^-1~phantom(), .(weight_type))) # atop + phantom() to create line break
units_for_land <- bquote(atop('m'^2*' t'^-1~phantom(), .(weight_type)))
units_for_nitrogen <- bquote(atop('kg N-eq t'^-1~phantom(), .(weight_type)))
units_for_phosphorus <- bquote(atop('kg P-eq t'^-1~phantom(), .(weight_type)))
units_for_water <- bquote(atop('m'^3*' t'^-1~phantom(), .(weight_type)))

# df_plot$taxa <- factor(df_plot$taxa, levels = rev(c("misc_diad", 
#                                                     "misc_marine", 
#                                                     "tilapia", 
#                                                     "shrimp", 
#                                                     "milkfish", 
#                                                     "catfish", 
#                                                     "oth_carp", 
#                                                     "trout", 
#                                                     "salmon", 
#                                                     "hypoph_carp", 
#                                                     "bivalves", 
#                                                     "plants")))

df_plot$plot_taxa_name <- factor(df_plot$plot_taxa_name, levels = plot_taxa_name_order)

base_size <- 12
base_family <- "sans"

# GHG plot
ghg_plot <- ggplot(df_plot %>% filter(stressor == "GHG"), 
                   aes(x = median, y = plot_taxa_name, fill = source)) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(values = c("#70468C", "#57D182")) +
  labs(title = "GHG", x = units_for_ghg, y = "") +
  theme(axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"), 
        axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"), 
        axis.text = element_text(size = ceiling(base_size*0.7), colour = "black"),
        axis.title = element_text(size = ceiling(base_size*0.8)), 
        axis.text.x = element_text(angle = 90),
        panel.grid.minor = element_blank(), 
        panel.grid.major.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        panel.background = element_blank(), panel.border = element_blank(),
        strip.background = element_rect(linetype = 1, fill = "white"), 
        strip.text = element_text(), 
        strip.text.x = element_text(vjust = 0.5), 
        strip.text.y = element_text(angle = -90), 
        strip.placement.y = "outside",
        legend.position="bottom",
        plot.title = element_text(size = ceiling(base_size*1.1), face = "bold"), 
        plot.subtitle = element_text(size = ceiling(base_size*1.05)))

# N plot
N_plot <- ggplot(df_plot %>% filter(stressor == "N"), 
                 aes(x = median, y = plot_taxa_name, fill = source)) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(values = c("#70468C", "#57D182")) +
  labs(title = "Nitrogen", x = units_for_nitrogen, y = "") +
  theme(axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"), 
        axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"), 
        axis.text = element_blank(),
        axis.title = element_text(size = ceiling(base_size*0.8)), 
        axis.text.x = element_text(angle = 90),
        panel.grid.minor = element_blank(), 
        panel.grid.major.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        panel.background = element_blank(), panel.border = element_blank(),
        strip.background = element_rect(linetype = 1, fill = "white"), 
        strip.text = element_text(), 
        strip.text.x = element_text(vjust = 0.5), 
        strip.text.y = element_text(angle = -90), 
        strip.placement.y = "outside",
        legend.position="bottom",
        plot.title = element_text(size = ceiling(base_size*1.1), face = "bold"), 
        plot.subtitle = element_text(size = ceiling(base_size*1.05)))

# P plot
P_plot <- ggplot(df_plot %>% filter(stressor == "P"), 
                 aes(x = median, y = plot_taxa_name, fill = source)) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(values = c("#70468C", "#57D182")) +
  labs(title = "Phosphorus", x = units_for_phosphorus, y = "") +
  theme(axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"), 
        axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"), 
        axis.text = element_blank(),
        axis.title = element_text(size = ceiling(base_size*0.8)), 
        axis.text.x = element_text(angle = 90),
        panel.grid.minor = element_blank(), 
        panel.grid.major.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        panel.background = element_blank(), panel.border = element_blank(),
        strip.background = element_rect(linetype = 1, fill = "white"), 
        strip.text = element_text(), 
        strip.text.x = element_text(vjust = 0.5), 
        strip.text.y = element_text(angle = -90), 
        strip.placement.y = "outside",
        legend.position="bottom",
        plot.title = element_text(size = ceiling(base_size*1.1), face = "bold"), 
        plot.subtitle = element_text(size = ceiling(base_size*1.05)))

# Land plot
land_plot <- ggplot(df_plot %>% filter(stressor == "Land"), 
                    aes(x = median, y = plot_taxa_name, fill = source)) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(values = c("#70468C", "#57D182")) +
  labs(title = "Land", x = units_for_land, y = "") +
  theme(axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"), 
        axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"), 
        axis.text = element_blank(),
        axis.title = element_text(size = ceiling(base_size*0.8)), 
        axis.text.x = element_text(angle = 90),
        panel.grid.minor = element_blank(), 
        panel.grid.major.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        panel.background = element_blank(), panel.border = element_blank(),
        strip.background = element_rect(linetype = 1, fill = "white"), 
        strip.text = element_text(), 
        strip.text.x = element_text(vjust = 0.5), 
        strip.text.y = element_text(angle = -90), 
        strip.placement.y = "outside",
        legend.position="bottom",
        plot.title = element_text(size = ceiling(base_size*1.1), face = "bold"), 
        plot.subtitle = element_text(size = ceiling(base_size*1.05)))

# Water plot
water_plot <- ggplot(df_plot %>% filter(stressor == "Water"), 
                     aes(x = median, y = plot_taxa_name, fill = source)) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(values = c("#70468C", "#57D182")) +
  labs(title = "Water", x = units_for_water, y = "") +
  theme(axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"), 
        axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"), 
        axis.text = element_blank(),
        axis.title = element_text(size = ceiling(base_size*0.8)), 
        axis.text.x = element_text(angle = 90),
        panel.grid.minor = element_blank(), 
        panel.grid.major.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        panel.background = element_blank(), panel.border = element_blank(),
        strip.background = element_rect(linetype = 1, fill = "white"), 
        strip.text = element_text(), 
        strip.text.x = element_text(vjust = 0.5), 
        strip.text.y = element_text(angle = -90), 
        strip.placement.y = "outside",
        legend.position="bottom",
        plot.title = element_text(size = ceiling(base_size*1.1), face = "bold"), 
        plot.subtitle = element_text(size = ceiling(base_size*1.05)))

png(file.path(outdir, "plot_stressor_by_source.png"), width = 8.5, height = 5, units = "in", res = 300)
ggarrange(ghg_plot, N_plot, P_plot, water_plot, land_plot, nrow = 1,
          common.legend = TRUE, legend = "bottom", align = "h", widths = c(1.1, 0.7, 0.7, 0.7, 0.7))
dev.off()
