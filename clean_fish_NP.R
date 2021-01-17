# Clean fish N and P data
# Author: Alon Shepon and Jessica Gephart

#------------------------------------------------------------------------------------
#--------------------------------------load packages---------------------------------
#------------------------------------------------------------------------------------
library("rfishbase")  #package to use fishbase in r
library("openxlsx")
library("rlang")
library("ggplot2")
library("dplyr")
library("tidyverse")
library("data.table")
library("stringr")
library("openxlsx")

source("Functions.R")

# Set data directory
datadir <- "/Volumes/jgephart/BFA Environment 2/Data"

# If not loaded, load and clean lca data
lca <- read.csv(file.path(datadir, "lca_clean_with_groups.csv"))
# Select species in lca database
lca_species <- lca %>% 
  select(clean_sci_name, Common.Name) %>%
  distinct()

#------------------------------------------------------------------------------------
#------------------------------------------functions---------------------------------
#------------------------------------------------------------------------------------

#first method: estimate from C content, Fig 1 of the above reference, linear regression values. Assuming metabolizable energy density is equal throughout whole fish

fishN_via_C <- function(Water,Energy.total.metabolizable.calculated.from.the.energy.producing.food.components.original.as.from.source.kcal){
  DM = (100 - Water) / 100     #convert to dry matter fraction
  
  #In aquacultured fish, whole-body C content was estimated from energy content using a converting factor of 1 g C = 11.4 kcal
  fish_C_percent <- 1/DM*Energy.total.metabolizable.calculated.from.the.energy.producing.food.components.original.as.from.source.kcal/11.4
  N_slope <- -0.175   #from fig 1
  N_intercept <- 18.506  #intercept (personal communication with authors)
  fish_N <- fish_C_percent * N_slope + N_intercept #in % of DM
  return(fish_N)}

fishP_via_C <- function(Water,Energy.total.metabolizable.calculated.from.the.energy.producing.food.components.original.as.from.source.kcal){
  
  DM = (100 - Water) / 100     #convert to dry matter fraction
  
  #In aquacultured fish, whole-body C content was estimated from energy content using a converting factor of 1 g C = 11.4 kcal
  fish_C_percent <- 1/DM*Energy.total.metabolizable.calculated.from.the.energy.producing.food.components.original.as.from.source.kcal/11.4
  P_slope <- -0.083  #from fig 1
  P_intercept <- 6.311  #intercept (personal communication with authors)
  fish_P <- fish_C_percent * P_slope + P_intercept #in % of DM
  return(fish_P)}

#second method: estimate from total lipids, figure 2, assuming lipid density is equal throughout whole body

fishN_viaFat <- function(Water,Fat.total){
  DM = (100 - Water) / 100     #Dry Matter percentage
  #Fat.total in units of g to 100 g
  fat = Fat.total * 1  / DM       #conversion to DM
  
  #fat percentage in whole body in % of DM. 
  fish_C_via_fat = fat * 0.31 + 38  #linear regression from fig 2
  fish_N_via_fat = fat * -0.08 + 12 #linear regression from fig 2
  return(fish_N_via_fat)}

fishP_viaFat <- function(Water,Fat.total){
  DM = (100 - Water) / 100     #Dry Matter percentage
  
  #Fat.total in units of g to 100 g
  fat = Fat.total * 1  / DM      #conversion to DM
  
  #fat percentage in whole body in % of DM. Assuming Fat is contained only (mostly) in the edible portion.
  fish_C_via_fat = fat * 0.31 + 38  #linear regression from fig 2
  fish_P_via_fat = fat * -0.04 + 3.2 #linear regression from fig 2
  return(fish_P_via_fat)}

#------------------------------------------------------------------------------------
#------------------------------------------load data---------------------------------
#------------------------------------------------------------------------------------

#calculate C,N,P content of fish (based on energy and fat content ) and compute discharge based on Czamanski et al 2011, Marine Biology,
#Carbon, nitrogen and phosphorus elemental stoichiometry in aquacultured and wild-caught fish and consequences
#for pelagic nutrient dynamics
#And also a simplified NPZ method approach

#N and P content of fish in % of DM
#read Zach's data and make calculations using built functions

fishNutrition <- read.csv(file.path(datadir, "AFCD_live.csv"), stringsAsFactors = FALSE)

#build taxa data 
fish_Nutrition_genus<-fishNutrition[,c("species","order","family","genus")] %>% 
  distinct()

fishNutrition1<-fishNutrition %>% mutate(N_C=fishN_via_C(Water,Energy.total.metabolizable.calculated.from.the.energy.producing.food.components.original.as.from.source.kcal),
                                         P_C=fishP_via_C(Water,Energy.total.metabolizable.calculated.from.the.energy.producing.food.components.original.as.from.source.kcal),
                                         N_fat=fishN_viaFat(Water,Fat.total),
                                         P_fat=fishP_viaFat(Water,Fat.total),
                                         N_built_in=Nitrogen.total*100/(100-Water),# Conversion to units of percentages in DM. Nitrogen.total is in units of g/100 edible gram; 
                                         P_built_in=Phosphorus/1000*100/(100-Water)) # Conversion to units of percentages in DM. Phosphorous in the dataset is in mg/100 edible gram
fishNutrition1$N_built_in[which(fishNutrition1$N_built_in>5)]=NaN   #remove large unrealstic numbers

fishNutrition<-fishNutrition1 %>% group_by(species, genus, family, order) %>% 
  summarise(P_byC=mean(P_C,na.rm = TRUE),N_byC=mean(N_C,na.rm = TRUE),
            P_fat_1=mean(P_fat,na.rm = TRUE),N_fat_1=mean(N_fat,na.rm = TRUE),
            P_built_in1=mean(P_built_in,na.rm = TRUE),N_built_in1=mean(N_built_in,na.rm = TRUE), 
            Water_content = mean(Water, na.rm = TRUE)) %>%                                 
  rowwise() %>%                                                    
  mutate(N_avg=mean(c(N_byC,N_fat_1,N_built_in1),na.rm = TRUE)) %>%
  mutate(P_avg=mean(c(P_byC,P_fat_1,P_built_in1),na.rm = TRUE)) %>%
  select(c("Scientific.Name" =  "species", "genus", "family", "order", "N_avg","P_avg", "Water_content")) %>%
  filter(!is.na(N_avg)) %>%
  filter(!is.na(P_avg)) %>%
  # Remove outlier value
  filter(Scientific.Name != "Macrobrachium lanchesteri")


#------------------------------------------------------------------------------------
# Add rows for missing taxa to complete merge
#------------------------------------------------------------------------------------
# FIX IT: Currently no matches for seaweeds - "Gracilaria chilensis", "Laminaria digitata", "Macrocystis pyrifera", "Saccharina latissima"

## No species, genus, family or order info for Chanos chanos. Currently using average across all species. FIX IT :Ask Chris
add_fishNutrition <- fishNutrition %>%
  ungroup() %>%
  summarise_all(mean, na.rm = TRUE)
colnames(add_fishNutrition) <- gsub(colnames(add_fishNutrition), pattern = "_1", replacement = "")
add_fishNutrition$Scientific.Name <- "Chanos chanos"
add_fishNutrition$genus <- ""
add_fishNutrition$order <- ""
add_fishNutrition$family <- ""
add_fishNutrition$Common.Name <- ""

fishNutrition <- fishNutrition %>%
  bind_rows(add_fishNutrition)

## Use Clarias genus average for "Clarias gariepinus"
add_fishNutrition <- fishNutrition %>%
  filter(genus == "Clarias") %>%
  ungroup() %>%
  summarise_all(mean, na.rm = TRUE)
colnames(add_fishNutrition) <- gsub(colnames(add_fishNutrition), pattern = "_1", replacement = "")
add_fishNutrition$Scientific.Name <- "Clarias gariepinus"
add_fishNutrition$genus <- ""
add_fishNutrition$order <- ""
add_fishNutrition$family <- ""
add_fishNutrition$Common.Name <- ""

fishNutrition <- fishNutrition %>%
  bind_rows(add_fishNutrition)

## Use Cyprinidae family average for "Ctenopharyngodon idella"
add_fishNutrition <- fishNutrition %>%
  filter(family == "Cyprinidae") %>%
  ungroup() %>%
  summarise_all(mean, na.rm = TRUE)
colnames(add_fishNutrition) <- gsub(colnames(add_fishNutrition), pattern = "_1", replacement = "")
add_fishNutrition$Scientific.Name <- "Ctenopharyngodon idella"
add_fishNutrition$genus <- ""
add_fishNutrition$order <- ""
add_fishNutrition$family <- ""
add_fishNutrition$Common.Name <- ""

fishNutrition <- fishNutrition %>%
  bind_rows(add_fishNutrition)

## Use Perciformes order average for "Cynoscion spp" FIX IT: Ask Chris if there is something else we can do here
add_fishNutrition <- fishNutrition %>%
  filter(order == "Perciformes") %>%
  ungroup() %>%
  summarise_all(mean, na.rm = TRUE)
colnames(add_fishNutrition) <- gsub(colnames(add_fishNutrition), pattern = "_1", replacement = "")
add_fishNutrition$Scientific.Name <- "Cynoscion spp"
add_fishNutrition$genus <- ""
add_fishNutrition$order <- ""
add_fishNutrition$family <- ""
add_fishNutrition$Common.Name <- ""

fishNutrition <- fishNutrition %>%
  bind_rows(add_fishNutrition)

## Use Cyprinidae family average for "Cyprinidae"
add_fishNutrition <- fishNutrition %>%
  filter(family == "Cyprinidae") %>%
  ungroup() %>%
  summarise_all(mean, na.rm = TRUE)
colnames(add_fishNutrition) <- gsub(colnames(add_fishNutrition), pattern = "_1", replacement = "")
add_fishNutrition$Scientific.Name <- "Cyprinidae"
add_fishNutrition$genus <- ""
add_fishNutrition$order <- ""
add_fishNutrition$family <- ""
add_fishNutrition$Common.Name <- ""

fishNutrition <- fishNutrition %>%
  bind_rows(add_fishNutrition)

## Use Perciformes order average for "Dicentrarchus labrax"
add_fishNutrition <- fishNutrition %>%
  filter(order == "Perciformes") %>%
  ungroup() %>%
  summarise_all(mean, na.rm = TRUE)
colnames(add_fishNutrition) <- gsub(colnames(add_fishNutrition), pattern = "_1", replacement = "")
add_fishNutrition$Scientific.Name <- "Dicentrarchus labrax"
add_fishNutrition$genus <- ""
add_fishNutrition$order <- ""
add_fishNutrition$family <- ""
add_fishNutrition$Common.Name <- ""

fishNutrition <- fishNutrition %>%
  bind_rows(add_fishNutrition)

## Use Epinephelus genus average for "Epinephelus spp"
add_fishNutrition <- fishNutrition %>%
  filter(genus == "Epinephelus") %>%
  ungroup() %>%
  summarise_all(mean, na.rm = TRUE)
colnames(add_fishNutrition) <- gsub(colnames(add_fishNutrition), pattern = "_1", replacement = "")
add_fishNutrition$Scientific.Name <- "Epinephelus spp"
add_fishNutrition$genus <- ""
add_fishNutrition$order <- ""
add_fishNutrition$family <- ""
add_fishNutrition$Common.Name <- ""

fishNutrition <- fishNutrition %>%
  bind_rows(add_fishNutrition)

## Use Perciformes order average for "Lates calcarifer"
add_fishNutrition <- fishNutrition %>%
  filter(order == "Perciformes") %>%
  ungroup() %>%
  summarise_all(mean, na.rm = TRUE)
colnames(add_fishNutrition) <- gsub(colnames(add_fishNutrition), pattern = "_1", replacement = "")
add_fishNutrition$Scientific.Name <- "Lates calcarifer"
add_fishNutrition$genus <- ""
add_fishNutrition$order <- ""
add_fishNutrition$family <- ""
add_fishNutrition$Common.Name <- ""

fishNutrition <- fishNutrition %>%
  bind_rows(add_fishNutrition)

## Use Macrobrachium species average for "Macrobrachium amazonicum"
add_fishNutrition <- fishNutrition %>%
  filter(str_detect(Scientific.Name, "Macrobrachium")) %>%
  ungroup() %>%
  summarise_all(mean, na.rm = TRUE)
colnames(add_fishNutrition) <- gsub(colnames(add_fishNutrition), pattern = "_1", replacement = "")
add_fishNutrition$Scientific.Name <- "Macrobrachium amazonicum"
add_fishNutrition$genus <- ""
add_fishNutrition$order <- ""
add_fishNutrition$family <- ""
add_fishNutrition$Common.Name <- ""

fishNutrition <- fishNutrition %>%
  bind_rows(add_fishNutrition)

## Use H. molitrix for "Mixed H. molitrix and H. nobilis"
add_fishNutrition <- fishNutrition %>%
  filter(str_detect(Scientific.Name, "molitrix")) %>%
  ungroup() %>%
  summarise_all(mean, na.rm = TRUE)
colnames(add_fishNutrition) <- gsub(colnames(add_fishNutrition), pattern = "_1", replacement = "")
add_fishNutrition$Scientific.Name <- "Mixed H. molitrix and H. nobilis"
add_fishNutrition$genus <- ""
add_fishNutrition$order <- ""
add_fishNutrition$family <- ""
add_fishNutrition$Common.Name <- ""

fishNutrition <- fishNutrition %>%
  bind_rows(add_fishNutrition)

## Use Pangasius average for "Pangasius spp"
add_fishNutrition <- fishNutrition %>%
  filter(str_detect(Scientific.Name, "Pangasius")) %>%
  ungroup() %>%
  summarise_all(mean, na.rm = TRUE)
colnames(add_fishNutrition) <- gsub(colnames(add_fishNutrition), pattern = "_1", replacement = "")
add_fishNutrition$Scientific.Name <- "Pangasius spp"
add_fishNutrition$genus <- ""
add_fishNutrition$order <- ""
add_fishNutrition$family <- ""
add_fishNutrition$Common.Name <- ""

fishNutrition <- fishNutrition %>%
  bind_rows(add_fishNutrition)

## Use Pangasius average for "Pangasianodon hypophthalmus"
add_fishNutrition$Scientific.Name <- "Pangasianodon hypophthalmus"
add_fishNutrition$genus <- ""
add_fishNutrition$order <- ""
add_fishNutrition$family <- ""
add_fishNutrition$Common.Name <- ""

fishNutrition <- fishNutrition %>%
  bind_rows(add_fishNutrition)

## Use Salvelinus, Salmo and Oncorhynchus species average for "Salmonidae"
add_fishNutrition <- fishNutrition %>%
  filter(family == "Salmonidae") %>%
  ungroup() %>%
  summarise_all(mean, na.rm = TRUE)
colnames(add_fishNutrition) <- gsub(colnames(add_fishNutrition), pattern = "_1", replacement = "")
add_fishNutrition$Scientific.Name <- "Salmonidae"
add_fishNutrition$genus <- ""
add_fishNutrition$order <- ""
add_fishNutrition$family <- ""
add_fishNutrition$Common.Name <- ""

fishNutrition <- fishNutrition %>%
  bind_rows(add_fishNutrition)

## Use Salvelinus species average for "Salvelinus alpinus"
add_fishNutrition <- fishNutrition %>%
  filter(str_detect(Scientific.Name, "Salvelinus")) %>%
  ungroup() %>%
  summarise_all(mean, na.rm = TRUE)
colnames(add_fishNutrition) <- gsub(colnames(add_fishNutrition), pattern = "_1", replacement = "")
add_fishNutrition$Scientific.Name <- "Salvelinus alpinus"
add_fishNutrition$genus <- ""
add_fishNutrition$order <- ""
add_fishNutrition$family <- ""
add_fishNutrition$Common.Name <- ""

fishNutrition <- fishNutrition %>%
  bind_rows(add_fishNutrition)

## Use Pleuronectiformes order average for "Scophthalmidae" 
add_fishNutrition <- fishNutrition %>%
  filter(order == "Pleuronectiformes") %>%
  ungroup() %>%
  summarise_all(mean, na.rm = TRUE)
colnames(add_fishNutrition) <- gsub(colnames(add_fishNutrition), pattern = "_1", replacement = "")
add_fishNutrition$Scientific.Name <- "Scophthalmidae"
add_fishNutrition$genus <- ""
add_fishNutrition$order <- ""
add_fishNutrition$family <- ""
add_fishNutrition$Common.Name <- ""

fishNutrition <- fishNutrition %>%
  bind_rows(add_fishNutrition)

## Use Seriola genus average for "Seriola rivoliana"
add_fishNutrition <- fishNutrition %>%
  filter(genus == "Seriola") %>%
  ungroup() %>%
  summarise_all(mean, na.rm = TRUE)
colnames(add_fishNutrition) <- gsub(colnames(add_fishNutrition), pattern = "_1", replacement = "")
add_fishNutrition$Scientific.Name <- "Seriola rivoliana"
add_fishNutrition$genus <- ""
add_fishNutrition$order <- ""
add_fishNutrition$family <- ""
add_fishNutrition$Common.Name <- ""

fishNutrition <- fishNutrition %>%
  bind_rows(add_fishNutrition)

## Use Sparidae family average for "Sparus aurata" 
add_fishNutrition <- fishNutrition %>%
  filter(family == "Sparidae") %>%
  ungroup() %>%
  summarise_all(mean, na.rm = TRUE)
colnames(add_fishNutrition) <- gsub(colnames(add_fishNutrition), pattern = "_1", replacement = "")
add_fishNutrition$Scientific.Name <- "Sparus aurata"
add_fishNutrition$genus <- ""
add_fishNutrition$order <- ""
add_fishNutrition$family <- ""
add_fishNutrition$Common.Name <- ""

fishNutrition <- fishNutrition %>%
  bind_rows(add_fishNutrition) %>%
  select(-Common.Name)

# Final check that all species are included
fishNutrition <- lca %>%
  select(clean_sci_name) %>%
  distinct() %>%
  left_join(fishNutrition, by = c("clean_sci_name" = "Scientific.Name"))

# Convert from % N and P dry weight to wet weight
fishNutrition <- fishNutrition %>%
  mutate(N_t_liveweight_t = (N_avg/100)*(1-(Water_content/100)),
         P_t_liveweight_t = (P_avg/100)*(1-(Water_content/100)))

# Manually add aquatic plants from Alon's calculation (use dry weight - from aquatic_plants_N_P.xls)
fishNutrition$N_t_liveweight_t[fishNutrition$clean_sci_name == "Laminaria digitata"] <- 6.85/100
fishNutrition$P_t_liveweight_t[fishNutrition$clean_sci_name == "Laminaria digitata"] <- 1.16/100

fishNutrition$N_t_liveweight_t[fishNutrition$clean_sci_name == "Saccharina latissima"] <- 3.46/100
fishNutrition$P_t_liveweight_t[fishNutrition$clean_sci_name == "Saccharina latissima"] <- 0.44/100


fishNutrition$N_t_liveweight_t[fishNutrition$clean_sci_name == "Macrocystis pyrifera"] <- 2.87/100
fishNutrition$P_t_liveweight_t[fishNutrition$clean_sci_name == "Macrocystis pyrifera"] <- 0.78/100

fishNutrition$N_t_liveweight_t[fishNutrition$clean_sci_name == "Gracilaria chilensis"] <- 2.69/100
fishNutrition$P_t_liveweight_t[fishNutrition$clean_sci_name == "Gracilaria chilensis"] <- 1.82/100

fishNutrition <- fishNutrition %>%
  select(clean_sci_name, N_t_liveweight_t, P_t_liveweight_t)

write.csv(fishNutrition, file.path(datadir, "fish_NP_clean.csv"), row.names = FALSE)
