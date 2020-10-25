# Functions

#_____________________________________________________________________________________________________#
# Clean LCA data
#_____________________________________________________________________________________________________#
clean.lca <- function(LCA_data){
  # Change columns to numeric where applicable
  LCA_data$Feed_soy_percent <- as.numeric(LCA_data$Feed_soy_percent)
  LCA_data$Electricity_kwh <- as.numeric(LCA_data$Electricity_kwh)
  LCA_data$Diesel_L <- as.numeric(LCA_data$Diesel_L)
  LCA_data$Petrol_L <- as.numeric(LCA_data$Petrol_L)
  LCA_data$NaturalGas_L <- as.numeric(LCA_data$NaturalGas_L)
  LCA_data$Yield_t_per_Ha <- as.numeric(LCA_data$Yield_t_per_Ha)
  LCA_data$Grow_out_period_days <- as.numeric(LCA_data$Grow_out_period_days)
  
  # Add country codes
  LCA_data$iso3c <- countrycode(LCA_data$Country, origin = "country.name", destination = "iso3c")
  
  # Scale feed percents to sum to 100%
  is.na(LCA_data$Feed_soy_percent) <- 0
  is.na(LCA_data$Feed_othercrops_percent) <- 0
  is.na(LCA_data$Feed_FMFO_percent) <- 0
  is.na(LCA_data$Feed_animal_percent) <- 0
  
  LCA_data <- LCA_data %>%
    mutate(sum_percent = Feed_soy_percent+Feed_othercrops_percent+Feed_FMFO_percent+Feed_animal_percent) %>%
    mutate(
    Feed_soy_percent = Feed_soy_percent/sum_percent,
    Feed_othercrops_percent = Feed_othercrops_percent/sum_percent,
    Feed_FMFO_percent = Feed_FMFO_percent/sum_percent,
    Feed_animal_percent = Feed_animal_percent/sum_percent
    ) %>%
    select(Year, Country, iso3c, Scientific.Name = Species.scientific.name, Common.Name = Species.common.name, 
           Production_system, Sample_size,
           Environment, Intensity, Yield_t_per_Ha, Grow_out_period_days, Mortality_rate, FCR, 
           Feed_type, Feed_soy_percent, Feed_othercrops_percent, Feed_FMFO_percent, Feed_animal_percent, Feed_method,
           Electricity_kwh, Diesel_L, Petrol_L, NaturalGas_L) %>%
    mutate(
      Production_system_group = case_when(
        (Production_system %in% c("Extensive raft culture", "Intensive lake net-pen", "Marine cages", "Marine net-pen",
                                  "Marine floating bag", "Lake-based net cages", "Integrated marine rafts", "Longline",
                                  "Bouchot culture", "Offshore cages", "Floating cages", "Net-pens", "Net pen",
                                  "Freshwater net pen (?)", "Saltwater net pen", "Semi-intensive cages")) ~ "open",
        (Production_system %in% c("Intensive pond", "Extensive pond polyculture", "Semi-intensive pond", "Extensive pond",
                                  "Earthen pond aquaculture", "Integrated pond, high input", "Integrated pond, medium inputs",
                                  "Solid-walled aquaculture system", "Reservoirs", "Earthen pond aquaculture integrated with pigs",
                                  "Earthern ponds", "Earthern/concrete ponds", "Ponds", "Silvo pond")) ~ "semi-open",
        (Production_system %in% c("Indoor recirculating", "Flow-through", "Land-based recirculating", "Onshore tanks",
                                  "Saltwater flow-through", "Freshwater flow-through", "Recirculating system", 
                                  "Land-based recirculating system", "Raceway", "Tanks / raceway", 
                                  "Semi-closed recirculating system")) ~ "closed",
        (Production_system %in% c("Integrated agri-aquaculture", "Unspecified", "Ecological farm", 
                                  "Floating net-cage in polyculture pond", "", "Ponds / recirculating", 
                                  "Ponds / pens")) ~ "not specified"
      ),  # Check groupings (some recirculating not really closed...)
      Intensity = case_when(
        (Intensity %in% c("Intensive")) ~ "Intensive",
        (Intensity %in% c("Semi-intensive", "Improved extensive", "Imp. extensive")) ~ "Semi-intensive",
        (Intensity %in% c("Extensive")) ~ "Extensive"
      ) # Many others can be identified based on the system description
      # Add fed/un-fed categories and replace FCR with 0 for all unfed 
      )
  
  # Note: Create column clean_sci_name - use this as the "official" scientific name column
  # Manually fill in blank scientific names based on Common.Name
  # Change osteichthyes (technically includes all terapods) to actinopterygii (bony fishes)
  # Simplify hybrid M. chrysops x M. saxatilis to its genus
  # Change outdated names (P. vannamei and P hypophthalmus)
  LCA_data <- LCA_data %>%
    mutate(clean_sci_name = case_when(Common.Name == "Freshwater prawn" ~ "Dendrobranchiata",
                                      Common.Name == "Indo-Pacific swamp crab" ~ "Brachyura",
                                      Common.Name == "Red crayfish" ~ "Astacidea", # crayfish are split into two superfamilies, so go to the next higher-classification, infraorder = Astacidea
                                      Common.Name == "Salmonids nei" ~ "Salmonidae",
                                      Common.Name == "Striped bass" ~ "Morone saxatilis",
                                      Common.Name == "Yellowtail_Seriola_Almaco jack" ~ "Seriola rivoliana",
                                      TRUE ~ Scientific.Name)) %>%
    mutate(clean_sci_name = case_when(str_detect(Scientific.Name, "spp") ~ str_replace(Scientific.Name, pattern = " spp\\.| spp", replacement = ""),
                                      Scientific.Name == "Morone chrysops x M. saxatilis" ~ "Morone",
                                      Scientific.Name == "Osteichthyes" ~ "Actinopterygii",
                                      Scientific.Name == "Penaeus vannamei" ~ "Litopenaeus vannamei",
                                      Scientific.Name == "Pangasius hypophthalmus" ~ "Pangasianodon hypophthalmus", 
                                      TRUE ~ clean_sci_name)) %>%
    mutate(FCR = case_when(str_detect(Scientific.Name, "Thunnus") ~ FCR/5,
                           TRUE ~ FCR)) 
  
  # sort(unique(LCA_data$clean_sci_name))
  # [1] "Acipenseridae"               "Actinopterygii"              "Anguilla"                    "Anoplopoma fimbria"          "Astacidea"                   "Brachyura"                   "Chanos chanos"              
  # [8] "Clarias batrachus"           "Clarias gariepinus"          "Cynoscion"                   "Cyprinus carpio"             "Dendrobranchiata"            "Dicentrarchus labrax"        "Epinephelus"                
  # [15] "Gadus morhua"                "Lates calcarifer"            "Litopenaeus vannamei"        "Macrobrachium"               "Macrobrachium amazonicum"    "Macrobrachium rosenbergii"   "Morone"                     
  # [22] "Mytilus edulis"              "Mytilus galloprovincialis"   "Oncorhynchus kisutch"        "Oncorhynchus mykiss"         "Oncorhynchus tshawytscha"    "Oreochromis niloticus"       "Pangasianodon hypophthalmus"
  # [29] "Pangasius"                   "Penaeus"                     "Penaeus monodon"             "Rachycentron canadum"        "Salmo salar"                 "Salmonidae"                  "Salvelinus alpinus"         
  # [36] "Sciaenops ocellatus"         "Scophthalmidae"              "Seriola rivoliana"           "Sparus aurata"               "Thunnus orientalis"          "Thunnus thynnus"   
  
  # Use Column clean_sci_name and common.name to create taxa groupings
  # CREATE "unassigned" category for: things that are not species-level Acipenseridae, Actinopterygii, Brachyura, Cynoscion spp, "Penaeus" can be fresh or marine
  # "Dicentrarchus labrax", "Lates calcarifer", "Morone" are migratory - i.e., including oceans, estuaries, and rivers - all categorized as other non-herbivore fin fish for now
  # Chanos chanos - mostly algae (but also inverts) - categorized as herbivore fin fish for now
  LCA_data <- LCA_data %>%
    mutate(taxa_group_name = case_when(#clean_sci_name %in% c("") ~ "algae", # none
                                  # clean_sci_name %in% c("") ~ "gastropods", # none
                                  clean_sci_name %in% c("Mytilus galloprovincialis", "Mytilus edulis") ~ "bivalves",
                                  clean_sci_name %in% c("Chanos chanos") ~ "herbivore finfish",
                                  clean_sci_name %in% c("Litopenaeus vannamei", "Penaeus monodon") ~ "marine shrimp",
                                  # clean_sci_name %in% c("") ~ "other marine crustacean",
                                  Common.Name %in% c("Red crayfish", "Freshwater prawn") ~ "freshwater crustacean",
                                  clean_sci_name %in% c("Macrobrachium", "Macrobrachium amazonicum", "Macrobrachium rosenbergii") ~ "freshwater crustacean", # note: freshwater crustacean assigned by sciname and commonnames
                                  clean_sci_name %in% c("Thunnus orientalis", "Thunnus thynnus") ~ "tuna",
                                  clean_sci_name %in% c("Oncorhynchus kisutch", "Oncorhynchus tshawytscha", "Salmo salar", "Salmonidae", "Salvelinus alpinus") ~ "salmon/char",
                                  clean_sci_name %in% c("Anoplopoma fimbria", "Dicentrarchus labrax", "Epinephelus", "Gadus morhua", "Lates calcarifer", "Morone", "Rachycentron canadum", "Sciaenops ocellatus", "Scophthalmidae", "Seriola rivoliana", "Sparus aurata") ~ "other non-herbivore marine finfish",
                                  clean_sci_name %in% c("Cyprinus carpio") ~ "carp",
                                  clean_sci_name %in% c("Oreochromis niloticus") ~ "tilapia",
                                  clean_sci_name %in% c("Clarias batrachus", "Clarias gariepinus", "Pangasianodon hypophthalmus", "Pangasius") ~ "catfish",
                                  clean_sci_name %in% c("Anguilla") ~ "eel",
                                  clean_sci_name %in% c("Oncorhynchus mykiss") ~ "trout",
                                  clean_sci_name %in% c("Pangasianodon hypophthalmus") ~ "other freshwater finfish",
                                  # clean_sci_name %in% c("") ~ "amphibians and reptiles", # none
                                  TRUE ~ "unassigned"
                                  )) %>%
    arrange(clean_sci_name) # LAST STEP arrange by sciname
  
}

#_____________________________________________________________________________________________________#
# Clean feed footprint data
#_____________________________________________________________________________________________________#
clean.feedFP <- function(feedFP_data){
  feedFP_data <- feedFP_data %>%
    filter(Unit != "kg PO4-eq") %>% # Not currently using
    group_by(Category, Unit) %>%
    # Revise this to change from arithmetic mean to a weighted mean
    summarise(FP_val = mean(Value, na.rm = TRUE), SD = sd(Value, na.rm = TRUE), .groups = 'drop') %>% 
    mutate(FP = case_when(
      (Unit == "kg CO2-eq") ~ "Carbon",
      (Unit == "m3") ~ "Water",
      (Unit == "m2a") ~ "Land",
      (Unit == "kg N-eq") ~ "Nitrogen",
      (Unit == "kg P-eq") ~ "Phosphorus"
    )) %>%
    select(FP, Category, FP_val, SD, Unit)
}

#_____________________________________________________________________________________________________#
# Estimate off farm, feed-associated footprint
#_____________________________________________________________________________________________________#
estimate.feedFP <- function(LCA_data, Feed_data, FP_option){
  # FP options are "Carbon", "Water", "Nitrogen", "Phosphorus", and "Land" 
  # Filter to FP option
  Feed_data <- Feed_data %>% 
    filter(FP == FP_option)
  
  soy <- Feed_data %>% filter(Category == "Soy") 
  soy <- soy$FP_val
  othercrops <- Feed_data %>% filter(Category == "Crop") %>% select(FP_val)
  othercrops <- othercrops$FP_val
  FMFO <- Feed_data %>% filter(Category == "Fishery") %>% select(FP_val)
  FMFO <- FMFO$FP_val
  animal <- Feed_data %>% filter(Category == "Animal by-products") %>% select(FP_val)
  animal <- animal$FP_val
  
  # Calculate estimated FP
  estimate_FP <- LCA_data$FCR * (soy*(LCA_data$Feed_soy_percent/100) + 
                                 othercrops*(LCA_data$Feed_othercrops_percent/100) + 
                                 FMFO * (LCA_data$Feed_FMFO_percent/100) + 
                                 animal * (LCA_data$Feed_animal_percent/100))
  return(estimate_FP)
}

#_____________________________________________________________________________________________________#
# Estimate on farm GHG footprint
#_____________________________________________________________________________________________________#
estimate.onfarm.GHG <- function(LCA_data, energy_data){
  diesel_perL_CO2eq <- 1 # Placeholder
  petrol_perL_CO2eq <- 1 # Placeholder
  naturalgas_perL_CO2eq <- 1 # Placeholder
  
  LCA_data <- LCA_data %>%
    left_join(energy_data, by = "iso3c")
  estimate_FP <- rowSums(cbind(c(LCA_data$Electricity_kwh * LCA_data$GWP_perkWh_CO2eq), 
                     c(LCA_data$Diesel_L * diesel_perL_CO2eq), 
                     c(LCA_data$Petrol_L * petrol_perL_CO2eq),
                     c(LCA_data$NaturalGas_L * naturalgas_perL_CO2eq)), 
                     na.rm = TRUE)
  
  return(estimate_FP)
}

#_____________________________________________________________________________________________________#
# Estimate on farm land footprint
#_____________________________________________________________________________________________________#
estimate.onfarm.land <- function(LCA_data){
  estimate_FP <- LCA_data$Yield_t_per_HA/LCA_data$harvest
  
  return(estimate_FP)
}

#_____________________________________________________________________________________________________#
# Estimate on farm water footprint
#_____________________________________________________________________________________________________#
estimate.onfarm.water <- function(LCA_data, evap_data, aerated = FALSE){
  # Requires calculating land footprint first, saved in a column "onfarm.land"
  
  if(aerated == TRUE){
    aeration.factor <- 1.25 # Placeholder of 25% increase in evaporation rate
  } else{
    aeration.factor <- 1
  }
  # Join lca and evap data by country
  LCA_data <- LCA_data %>%
    left_join(evap_data, by = "iso3c")
  
  # Need land FP in HA, evap rate in volume/HA/day and grow out period in days
  estimate_FP <- LCA_data$onfarm.land * aeration.factor*LCA_data$evap_rate *  LCA_data$Grow_out_period_days
  
  return(estimate_FP)
}

#_____________________________________________________________________________________________________#
# Estimate on farm N and P footprint
#_____________________________________________________________________________________________________#
estimate.feed.NP <- function(LCA_data, Feed_data, FP){
  # FP set to N or P
  if(FP == "N"){
    Feed_data <- Feed_data %>%
      filter(element == "N")
  }else if(FP == "P"){
    Feed_data <- Feed_data %>%
      filter(element == "P")
  }else print("error: FP must equal N or P")
  
  soy <- Feed_data %>% filter(ingredient == "Soy") 
  soy <- soy$value
  othercrops <- Feed_data %>% filter(ingredient == "Crop")
  othercrops <- othercrops$value
  FMFO <- Feed_data %>% filter(ingredient == "Fishery")
  FMFO <- FMFO$value
  animal <- Feed_data %>% filter(ingredient == "Animal by-products") 
  animal <- animal$value
  
  # Calculate estimated FP
  estimate_FP <- soy*(LCA_data$Feed_soy_percent/100) + 
                                   othercrops*(LCA_data$Feed_othercrops_percent/100) + 
                                   FMFO * (LCA_data$Feed_FMFO_percent/100) + 
                                   animal * (LCA_data$Feed_animal_percent/100)
  
  return(estimate_FP)
}

#_____________________________________________________________________________________________________#
# Estimate on-farm N and P emissions
#_____________________________________________________________________________________________________#
clean_feedNutrition <- function(feedNutrition_data){
  # Clean columns
  feedNutrition_data$Phosphorus <- as.numeric(feedNutrition_data$Phosphorus....)
  feedNutrition_data <- feedNutrition_data %>% 
    mutate(Nitrogen = Crude.protein..../ 6.25) 
  feedNutrition_data$Nitrogen <- as.numeric(feedNutrition_data$Nitrogen)   #N in % of DM
  
  #soy: Entry numbers: 601-620
  s=c(602,604,608,610,612,614,616,618,620)
  
  #animal by products, Entry number: meat byproducts 385-392, poultry 479-481 (take only 100% DM)
  a=c(386,388,390,392,480,482,484)
  
  #other crops: cassava 141-142, peanut extr. 461-464, linseed extr. 345-348 , 
  #corn gluten meal 223-228, pea (not protein concentrate) 455-456, rape 493-496, sunflower oil and meal 633-636, 
  #wheat byproduct 683-688, Sorghum grain 573-580, triticale 667-668, Navy beans (instead of faba) 57-58,
  cr=c(142,462,464,346,348,224,226,228,456,494,496,634,636,684,686,688,578,580,668,58)
  
  #fishmeal and oil
  # fish 301-304, fish alewife meal 305-306,anchovy meal 307-308, catfish meal 313-314, herring meal 317-318, 
  # mackerel meal 322, manhaden meal 323-324, redfish meal 325-326, salmon meal 329-330, sardine meal 333-334
  #tuna meal 337-338, white fish 341-342
  f=c(302,304,306,308,314,318,322,324,326,330,334,338,342)
  
  # Create dataframe with crop N and P values
  out <- data.frame(ingredient = c("Animal by-products", "Animal by-products", "Crop", "Crop", "Fishery", "Fishery", "Soy", "Soy"), 
                    element = rep(c("N", "P"), 4), value = NA, sd = NA)
  
  # Calculate the mean and standard deviation by group
  # Change this to calculated weighted mean
  out$value[out$ingredient=="Animal by-products" & out$element == "N"] <- 
    mean(feedNutrition_data$Nitrogen[which(feedNutrition_data$Entry.Number %in% a)],na.rm=TRUE)
  out$sd[out$ingredient=="Animal by-products" & out$element == "N"] <- 
    sd(feedNutrition_data$Nitrogen[which(feedNutrition_data$Entry.Number %in% a)],na.rm=TRUE)
  
  out$value[out$ingredient=="Animal by-products" & out$element == "P"] <- 
    mean(feedNutrition_data$Phosphorus[which(feedNutrition_data$Entry.Number %in% a)],na.rm=TRUE)
  out$sd[out$ingredient=="Animal by-products" & out$element == "P"] <- 
    sd(feedNutrition_data$Phosphorus[which(feedNutrition_data$Entry.Number %in% a)],na.rm=TRUE)
  
  out$value[out$ingredient=="Crop" & out$element == "N"] <- 
    mean(feedNutrition_data$Nitrogen[which(feedNutrition_data$Entry.Number %in% cr)],na.rm=TRUE)
  out$sd[out$ingredient=="Crop" & out$element == "N"] <- 
    sd(feedNutrition_data$Nitrogen[which(feedNutrition_data$Entry.Number %in% cr)],na.rm=TRUE)
  
  out$value[out$ingredient=="Crop" & out$element == "P"] <- 
    mean(feedNutrition_data$Phosphorus[which(feedNutrition_data$Entry.Number %in% cr)],na.rm=TRUE)
  out$sd[out$ingredient=="Crop" & out$element == "P"] <- 
    sd(feedNutrition_data$Phosphorus[which(feedNutrition_data$Entry.Number %in% cr)],na.rm=TRUE)
  
  out$value[out$ingredient=="Fishery" & out$element == "N"] <- 
    mean(feedNutrition_data$Nitrogen[which(feedNutrition_data$Entry.Number %in% f)],na.rm=TRUE)
  out$sd[out$ingredient=="Fishery" & out$element == "N"] <- 
    sd(feedNutrition_data$Nitrogen[which(feedNutrition_data$Entry.Number %in% f)],na.rm=TRUE)
  
  out$value[out$ingredient=="Fishery" & out$element == "P"] <- 
    mean(feedNutrition_data$Phosphorus[which(feedNutrition_data$Entry.Number %in% f)],na.rm=TRUE)
  out$sd[out$ingredient=="Fishery" & out$element == "P"] <- 
    sd(feedNutrition_data$Phosphorus[which(feedNutrition_data$Entry.Number %in% f)],na.rm=TRUE)
  
  out$value[out$ingredient=="Soy" & out$element == "N"] <- 
    mean(feedNutrition_data$Nitrogen[which(feedNutrition_data$Entry.Number %in% s)],na.rm=TRUE)
  out$sd[out$ingredient=="Soy" & out$element == "N"] <- 
    sd(feedNutrition_data$Nitrogen[which(feedNutrition_data$Entry.Number %in% s)],na.rm=TRUE)
  
  out$value[out$ingredient=="Soy" & out$element == "P"] <- 
    mean(feedNutrition_data$Phosphorus[which(feedNutrition_data$Entry.Number %in% s)],na.rm=TRUE)
  out$sd[out$ingredient=="Soy" & out$element == "P"] <- 
    sd(feedNutrition_data$Phosphorus[which(feedNutrition_data$Entry.Number %in% s)],na.rm=TRUE)
  
  return(out)
}

# Method 1: estimate from C content, Fig 1 of the above reference, linear regression values. 
# Assuming metabolizable energy density is equal throughout whole fish
fishN <- function(Water,
           Energy.total.metabolizable.calculated.from.the.energy.producing.food.components.original.as.from.source.kcal) {
  # Convert to dry matter fraction
    DM = (100 - Water) / 100     
    
    #In aquacultured fish, whole-body C content was estimated from energy content using a converting factor of 1 g C = 11.4 kcal
    fish_C_percent <-
      1 / DM * Energy.total.metabolizable.calculated.from.the.energy.producing.food.components.original.as.from.source.kcal / 11.4
    N_slope <- -0.175   # from Czamanski et al 2011 Fig 1
    N_intercept <- 18.506  #intercept (personal communication with authors)
    
    fish_N <- fish_C_percent * N_slope + N_intercept #in % of DM
    
    return(fish_N)
  }

fishP <- function(Water,
           Energy.total.metabolizable.calculated.from.the.energy.producing.food.components.original.as.from.source.kcal) {
  # Convert to dry matter fraction
    DM = (100 - Water) / 100     
    
    #In aquacultured fish, whole-body C content was estimated from energy content using a converting factor of 1 g C = 11.4 kcal
    fish_C_percent <-
      1 / DM * Energy.total.metabolizable.calculated.from.the.energy.producing.food.components.original.as.from.source.kcal / 11.4
    
    P_slope <- -0.083  # from Czamanski et al 2011 Fig 1
    P_intercept <- 6.311  #intercept (personal communication with authors)
    
    fish_P <- fish_C_percent * P_slope + P_intercept #in % of DM
    
    return(fish_P)
  }

# Method 2: estimate from total lipids, figure 2, assuming lipid density is equal throughout whole body
fishN_viaFat <- function(Water, Fat.total) {
  # Dry Matter percentage
  DM = (100 - Water) / 100     
  
  # Fat.total in units of g to 100 g
  # Convert to DM
  fat = Fat.total * 1  / DM  
    
  
  # Fat percentage in whole body in % of DM.
  fish_C_via_fat = fat * 0.31 + 38  # linear regression from Czamanski et al 2011 Fig 2
  fish_N_via_fat = fat * -0.08 + 12 # linear regression from Czamanski et al 2011 Fig 2
  
  return(fish_N_via_fat)
}

fishP_viaFat <- function(Water, Fat.total) {
  # Dry Matter percentage
  DM = (100 - Water) / 100     
  
  # Fat.total in units of g to 100 g
  # Convert to DM
  fat = Fat.total * 1  / DM 
    
  # Fat percentage in whole body in % of DM. Assuming Fat is contained only (mostly) in the edible portion.
  fish_C_via_fat = fat * 0.31 + 38  #linear regression from fig 2
  fish_P_via_fat = fat * -0.04 + 3.2 #linear regression from fig 2
  
  return(fish_P_via_fat)
}

# N and P discharge model
estimate.onfarm.NP <- function(LCA_data){
  # From Alon's branch N_P_discharge_simple_model
  # Simplified model of discharge, similar to many NPZ models 
  
  LCA_data <- LCA_data %>%
    #conservation of nutrient: feed-fish_content=discharge to enviroment. units: weights as FCR
    mutate(onfarm.N = (feed.N.percent/100*(FCR)-feed.N.percent/100*1),
           onfarm.P = (feed.P.percent/100*(FCR)-feed.P.percent/100*1))
  
    # To avoid negative emissions
  LCA_data$onfarm.N[LCA_data$onfarm.N < 0] <- 0
  LCA_data$onfarm.P[LCA_data$onfarm.P < 0] <- 0

  #in units of FCR e.g. kg N and P for 1 kg of live weight fish
  
  return(LCA_data) 
}

