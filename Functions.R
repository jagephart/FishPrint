# Functions

#_____________________________________________________________________________________________________#
# Estimate feed-associated footprint
#_____________________________________________________________________________________________________#
estimate.feedFP <- function(LCA_data, Feed_data, FP){
  # FP options are "GHG", "water", "N", "P", "land" 
  # Filter to FP option
  Feed_data <- Feed_data %>% 
    filter(FP_name == FP)
  
  soy <- Feed_data %>% filter(feed_component == "soy") 
  soy <- soy$FP_val
  othercrops <- Feed_data %>% filter(feed_component == "othercrops") %>% select(FP_val)
  othercrops <- othercrops$FP_val
  FMFO <- Feed_data %>% filter(feed_component == "FMFO") %>% select(FP_val)
  FMFO <- FMFO$FP_val
  animal <- Feed_data %>% filter(feed_component == "animal") %>% select(FP_val)
  animal <- animal$FP_val
  
  # Calculate estimated FP
  estimate_FP <- LCA_data$FCR * (soy*(LCA_data$Feed_soy_percent/100) + 
                                 othercrops*(LCA_data$Feed_othercrops_percent/100) + 
                                 FMFO * (LCA_data$Feed_FMFO_percent/100) + 
                                 animal * (LCA_data$Feed_animal_percent/100))
  return(estimate_FP)
}

#_____________________________________________________________________________________________________#
# Estimate non-feed GHG footprint
#_____________________________________________________________________________________________________#
estimate.nonfeed.GHG <- function(LCA_data, energy_data){
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
# Estimate non-feed N and P footprint
#_____________________________________________________________________________________________________#
estimate.nonfeed.NP <- function(LCA_data, Feed_data, FP){
  # FP set to N or P
  if(FP == "N"){
    protein_factor <- 2.5 # Placeholder
  }else if(FP == "P"){
    protein_factor <- 3 # Placeholder
  }else print("error: FP must equal N or P")
  
  Feed_data <- Feed_data %>% 
    filter(FP_name == "protein_content")
  
  soy <- Feed_data %>% filter(feed_component == "soy") 
  soy <- soy$FP_val
  othercrops <- Feed_data %>% filter(feed_component == "othercrops") %>% select(FP_val)
  othercrops <- othercrops$FP_val
  FMFO <- Feed_data %>% filter(feed_component == "FMFO") %>% select(FP_val)
  FMFO <- FMFO$FP_val
  animal <- Feed_data %>% filter(feed_component == "animal") %>% select(FP_val)
  animal <- animal$FP_val
  
  # Calculate estimated FP
  estimate_FP <- LCA_data$FCR * protein_factor*(soy*(LCA_data$Feed_soy_percent/100) + 
                                   othercrops*(LCA_data$Feed_othercrops_percent/100) + 
                                   FMFO * (LCA_data$Feed_FMFO_percent/100) + 
                                   animal * (LCA_data$Feed_animal_percent/100)) -
                LCA_data$fish_protein_content
  return(estimate_FP)
}


#_____________________________________________________________________________________________________#
# Estimate non-feed land footprint
#_____________________________________________________________________________________________________#
estimate.nonfeed.land <- function(LCA_data){
  estimate_FP <- LCA_data$Yield_per_HA/LCA_data$harvest
  
  return(estimate_FP)
}

#_____________________________________________________________________________________________________#
# Estimate non-feed water footprint
#_____________________________________________________________________________________________________#
estimate.nonfeed.water <- function(LCA_data, evap_data){
  
}