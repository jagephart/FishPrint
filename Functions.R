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
  LCA_data$Yield_t_per_HA <- as.numeric(LCA_data$Yield_t_per_HA)
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
    select(Year, Country, iso3c, Scientific.Name = Species.scientific.name, Production_system, Sample_size,
           Environment, Intensity, Yield_t_per_HA, Grow_out_period_days, Mortality_rate, FCR, 
           Feed_type, Feed_soy_percent, Feed_othercrops_percent, Feed_FMFO_percent, Feed_animal_percent, Feed_method,
           Electricity_kwh, Diesel_L, Petrol_L, NaturalGas_L)
  
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
  estimate_FP <- LCA_data$Yield_per_HA/LCA_data$harvest
  
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
      filter(element == "N")
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
  
  #define categorial feeds: soy, other products, animal and FM&O 
  # FIX IT - Alon, these don't look correct, please check (e.g., I see hay in the soy range)
  #soy: Entry numbers: 601-620
  s=seq(601,620)
  
  #animal by products, Entry number: meat byproducts 385-392, poultry 479-481
  a=c(seq(385,392),seq(479,481))
  
  #other crops: cassava 141-142, peanut extr. 461-464, linseed extr. 345-348 , 
  #corn gluten meal 223-228, pea (not protein conentrate) 455-456, rape 493-496, sunflower oil and meal 633-636, 
  #wheat byproduct 683-688, Sorghum grain 573-580, triticale 667-668, Navy beans (instead of faba) 57-58,
  cr=c(141,142,seq(461,464),seq(345,348),seq(223,228),seq(455,456),seq(493,496),seq(633,636),seq(683,688),seq(573,580),seq(667,668),seq(57,58))
  
  #fishmeal and oil
  # fish alewife meal 305-306,anchovy meal 307-308, catfish meal 313-314, herring meal 317-318, 
  # mackerel meal 322, manhaden meal 323-324, redfish meal 325-326, salmon meal 329-330, sardine meal 333-334
  #tuna meal 337-338, white fish 341-342
  f=c(305,306,307,308,313,314,317,318,322,323,324,325,326,329,330,333,334,337,338,341,342)
  
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
# FIX IT: Alon - which reference and please write up a narrative description of the methods
fishN <- function(Water,
           Energy.total.metabolizable.calculated.from.the.energy.producing.food.components.original.as.from.source.kcal) {
  # Convert to dry matter fraction
    DM = (100 - Water) / 100     
    
    #In aquacultured fish, whole-body C content was estimated from energy content using a converting factor of 1 g C = 11.4 kcal
    fish_C_percent <-
      1 / DM * Energy.total.metabolizable.calculated.from.the.energy.producing.food.components.original.as.from.source.kcal / 11.4
    N_slope <- -0.175   #from fig 1 # FIX IT: Alon add ref here
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
    
    P_slope <- -0.083  #from fig 1
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
  fat = Fat.total * 1  /  1 / DM # FIX IT: Alon - check for any need for parentheses here
    
  
  # Fat percentage in whole body in % of DM.
  fish_C_via_fat = fat * 0.31 + 38  #linear regression from fig 2
  fish_N_via_fat = fat * -0.08 + 12 #linear regression from fig 2
  
  return(fish_N_via_fat)
}

fishP_viaFat <- function(Water, Fat.total) {
  # Dry Matter percentage
  DM = (100 - Water) / 100     
  
  # Fat.total in units of g to 100 g
  # Convert to DM
  fat = Fat.total * 1  / 1 / DM
    
  # Fat percentage in whole body in % of DM. Assuming Fat is contained only (mostly) in the edible portion.
  fish_C_via_fat = fat * 0.31 + 38  #linear regression from fig 2
  fish_P_via_fat = fat * -0.04 + 3.2 #linear regression from fig 2
  
  return(fish_P_via_fat)
}

# N and P discharge models
N_P_discharge_model <- function(fish_N_percent,fish_P_percent,feed_N_percent,feed_P_percent,FCR){
  #G=I-[E+F]   from the above reference "Methods", where G is growth, I ingestion, E emission to enviroment (resupply) and F is feces
  #FCR = I/G  FCR assumes I=1 unit weight -->E+F=FCR-1
  #in terms of nutrient conservation: E_N+F_N=feed_N*(FCR-1)-fish_N*1; E_P+F_P=feed_P(FCR-1)-fish_P*1
  #convert N and P to molar weight 
  Lm=0.75;
  #browser()
  fish_N=fish_N_percent/14;feed_N=feed_N_percent/14; #stoichometric molar ratio
  fish_P=fish_P_percent/31;feed_P=feed_P_percent/31;
  beta_N<-0.91
  R_fish<-fish_N/fish_P;
  R_feed<-feed_N/feed_P;
  beta_P<- -0.19*R_fish/R_feed+0.8   #from fig 4 of reference
  if (R_feed*beta_N/beta_P<R_fish){
    alpha_N<-Lm;
    alpha_P<-Lm*R_feed*beta_N/beta_P/R_fish}
  else if (R_feed*beta_N/beta_P>R_fish){
    alpha_P<-Lm;
    alpha_N<-Lm*(R_fish)/(R_feed*beta_N/beta_P)}
  
  R_resupply<-((1-alpha_N)*beta_N)/((1-alpha_P)*beta_P)*R_feed #equation 7 from reference # FIX IT: Alon - add reference
  R_feces<-R_feed*(1-beta_N)/(1-beta_P)
  constant_N<-(feed_N_percent/100*(FCR)-fish_N_percent/100*1)/0.014; #conservation of nutrient: feed-fish_content=feces+resupply: units of molar weight: 0.014 kg/mole
  constant_P<-(feed_P_percent/100*(FCR)-fish_P_percent/100*1)/0.031; #conservation of nutrient: feed-fish_content=feces+resupply: units of molar weight: 0.031 kg/mole
  E_P<-(constant_N-R_feces*constant_P)/(R_resupply-R_feces);
  E_N<-R_resupply*E_P
  F_N<-constant_N-E_N;
  F_P<-constant_P-E_P;
  output<-c(E_N*0.014,F_N*0.014,E_P*0.031,F_P*0.031)
  return(output)
}


N_P_discharge_simple_model <- function(fish_N_percent,fish_P_percent,feed_N_percent,feed_P_percent,FCR){
  #simplified model of discharge, similar to many NPZ models 
  discharge_N<-(feed_N_percent/100*(FCR)-fish_N_percent/100*1); #conservation of nutrient: feed-fish_content=discharge to enviroment. units: weights as FCR
  discharge_P<-(feed_P_percent/100*(FCR)-fish_P_percent/100*1); #conservation of nutrient: feed-fish_content=discharge to enviromen.  weights as FCR 
  if(discharge_N<=0){discharge_N=0}else{discharge_N=discharge_N} #to avoid negative emissions
  if(discharge_P<=0){discharge_P=0}else{discharge_P=discharge_P}
  RT=c(discharge_N,discharge_P)  #in units of FCR e.g. kg N and P for 1 kg of live weight fish
  
  return(RT) 
}

