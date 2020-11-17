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
  LCA_data$Yield_kg_per_m3 <- as.numeric(LCA_data$Yield_kg_per_m3)
  LCA_data$Grow_out_period_days <- as.numeric(LCA_data$Grow_out_period_days)
  
  # Add country codes
  LCA_data$iso3c <- countrycode(LCA_data$Country, origin = "country.name", destination = "iso3c")
  
  # Scale feed percents to sum to 100%
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
           Environment, Intensity, Yield_t_per_Ha, Yield_kg_per_m3, Grow_out_period_days, Mortality_rate, FCR, 
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
                                  "Earthern ponds", "Earthern/concrete ponds", "Ponds", "Silvo pond")) ~ "semi",
        (Production_system %in% c("Indoor recirculating", "Flow-through", "Land-based recirculating", "Onshore tanks",
                                  "Saltwater flow-through", "Freshwater flow-through", "Recirculating system", 
                                  "Land-based recirculating system", "Raceway", "Tanks / raceway", 
                                  "Semi-closed recirculating system")) ~ "closed",
        (Production_system %in% c("Integrated agri-aquaculture", "Unspecified", "Ecological farm", 
                                  "Floating net-cage in polyculture pond", "", "Ponds / recirculating", 
                                  "Ponds / pens")) ~ "not specified"
      ),  # Check groupings (some recirculating not really closed...)
      Intensity = case_when(
        (Intensity %in% c("Intensive")) ~ "intensive",
        (Intensity %in% c("Semi-intensive", "Improved extensive", "Imp. extensive")) ~ "semi",
        (Intensity %in% c("Extensive")) ~ "extensive"
      ),
      # Many others can be identified based on the system description
      # Add fed/un-fed categories and replace FCR with 0 for all unfed 
      ) %>%
    # Convert "not specified" to NA
    mutate(Production_system_group = na_if(Production_system_group, "not specified"))
  
  # Make adjustments to data:
  # Divide tuna FCR by 5
  LCA_data <- LCA_data %>%
    mutate(FCR = case_when(str_detect(Scientific.Name, "Thunnus") ~ FCR/5,
                           TRUE ~ FCR)) 
    
  
  # Create column clean_sci_name - use this as the "official" scientific name column
  # When no valid sci name exists just use common name
  # Simplify hybrid M. chrysops x M. saxatilis to its genus
  # Change outdated names (P. vannamei and P hypophthalmus)
  LCA_data <- LCA_data %>%
    mutate(Scientific.Name = case_when(str_detect(Scientific.Name, "spp") ~ str_replace(Scientific.Name, pattern = " spp\\.| spp", replacement = ""),
                                       TRUE ~ Scientific.Name)) %>%
    mutate(clean_sci_name = case_when(Common.Name == "Freshwater prawn" ~ "Freshwater prawn",
                                      Common.Name == "Indo-Pacific swamp crab" ~ "Brachyura",
                                      Common.Name == "Red crayfish" ~ "Red crayfish", # crayfish are split into two superfamilies, so go to the next higher-classification, infraorder = Astacidea
                                      Common.Name == "River eels nei" ~ "Freshwater eels",
                                      Common.Name == "Salmonids nei" ~ "Salmonidae",
                                      Common.Name == "Striped bass" ~ "Morone saxatilis",
                                      Common.Name == "Yellowtail_Seriola_Almaco jack" ~ "Seriola rivoliana",
                                      TRUE ~ Scientific.Name)) %>%
    mutate(clean_sci_name = case_when(Scientific.Name == "Morone chrysops x M. saxatilis" ~ "Morone",
                                      Scientific.Name == "Labeo rohita and Catla Catla" ~ "Mixed carps",
                                      Scientific.Name == "Osteichthyes" ~ "Freshwater fishes",
                                      Scientific.Name == "Penaeus vannamei" ~ "Litopenaeus vannamei",
                                      Scientific.Name == "Pangasius hypophthalmus" ~ "Pangasianodon hypophthalmus",
                                      TRUE ~ clean_sci_name))
  
  # Remove Sturgeons since these are farmed for caviar
  # Remove eels
  LCA_data <- LCA_data %>%
    filter(Scientific.Name != "Acipenseridae") %>%
    filter(Scientific.Name != "Anguilla")
  

}

#_____________________________________________________________________________________________________#
# Rebuild FAO fish from zip file
#_____________________________________________________________________________________________________#

rebuild_fish <- function(path_to_zipfile) {
  require(tools) # needed for file_path_sans_ext
  require(dplyr)
  require(purrr)
  require(readxl) # part of tidyverse but still need to load readxl explicitly, because it is not a core tidyverse package
  
  # The following ensures unzipped folder is created in the same directory as the zip file (can be different from the working directory)
  # set outdir
  if (file.exists(basename(path_to_zipfile))) { # if file is in current directory and only file name was given
    outdir <- getwd()
  } else if (file.exists(path_to_zipfile)) { # if file path was given
    outdir <- dirname(path_to_zipfile)
  } else {
    stop("Check path_to_zipfile")
  }
  
  foldername <- file_path_sans_ext(basename(path_to_zipfile))
  outfolder <- paste(outdir, foldername, sep = "/")
  unzip(path_to_zipfile, exdir = outfolder) # Problem: if unable to unzip folder, still creates outfolder how to supress this?
  # setwd(outfolder)
  # list files
  fish_files <- list.files(outfolder)
  
  # read .xlsx file (explains data structure of time series)
  # IMPORTANT: column ORDER (ABCDEF) in DS file should match columns ABCDEF in time series for looping to work below
  # each row gives info for how this time series column should be merged with a code list (CL) file
  ds_file <- fish_files[grep("DSD", fish_files)]
  path_to_ds <- paste(outfolder, ds_file, sep = "/")
  
  # skip removes title row
  ds <- read_excel(path_to_ds, skip=1)
  
  # manually correct ds file's codelist ID column:
  ds <- ds %>%
    mutate(Codelist_Code_id = case_when(
      Concept_id == "SOURCE" ~ "IDENTIFIER",
      Concept_id == "SYMBOL" ~ "SYMBOL",
      Concept_id != "SYMBOL|SOURCE" ~ Codelist_Code_id
    ))
  
  # Multiple CL files have the following column names in common: "Identifier" and "Code"
  # Which means after merge, below, you get "Identifier.x" and "Identifier.y", etc.
  # To disambiguate, Append Codelist with Concept_id
  code_ids_to_change<-ds$Codelist_Code_id[grep("IDENTIFIER|CODE", ds$Codelist_Code_id)]
  concept_ids_to_append<-ds$Concept_id[grep("IDENTIFIER|CODE", ds$Codelist_Code_id)]
  new_code_ids <- paste(concept_ids_to_append, code_ids_to_change, sep = "_")
  ds$Codelist_Code_id[grep("IDENTIFIER|CODE", ds$Codelist_Code_id)]<-new_code_ids
  
  # remove non CSVs (do this to ignore "CL_History.txt" file)
  fish_files <- fish_files[grep(".csv", fish_files)]
  
  # read in time series.csv
  time_files <- fish_files[grep("TS", fish_files)]
  path_to_ts <- paste(outfolder, time_files, sep = "/")
  time_series <- read.csv(path_to_ts)
  names(time_series) <- tolower(names(time_series))
  time_series_join <- time_series
  
  for (i in 1:nrow(ds)) {
    # TRUE/FALSE: is there a filename listed in Codelist_id?
    if (!is.na(ds$Codelist_id[i])) {
      # Use ds file to generate path_to_cl individually
      code_file_i <- paste(ds$Codelist_id[i], ".csv", sep = "")
      path_to_cl <- paste(outfolder, code_file_i, sep = "/")
      cl_i <- read.csv(path_to_cl, check.names = FALSE) # check.names = FALSE to prevent R from adding "X" in front of column "3Alpha_Code" - creates problems because this is the matching column for merging with time series
      
      # Many CL files have "Name" as a column, also Name_En, Name_Fr, Name_Es, etc
      # Also, "Identifier", "Major Group", and "Code" are common across some CL files
      # To disambiguate, append "Concept_ID" from DS file to all columns in CL that contain these terms
      concept_names <- paste(ds$Concept_id[i], names(cl_i)[grep("Name|Major_Group|Identifier|Code", names(cl_i))], sep = "_")
      names(cl_i)[grep("Name|Major_Group|Identifier|Code", names(cl_i))] <- concept_names
      
      
      names(cl_i) <- tolower(names(cl_i)) # convert all cl headers to lowercase
      merge_col <- tolower(ds$Codelist_Code_id[i]) # do the same to DS file's code ID so it matches with cl
      
      
      # If factor...
      #if (is.factor(cl_i[[merge_col]])) {
      # ...Test if factor levels need to be merged?
      #if (!nlevels(cl_i[[merge_col]]) == nlevels(time_series_join[[names(time_series_join)[i]]])) {
      # combined <- sort(union(time_series_join[[names(time_series_join)[i]]], levels(cl_i[[merge_col]])))
      #    levels(time_series_join[[names(time_series_join)[i]]]) <- levels(cl_i[[merge_col]])
      #  }
      #}
      # This avoids warnings about unequal factor levels below
      
      # Try converting to character first instead
      if (is.factor(cl_i[[merge_col]])){
        cl_i[[merge_col]]<-as.character(cl_i[[merge_col]])
        time_series_join[[names(time_series_join)[i]]]<-as.character(time_series_join[[names(time_series_join)[i]]])
      }
      
      
      # Can't just merge by column number:
      # In Time Series, column COUNTRY, AREA, SOURCE, SPECIES, and UNIT correspond to column 1 in their respective CL files
      # but in Time Series, column SYMBOL corresponds to column 2
      
      # Note: the following code does not work: #time_series_join<-left_join(time_series, cl_i, by = c(names(time_series)[i] = merge_col))
      # the argument "by" needs to take on the form of join_cols as shown below
      firstname <- names(time_series_join)[i]
      join_cols <- merge_col
      names(join_cols) <- firstname
      
      
      time_series_join <- left_join(time_series_join, cl_i, by = join_cols)
      
      # Convert back to factor
      if (is.character(time_series_join[[names(time_series_join)[i]]])){
        time_series_join[[names(time_series_join)[i]]]<-as.factor(time_series_join[[names(time_series_join)[i]]])
      }
    }
    # Expected warning: Coerces from factor to character because time_series$SPECIES (nlevels=2341) and CL_FI_SPECIES_GROUPS.csv column "3alpha_code" (nlevels = 12751) have different number of factor levels
    # Expected warning: Coerces from factor to chracter because time_series$UNIT and CL_FILE_UNIT.csv column "code" have different number of factor levels
    # Expected warning: Coerces from factor to character because time_series$SYMBOL and CL_FI_SYMBOL.csv column "symbol" have diff number of factors
  }
  
  return(time_series_join)
}

#_____________________________________________________________________________________________________#
# Add taxa grouping
#_____________________________________________________________________________________________________#

add_taxa_group <- function(lca_dat_clean, fishstat_dat){
  isscaap_lookup <- fishstat_dat %>%
    select(species_scientific_name, isscaap_group) %>%
    unique()
  
  lca_dat_out <- lca_dat_clean %>%
    left_join(isscaap_lookup, by = c("clean_sci_name" = "species_scientific_name")) 
  
  #setdiff(lca_dat_clean$clean_sci_name, unique(fishstat_dat$species_scientific_name))
  #sort(unique(fishstat_dat$isscaap_group))
  #sort(unique(lca_dat_out$isscaap_group))

  lca_dat_out <- lca_dat_out %>% 
    # First pass is to deal with NAs in ISSCAAP group
    mutate(isscaap_group = case_when(clean_sci_name %in% c("Cynoscion", "Epinephelus") ~ "Miscellaneous marine fishes",
                                     clean_sci_name %in% c("Freshwater fishes", "Pangasius") ~ "Miscellaneous freshwater fishes",
                                     clean_sci_name %in% c("Red crayfish") ~ "Freshwater crustaceans",
                                     clean_sci_name %in% c("Freshwater prawn") ~ "Shrimps, prawns",
                                     clean_sci_name %in% c("Litopenaeus vannamei", "Macrobrachium", "Macrobrachium amazonicum", "Macrobrachium rosenbergii", "Penaeus") ~ "Shrimps, prawns",
                                     clean_sci_name %in% c("Mixed carps") ~ "Carps, barbels and other cyprinids",
                                     TRUE ~ isscaap_group)) 
  
  # Create taxa groups
  lca_dat_out <- lca_dat_out %>%
    # Split up carps
    mutate(taxa_group_name = case_when(clean_sci_name == "Cyprinus carpio" ~ "Common carp",
                                       clean_sci_name == "Mixed carps" ~ "Other Carps, barbels and cyprinids",
                                       # Split salmons and trouts
                                       str_detect(Common.Name, "salmon|Salmonids") ~ "Salmon",
                                       str_detect(Common.Name, "trout") ~ "Trout",
                                       str_detect(Common.Name, "char") ~ "Miscellaneous diadromous fishes",
                                       # Remove milkfish from misc diadromous fishes
                                       clean_sci_name == "Chanos chanos" ~ "Miscellaneous diadromous fishes",
                                       # Combine groups into misc marine fishes
                                       isscaap_group %in% c("Miscellaneous coastal fishes", "Miscellaneous demersal fishes", "Miscellaneous pelagic fishes", "Flounders, halibuts, soles") ~ "Miscellaneous marine fishes",
                                       TRUE ~ isscaap_group))
                                                                          
  # TEST:
  # lca_dat_out%>%
  #   #filter(isscaap_group == "Salmons, trouts, smelts") %>%
  #   #filter(str_detect(Common.Name, "char")) %>%
  #   #filter(isscaap_group == "Miscellaneous diadromous fishes") %>%
  #   select(clean_sci_name, isscaap_group, taxa_group_name) %>%
  #   unique() %>%
  #   arrange(taxa_group_name)
  
}


#_____________________________________________________________________________________________________#
# Replicate data based on Sample_size column
#_____________________________________________________________________________________________________#
rep_data <- function(lca_dat_clean){
  lca_dat_clean_rep <- lca_dat_clean %>%
    # Clean up sample size column
    # First ignore numbers that are percentages, then find and extract any numbers, then fill the rest in with 1s
    mutate(clean_sample_size = case_when(str_detect(Sample_size, "%") ~ 1,
                                         str_detect(Sample_size, "[0-9]+") ~ as.numeric(str_extract(Sample_size, pattern = "[0-9]+")),
                                         TRUE ~ 1)) 
  
  lca_dat_clean_rep <- as.data.frame(lapply(lca_dat_clean_rep, rep, lca_dat_clean_rep$clean_sample_size))
}

#_____________________________________________________________________________________________________#
# Clean feed footprint data
#_____________________________________________________________________________________________________#
clean.feedFP <- function(feedFP_data){
  feedFP_data <- feedFP_data %>%
    filter(Unit != "kg PO4-eq") %>% # Not currently using
    group_by(Category, Unit, Allocation.method) %>%
    # Revise this to change from arithmetic mean to a weighted mean
    summarise(FP_val = mean(Value, na.rm = TRUE), SD = sd(Value, na.rm = TRUE), .groups = 'drop') %>% 
    mutate(FP = case_when(
      (Unit == "kg CO2-eq") ~ "Carbon",
      (Unit == "m3") ~ "Water",
      (Unit == "m2a") ~ "Land",
      (Unit == "kg N-eq") ~ "Nitrogen",
      (Unit == "kg P-eq") ~ "Phosphorus"
    )) %>%
    select(FP, Category, Allocation.method, FP_val, SD, Unit)
}

#_____________________________________________________________________________________________________#
# Estimate off farm, feed-associated footprint
#_____________________________________________________________________________________________________#
estimate.feedFP <- function(LCA_data, Feed_data, FP_option, allocation){
  # FP options are "Carbon", "Water", "Nitrogen", "Phosphorus", and "Land" 
  # Filter to FP option
  Feed_data <- Feed_data %>% 
    filter(FP == FP_option) %>%
    filter(Allocation.method == allocation)
  
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
  estimate_FP <- LCA_data$Yield_t_per_Ha/LCA_data$harvest
  
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

