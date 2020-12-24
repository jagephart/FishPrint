# Functions

#_____________________________________________________________________________________________________#
# Clean LCA data
#_____________________________________________________________________________________________________#
clean.lca <- function(LCA_data){
  # Change columns to numeric where applicable
  LCA_data$Feed_soy_percent <- as.numeric(LCA_data$Feed_soy_percent)
  LCA_data$Feed_othercrops_percent <- as.numeric(LCA_data$Feed_othercrops_percent)
  LCA_data$Feed_FMFO_percent <- as.numeric(LCA_data$Feed_FMFO_percent)
  LCA_data$Feed_animal_percent <- as.numeric(LCA_data$Feed_animal_percent)
  LCA_data$Electricity_kwh <- as.numeric(LCA_data$Electricity_kwh)
  LCA_data$Diesel_L <- as.numeric(LCA_data$Diesel_L)
  LCA_data$Petrol_L <- as.numeric(LCA_data$Petrol_L)
  LCA_data$NaturalGas_L <- as.numeric(LCA_data$NaturalGas_L)
  LCA_data$Yield_t_per_Ha <- as.numeric(LCA_data$Yield_t_per_Ha)
  LCA_data$Yield_kg_per_m3 <- as.numeric(LCA_data$Yield_kg_per_m3)
  LCA_data$Grow_out_period_days <- as.numeric(LCA_data$Grow_out_period_days)
  LCA_data$FCR_overall <- as.numeric(LCA_data$FCR_overall)
  
  # Add country codes
  LCA_data$iso3c <- countrycode(LCA_data$Country, origin = "country.name", destination = "iso3c")
  
  # Remove experimental, polyculture, IAA, tuna, eel, and caviar studies
  LCA_data <- LCA_data %>%
    filter(Drop_study_flag == "") # ie, filter out Drop_study_flag %in% c(Experimental, Hypothetical, IAA, Polyculture, Tuna)
  
  LCA_data <- LCA_data %>%
    # Scale feed percents to sum to 100%
    # mutate(sum_percent = Feed_soy_percent+Feed_othercrops_percent+Feed_FMFO_percent+Feed_animal_percent) %>%
    # mutate(
    # Feed_soy_percent = Feed_soy_percent/sum_percent,
    # Feed_othercrops_percent = Feed_othercrops_percent/sum_percent,
    # Feed_FMFO_percent = Feed_FMFO_percent/sum_percent,
    # Feed_animal_percent = Feed_animal_percent/sum_percent
    # ) %>%
    # Convert all yield to the same units # FIX IT: Check this is done correctly
    mutate(Yield_m2_per_t = 
             ifelse(m2a_t > 0, m2a_t,
           ifelse(m3a_t  > 0, m3a_t/1, # Add in average depth when we have this; 1 is a placeholder
           ifelse(Yield_t_per_Ha > 0, 1/(Yield_t_per_Ha*1/10000),
           ifelse(Yield_kg_per_m3  > 0, 1/(Yield_kg_per_m3*1*(1/1000)), # Add in average depth when we have this; 1 is a placeholder
           NA))))) %>%
    # Create system group
    mutate(
      Production_system_group = case_when(
        (Production_system %in% c("Extensive raft culture", "Marine floating bag", "Integrated marine rafts", "Longline",
                                  "Bouchot culture", "Wooden stakes", "Long-lines", "Suspended baskets", 
                                  "Bottom planted", "On or off bottom")) ~ "On- and off-bottom",
        
        (Production_system %in% c("Intensive lake net-pen", "Marine cages", "Marine net-pen", "Reservoirs", "Lakes", 
                                  "Ponds, lakes, and reservoirs", "Lake-based net cages", "Offshore cages", "Floating cages", 
                                  "Net-pens", "Net pen", "Freshwater net pen (?)", "Saltwater net pen", 
                                  "Semi-intensive cages")) ~ "Cages & pens",
        
        (Production_system %in% c("Intensive pond", "Extensive pond polyculture", "Semi-intensive pond", "Extensive pond",
                                  "Earthen pond aquaculture", "Integrated pond, high input", "Integrated pond, medium inputs",
                                  "Solid-walled aquaculture system",  "Earthen pond aquaculture integrated with pigs",
                                  "Earthern ponds", "Earthern/concrete ponds", "Ponds", "Silvo pond", "Lined pond aquaculture",
                                  "Pond aquaculture", "Earthen pond monoculture", "Pond")) ~ "Ponds",
        
        (Production_system %in% c("Indoor recirculating", "Flow-through", "Land-based recirculating", "Onshore tanks",
                                  "Saltwater flow-through", "Freshwater flow-through", "Recirculating system", 
                                  "Land-based recirculating system", "Raceway", "Tanks / raceway", 
                                  "Semi-closed recirculating system", "Concrete tanks")) ~ "Recirculating and tanks",
        
        (Production_system %in% c("Unspecified", "", "Ponds / recirculating", "Ponds / pens")) ~ "not specified"
      ) 
      ) %>%
    # Convert "not specified" to NA
    mutate(Production_system_group = na_if(Production_system_group, "not specified")) %>%
    mutate(Intensity = na_if(Intensity, ""))
  
  # Create column clean_sci_name - use this as the "official" scientific name column
  # When no valid sci name exists just use common name
  # Simplify hybrid M. chrysops x M. saxatilis to its genus
  # Change outdated names (P. vannamei and P hypophthalmus)
  LCA_data <- LCA_data %>%
    mutate(Species.scientific.name = case_when(str_detect(Species.scientific.name, " spp\\.") ~ str_replace(Species.scientific.name, pattern = " spp\\.", replacement = " spp"),
                                       TRUE ~ Species.scientific.name)) %>%
    mutate(clean_sci_name = case_when(Species.common.name == "Carps" ~ "Cyprinidae", 
                                      Species.common.name == "Common carp (33%), Grass carp (21%), Crucian carp (9%), silver carp (9%), bighead carp (7%), and other carp (21%)" ~ "Cyprinidae",
                                      Species.common.name == "Grass carp (76%), bighead carp (5%), silver carp (8%), and crucian carp (12%)" ~ "Cyprinidae",
                                      Species.common.name == "Freshwater prawn" ~ "Freshwater prawn",
                                      Species.common.name == "Indo-Pacific swamp crab" ~ "Brachyura",
                                      Species.common.name == "River eels nei" ~ "Freshwater eels",
                                      Species.common.name == "Salmonids nei" ~ "Salmonidae",
                                      Species.common.name == "Striped bass" ~ "Morone saxatilis",
                                      Species.common.name == "Yellowtail_Seriola_Almaco jack" ~ "Seriola rivoliana",
                                      TRUE ~ Species.scientific.name)) %>%
    mutate(clean_sci_name = case_when(Species.scientific.name == "Chinese fed carp" ~ "Cyprinidae",
                                      Species.scientific.name == "Ctenopharyngodon idella; Carassius carassius" ~ "Cyprinidae",
                                      # Instead of Cyprinidae, be specific that this next group are Hypophtalmichthys so that they can be pulled into their own taxa group in add_taxa_group
                                      Species.scientific.name == "Hypophthalmichthys molitrix and Hypophthalmichthys nobilis" ~ "Mixed H. molitrix and H. nobilis",
                                      Species.scientific.name == "Morone chrysops x M. saxatilis" ~ "Morone hybrid",
                                      Species.scientific.name == "Labeo rohita and Catla Catla" ~ "Cyprinidae",
                                      Species.scientific.name == "Labeo rohita and Catla catla" ~ "Cyprinidae",
                                      Species.scientific.name == "Osteichthyes" ~ "Freshwater fishes",
                                      Species.scientific.name == "Penaeus vannamei" ~ "Litopenaeus vannamei",
                                      Species.scientific.name == "Pangasius hypophthalmus" ~ "Pangasianodon hypophthalmus",
                                      TRUE ~ clean_sci_name))

  # Check if any have unassigned clean_sci_name
  LCA_data %>%
    filter(clean_sci_name == "") %>%
    select(Species.common.name, Species.scientific.name, clean_sci_name) %>%
    unique()

  LCA_data <- LCA_data %>%
    # Divide by 5 for moist pellets
    mutate(FCR = case_when(Feed_type == "Moist pellet" ~ FCR_overall/5,
                           TRUE ~ FCR_overall)) %>%
    # Add 0s to Feed percentages where there should actually be 0s (instead of NAs)
    mutate(Feed_soy_percent = if_else(is.na(Feed_soy_percent) & (!is.na(Feed_othercrops_percent) | !is.na(Feed_FMFO_percent) | !is.na(Feed_animal_percent)), true = 0, false = Feed_soy_percent),
           Feed_othercrops_percent = if_else(is.na(Feed_othercrops_percent) & (!is.na(Feed_soy_percent) | !is.na(Feed_FMFO_percent) | !is.na(Feed_animal_percent)), true = 0, false = Feed_othercrops_percent),
           Feed_FMFO_percent = if_else(is.na(Feed_FMFO_percent) & (!is.na(Feed_soy_percent) | !is.na(Feed_othercrops_percent) | !is.na(Feed_animal_percent)), true = 0, false = Feed_FMFO_percent),
           Feed_animal_percent = if_else(is.na(Feed_animal_percent) & (!is.na(Feed_soy_percent) | !is.na(Feed_othercrops_percent) | !is.na(Feed_FMFO_percent)), true = 0, false = Feed_animal_percent)) %>%
    # Normalize the FINAL feed proportion values to be greater than 0 and no less than 0.01
    mutate(feed_soy_new = if_else(Feed_soy_percent == 0, true = 0.0105, false = Feed_soy_percent),
           feed_crops_new = if_else(Feed_othercrops_percent == 0, true = 0.0105, false = Feed_othercrops_percent),
           feed_fmfo_new = if_else(Feed_FMFO_percent == 0, true = 0.0105, false = Feed_FMFO_percent),
           feed_animal_new = if_else(Feed_animal_percent == 0, true = 0.0105, false = Feed_animal_percent)) %>%
    # Renomoralize values so they sum to 1
    mutate(sum_for_rescaling_feed = rowSums(select(., contains("new")))) %>%
    mutate(feed_soy_new = feed_soy_new / sum_for_rescaling_feed,
           feed_crops_new = feed_crops_new / sum_for_rescaling_feed,
           feed_fmfo_new = feed_fmfo_new / sum_for_rescaling_feed,
           feed_animal_new = feed_animal_new / sum_for_rescaling_feed) %>%
    # Add 0s to Energy data where there should actually be 0s (instead of NAs)
    mutate(Electricity_kwh = if_else(is.na(Electricity_kwh) & (!is.na(Diesel_L) | !is.na(Petrol_L) | !is.na(NaturalGas_L)), true = 0, false = Electricity_kwh),
           Diesel_L = if_else(is.na(Diesel_L) & (!is.na(Electricity_kwh) | !is.na(Petrol_L) | !is.na(NaturalGas_L)), true = 0, false = Diesel_L),
           Petrol_L = if_else(is.na(Petrol_L) & (!is.na(Electricity_kwh) | !is.na(Diesel_L) | !is.na(NaturalGas_L)), true = 0, false = Petrol_L),
           NaturalGas_L = if_else(is.na(NaturalGas_L) & (!is.na(Electricity_kwh) | !is.na(Diesel_L) | !is.na(Petrol_L)), true = 0, false = NaturalGas_L))
  
  # Option 1:
  # Replicate data by sqrt(n_farms):
  # LCA_data <- LCA_data %>%
  #   mutate(Sample_replication = round(sqrt(Sample_size_n_farms)))
  # LCA_data <- as.data.frame(lapply(LCA_data, rep, LCA_data$Sample_replication))
  
  # Options 2: Replicate data by "n" farms
  LCA_data <- as.data.frame(lapply(LCA_data, rep, LCA_data$Sample_size_n_farms))
  
  # FINAL STEP: Create study_id column, use this in all analyses to bind predictions from multiple models back together
  LCA_data <- LCA_data %>%
    mutate(study_id = row_number()) %>%
    # Drop columns we no longer need
    # Subset to columns to keep
    select(SeaWEED.ID, Source, study_id, Country, iso3c, Common.Name = Species.common.name, Scientific.Name = Species.scientific.name, clean_sci_name, 
           Production_system_group, Intensity, Product, Yield_m2_per_t, Grow_out_period_days, FCR = FCR_overall, 
           Feed_type, feed_soy_new, feed_crops_new, feed_fmfo_new, feed_animal_new,
           Electricity_kwh, Diesel_L, Petrol_L, NaturalGas_L)
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
  
  # Which sci_names have NA for isscaap_group after joining
  lca_dat_out %>%
    filter(is.na(isscaap_group)) %>%
    select(clean_sci_name, isscaap_group) %>%
    unique()

  lca_dat_out <- lca_dat_out %>% 
    # First pass is to assign NAs to ISSCAAP group
    mutate(isscaap_group = case_when(clean_sci_name %in% c("Mixed Hypophthalmichthys molitrix and H. nobilis",
                                                           "Ctenopharyngodon idella") ~ "Carps, barbels and other cyprinids",
                                     clean_sci_name %in% c("Freshwater prawn", "Macrobrachium amazonicum") ~ "Freshwater crustaceans",
                                     clean_sci_name %in% c("Morone hybrid") ~ "Miscellaneous diadromous fishes",
                                     clean_sci_name %in% c("Gracilaria chilensis") ~ "Red seaweeds",
                                     clean_sci_name %in% c("Litopenaeus vannamei") ~ "Shrimps, prawns",
                                     TRUE ~ isscaap_group)) 
  
  # Inspect assignment to isscaap_group
  lca_dat_out %>%
    select(clean_sci_name, isscaap_group) %>%
    unique() %>%
    arrange(isscaap_group)
  
  
  # Once all groups have an isscaap group, Create taxa groups either manually or from isscaap_group
  lca_dat_out <- lca_dat_out %>%
    # Split up carps
    mutate(taxa_group_name = case_when(clean_sci_name %in% c("Mixed H. molitrix and H. nobilis", 
                                                             "Hypophthalmichthys molitrix") ~ "Silver and bighead carp",
                                       clean_sci_name %in% c("Carassius carassius", 
                                                             "Ctenopharyngodon idella", 
                                                             "Cyprinidae", 
                                                             "Cyprinus carpio", 
                                                             "Mixed Ctenopharyngodon idella and Carassius carassius",
                                                             "Mixed Labeo rohita and Catla catla") ~ "Other carps, barbels and cyprinids",
                                       # Remove catfish from Miscellaneous freshwater fishes as its own grouping
                                       clean_sci_name %in% c("Clarias batrachus", "Clarias gariepinus", "Pangasianodon hypophthalmus", "Pangasius spp") ~ "Catfish",
                                       # Split salmons and trouts
                                       str_detect(Common.Name, "salmon|Salmonids") ~ "Salmon",
                                       str_detect(Common.Name, "trout") ~ "Trout",
                                       str_detect(Common.Name, "char") ~ "Miscellaneous diadromous fishes",
                                       # Remove milkfish from misc diadromous fishes
                                       clean_sci_name == "Chanos chanos" ~ "Milkfish",
                                       # Combine groups into misc marine fishes
                                       isscaap_group %in% c("Miscellaneous coastal fishes", "Miscellaneous demersal fishes", "Miscellaneous pelagic fishes", "Flounders, halibuts, soles") ~ "Miscellaneous marine fishes",
                                       # Combine seaweeds into aquatic plants
                                       isscaap_group %in% c("Brown seaweeds", "Red seaweeds") ~ "Aquatic plants",
                                       # Combine bivalves
                                       isscaap_group %in% c("Mussels", "Oysters") ~ "Bivalves",
                                       TRUE ~ isscaap_group))
  
  # Inspect assignment of taxa_group_name
  lca_dat_out %>%
    select(clean_sci_name, isscaap_group, taxa_group_name) %>%
    unique() %>%
    arrange(taxa_group_name)
  
  # Clean up taxa_group_names
  lca_dat_out <- lca_dat_out %>%
    mutate(taxa = case_when(taxa_group_name == "Aquatic plants" ~ "plants",
                            taxa_group_name == "Bivalves" ~ "bivalves",
                            taxa_group_name == "Catfish" ~ "catfish",
                            taxa_group_name == "Crabs, sea-spiders" ~ "crabs",
                            taxa_group_name == "Freshwater crustaceans" ~ "fresh_crust",
                            taxa_group_name == "Milkfish" ~ "milkfish",
                            taxa_group_name == "Miscellaneous diadromous fishes" ~ "misc_diad",
                            taxa_group_name == "Miscellaneous freshwater fishes" ~ "misc_fresh",
                            taxa_group_name == "Miscellaneous marine fishes" ~ "misc_marine",
                            taxa_group_name == "Other carps, barbels and cyprinids" ~ "oth_carp",
                            taxa_group_name == "Salmon" ~ "salmon",
                            taxa_group_name == "Shrimps, prawns" ~ "shrimp",
                            taxa_group_name == "Silver and bighead carp" ~ "hypoph_carp",
                            taxa_group_name == "Tilapias and other cichlids" ~ "tilapia",
                            taxa_group_name == "Trout" ~ "trout",
                            TRUE ~ "unassigned"))

  
}


#_____________________________________________________________________________________________________#
# Replicate data based on Sample_size_n_farms column
#_____________________________________________________________________________________________________#

rep_data <- function(lca_dat_clean){
  # No longer need to clean sample size column - now has column Sample_size_n_farms
  # lca_dat_clean_rep <- lca_dat_clean %>%
  #   # Clean up sample size column
  #   # First ignore numbers that are percentages, then find and extract any numbers, then fill the rest in with 1s
  #   mutate(clean_sample_size = case_when(str_detect(Sample_size, "%") ~ 1,
  #                                        str_detect(Sample_size, "[0-9]+") ~ as.numeric(str_extract(Sample_size, pattern = "[0-9]+")),
  #                                        TRUE ~ 1)) 
  
  lca_dat_clean_rep <- as.data.frame(lapply(lca_dat_clean, rep, lca_dat_clean$Sample_size_n_farms))
}

#_____________________________________________________________________________________________________#
# Clean feed footprint data
#_____________________________________________________________________________________________________#
clean.feedFP <- function(feedFP_data){
  feedFP_data <- feedFP_data %>%
    #filter(Units != "kg PO4-eq") %>% # Not currently using
    group_by(Input.type, Impact.category, Units, Allocation) %>%
    # Revise this to change from arithmetic mean to a weighted mean
    summarise(impact_val = mean(Value, na.rm = TRUE), SD = sd(Value, na.rm = TRUE), .groups = 'drop') %>% 
    mutate(impact_factor = case_when(
      (Units == "kg CO2-eq / kg") ~ "Carbon",
      (Units == "m3 / kg") ~ "Water",
      (Units == "m2a / kg") ~ "Land",
      (Units == "kg N-eq / kg") ~ "Nitrogen",
      (Units == "kg P-eq / kg") ~ "Phosphorus"
    )) %>%
    select(impact_factor, Impact.category, Allocation, impact_val, SD, Units)
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

# Refer to clean_fish_NP.R for functions to estimate fish N and P content



#_____________________________________________________________________________________________________#
# Create visualizations of brms gamma (or dirichlet) regression outputs
#_____________________________________________________________________________________________________#

plot_for_si <- function(name_of_fit, name_of_data, name_of_var, regression_type = "gamma"){
  # Mac
  datadir <- "/Volumes/jgephart/BFA Environment 2/Data"
  outdir <- "/Volumes/jgephart/BFA Environment 2/Outputs"
  
  library(tidyverse)
  library(modelr)
  library(ggdist)
  library(tidybayes)
  library(ggplot2)
  library(cowplot)
  library(rstan)
  library(brms)
  library(ggrepel)
  library(RColorBrewer)
  
  # SET THEME
  plot_theme <- theme(title = element_text(size = 20),
                      axis.title.x = element_text(size = 20),
                      axis.text=element_text(size=20, color = "black"))
  
  compiled_dat_clean <- read.csv(file.path(datadir, "lca_clean_with_groups.csv"))
  brms_output <- get(name_of_fit)
  full_dat <- get(name_of_data)
  
  # FCR
  # name_of_fit <- "fit_fcr_no_na"
  # name_of_data <- "full_fcr_dat"
  # name_of_var <- "fcr"
  
  # Yield
  # name_of_fit <- "fit_yield_no_na"
  # name_of_data <- "full_yield_dat"
  # name_of_var <- "yield"
  
  # Electricity
  # name_of_fit <- "fit_electric_no_na"
  # name_of_data <- "full_electric_dat"
  # name_of_var <- "electric"
  
  # Diesel
  # name_of_fit <- "fit_diesel_no_na"
  # name_of_data <- "full_diesel_dat"
  # name_of_var <- "diesel"
  
  # Petrol
  # name_of_fit <- "fit_petrol_no_na"
  # name_of_data <- "full_petrol_dat"
  # name_of_var <- "petrol"
  
  # Natural gas
  # name_of_fit <- "fit_natgas_no_na"
  # name_of_data <- "full_natgas_dat"
  # name_of_var <- "natgas"
  
  if (regression_type == "gamma") {
    
    # Combine modeled data (both data and predictions) with the full clean LCA dataset and output this
    dat_for_si <- compiled_dat_clean %>%
      left_join(full_dat, by = c("study_id", "clean_sci_name", "taxa", "Intensity" = "intensity", "Production_system_group" = "system"))
    write.csv(dat_for_si, file.path(outdir, paste("lca_clean_with_model_predictions-", name_of_var, ".csv", sep = "")), row.names = FALSE)
    
    summary(brms_output)
    
    
    
    # Default posterior predictive check is a density plot:
    # Specify response variable in resp
    pp_check(brms_output, resp = "y", nsamples = 50) + 
      ggtitle(paste("Posterior predictive check: ", name_of_var, sep = ""))
    ggsave(filename = file.path(outdir, paste("plot_gamma-regression_", name_of_var, "_post-pred-checks_density.png", sep = "")), width  = 8, height = 11.5)
    pp_check(brms_output, type = "scatter_avg", nsamples = 1000, resp = "y") +
      ggtitle(paste("Posterior predictive check: ", name_of_var, sep = ""))
    ggsave(filename = file.path(outdir, paste("plot_gamma-regression_", name_of_var, "_post-pred-checks_scatter.png", sep = "")), width  = 8, height = 11.5)
    
    # Other posterior predictive checks
    # pp_check(brms_output, type = "error_hist", nsamples = 5, resp = "y")
    # pp_check(brms_output, type = "stat_2d", resp = "y")
    # pp_check(brms_output, type = "stat", resp = "y")
    
    get_variables(brms_output)
    
    # Plot coefficients
    p <- mcmc_plot(brms_output, pars = c("^b_"))
    
    coeff_data <- p$data %>% 
      mutate(effect = case_when(m > 0 ~ "positive",
                                m < 0 ~ "negative",
                                TRUE ~ "none"))
    
    ggplot(data = coeff_data) +
      geom_segment(data = coeff_data, mapping = aes(xend = hh, yend = parameter, x = ll, y = parameter)) +
      geom_segment(data = coeff_data, mapping = aes(xend = h, yend = parameter, x = l, y = parameter), size = 2) +
      geom_point(aes(x = m, y = parameter, color = effect), size = 3) +
      labs(x = "", y = "", title = paste("Coefficients for ", name_of_var, sep = "")) +
      theme_classic() +
      plot_theme
    
    ggsave(filename = file.path(outdir, paste("plot_gamma-regression_", name_of_var, "_coeffs.png", sep = "")), height = 8.5, width = 11)
    
    # Boxplot of data + predictions per taxa group
    ggplot(full_dat, aes(x = !!sym(name_of_var), y = taxa)) + 
      geom_boxplot(outlier.shape = NA) + geom_jitter(aes(color = data_type)) +
      theme_classic() +
      plot_theme
    
    ggsave(filename = file.path(outdir, paste("plot-data_", name_of_var, "-ghg_tx-level.png", sep = "")), width = 11, height = 8.5)
    
    # Make separate plot for each taxa group
    for (i in 1:length(unique(full_dat$taxa))){
      taxa_i <- unique(full_dat$taxa)[i]
      dat_taxa_i <- full_dat %>%
        filter(taxa == taxa_i) %>%
        mutate(rowname = row_number()) %>%
        replace_na(replace = list(taxa = "unknown", intensity = "unknown", system = "unknown"))
      p <- ggplot() +
        geom_pointinterval(aes(y = rowname, x = !!sym(name_of_var), xmin = .lower, xmax = .upper, shape = system, point_color = intensity), size = 2, data = dat_taxa_i) +
        geom_point(aes(y = rowname, x = !!sym(name_of_var), shape = system, color = intensity), data = dat_taxa_i, size = 3) +
        #scale_interval_shape(drop = FALSE) +
        scale_shape_discrete(drop = FALSE) +
        #coord_cartesian(xlim = c(0, 10)) +
        theme_classic() +
        labs(x = name_of_var, y = "", title = taxa_i) +
        theme(axis.text = element_text(size = 16),
              axis.title = element_text(size = 20),
              axis.text.y = element_blank())
      plot(p)
      file_i <- paste("plot_gamma-regression_", name_of_var, "_missing-dat-predictions_taxa-", taxa_i, ".png", sep = "")
      ggsave(filename = file.path(outdir, file_i), width = 8, height = 11.5)
    }
  }
  
  if (regression_type == "dirichlet"){
    full_dat_for_merge <- full_dat %>%
      select(-c(.row, rowname, .lower, .upper, .width, .point, .interval)) %>% 
      pivot_wider(names_from = .category, values_from = feed_proportion) %>%
      group_by(study_id, clean_sci_name, taxa, intensity, system, data_type) %>%
      # SUM to collapse into a single row
      mutate(feed_soy = sum(feed_soy, na.rm = TRUE),
             feed_crops = sum(feed_crops, na.rm = TRUE),
             feed_fmfo = sum(feed_fmfo, na.rm = TRUE),
             feed_animal = sum(feed_animal, na.rm = TRUE)) %>%
      distinct()
    
    # Combine modeled data (both data and predictions) with the full clean LCA dataset and output this
    dat_for_si <- compiled_dat_clean %>%
      left_join(full_dat_for_merge, by = c("study_id", "clean_sci_name", "taxa", "Intensity" = "intensity", "Production_system_group" = "system"))
    write.csv(dat_for_si, file.path(outdir, paste("lca_clean_with_model_predictions-", name_of_var, ".csv", sep = "")), row.names = FALSE)
    
    summary(brms_output)
    
    ## PP_CHECK not implemented for dirichlet regression models
    # pp_check(brms_output, nsamples = 50) + 
    #   ggtitle("Posterior predictive check")
    # ggsave(filename = file.path(outdir, "plot_gamma-regression_post-pred-checks_density.png"), width  = 11.5, height = 8)
    # pp_check(brms_output, type = "error_hist", nsamples = 5)
    # pp_check(brms_output, type = "scatter_avg", nsamples = 1000)
    # pp_check(brms_output, type = "stat_2d")
    # pp_check(brms_output, type = "stat")
    
    get_variables(brms_output)
    
    # SET THEME
    plot_theme <- theme(title = element_text(size = 20),
                        axis.title.x = element_text(size = 20),
                        axis.text=element_text(size=20, color = "black"))
    
    # Plot coefficients (separate for each feed component)
    feed_coeffs <- c("b_mufeedcrops", "b_mufeedfmfo", "b_mufeedanimal")
    for (i in 1:length(feed_coeffs)){
      p <- mcmc_plot(brms_output, pars = feed_coeffs[i])
      coeff_data <- p$data %>% 
        mutate(effect = case_when(m > 0 ~ "positive",
                                  m < 0 ~ "negative",
                                  TRUE ~ "none"))
      p_custom <- ggplot(data = coeff_data) +
        geom_segment(data = coeff_data, mapping = aes(xend = hh, yend = parameter, x = ll, y = parameter)) +
        geom_segment(data = coeff_data, mapping = aes(xend = h, yend = parameter, x = l, y = parameter), size = 2) +
        geom_point(aes(x = m, y = parameter, color = effect), size = 3) +
        labs(x = "", y = "") +
        theme_classic() +
        plot_theme
      plot(p_custom)
      
      ggsave(filename = file.path(outdir, paste("plot_dirichlet-regression_coeffs_", feed_coeffs[i], ".png", sep = "")), height = 8.5, width = 11)
    }
    
    # feed proportion all taxa
    feed_vars <- c("soy", "crops", "fmfo", "animal")
    for (i in 1:length(feed_vars)) {
      p <- ggplot(data = full_feed_dat %>% 
                    filter(str_detect(.category, feed_vars[i])), aes(y = taxa, x = feed_proportion)) +
        geom_boxplot(outlier.shape = NA) +
        #geom_violin(aes(color = taxa), scale = "width") +
        geom_jitter(aes(color = data_type), size = 3) +
        theme_classic() +
        plot_theme +
        labs(title = paste("Boxplots of ", feed_vars[i], " feed proportions", sep = ""),
             x = "",
             y = "")  +
        theme(axis.text.x = element_text(hjust = 1))
      print(p)
      ggsave(file.path(outdir, paste("boxplot_feed-prop_", feed_vars[i], ".png", sep = "")), height = 8, width = 11.5)
    }
    
    # Feed proportions with facet_wrap
    p <- ggplot(data = full_feed_dat, aes(y = taxa, x = feed_proportion)) +
      geom_boxplot(outlier.shape = NA) +
      #geom_violin(aes(color = taxa), scale = "width") +
      geom_jitter(aes(color = data_type), size = 3) +
      theme_classic() +
      plot_theme +
      labs(title = paste("Boxplots of all feed proportions", sep = ""),
           x = "",
           y = "")  +
      facet_wrap(~.category, nrow = 1)
    print(p)
    ggsave(file.path(outdir, "boxplot_feed-prop_all-feeds.png"), height = 8, width = 11.5)
  }
  
}
