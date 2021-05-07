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
  LCA_data <- LCA_data %>%
    mutate(Sample_replication = round(sqrt(Sample_size_n_farms)))
  LCA_data <- as.data.frame(lapply(LCA_data, rep, LCA_data$Sample_replication))
  
  # Options 2: Replicate data by "n" farms
  # LCA_data <- as.data.frame(lapply(LCA_data, rep, LCA_data$Sample_size_n_farms))
  
  # FINAL STEP: Create study_id column, use this in all analyses to bind predictions from multiple models back together
  LCA_data <- LCA_data %>%
    mutate(study_id = row_number()) %>%
    # Drop columns we no longer need
    # Subset to columns to keep
    select(SeaWEED.ID, Source, Sample_size_n_farms, study_id, Country, iso3c, Common.Name = Species.common.name, Scientific.Name = Species.scientific.name, clean_sci_name, 
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
# Clean aquaculture priors data, add taxa grouping
# Can ignore warning: NAs introduced by coercion (inserts NAs for blank cells)
#_____________________________________________________________________________________________________#
clean_priors <- function(priors_dat){
  priors_csv <- read.csv(file.path(datadir, priors_dat)) %>%
    select(Group.name, Mean.Annual.Yield.t.ha, Ave.FCR, Upper.FCR, Lower.FCR) %>%
    # Deal with "-" entries
    mutate(across(everything(), str_replace, pattern = "-", replacement = "")) %>%
    # Convert data cols to numeric
    mutate(across(contains(c("Yield", "FCR")), as.numeric)) %>%
    # Add taxa group name:
    mutate(taxa = case_when(Group.name %in% c("Mussels", "Oysters") ~ "bivalves",
                            Group.name %in% c("Catfish", "Miscellaneous freshwater fishes") ~ "catfish",
                            Group.name == "Silver and bighead carp" ~ "hypoph_carp",
                            Group.name == "Milkfish" ~ "milkfish",
                            Group.name == "Miscellaneous diadromous fishes" ~ "misc_diad",
                            Group.name == "Other miscellaneous marine fishes" ~ "misc_marine",
                            Group.name %in% c("Common carp", "Crucian carp", "Grass carp", "Other Carps, barbels and other cyprinids") ~ "oth_carp",
                            Group.name == "Aquatic plants" ~ "plants",
                            Group.name == "Salmon" ~ "salmon",
                            Group.name == "Shrimps, prawns" ~ "shrimp",
                            Group.name == "Tilapias and other cichlids" ~ "tilapia",
                            Group.name == "Trouts" ~ "trout",
                            TRUE ~ "not_grouped")) %>%
    group_by(taxa) %>%
    summarise(Mean.Annual.Yield.t.ha = mean(Mean.Annual.Yield.t.ha, na.rm = TRUE), 
              Ave.FCR = mean(Ave.FCR, na.rm = TRUE),
              Upper.FCR = mean(Upper.FCR, na.rm = TRUE),
              Lower.FCR = mean(Lower.FCR, na.rm = TRUE)) %>%
    ungroup() %>%
    # Manually add other freshwater fishes to catfish
    filter(taxa != "not_grouped") %>%
    arrange(taxa) %>%
    # Convert units for Yield
    mutate(Mean.Annual.Yield.m2.per.tonne = 1/(Mean.Annual.Yield.t.ha * 1/10000))
  
}

#_____________________________________________________________________________________________________#
# Clean WILD CAPTURE priors data, add taxa grouping
# Can ignore warning: NAs introduced by coercion (inserts NAs for blank cells)
#_____________________________________________________________________________________________________#
clean_wild_priors <- function(wild_priors_dat = "Priors - Capture.csv"){
  priors_csv <- read.csv(file.path(datadir, wild_priors_dat)) %>%
    select(Group.name, Mean.FUI, Median.FUI) %>%
    # Convert to GHG values (as in Parker et al. 2018, see Rob's note in CSV file)
    mutate(Mean.GHG = Mean.FUI*3.3*1.33,
           Median.GHG = Median.FUI*3.3*1.33)  %>%
    mutate(Group.name = case_when(str_detect(Group.name, "Lobsters") ~ "Lobsters",
                                  str_detect(Group.name, "Shrimps") ~ "Shrimps",
                                  TRUE ~ Group.name)) %>%
    drop_na()
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
# Calculate the weighted averages for the feed components
#_____________________________________________________________________________________________________#

# Calculate soy weightings
# Since the specific soy products don't map onto the FAO items, use all soy exports for all soy products
calc_soy_weights <- function(faostat = faostat, feed_fp = feed_fp, deforestation_free = FALSE){
  weightings <-  faostat %>% 
    filter(Unit == "tonnes") %>% 
    filter(Item %in% c("Soybeans", "Cake, soybeans", "Oil, soybean")) %>%
    {if (deforestation_free == TRUE) filter(., iso3c != "BRA")
      else .} %>%
    group_by(iso3c) %>%
    summarise(Exports = sum(Value, na.rm = TRUE)) %>%
    filter(Exports > 0) %>% 
    left_join(feed_fp %>% filter(Input.type == "Soy"), by = c("iso3c")) %>%
    filter(is.na(Input.type) == FALSE) %>%
    group_by(iso3c) %>%
    summarise(Exports = sum(Exports, na.rm = TRUE)) %>%
    mutate(weighting = Exports/sum(Exports)) %>%
    select(iso3c, weighting)
  
  weighted_soy <- feed_fp %>% 
    filter(Input.type == "Soy") %>%
    left_join(weightings, by = "iso3c") %>%
    filter(is.na(weighting) == FALSE) %>%
    group_by(Input.type, Input, Impact.category, Allocation, Units) %>%
    # Normalize weightings to sum to 1
    mutate(reweighting = weighting/sum(weighting, na.rm = TRUE)) %>%
    summarise(Value = sum(Value * reweighting)) %>%
    # If weighting ingredient types, do here along with country weightings
    group_by(Input.type, Impact.category, Allocation, Units) %>%
    summarise(ave_stressor = mean(Value, na.rm = TRUE))
  
  return(weighted_soy)
}
  
# Calculate crop weightings
calc_crop_weights <- function(faostat = faostat, feed_fp = feed_fp, deforestation_free = FALSE){
  weightings <- faostat %>% 
    filter(Unit == "tonnes") %>% 
    {if (deforestation_free == TRUE) filter(., iso3c != "ARG")
      else .} %>%
    group_by(iso3c, Item) %>%
    summarise(Exports = sum(Value, na.rm = TRUE)) %>%
    filter(Exports > 0) %>% 
    mutate(Input = case_when(
      (Item %in% c("Cassava Equivalent")) ~ "Cassava",
      (Item %in% c("Maize")) ~ "Maize",
      (Item %in% c("Cake, maize")) ~ "Corn gluten meal",
      (Item %in% c("Cake, groundnuts")) ~ "Peanut meal",
      (Item %in% c("Rape and Mustard Oils")) ~ "Rapeseed oil",
      (Item %in% c("Cake, rapeseed")) ~ "Rapeseed meal",
      (Item %in% c("Wheat")) ~ "Wheat",
      (Item %in% c("Bran, wheat")) ~ "Wheat bran",
      (Item %in% c("Rice")) ~ "Rice bran",
      (Item %in% c("Cake, sunflower")) ~ "Sunflower meal"
    )) %>%
    filter(is.na(Input) == FALSE) %>% 
    left_join(feed_fp %>% filter(Input.type == "Crop"), by = c("iso3c", "Input")) %>%
    filter(is.na(Input.type) == FALSE) %>%
    ungroup() %>%
    group_by(Input, iso3c) %>%
    summarise(Exports = sum(Exports, na.rm = TRUE)) %>%
    mutate(weighting = Exports/sum(Exports)) %>%
    select(iso3c, Input, weighting)
  
  weighted_crop <- feed_fp %>% 
    filter(Input.type == "Crop") %>%
    left_join(weightings, by = c("iso3c", "Input")) %>%
    filter(is.na(weighting) == FALSE) %>% 
    # If weighting soy ingredient types, do here along with country weightings
    group_by(Input.type, Input, Impact.category, Allocation, Units) %>%
    # Normalize weightings to sum to 1
    mutate(reweighting = weighting/sum(weighting, na.rm = TRUE)) %>%
    summarise(Value = sum(Value * reweighting)) %>%
    # If weighting ingredient types, do here along with country weightings
    group_by(Input.type, Impact.category, Allocation, Units) %>%
    summarise(ave_stressor = mean(Value, na.rm = TRUE))
  
  return(weighted_crop)
}

# Calculate chicken weightings
calc_chicken_weights <- function(faostat = faostat, feed_fp = feed_fp){
  
  weightings <-  faostat %>%
    filter(Unit == "tonnes") %>%
    filter(Item %in% c("Poultry Meat")) %>%
    mutate(iso3c = ifelse(iso3c %in% c("FRA", "ITA"), "EUR", iso3c)) %>%
    group_by(iso3c) %>%
    summarise(Exports = sum(Value, na.rm = TRUE)) %>%
    filter(Exports > 0) %>%
    left_join(feed_fp %>%
                filter(Input %in% c("Chicken by-product meal", "Chicken by-product oil")), by = c("iso3c")) %>%
    filter(is.na(Input.type) == FALSE) %>%
    group_by(iso3c) %>%
    summarise(Exports = sum(Exports, na.rm = TRUE)) %>%
    mutate(weighting = Exports/sum(Exports)) %>%
    select(iso3c, weighting)
  
  weighted_chicken <- feed_fp %>% 
    filter(Input %in% c("Chicken by-product meal", "Chicken by-product oil")) %>%
    left_join(weightings, by = "iso3c") %>%
    filter(is.na(weighting) == FALSE) %>% 
    # If weighting soy ingredient types, do here along with country weightings
    group_by(Input.type, Input, Impact.category, Allocation, Units) %>%  
    # Normalize weightings to sum to 1
    mutate(reweighting = weighting/sum(weighting, na.rm = TRUE)) %>%
    summarise(Value = sum(Value * reweighting)) %>%
    # If weighting ingredient types, do here along with country weightings
    group_by(Input.type, Impact.category, Allocation, Units) %>%
    summarise(ave_stressor = mean(Value, na.rm = TRUE))
  
  return(weighted_chicken)
}

# Calculate fishery weights
calc_fishery_weights <- function(fmfo_prod = fmfo_prod, feed_fp = feed_fp){
  weighted_fishery <- feed_fp %>%
    filter(Input.type == c("Fishery")) %>% 
    left_join(fmfo_prod, by = c("Input" = "Name")) %>%
    group_by(Input.type, Impact.category, Allocation, Units) %>%
    # Normalize weightings to sum to 1
    mutate(reweighting = Weighting/sum(Weighting, na.rm = TRUE)) %>%
    summarise(Value = sum(Value * reweighting)) %>%
    group_by(Impact.category, Allocation, Units) %>%
    summarise(ave_stressor = mean(Value, na.rm = TRUE))
  return(weighted_fishery)
}

# Calculate fish byproduct weights
calc_byproduct_weights <- function(fmfo_prod = fmfo_prod, feed_fp = feed_fp){
  weighted_fishbyproduct <- feed_fp %>%
    filter(Input.type == c("Fishery by-product")) %>% 
    left_join(fmfo_prod, by = c("Input" = "Name")) %>%
    group_by(Input.type,  Impact.category, Allocation, Units) %>%
    # Normalize weightings to sum to 1
    mutate(reweighting = Weighting/sum(Weighting, na.rm = TRUE)) %>%
    summarise(Value = sum(Value * reweighting)) %>%
    group_by(Impact.category, Allocation, Units) %>%
    summarise(ave_stressor = mean(Value, na.rm = TRUE))
  return(weighted_fishbyproduct)
}

# Combine fishery and byproduct weights
combine_fish_weights <- function(weighted_fishery = weighted_fishery, weighted_fishbyproduct = weighted_fishbyproduct){
  weighted_fish <- weighted_fishery %>% 
    left_join(weighted_fishbyproduct, by = c("Impact.category", "Allocation", "Units")) %>%
    mutate(ave_stressor = 0.675*ave_stressor.x + 0.325*ave_stressor.y) %>%
    select(-ave_stressor.x, -ave_stressor.y)
  return(weighted_fish)
}

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

#_______________________________________________________________________________________________________________________#
# Change in stressors with change in parameters function
#_______________________________________________________________________________________________________________________#

stressor_sensitivity <- function(data_lca, data_feed_stressors, data_feed_NP, data_energy,
                                 delta_FCR = 0, delta_soy = 0, delta_crops = 0, delta_animal = 0, delta_fmfo = 0,
                                 delta_electricity = 0, delta_diesel = 0, delta_petrol = 0, delta_natgas = 0, delta_yield = 0,
                                 # Optional changes to constants 
                                 land_change_free_soy = NULL, land_change_free_crops = NULL, 
                                 fmfo_with_land_change_free_crops = NULL,
                                 fishery_byproducts = NULL, low_impact_fishery_byproducts = NULL,
                                 fmfo_low_impact_fishery_byproducts = NULL, zero_emissions_electric = FALSE){
  # Set feed component stressor constants
  soy_stressor_ghg <- data_feed_stressors %>% filter(stressor == "ghg" & feed_type == "soy") %>% pull(ave_stressor)
  soy_stressor_N <- data_feed_stressors %>% filter(stressor == "N" & feed_type == "soy") %>% pull(ave_stressor)
  soy_stressor_P <- data_feed_stressors %>% filter(stressor == "P" & feed_type == "soy") %>% pull(ave_stressor)
  soy_stressor_land <- data_feed_stressors %>% filter(stressor == "land" & feed_type == "soy") %>% pull(ave_stressor)
  soy_stressor_water <- data_feed_stressors %>% filter(stressor == "water" & feed_type == "soy") %>% pull(ave_stressor)
  
  if(is.null(land_change_free_soy) == FALSE){
    soy_stressor_ghg <- land_change_free_soy %>% filter(stressor == "ghg" & feed_type == "soy") %>% pull(ave_stressor)
    soy_stressor_N <- land_change_free_soy %>% filter(stressor == "N" & feed_type == "soy") %>% pull(ave_stressor)
    soy_stressor_P <- land_change_free_soy %>% filter(stressor == "P" & feed_type == "soy") %>% pull(ave_stressor)
    soy_stressor_land <- land_change_free_soy %>% filter(stressor == "land" & feed_type == "soy") %>% pull(ave_stressor)
    soy_stressor_water <- land_change_free_soy %>% filter(stressor == "water" & feed_type == "soy") %>% pull(ave_stressor)
  }else{}
  
  crops_stressor_ghg <- data_feed_stressors %>% filter(stressor == "ghg" & feed_type == "crops") %>% pull(ave_stressor)
  crops_stressor_N <- data_feed_stressors %>% filter(stressor == "N" & feed_type == "crops") %>% pull(ave_stressor)
  crops_stressor_P <- data_feed_stressors %>% filter(stressor == "P" & feed_type == "crops") %>% pull(ave_stressor)
  crops_stressor_land <- data_feed_stressors %>% filter(stressor == "land" & feed_type == "crops") %>% pull(ave_stressor)
  crops_stressor_water <- data_feed_stressors %>% filter(stressor == "water" & feed_type == "crops") %>% pull(ave_stressor)
  
  if(is.null(land_change_free_crops) == FALSE){
    crops_stressor_ghg <- land_change_free_crops %>% filter(stressor == "ghg" & feed_type == "crops") %>% pull(ave_stressor)
    crops_stressor_N <- land_change_free_crops %>% filter(stressor == "N" & feed_type == "crops") %>% pull(ave_stressor)
    crops_stressor_P <- land_change_free_crops %>% filter(stressor == "P" & feed_type == "crops") %>% pull(ave_stressor)
    crops_stressor_land <- land_change_free_crops %>% filter(stressor == "land" & feed_type == "crops") %>% pull(ave_stressor)
    crops_stressor_water <- land_change_free_crops %>% filter(stressor == "water" & feed_type == "crops") %>% pull(ave_stressor)
  }else{}
  
  animal_stressor_ghg <- data_feed_stressors %>% filter(stressor == "ghg" & feed_type == "animal") %>% pull(ave_stressor)
  animal_stressor_N <- data_feed_stressors %>% filter(stressor == "N" & feed_type == "animal") %>% pull(ave_stressor)
  animal_stressor_P <- data_feed_stressors %>% filter(stressor == "P" & feed_type == "animal") %>% pull(ave_stressor)
  animal_stressor_land <- data_feed_stressors %>% filter(stressor == "land" & feed_type == "animal") %>% pull(ave_stressor)
  animal_stressor_water <- data_feed_stressors %>% filter(stressor == "water" & feed_type == "animal") %>% pull(ave_stressor)
  
  fmfo_stressor_ghg <- data_feed_stressors %>% filter(stressor == "ghg" & feed_type == "fmfo") %>% pull(ave_stressor)
  fmfo_stressor_N <- data_feed_stressors %>% filter(stressor == "N" & feed_type == "fmfo") %>% pull(ave_stressor)
  fmfo_stressor_P <- data_feed_stressors %>% filter(stressor == "P" & feed_type == "fmfo") %>% pull(ave_stressor)
  fmfo_stressor_land <- data_feed_stressors %>% filter(stressor == "land" & feed_type == "fmfo") %>% pull(ave_stressor)
  fmfo_stressor_water <- data_feed_stressors %>% filter(stressor == "water" & feed_type == "fmfo") %>% pull(ave_stressor)
  
  if(is.null(fishery_byproducts) == FALSE){
    fmfo_stressor_ghg <- fishery_byproducts %>% filter(stressor == "ghg" & feed_type == "fmfo") %>% pull(ave_stressor)
    fmfo_stressor_N <- fishery_byproducts %>% filter(stressor == "N" & feed_type == "fmfo") %>% pull(ave_stressor)
    fmfo_stressor_P <- fishery_byproducts %>% filter(stressor == "P" & feed_type == "fmfo") %>% pull(ave_stressor)
    fmfo_stressor_land <- fishery_byproducts %>% filter(stressor == "land" & feed_type == "fmfo") %>% pull(ave_stressor)
    fmfo_stressor_water <- fishery_byproducts %>% filter(stressor == "water" & feed_type == "fmfo") %>% pull(ave_stressor)
  }else{}
  
  if(is.null(low_impact_fishery_byproducts) == FALSE){
    fmfo_stressor_ghg <- low_impact_fishery_byproducts %>% filter(stressor == "ghg" & feed_type == "fmfo") %>% pull(ave_stressor)
    fmfo_stressor_N <- low_impact_fishery_byproducts %>% filter(stressor == "N" & feed_type == "fmfo") %>% pull(ave_stressor)
    fmfo_stressor_P <- low_impact_fishery_byproducts %>% filter(stressor == "P" & feed_type == "fmfo") %>% pull(ave_stressor)
    fmfo_stressor_land <- low_impact_fishery_byproducts %>% filter(stressor == "land" & feed_type == "fmfo") %>% pull(ave_stressor)
    fmfo_stressor_water <- low_impact_fishery_byproducts %>% filter(stressor == "water" & feed_type == "fmfo") %>% pull(ave_stressor)
  }else{}
  
  if(is.null(fmfo_with_land_change_free_crops) == FALSE){
    fmfo_stressor_ghg <- fmfo_with_land_change_free_crops %>% filter(stressor == "ghg" & feed_type == "fmfo") %>% pull(ave_stressor)
    fmfo_stressor_N <- fmfo_with_land_change_free_crops %>% filter(stressor == "N" & feed_type == "fmfo") %>% pull(ave_stressor)
    fmfo_stressor_P <- fmfo_with_land_change_free_crops %>% filter(stressor == "P" & feed_type == "fmfo") %>% pull(ave_stressor)
    fmfo_stressor_land <- fmfo_with_land_change_free_crops %>% filter(stressor == "land" & feed_type == "fmfo") %>% pull(ave_stressor)
    fmfo_stressor_water <- fmfo_with_land_change_free_crops %>% filter(stressor == "water" & feed_type == "fmfo") %>% pull(ave_stressor)
  }else{}
  
  if(is.null(fmfo_low_impact_fishery_byproducts) == FALSE){
    fmfo_stressor_ghg <- fmfo_low_impact_fishery_byproducts %>% filter(stressor == "ghg" & feed_type == "fmfo") %>% pull(ave_stressor)
    fmfo_stressor_N <- fmfo_low_impact_fishery_byproducts %>% filter(stressor == "N" & feed_type == "fmfo") %>% pull(ave_stressor)
    fmfo_stressor_P <- fmfo_low_impact_fishery_byproducts %>% filter(stressor == "P" & feed_type == "fmfo") %>% pull(ave_stressor)
    fmfo_stressor_land <- fmfo_low_impact_fishery_byproducts %>% filter(stressor == "land" & feed_type == "fmfo") %>% pull(ave_stressor)
    fmfo_stressor_water <- fmfo_low_impact_fishery_byproducts %>% filter(stressor == "water" & feed_type == "fmfo") %>% pull(ave_stressor)
  }else{}
  
  # Set energy constants
  diesel_ghg <- data_energy %>% filter(Input == "Diesel_L") %>% pull(Value)
  petrol_ghg <- data_energy %>% filter(Input == "Petrol_L") %>% pull(Value)
  natgas_ghg <- data_energy %>% filter(Input == "NaturalGas_L") %>% pull(Value)
  
  if(zero_emissions_electric == TRUE){
    data_lca$electricity_ghg_weighted_mean <- 0
  }else{}
  
  # Set feed N and P constants
  N_content_soy <- data_feed_NP %>% filter(feed_type == "soy") %>% pull(N)
  P_content_soy <- data_feed_NP %>% filter(feed_type == "soy") %>% pull(P)
  
  N_content_crops <- data_feed_NP %>% filter(feed_type == "crops") %>% pull(N)
  P_content_crops <- data_feed_NP %>% filter(feed_type == "crops") %>% pull(P)
  
  N_content_animal <- data_feed_NP %>% filter(feed_type == "animal") %>% pull(N)
  P_content_animal <- data_feed_NP %>% filter(feed_type == "animal") %>% pull(P)
  
  N_content_fmfo <- data_feed_NP %>% filter(feed_type == "fmfo") %>% pull(N)
  P_content_fmfo <- data_feed_NP %>% filter(feed_type == "fmfo") %>% pull(P)
  
  if(is.null(fmfo_with_land_change_free_crops) == FALSE){
    N_content_fmfo <- data_feed_NP %>% filter(feed_type == "soy") %>% pull(N)
    P_content_fmfo <- data_feed_NP %>% filter(feed_type == "soy") %>% pull(P)
  }else{}
  
  data_lca <- data_lca %>%
    # Add deltas to calculate perturbation
    mutate(fcr_weighted_mean = fcr_weighted_mean + delta_FCR,
           feed_soy_weighted_mean = feed_soy_weighted_mean + delta_soy, 
           feed_crops_weighted_mean = feed_crops_weighted_mean + delta_crops,
           feed_animal_weighted_mean = feed_animal_weighted_mean + delta_animal,
           feed_fmfo_weighted_mean = feed_fmfo_weighted_mean + delta_fmfo,
           electricity_weighted_mean = electricity_weighted_mean + delta_electricity,
           diesel_weighted_mean = diesel_weighted_mean + delta_diesel, 
           petrol_weighted_mean = petrol_weighted_mean + delta_petrol, 
           naturalgas_weighted_mean = naturalgas_weighted_mean + delta_natgas, 
           yield_weighted_mean = yield_weighted_mean + delta_yield) %>% 
    # Reweight feed components to ensure they still sum to 1
    mutate(feed_soy_weighted_mean_new = feed_soy_weighted_mean/(feed_soy_weighted_mean + feed_crops_weighted_mean +
                                                                  feed_animal_weighted_mean + feed_fmfo_weighted_mean),
           feed_crops_weighted_mean_new = feed_crops_weighted_mean/(feed_soy_weighted_mean + feed_crops_weighted_mean +
                                                                      feed_animal_weighted_mean + feed_fmfo_weighted_mean),
           feed_animal_weighted_mean_new = feed_animal_weighted_mean/(feed_soy_weighted_mean + feed_crops_weighted_mean +
                                                                        feed_animal_weighted_mean + feed_fmfo_weighted_mean),
           feed_fmfo_weighted_mean_new = feed_fmfo_weighted_mean/(feed_soy_weighted_mean + feed_crops_weighted_mean +
                                                                    feed_animal_weighted_mean + feed_fmfo_weighted_mean)) %>%
    mutate(feed_soy_weighted_mean_new = ifelse(is.na(feed_soy_weighted_mean_new), 0, feed_soy_weighted_mean_new),
           feed_crops_weighted_mean_new = ifelse(is.na(feed_crops_weighted_mean_new), 0, feed_crops_weighted_mean_new),
           feed_animal_weighted_mean_new = ifelse(is.na(feed_animal_weighted_mean_new), 0, feed_animal_weighted_mean_new),
           feed_fmfo_weighted_mean_new = ifelse(is.na(feed_fmfo_weighted_mean_new), 0, feed_fmfo_weighted_mean_new)) %>%
    # Add grow out period constants
    mutate(grow_out_yr_prop = case_when(
      taxa == "oth_carp" ~ 300/365,
      taxa == "hypoph_carp" ~ 300/365,
      taxa == "catfish" ~ 210/365,
      taxa == "tilapia" ~ 200/365,
      taxa == "trouts" ~ 365/365,
      taxa == "fresh_crust" ~ 240/365
    )) %>%
    # Calculate all stressors
    mutate(
      # Calculate feed-associate components
      feed_ghg = 1000*fcr_weighted_mean*(feed_soy_weighted_mean_new*soy_stressor_ghg + 
                                           feed_crops_weighted_mean_new*crops_stressor_ghg + 
                                           feed_animal_weighted_mean_new*animal_stressor_ghg + 
                                           feed_fmfo_weighted_mean_new*fmfo_stressor_ghg),
      
      feed_N = 1000*fcr_weighted_mean*(feed_soy_weighted_mean_new*soy_stressor_N + 
                                         feed_crops_weighted_mean_new*crops_stressor_N + 
                                         feed_animal_weighted_mean_new*animal_stressor_N + 
                                         feed_fmfo_weighted_mean_new*fmfo_stressor_N),
      
      feed_P = 1000*fcr_weighted_mean*(feed_soy_weighted_mean_new*soy_stressor_P + 
                                         feed_crops_weighted_mean_new*crops_stressor_P + 
                                         feed_animal_weighted_mean_new*animal_stressor_P + 
                                         feed_fmfo_weighted_mean_new*fmfo_stressor_P),
      
      feed_land = 1000*fcr_weighted_mean*(feed_soy_weighted_mean_new*soy_stressor_land + 
                                            feed_crops_weighted_mean_new*crops_stressor_land + 
                                            feed_animal_weighted_mean_new*animal_stressor_land + 
                                            feed_fmfo_weighted_mean_new*fmfo_stressor_land),
      
      feed_water = 1000*fcr_weighted_mean*(feed_soy_weighted_mean_new*soy_stressor_water + 
                                             feed_crops_weighted_mean_new*crops_stressor_water + 
                                             feed_animal_weighted_mean_new*animal_stressor_water + 
                                             feed_fmfo_weighted_mean_new*fmfo_stressor_water), 
      
      # Calculate on farm GHG
      onfarm_ghg = electricity_weighted_mean*electricity_ghg_weighted_mean +
        diesel_weighted_mean*diesel_ghg + petrol_weighted_mean*petrol_ghg + naturalgas_weighted_mean*natgas_ghg,
      
      # Calculate on farm N and P
      onfarm_N = 1000*((fcr_weighted_mean*(feed_soy_weighted_mean_new*N_content_soy + 
                                             feed_crops_weighted_mean_new*N_content_crops + 
                                             feed_animal_weighted_mean_new*N_content_animal + 
                                             feed_fmfo_weighted_mean_new*N_content_fmfo)) - fish_N_weighted_mean),
      
      onfarm_P = 1000*((fcr_weighted_mean*(feed_soy_weighted_mean_new*P_content_soy + 
                                             feed_crops_weighted_mean_new*P_content_crops + 
                                             feed_animal_weighted_mean_new*P_content_animal + 
                                             feed_fmfo_weighted_mean_new*P_content_fmfo)) - fish_P_weighted_mean),
      
      # Calculate on farm land
      onfarm_land = yield_weighted_mean,
      
      # Calculate on farm water (for freshwater species only)
      onfarm_water =  ifelse(taxa %in% c("oth_carp", "catfish", "hypoph_carp", "tilapia", "trouts", "fresh_crust"), 
                             onfarm_land*evap_weighted_mean*grow_out_yr_prop, 0),
      
      # Calculate totals
      total_ghg = feed_ghg + onfarm_ghg,
      total_N = feed_N + onfarm_N,
      total_P = feed_P + onfarm_P,
      total_land = feed_land + onfarm_land,
      total_water = feed_water + onfarm_water
    ) %>%
    select(taxa, feed_ghg:total_water)
  
  return(data_lca)
}

#_______________________________________________________________________________________________________________________#
# Function to compare change in variable with change (relative to sd)
#_______________________________________________________________________________________________________________________#
compare_perturbations_sd <- function(n.sd){
  # Baseline data
  stressor_0 <- stressor_sensitivity(data_lca = df_taxa, data_feed_stressors = feed_fp, data_feed_NP = feed_NP, data_energy = energy_gwp,
                                     delta_FCR = 0, delta_soy = 0, delta_crops = 0, delta_animal = 0, delta_fmfo = 0,
                                     delta_electricity = 0, delta_diesel = 0, delta_petrol = 0, delta_natgas = 0, delta_yield = 0)
  stressor_0 <- stressor_0 %>%
    select(taxa, "baseline_total_ghg" = "total_ghg", "baseline_total_N" = "total_N", "baseline_total_P" = "total_P", 
           "baseline_total_land" = "total_land", "baseline_total_water" = "total_water")
  
  # FCR perturbation
  stressor_1 <- stressor_sensitivity(data_lca = df_taxa, data_feed_stressors = feed_fp, data_feed_NP = feed_NP, data_energy = energy_gwp,
                                     delta_FCR = n.sd*df_taxa$fcr_weighted_sd, delta_soy = 0, delta_crops = 0, delta_animal = 0, delta_fmfo = 0,
                                     delta_electricity = 0, delta_diesel = 0, delta_petrol = 0, delta_natgas = 0, delta_yield = 0)
  stressor_1 <- stressor_1 %>%
    select(taxa, "fcr_total_ghg" = "total_ghg", "fcr_total_N" = "total_N", "fcr_total_P" = "total_P", 
           "fcr_total_land" = "total_land", "fcr_total_water" = "total_water")
  
  stressor_0 <- stressor_0 %>%
    left_join(stressor_1, by = "taxa") %>%
    mutate(fcr_ghg_percent_change = 100*(fcr_total_ghg - baseline_total_ghg)/baseline_total_ghg,
           fcr_N_percent_change = 100*(fcr_total_N - baseline_total_N)/baseline_total_N,
           fcr_P_percent_change = 100*(fcr_total_P - baseline_total_P)/baseline_total_P,
           fcr_land_percent_change = 100*(fcr_total_land - baseline_total_land)/baseline_total_land,
           fcr_water_percent_change = 100*(fcr_total_water - baseline_total_water)/baseline_total_water
    ) 
  
  # Soy perturbation
  stressor_1 <- stressor_sensitivity(data_lca = df_taxa, data_feed_stressors = feed_fp, data_feed_NP = feed_NP, data_energy = energy_gwp,
                                     delta_FCR = 0, delta_soy = n.sd*df_taxa$feed_soy_weighted_sd, delta_crops = 0, delta_animal = 0, delta_fmfo = 0,
                                     delta_electricity = 0, delta_diesel = 0, delta_petrol = 0, delta_natgas = 0, delta_yield = 0)
  stressor_1 <- stressor_1 %>%
    select(taxa, "soy_total_ghg" = "total_ghg", "soy_total_N" = "total_N", "soy_total_P" = "total_P", 
           "soy_total_land" = "total_land", "soy_total_water" = "total_water")
  
  stressor_0 <- stressor_0 %>%
    left_join(stressor_1, by = "taxa") %>%
    mutate(soy_ghg_percent_change = 100*(soy_total_ghg - baseline_total_ghg)/baseline_total_ghg,
           soy_N_percent_change = 100*(soy_total_N - baseline_total_N)/baseline_total_N,
           soy_P_percent_change = 100*(soy_total_P - baseline_total_P)/baseline_total_P,
           soy_land_percent_change = 100*(soy_total_land - baseline_total_land)/baseline_total_land,
           soy_water_percent_change = 100*(soy_total_water - baseline_total_water)/baseline_total_water
    )
  
  # Crop perturbation
  stressor_1 <- stressor_sensitivity(data_lca = df_taxa, data_feed_stressors = feed_fp, data_feed_NP = feed_NP, data_energy = energy_gwp,
                                     delta_FCR = 0, delta_soy = 0, delta_crops = n.sd*df_taxa$feed_crops_weighted_sd, delta_animal = 0, delta_fmfo = 0,
                                     delta_electricity = 0, delta_diesel = 0, delta_petrol = 0, delta_natgas = 0, delta_yield = 0)
  stressor_1 <- stressor_1 %>%
    select(taxa, "crops_total_ghg" = "total_ghg", "crops_total_N" = "total_N", "crops_total_P" = "total_P", 
           "crops_total_land" = "total_land", "crops_total_water" = "total_water")
  
  stressor_0 <- stressor_0 %>%
    left_join(stressor_1, by = "taxa") %>%
    mutate(crops_ghg_percent_change = 100*(crops_total_ghg - baseline_total_ghg)/baseline_total_ghg,
           crops_N_percent_change = 100*(crops_total_N - baseline_total_N)/baseline_total_N,
           crops_P_percent_change = 100*(crops_total_P - baseline_total_P)/baseline_total_P,
           crops_land_percent_change = 100*(crops_total_land - baseline_total_land)/baseline_total_land,
           crops_water_percent_change = 100*(crops_total_water - baseline_total_water)/baseline_total_water
    )
  
  # Animal perturbation
  stressor_1 <- stressor_sensitivity(data_lca = df_taxa, data_feed_stressors = feed_fp, data_feed_NP = feed_NP, data_energy = energy_gwp,
                                     delta_FCR = 0, delta_soy = 0, delta_crops = 0, delta_animal = n.sd*df_taxa$feed_animal_weighted_sd, delta_fmfo = 0,
                                     delta_electricity = 0, delta_diesel = 0, delta_petrol = 0, delta_natgas = 0, delta_yield = 0)
  stressor_1 <- stressor_1 %>%
    select(taxa, "animal_total_ghg" = "total_ghg", "animal_total_N" = "total_N", "animal_total_P" = "total_P", 
           "animal_total_land" = "total_land", "animal_total_water" = "total_water")
  
  stressor_0 <- stressor_0 %>%
    left_join(stressor_1, by = "taxa") %>%
    mutate(animal_ghg_percent_change = 100*(animal_total_ghg - baseline_total_ghg)/baseline_total_ghg,
           animal_N_percent_change = 100*(animal_total_N - baseline_total_N)/baseline_total_N,
           animal_P_percent_change = 100*(animal_total_P - baseline_total_P)/baseline_total_P,
           animal_land_percent_change = 100*(animal_total_land - baseline_total_land)/baseline_total_land,
           animal_water_percent_change = 100*(animal_total_water - baseline_total_water)/baseline_total_water
    )
  
  # FMFO perturbation
  stressor_1 <- stressor_sensitivity(data_lca = df_taxa, data_feed_stressors = feed_fp, data_feed_NP = feed_NP, data_energy = energy_gwp,
                                     delta_FCR = 0, delta_soy = 0, delta_crops = 0, delta_animal = 0, delta_fmfo = n.sd*df_taxa$feed_fmfo_weighted_sd,
                                     delta_electricity = 0, delta_diesel = 0, delta_petrol = 0, delta_natgas = 0, delta_yield = 0)
  stressor_1 <- stressor_1 %>%
    select(taxa, "fmfo_total_ghg" = "total_ghg", "fmfo_total_N" = "total_N", "fmfo_total_P" = "total_P", 
           "fmfo_total_land" = "total_land", "fmfo_total_water" = "total_water")
  
  stressor_0 <- stressor_0 %>%
    left_join(stressor_1, by = "taxa") %>%
    mutate(fmfo_ghg_percent_change = 100*(fmfo_total_ghg - baseline_total_ghg)/baseline_total_ghg,
           fmfo_N_percent_change = 100*(fmfo_total_N - baseline_total_N)/baseline_total_N,
           fmfo_P_percent_change = 100*(fmfo_total_P - baseline_total_P)/baseline_total_P,
           fmfo_land_percent_change = 100*(fmfo_total_land - baseline_total_land)/baseline_total_land,
           fmfo_water_percent_change = 100*(fmfo_total_water - baseline_total_water)/baseline_total_water
    )
  
  # Energy perturbation
  stressor_1 <- stressor_sensitivity(data_lca = df_taxa, data_feed_stressors = feed_fp, data_feed_NP = feed_NP, data_energy = energy_gwp,
                                     delta_FCR = 0, delta_soy = 0, delta_crops = 0, delta_animal = 0, delta_fmfo = 0,
                                     delta_electricity = n.sd*df_taxa$electricity_weighted_sd, 
                                     delta_diesel = n.sd*df_taxa$diesel_weighted_sd, 
                                     delta_petrol = n.sd*df_taxa$petrol_weighted_sd, 
                                     delta_natgas = n.sd*df_taxa$naturalgas_weighted_sd, delta_yield = 0)
  stressor_1 <- stressor_1 %>%
    select(taxa, "energy_total_ghg" = "total_ghg", "energy_total_N" = "total_N", "energy_total_P" = "total_P", 
           "energy_total_land" = "total_land", "energy_total_water" = "total_water")
  
  stressor_0 <- stressor_0 %>%
    left_join(stressor_1, by = "taxa") %>%
    mutate(energy_ghg_percent_change = 100*(energy_total_ghg - baseline_total_ghg)/baseline_total_ghg,
           energy_N_percent_change = 100*(energy_total_N - baseline_total_N)/baseline_total_N,
           energy_P_percent_change = 100*(energy_total_P - baseline_total_P)/baseline_total_P,
           energy_land_percent_change = 100*(energy_total_land - baseline_total_land)/baseline_total_land,
           energy_water_percent_change = 100*(energy_total_water - baseline_total_water)/baseline_total_water
    )
  
  # Yield perturbation
  stressor_1 <- stressor_sensitivity(data_lca = df_taxa, data_feed_stressors = feed_fp, data_feed_NP = feed_NP, data_energy = energy_gwp,
                                     delta_FCR = 0, delta_soy = 0, delta_crops = 0, delta_animal = 0, delta_fmfo = 0,
                                     delta_electricity = 0, delta_diesel = 0, delta_petrol = 0, delta_natgas = 0, 
                                     delta_yield = n.sd*df_taxa$yield_weighted_sd)
  stressor_1 <- stressor_1 %>%
    select(taxa, "yield_total_ghg" = "total_ghg", "yield_total_N" = "total_N", "yield_total_P" = "total_P", 
           "yield_total_land" = "total_land", "yield_total_water" = "total_water")
  
  stressor_0 <- stressor_0 %>%
    left_join(stressor_1, by = "taxa") %>%
    mutate(yield_ghg_percent_change = 100*(yield_total_ghg - baseline_total_ghg)/baseline_total_ghg,
           yield_N_percent_change = 100*(yield_total_N - baseline_total_N)/baseline_total_N,
           yield_P_percent_change = 100*(yield_total_P - baseline_total_P)/baseline_total_P,
           yield_land_percent_change = 100*(yield_total_land - baseline_total_land)/baseline_total_land,
           yield_water_percent_change = 100*(yield_total_water - baseline_total_water)/baseline_total_water
    )
  
  # Only keep percent change columns
  stressor_0 <- stressor_0 %>%
    select(taxa, contains("percent_change"))
  
  return(stressor_0)
}



#_______________________________________________________________________________________________________________________#
# Function to compare change in variable with change (relative to mean)
#_______________________________________________________________________________________________________________________#
compare_perturbations_mean <- function(n){
  # Baseline data
  stressor_0 <- stressor_sensitivity(data_lca = df_taxa, data_feed_stressors = feed_fp, data_feed_NP = feed_NP, data_energy = energy_gwp,
                                     delta_FCR = 0, delta_soy = 0, delta_crops = 0, delta_animal = 0, delta_fmfo = 0,
                                     delta_electricity = 0, delta_diesel = 0, delta_petrol = 0, delta_natgas = 0, delta_yield = 0)
  stressor_0 <- stressor_0 %>%
    select(taxa, "baseline_total_ghg" = "total_ghg", "baseline_total_N" = "total_N", "baseline_total_P" = "total_P", 
           "baseline_total_land" = "total_land", "baseline_total_water" = "total_water")
  
  # FCR perturbation
  stressor_1 <- stressor_sensitivity(data_lca = df_taxa, data_feed_stressors = feed_fp, data_feed_NP = feed_NP, data_energy = energy_gwp,
                                     delta_FCR = n*df_taxa$fcr_weighted_mean, delta_soy = 0, delta_crops = 0, delta_animal = 0, delta_fmfo = 0,
                                     delta_electricity = 0, delta_diesel = 0, delta_petrol = 0, delta_natgas = 0, delta_yield = 0)
  stressor_1 <- stressor_1 %>%
    select(taxa, "fcr_total_ghg" = "total_ghg", "fcr_total_N" = "total_N", "fcr_total_P" = "total_P", 
           "fcr_total_land" = "total_land", "fcr_total_water" = "total_water")
  
  stressor_0 <- stressor_0 %>%
    left_join(stressor_1, by = "taxa") %>%
    mutate(fcr_ghg_percent_change = 100*(fcr_total_ghg - baseline_total_ghg)/baseline_total_ghg,
           fcr_N_percent_change = 100*(fcr_total_N - baseline_total_N)/baseline_total_N,
           fcr_P_percent_change = 100*(fcr_total_P - baseline_total_P)/baseline_total_P,
           fcr_land_percent_change = 100*(fcr_total_land - baseline_total_land)/baseline_total_land,
           fcr_water_percent_change = 100*(fcr_total_water - baseline_total_water)/baseline_total_water
    ) 
  
  # Soy perturbation
  stressor_1 <- stressor_sensitivity(data_lca = df_taxa, data_feed_stressors = feed_fp, data_feed_NP = feed_NP, data_energy = energy_gwp,
                                     delta_FCR = 0, delta_soy = n*df_taxa$feed_soy_weighted_mean, 
                                     delta_crops = (-1*n*df_taxa$feed_soy_weighted_mean)/3, # change other feeds 
                                     delta_animal = (-1*n*df_taxa$feed_soy_weighted_mean)/3, 
                                     delta_fmfo = (-1*n*df_taxa$feed_soy_weighted_mean)/3,
                                     delta_electricity = 0, delta_diesel = 0, delta_petrol = 0, delta_natgas = 0, delta_yield = 0)
  stressor_1 <- stressor_1 %>%
    select(taxa, "soy_total_ghg" = "total_ghg", "soy_total_N" = "total_N", "soy_total_P" = "total_P", 
           "soy_total_land" = "total_land", "soy_total_water" = "total_water")
  
  stressor_0 <- stressor_0 %>%
    left_join(stressor_1, by = "taxa") %>%
    mutate(soy_ghg_percent_change = 100*(soy_total_ghg - baseline_total_ghg)/baseline_total_ghg,
           soy_N_percent_change = 100*(soy_total_N - baseline_total_N)/baseline_total_N,
           soy_P_percent_change = 100*(soy_total_P - baseline_total_P)/baseline_total_P,
           soy_land_percent_change = 100*(soy_total_land - baseline_total_land)/baseline_total_land,
           soy_water_percent_change = 100*(soy_total_water - baseline_total_water)/baseline_total_water
    )
  
  # Crop perturbation
  stressor_1 <- stressor_sensitivity(data_lca = df_taxa, data_feed_stressors = feed_fp, data_feed_NP = feed_NP, data_energy = energy_gwp,
                                     delta_FCR = 0, delta_soy = (-1*n*df_taxa$feed_crops_weighted_mean)/3, 
                                     delta_crops = n*df_taxa$feed_crops_weighted_mean, 
                                     delta_animal = (-1*n*df_taxa$feed_crops_weighted_mean)/3, 
                                     delta_fmfo = (-1*n*df_taxa$feed_crops_weighted_mean)/3,
                                     delta_electricity = 0, delta_diesel = 0, delta_petrol = 0, delta_natgas = 0, delta_yield = 0)
  stressor_1 <- stressor_1 %>%
    select(taxa, "crops_total_ghg" = "total_ghg", "crops_total_N" = "total_N", "crops_total_P" = "total_P", 
           "crops_total_land" = "total_land", "crops_total_water" = "total_water")
  
  stressor_0 <- stressor_0 %>%
    left_join(stressor_1, by = "taxa") %>%
    mutate(crops_ghg_percent_change = 100*(crops_total_ghg - baseline_total_ghg)/baseline_total_ghg,
           crops_N_percent_change = 100*(crops_total_N - baseline_total_N)/baseline_total_N,
           crops_P_percent_change = 100*(crops_total_P - baseline_total_P)/baseline_total_P,
           crops_land_percent_change = 100*(crops_total_land - baseline_total_land)/baseline_total_land,
           crops_water_percent_change = 100*(crops_total_water - baseline_total_water)/baseline_total_water
    )
  
  # Animal perturbation
  stressor_1 <- stressor_sensitivity(data_lca = df_taxa, data_feed_stressors = feed_fp, data_feed_NP = feed_NP, data_energy = energy_gwp,
                                     delta_FCR = 0, delta_soy = (-1*n*df_taxa$feed_animal_weighted_mean)/3, 
                                     delta_crops = (-1*n*df_taxa$feed_animal_weighted_mean)/3,
                                     delta_animal = n*df_taxa$feed_animal_weighted_mean, 
                                     delta_fmfo = (-1*n*df_taxa$feed_animal_weighted_mean)/3,
                                     delta_electricity = 0, delta_diesel = 0, delta_petrol = 0, delta_natgas = 0, delta_yield = 0)
  stressor_1 <- stressor_1 %>%
    select(taxa, "animal_total_ghg" = "total_ghg", "animal_total_N" = "total_N", "animal_total_P" = "total_P", 
           "animal_total_land" = "total_land", "animal_total_water" = "total_water")
  
  stressor_0 <- stressor_0 %>%
    left_join(stressor_1, by = "taxa") %>%
    mutate(animal_ghg_percent_change = 100*(animal_total_ghg - baseline_total_ghg)/baseline_total_ghg,
           animal_N_percent_change = 100*(animal_total_N - baseline_total_N)/baseline_total_N,
           animal_P_percent_change = 100*(animal_total_P - baseline_total_P)/baseline_total_P,
           animal_land_percent_change = 100*(animal_total_land - baseline_total_land)/baseline_total_land,
           animal_water_percent_change = 100*(animal_total_water - baseline_total_water)/baseline_total_water
    )
  
  # FMFO perturbation
  stressor_1 <- stressor_sensitivity(data_lca = df_taxa, data_feed_stressors = feed_fp, data_feed_NP = feed_NP, data_energy = energy_gwp,
                                     delta_FCR = 0, delta_soy = (-1*n*df_taxa$feed_fmfo_weighted_mean)/3, 
                                     delta_crops = (-1*n*df_taxa$feed_fmfo_weighted_mean)/3, 
                                     delta_animal = (-1*n*df_taxa$feed_fmfo_weighted_mean)/3, 
                                     delta_fmfo = n*df_taxa$feed_fmfo_weighted_mean,
                                     delta_electricity = 0, delta_diesel = 0, delta_petrol = 0, delta_natgas = 0, delta_yield = 0)
  stressor_1 <- stressor_1 %>%
    select(taxa, "fmfo_total_ghg" = "total_ghg", "fmfo_total_N" = "total_N", "fmfo_total_P" = "total_P", 
           "fmfo_total_land" = "total_land", "fmfo_total_water" = "total_water")
  
  stressor_0 <- stressor_0 %>%
    left_join(stressor_1, by = "taxa") %>%
    mutate(fmfo_ghg_percent_change = 100*(fmfo_total_ghg - baseline_total_ghg)/baseline_total_ghg,
           fmfo_N_percent_change = 100*(fmfo_total_N - baseline_total_N)/baseline_total_N,
           fmfo_P_percent_change = 100*(fmfo_total_P - baseline_total_P)/baseline_total_P,
           fmfo_land_percent_change = 100*(fmfo_total_land - baseline_total_land)/baseline_total_land,
           fmfo_water_percent_change = 100*(fmfo_total_water - baseline_total_water)/baseline_total_water
    )
  
  # Energy perturbation
  stressor_1 <- stressor_sensitivity(data_lca = df_taxa, data_feed_stressors = feed_fp, data_feed_NP = feed_NP, data_energy = energy_gwp,
                                     delta_FCR = 0, delta_soy = 0, delta_crops = 0, delta_animal = 0, delta_fmfo = 0,
                                     delta_electricity = n*df_taxa$electricity_weighted_mean, 
                                     delta_diesel = n*df_taxa$diesel_weighted_mean, 
                                     delta_petrol = n*df_taxa$petrol_weighted_mean, 
                                     delta_natgas = n*df_taxa$naturalgas_weighted_mean, delta_yield = 0)
  stressor_1 <- stressor_1 %>%
    select(taxa, "energy_total_ghg" = "total_ghg", "energy_total_N" = "total_N", "energy_total_P" = "total_P", 
           "energy_total_land" = "total_land", "energy_total_water" = "total_water")
  
  stressor_0 <- stressor_0 %>%
    left_join(stressor_1, by = "taxa") %>%
    mutate(energy_ghg_percent_change = 100*(energy_total_ghg - baseline_total_ghg)/baseline_total_ghg,
           energy_N_percent_change = 100*(energy_total_N - baseline_total_N)/baseline_total_N,
           energy_P_percent_change = 100*(energy_total_P - baseline_total_P)/baseline_total_P,
           energy_land_percent_change = 100*(energy_total_land - baseline_total_land)/baseline_total_land,
           energy_water_percent_change = 100*(energy_total_water - baseline_total_water)/baseline_total_water
    )
  
  # Yield perturbation
  stressor_1 <- stressor_sensitivity(data_lca = df_taxa, data_feed_stressors = feed_fp, data_feed_NP = feed_NP, data_energy = energy_gwp,
                                     delta_FCR = 0, delta_soy = 0, delta_crops = 0, delta_animal = 0, delta_fmfo = 0,
                                     delta_electricity = 0, delta_diesel = 0, delta_petrol = 0, delta_natgas = 0, 
                                     delta_yield = n*df_taxa$yield_weighted_mean)
  stressor_1 <- stressor_1 %>%
    select(taxa, "yield_total_ghg" = "total_ghg", "yield_total_N" = "total_N", "yield_total_P" = "total_P", 
           "yield_total_land" = "total_land", "yield_total_water" = "total_water")
  
  stressor_0 <- stressor_0 %>%
    left_join(stressor_1, by = "taxa") %>%
    mutate(yield_ghg_percent_change = 100*(yield_total_ghg - baseline_total_ghg)/baseline_total_ghg,
           yield_N_percent_change = 100*(yield_total_N - baseline_total_N)/baseline_total_N,
           yield_P_percent_change = 100*(yield_total_P - baseline_total_P)/baseline_total_P,
           yield_land_percent_change = 100*(yield_total_land - baseline_total_land)/baseline_total_land,
           yield_water_percent_change = 100*(yield_total_water - baseline_total_water)/baseline_total_water
    )
  
  # Only keep percent change columns
  stressor_0 <- stressor_0 %>%
    select(taxa, contains("percent_change"))
}


