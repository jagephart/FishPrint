# Combine off-farm (feed) and on-farm (non-feed associated) WATER impacts 

datadir <- "/Volumes/jgephart/BFA Environment 2/Data"
outdir <- "/Volumes/jgephart/BFA Environment 2/Outputs"

# Load off-farm (feed-associated) variables
load(file.path(outdir, "2020-12-20_off-farm-all-impacts-all-data-prior-to-aggregation.RData"))

# Organize off-farm feed data (FCR and feed proportions)
# FIX IT - change clean.lca in Functions.R so no less than 0.01
# Adjust feed proportions to be no less than 0.01
# TEMPORARY FIX:
# Normalize the FINAL feed proportion values to be greater than 0 and no less than 0.01
full_feed_dat <- full_feed_dat %>%
  mutate(feed_proportion = if_else(feed_proportion < 0.01, true = 0.0105, false = feed_proportion))

# Format feed_dat
# Standardize number of significant digits (round to digist = 3)
# Re-scale so rowsums equals 1
feed_dat_merge <- full_feed_dat %>%
  select(study_id, clean_sci_name, taxa, intensity, system, feed_data_type = data_type, .category, contains("feed")) %>%
  pivot_wider(names_from = .category, values_from = feed_proportion) %>%
  mutate(feed_soy = round(feed_soy, digits = 4),
         feed_crops = round(feed_crops, digits = 4),
         feed_fmfo = round(feed_fmfo, digits = 4),
         feed_animal = round(feed_animal, digits = 4)) %>%
  mutate(new_scale = feed_soy + feed_crops + feed_fmfo + feed_animal) %>%
  mutate(feed_soy = feed_soy / new_scale,
         feed_crops = feed_crops / new_scale,
         feed_fmfo = feed_fmfo / new_scale,
         feed_animal = feed_animal / new_scale)

# Format FCR dat
# Rename data_type to keep track of both feed vs fcr data types
fcr_dat_merge <- full_fcr_dat %>%
  select(study_id, clean_sci_name, taxa, intensity, system, fcr_data_type = data_type, fcr)

# Load on-farm evaporative water loss (NOAA data - country-level mean of monthly climatological means 1981-2010)
evap_clim <- read.csv(file.path(datadir, "20201222_clim_summarise_by_country.csv")) %>%
  mutate(iso3c = countrycode(admin, origin = "country.name", destination = "iso3c")) %>%
  select(-X) %>%
  drop_na()

# Merge evap dat with Countries in lca data, only need mean_evap_mm and ID columns: study_id, Country, iso3c
# Then merge with feed data to get full model dataset
evap_dat <- lca_dat_clean_groups %>%
  left_join(evap_clim, by = "iso3c") %>%
  # Deal with some countries manually: Calculate mean for Germany, Denmark; Assign global mean to N/A
  mutate(mean_evap_mm = case_when(Country == "Germany, Denmark" ~ evap_clim %>% filter(iso3c %in% c("DEU", "DNK")) %>% pull(mean_evap_mm) %>% mean(),
                                  Country == "Scotland" ~ evap_clim %>% filter(iso3c == "GBR") %>% pull(mean_evap_mm),
                                  Country == "N/A" ~ mean(evap_clim$mean_evap_mm),
                                  TRUE ~ mean_evap_mm)) %>%
  select(study_id, Country, iso3c, clean_sci_name, taxa, intensity = Intensity, system = Production_system_group, mean_evap_mm) %>%
  # Merge with feed data for full model dataset
  full_join(feed_dat_merge, by = intersect(names(feed_dat_merge), names(.))) %>%
  full_join(fcr_dat_merge, by = intersect(names(feed_dat_merge), names(fcr_dat_merge)))
  
# LEFT OFF HERE: 
# Notice there are NAs in evap_dat, sill need to place conditionals in STAN (example for bivalves, only do on-farm water and skip feed-associated water)

# Check which studies have no evap data:
#evap_dat %>% filter(is.na(iso3c)) %>% select(Country, iso3c, mean_evap_mm)
