# Kelvin Gorospe
# Process NOAA global evaporation dataset 

rm(list=ls())
library(tidyverse)
library(stars) # newer package for working with raster data; better integration with sf package
library(sf)
library(raster) # needed for rotate to reset climate data from 0 to 360 to -180 to 180
library(rnaturalearth)
# library(rgeos) # required for rnaturalearth
# OLD spatial packages:
#library(raster) # for raster()
#library(rgdal) # for GDALinfo()

####################################################### READ IN monthly climatology files
datadir <- "/Volumes/jgephart/BFA Environment 2/Data"
datadir_evap <- "/Volumes/jgephart/BFA Environment 2/Data/Evaporation"
outdir_evap <- "/Volumes/jgephart/BFA Environment 2/Outputs/Evaporation"

# Use rotate() to transform coordinate system to have standard x coordinates: xmin = -180 to xmax = 180
# Convert longitude from 0 - 360 to -180 to 180 (Standard transformation for climate data)
# https://stackoverflow.com/questions/25730625/how-to-convert-longitude-from-0-360-to-180-180

clim_files <- list.files(file.path(datadir_evap, "Monthly Climatology 1981-2010"))
clim_raster <- lapply(clim_files, function(i){raster(file.path(datadir_evap, "Monthly Climatology 1981-2010", i))}) # read in with raster() function so that rotate() function works
clim_rotate <- lapply(clim_raster, function(i){rotate(i)})
clim_stars_list <- lapply(clim_rotate, function(i){st_as_stars(i)})

# Generate expression for c() based on length(clim_stars) and evaluate with eval()
stars_objects <- paste("clim_stars_list[[", 1:length(clim_stars_list), "]]", sep = "", collapse = ", ")
clim_stars <- eval(parse(text = paste("c(", stars_objects, ", along = list(month = clim_files))", sep = "")))

# Plot with ggplot:
ggplot() +
  geom_stars(data = clim_stars) +
  facet_wrap(~month)

ggsave(file.path(outdir_evap, "clim_panels.png"), width = 6, height = 4, unit = "in")

# Confirm that month labels in clim_dat are the same as clim_dat_list:
clim_dat_list <- read_stars(file.path(datadir_evap, "Monthly Climatology 1981-2010", clim_files))
for (i in 1:length(clim_dat_list)){
  ggplot() +
    geom_stars(data = clim_dat_list[i])
  
  ggsave(file.path(outdir_evap, paste("clim_list_", i, ".png", sep = "")), width = 6, height = 4, unit = "in")
}

####################################################### Process raster data
# <calculate mean across climatology months>
# Compute mean for each pixel (x,y) across all months

# use sum(is.na(clim_dat_list[[i]])) to check for NA's; no NAs found but keep na.rm = TRUE anyway
clim_mean_stars <- st_apply(clim_stars, c("x", "y"), mean, na.rm = TRUE) 

# Plot with ggplot:
ggplot() +
  geom_stars(data = clim_mean_stars)
ggsave(file.path(outdir_evap, "clim_mean.png"), width = 6, height = 4, unit = "in")

####################################################### Add country polygons
# Get countries from naturalearth package and set to same crs as clim_mean_stars
# datum = WGS84 i.e., standard lat/long and GPS coordinate system
# projection (instructions on how to transform datum onto a flat surface)
world_sf<- ne_countries(returnclass='sf') %>%
  dplyr::select(admin) %>%
  #st_transform("+proj=moll")
  st_transform(st_crs(clim_mean_stars))

st_crs(world_sf) == st_crs(clim_mean_stars)

ggplot() +
  geom_stars(data = clim_mean_stars) +
  #geom_stars(data = clim_mean_stars, mapping = aes(x = x, y = y, fill = mean)) +
  geom_sf(data = world_sf, alpha = 0.1, color = "black")

ggsave(file.path(outdir_evap, "clim_mean_world_map.png"), width = 6, height = 4, unit = "in")

# CROP raster and replot
clim_mean_stars <- st_crop(clim_mean_stars, world_sf)
ggplot() +
  geom_stars(data = clim_mean_stars) +
  geom_sf(data = world_sf, alpha = 0.1, color = "black")
ggsave(file.path(outdir_evap, "clim_mean_world_map_2.png"), width = 6, height = 4, unit = "in")

## although coordinates are longitude/latitude, st_intersects assumes that they are planar
# See: https://www.r-spatial.org/r/2020/06/17/s2.html (The Earth is no longer flat in r-spatial)

####################################################### Summarize data at the country-level
clim_mean_sf <- st_as_sf(clim_mean_stars)
clim_world_sf <- st_intersection(clim_mean_sf, world_sf)

clim_world_sf <- set_units(clim_world_sf$mean, mm)

# FIX IT - do we want the mean of all pixels or the sum of all pixels weighted by pixel area? Note: pixel's whose centroid are not within the polygon are assigned NA
# Calculate mean within each country
clim_by_country <- clim_world_sf %>%
  group_by(admin) %>%
  summarise(mean_evap_mm = mean(mean, na.rm = TRUE)) %>% # units are in "mm"/ month
  ungroup() 

ggplot(clim_by_country) +
  geom_sf(mapping = aes(fill = mean_evap_mm)) +
  labs(fill = "Country-level evaporation")
ggsave(file.path(outdir_evap, "clim_summarise_by_country_mean_of_pixels.png"), width = 6, height = 4, unit = "in")
 
write.csv(clim_by_country %>% st_set_geometry(NULL), file = file.path(datadir, "clim_summarise_by_country.csv"), quote = FALSE)

# BAR GRAPH:
#ggplot(clim_by_country) +
#  geom_col(mapping = aes(x = reorder(admin, desc(country_level)), y = country_level))

#country_level_clim <- aggregate(clim_world_sf, by = admin, FUN = mean)
