# Kelvin Gorospe kdgorospe@gmail.com
# Download NOAA global evaporation dataset 

rm(list=ls())
library(RCurl) # for getURL 

# GET FILES from NOAA FTP:
# FTP: ftp://ftp.cpc.ncep.noaa.gov/wd51yf/global_monthly

# Historical Data:
datadir_evap <- "/Volumes/jgephart/BFA Environment 2/Data/Evaporation"
url <- "ftp://ftp.cpc.ncep.noaa.gov/wd51yf/global_monthly/GeoTIFF/monthly/"
filenames <- getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
filenames <- strsplit(filenames, "\n")
filenames <- unlist(filenames)

for (filename in filenames) {
  download.file(url = paste(url, filename, sep = ""), destfile = file.path(datadir_evap, "Historical", filename))
}


# Cliamtology Data:
url <- "ftp://ftp.cpc.ncep.noaa.gov/wd51yf/global_monthly/GeoTIFF/clim/"
filenames <- getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
filenames <- strsplit(filenames, "\n")
filenames <- unlist(filenames)

for (filename in filenames) {
  download.file(url = paste(url, filename, sep = ""), destfile = file.path(datadir_evap, "Monthly Climatology 1981-2010", filename))
}
