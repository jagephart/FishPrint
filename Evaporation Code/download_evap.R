# Kelvin Gorospe kdgorospe@gmail.com
# Download NOAA global evaporation dataset 

rm(list=ls())
library(RCurl) # for getURL 

# GET FILES from NOAA FTP:
# FTP: ftp://ftp.cpc.ncep.noaa.gov/wd51yf/global_monthly

# Historical Data:
datadir <- "Data/Historical Monthly"
url <- "ftp://ftp.cpc.ncep.noaa.gov/wd51yf/global_monthly/GeoTIFF/monthly/"
filenames <- getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
filenames <- strsplit(filenames, "\n")
filenames <- unlist(filenames)

for (filename in filenames) {
  download.file(url = paste(url, filename, sep = ""), destfile = file.path(datadir, filename))
}


# Cliamtology Data:
datadir <- "Data/Monthly Climatology 1981-2010"
url <- "ftp://ftp.cpc.ncep.noaa.gov/wd51yf/global_monthly/GeoTIFF/clim/"
filenames <- getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
filenames <- strsplit(filenames, "\n")
filenames <- unlist(filenames)

for (filename in filenames) {
  download.file(url = paste(url, filename, sep = ""), destfile = file.path(datadir, filename))
}
