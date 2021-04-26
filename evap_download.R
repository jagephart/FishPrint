# Kelvin Gorospe kdgorospe@gmail.com
# Download NOAA global evaporation dataset 

rm(list=ls())
library(RCurl) # for getURL 

# GET FILES from NOAA FTP:
# https://ftp.cpc.ncep.noaa.gov/wd51yf/global_monthly

# Historical Data:
# REFERENCE paper for data: https://www.cpc.ncep.noaa.gov/products/Soilmst_Monitoring/Papers/2003JD004345.pdf
datadir_evap <- "/Volumes/jgephart/BFA Environment 2/Data/Evaporation"
url <- "ftp://ftp.cpc.ncep.noaa.gov/wd51yf/global_monthly/GeoTIFF/monthly/"
# Website: https://ftp.cpc.ncep.noaa.gov/wd51yf/global_monthly/GeoTIFF/monthly/
filenames <- getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
filenames <- strsplit(filenames, "\n")
filenames <- unlist(filenames)

for (filename in filenames) {
  download.file(url = paste(url, filename, sep = ""), destfile = file.path(datadir_evap, "Historical", filename))
}


# Cliamtology Data:
url <- "ftp://ftp.cpc.ncep.noaa.gov/wd51yf/global_monthly/GeoTIFF/clim/"
# Website: https://ftp.cpc.ncep.noaa.gov/wd51yf/global_monthly/GeoTIFF/clim/
# The data files in GeoTIFF/clim folder are monthly climatology (mean values for period 1981-2010 from data in GeoTiff/monthly folder). 
# The web page for GeoTiff/monthly data is here: https://www.cpc.ncep.noaa.gov/soilmst/leaky_glb.htm  (forced with a better version of CPC precipitation). 
filenames <- getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
filenames <- strsplit(filenames, "\n")
filenames <- unlist(filenames)

for (filename in filenames) {
  download.file(url = paste(url, filename, sep = ""), destfile = file.path(datadir_evap, "Monthly Climatology 1981-2010", filename))
}
