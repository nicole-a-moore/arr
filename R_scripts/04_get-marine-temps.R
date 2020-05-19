## extracting marine temperature data from NOAA OI SST V2 High Resolution Dataset 
## SST, Daily Optimum Interpolation (OI), AVHRR Only, Version 2, Final+Preliminary, 1981-present
library(tidyverse)
library(rerddap)
library(ncdf4)

## bring in population data
cadillac <- read.csv("./data-processed/intratherm-may-2020-squeaky-clean.csv")

## filter out rows of data we cannot use and that are not marine
marine <- cadillac %>%
	subset(subset = !is.na(latitude)) %>%
	subset(subset = !is.na(longitude)) %>%
	filter(realm_general2 == "Marine")

marine <- droplevels(marine)

## get rid of populations with same latitude, longitude and elevation since all will have the same temp data
unique_pairs <- marine[!duplicated(marine[,c("latitude", "longitude", "elevation_of_collection")]),]


## load info about NOAA data:
info <- info("ncdcOisst2Agg_LonPM180")


## make latitude and longitude vectors based on NOAA format
## longitude: Uniform grid with centers from -179.875 to 179.875 by 0.25 degrees.
lon <- rep(-179.875, times = 1439)
n = 2
while (n < 1441) {
	lon[n] <- lon[n -1] + 0.25
	n = n+1
}
## latitude: Uniform grid with centers from -89.875 to 89.875 by 0.25 degrees.
lat <- rep(-89.875, times = 719)
n = 2
while (n < 721) {
	lat[n] <- lat[n -1] + 0.25
	n = n+1
}


unique_pairs$grid_lat <- c()
unique_pairs$grid_lon <- c()

## find closest lat lon grid cell to each population collection location 
num_unique <- 1
while (num_unique < length(unique_pairs$population_id) + 1) {
	loc_lon_index <- which.min(abs(lon - unique_pairs$longitude[num_unique]))
	loc_lat_index <- which.min(abs(lat - unique_pairs$latitude[num_unique]))
	
	unique_pairs$grid_lon[num_unique] <- lon[loc_lon_index]
	unique_pairs$grid_lat[num_unique] <- lat[loc_lat_index]
	
	num_unique = num_unique + 1
}


## create dataframe for temp data
## nValues for time attribute = 14096
temperature_data <- data.frame(matrix(nrow = 14096))
colnames(temperature_data) = c("date")

## loop through each population getting temp data for its grid cell and adding to temp data
num_unique <- 1
while (num_unique < length(unique_pairs$population_id) + 1) {
	print(paste("On population number", num_unique))
	time_series <- griddap(info,
						   time = c("1981-09-01", "2020-04-04"), 
						   latitude = c(unique_pairs$grid_lat[num_unique],unique_pairs$grid_lat[num_unique]),
						   longitude = c(unique_pairs$grid_lon[num_unique], unique_pairs$grid_lon[num_unique]),
						   url = "https://upwell.pfeg.noaa.gov/erddap/")
	temps <- time_series$data$sst
	print(paste("Successfully got time series of length", length(temps)))
	
	if (num_unique == 1) {
		times <- time_series$data$time
		temperature_data$date <- times
	}
	
	temperature_data$temp <- temps
	pop_id <- paste(unique_pairs$population_id[num_unique], unique_pairs$longitude[num_unique], sep = "_")
	colnames(temperature_data)[num_unique+1]<- pop_id
	
	print("Stored data in temperature_data and moving on to next population!")
	
	num_unique <- num_unique + 1
}

precious_marine_temps <- temperature_data
saveRDS(precious_marine_temps, "~/Documents/SUNDAY LAB/Intratherm/Data sheets/precious_marine_temps.rds")

temperature_data <- readRDS("~/Documents/SUNDAY LAB/Intratherm/Data sheets/precious_marine_temps.rds")

