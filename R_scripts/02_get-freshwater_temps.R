## getting modelled monthly freshwater temperature time series' from 1981-2014 for each freshwater population 
## uses data from https://zenodo.org/record/3337659#.Xr3YcxNKjOR

library(ncdf4)
library(tidyverse)

## open nc file and get lat lon and time vectors
filename <- paste("/Volumes/TimeMachine/freshwater data/waterTemperature_monthly_1981-2014.nc", sep = "")
ncfile <- nc_open(filename)


lon <- ncvar_get(ncfile, "lon") ## units: degrees - intervals of 0.0833 (probably 0.5/60)
lat <- ncvar_get(ncfile, "lat") ## units: degrees - intervals of 0.08333 (probably 0.5/60)
time <- ncvar_get(ncfile, "time") ## units: days since 1901-01-01

## close the file
nc_close(ncfile)

## bring in population data
cadillac <- read.csv("./data-processed/intratherm-may-2020-squeaky-clean.csv")

## filter out rows of data we cannot use and that are not freshwater
freshwater <- cadillac %>%
	subset(subset = !is.na(latitude)) %>%
	subset(subset = !is.na(longitude)) %>%
	filter(realm_general2 == "Freshwater")

freshwater <- droplevels(freshwater)

## get rid of populations with same latitude, longitude and elevation since all will have the same temp data
unique_pairs <- freshwater[!duplicated(freshwater[,c("latitude", "longitude", "elevation_of_collection")]),]

## temps:
temperature_data <- data.frame(matrix(nrow = length(time)))
colnames(temperature_data) = c("date")
temperature_data$date <- time

## for each population:
num_unique <- 1
while (num_unique < length(unique_pairs$intratherm_id) + 1) {
	
	## find closest lat lon coordinates to population collection location
	loc_lon_index <- which.min(abs(lon - unique_pairs$longitude[num_unique]))
	loc_lat_index <- which.min(abs(lat - unique_pairs$latitude[num_unique]))
	
	## get waterTemp time series for closest lat lon coordinates 
	ncfile <- nc_open(filename)
	waterTemp <- ncvar_get(ncfile, "waterTemp", start = c(loc_lon_index, loc_lat_index, 1), count = c(1, 1, -1))
	nc_close(ncfile)
	
	## add to column in temperature_data and rename after column's population_id with longitude added onto the end
	temperature_data$temp <- waterTemp
	pop_id <- paste(unique_pairs$population_id[num_unique], unique_pairs$longitude[num_unique], sep = "_")
	colnames(temperature_data)[num_unique+1]<- pop_id

	num_unique <- num_unique + 1
}

## save to RDS
precious_temps <- temperature_data

temperature_data <- readRDS("~/Documents/SUNDAY LAB/Intratherm/Data sheets/precious_temps.rds")

## convert from degrees K to degrees C
converted <- temperature_data
converted[, 2:192] <- converted[, 2:192] - 273.15

## convert date to years and month fractions 
## starts at Jan 1981

months_num <- c((1 - 0.5)/12, (2 - 0.5)/12, (3 - 0.5)/12, (4 - 0.5)/12, (5 - 0.5)/12, (6 - 0.5)/12, (7 - 0.5)/12, (8 - 0.5)/12, (9 - 0.5)/12, (10 - 0.5)/12, (11 - 0.5)/12, (12-0.5)/12)

months_name <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

date_num <- c()
date_name <- c()
year <- 1981

x <- 1
while (x < 408/12 + 1) {
		
	date_num <- append(date_num, year + months_num)
	date_name <- append(date_name, paste(year, months_name))

	year = year + 1
	x = x + 1
}


## add columns
converted$date <- as.vector(date_name)
converted$date_number <- as.vector(date_num)

converted <- converted %>%
	select(date, date_number, everything())


## write temperature file to processed data 
write.csv(converted, "./data-processed/intratherm-freshwater-temp-data.csv")
