## getting daily data that Rens gave me:
########################################
## this data has a 30' spatial resolution 
## he says he will try and get me 5' resolution, but in the meantime we can use this 
library(ncdf4)
library(tidyverse)

## open nc file and get lat lon and time vectors
filename <- paste("watertemperature_wfd_historical_1958-2001.nc", sep = "")
ncfile <- nc_open(filename)

lon <- ncvar_get(ncfile, "longitude") ## units: degrees - intervals of 0.5 (30')
lat <- ncvar_get(ncfile, "latitude") ## units: degrees - intervals of 0.5 (30')
time <- ncvar_get(ncfile, "time") ## units: hours since 1901-01-01 (first time is 1958-01-01) 

## close the file
nc_close(ncfile)


## bring in population data
cadillac <- read.csv("./data-raw/intratherm-may-2020-squeaky-clean.csv")

## filter out rows of data we cannot use  
cadillac <- cadillac %>%
  filter(!is.na(latitude)) %>%
  filter(!is.na(longitude)) %>%
  filter(!is.na(acclim_temp)) %>%
  filter(!is.na(lifespan_days)) %>%
  filter(lifespan_days != "unk" || "kunk") %>%
  filter(parameter_tmax_or_tmin == "tmax") 

cadillac <- droplevels(cadillac)


## filter out rows of data that are not freshwater
freshwater <- cadillac %>%
  filter(realm_general2 == "Freshwater")

freshwater <- droplevels(freshwater)

## get rid of duplicate population rows since all will have the same temp data
unique_pairs <- freshwater[!duplicated(freshwater[,c("latitude", "longitude")]),]

## temps:
temperature_data <- data.frame(matrix(nrow = length(time)))
colnames(temperature_data) = c("date")
temperature_data$date <- time

## for each population:
num_unique <- 1
while (num_unique < nrow(unique_pairs) + 1) {
  
  ## find closest lat lon coordinates to population collection location
  loc_lon_index <- which.min(abs(lon - unique_pairs$longitude[num_unique]))
  loc_lat_index <- which.min(abs(lat - unique_pairs$latitude[num_unique]))
  
  ## get waterTemp time series for closest lat lon coordinates 
  ncfile <- nc_open(filename)
  waterTemp <- ncvar_get(ncfile, "waterTemperature", start = c(loc_lon_index, loc_lat_index, 1),
                         count = c(1, 1, -1))
  nc_close(ncfile)
  
  ## add to column in temperature_data and rename after column's population_id with longitude added onto the end
  temperature_data$temp <- waterTemp
  pop_id <- as.character(paste(unique_pairs$genus_species[num_unique], unique_pairs$latitude[num_unique], unique_pairs$longitude[num_unique], sep = "_"))
  colnames(temperature_data)[num_unique+1]<- pop_id
  
  print(paste("Done with ", pop_id, " (population number ", num_unique, "/188)", sep = ""))
  num_unique <- num_unique + 1
}

## save to RDS
precious_temps_freshdaily_arr <- temperature_data

temperature_data <- readRDS("./data-processed/precious_temps_freshdaily_arr.rds")

## convert from degrees K to degrees C
converted <- temperature_data
converted[, 2:189] <- converted[, 2:189] - 273.15

## convert date 
## starts at 1958-01-01
year <- c(round(0.5/365, digits = 3))
leap_year <- c(round(0.5/366, digits = 3))

i = 1
while (i < 366) {
  if (i < 365) {
    year = append(year, round((i+0.5)/365, digits = 3))
  }
  leap_year = append(leap_year, round((i+0.5)/366, digits = 3))
  i = i+1
}

rep = 1958
last_year = 2002
date <- c()

while (rep < last_year) {
  if (rep %% 4 == 0){
    if (rep == 1900) { 
      date <- append(date, rep+year, after = length(date))
    }
    else {
      date <- append(date, rep+leap_year, after = length(date))
    }
  }
  else {
    date <- append(date, rep+year, after = length(date))
  }
  rep = rep + 1
}

## replace column for date
converted$date <- as.vector(date)


temperature_data <- converted

## add rows for populations with the same lat and long but different species 
## will have same temperature data 
populations <- freshwater[!duplicated(freshwater[,c("genus_species", "latitude", "longitude")]),]

z = 1
while (z < nrow(populations) + 1) {
  i = 1
  while (i < nrow(unique_pairs) + 1) {
    
    same_pop <- (unique_pairs$latitude[i] == populations$latitude[z] & 
                   unique_pairs$longitude[i] == populations$longitude[z] & 
                   unique_pairs$genus_species[i] == populations$genus_species[z])
    
    same_loc <- (unique_pairs$latitude[i] == populations$latitude[z] & 
                   unique_pairs$longitude[i] == populations$longitude[z])
    
    if (!same_pop & same_loc) {
      pop_id <- as.character(paste(populations$genus_species[z], populations$latitude[z],populations$longitude[z], sep = "_"))
      temperature_data$temps <- temperature_data[,i+1]
      colnames(temperature_data)[length(temperature_data)]<- pop_id
      same_pop = FALSE
      same_loc = FALSE
      i = i+1
    }
    else {
      i = i + 1
    }
    
  }
  z = z+1
}



## write temperature file to processed data 
write.csv(temperature_data, "./data-processed/arr_freshwater-temp-data.csv", row.names = FALSE)
