## script used to get temperature data for terrestrial species from Berkeley Earth gridded TMAX data 
library(ncdf4)
library(janitor)
library(tidyverse)

## read in gridded data in nc file for the file_index from berkeley earth and store data in R workspace 
filename <- paste("./Gridded_daily_TMAX/Complete_TMAX_Daily_LatLong1_1880.nc", sep = "")
ncfile <- nc_open(filename)

## create variables for things needed to use data
lat <- ncvar_get(ncfile, "latitude")
long <- ncvar_get(ncfile, "longitude")

## close the file
nc_close(ncfile)

## bring in data
cadillac <- read.csv("data-processed/intratherm_version-for-arr.csv")

## filter out rows of data we cannot use  
cadillac <- cadillac %>%
  filter(!is.na(latitude)) %>%
  filter(!is.na(longitude)) %>%
  filter(!is.na(acclim_temp)) %>%
  filter(!is.na(lifespan_days)) %>%
  filter(lifespan_days != "unk" || "kunk") %>%
  filter(parameter_tmax_or_tmin == "tmax") 

cadillac <- droplevels(cadillac)

## get rid of duplicate population rows since all will have the same temp data
unique_pairs <- cadillac[!duplicated(cadillac[,c("latitude", "longitude")]),]

## temps:
temperature_data <- data.frame(matrix(nrow = 51043))

## for each population:
num_unique <- 1
while (num_unique < nrow(unique_pairs)+1) {
  
  loc_cumulative <- data.frame(matrix(ncol=2))
  colnames(loc_cumulative) <- c("date", "temp_value")
  rep = 1880
  
  while (rep < 2020) {
    print(paste("On population number ", num_unique, ", getting temp data from ", rep, sep = ""))
    ## read in gridded data in nc file for the file_index from berkeley earth and store data in R workspace 
    filename <- paste("./Gridded_daily_TMAX/Complete_TMAX_Daily_LatLong1_", rep, ".nc", sep = "")
    ncfile <- nc_open(filename)
    
    ## create variables for things needed to use data
    date <- ncvar_get(ncfile, "date_number")
    arr.anom <-ncvar_get(ncfile, "temperature")
    arr.clim <- ncvar_get(ncfile, "climatology")
    
    nc_close(ncfile)
    
    ## get clim and anom data for collection location of species 
    ## NaN here if location does not have data 
    loc_long_index <- which.min(abs(long - unique_pairs$longitude[num_unique]))
    loc_lat_index <- which.min(abs(lat - unique_pairs$latitude[num_unique]))
    loc.anom <- arr.anom[loc_long_index,loc_lat_index,]
    loc.clim.365d <- arr.clim[loc_long_index,loc_lat_index,]
    
    ## account for leap year - duplicate day index added on feb 28 in clim array (this seems to be how they dealt with it when calculating anom)
    index_59 <- loc.clim.365d[59]
    loc.clim.366d <- append(loc.clim.365d, index_59, after = 59)
    
    ## repeat day list loc.clim.365d on normal years + loc.clim.366d on leap years 
    last_year = (rep + 10)
    loc.clim <- c()
    while (rep < last_year) {
      if (rep == 2019) {
        loc.clim <- append(loc.clim, loc.clim.365d[1:273], after = length(loc.clim))
      }
      else if (rep %% 4 == 0){
        if (rep == 1900) { 
          loc.clim <- append(loc.clim, loc.clim.365d, after = length(loc.clim))
        }
        else {
          loc.clim <- append(loc.clim, loc.clim.366d, after = length(loc.clim))
        }
      }
      else {
        loc.clim <- append(loc.clim, loc.clim.365d, after = length(loc.clim))
      }
      rep = rep + 1
    }
    
    ## create dataframe of actual temp values at location by adding anomaly and climatology values over the 10 years 
    temp_list <- c()
    max <- rep
    rep <- rep - 10
    d <- 1
    
    while (rep < max) {
      if (rep == 2019) {
        temps <- loc.anom[d:(d+272)] + loc.clim[d:(d+272)]
        temp_list <- append(temp_list, temps, after = length(temp_list))
        d = d + 273
      }
      else if (rep %% 4 == 0) {
        if (rep == 1900) {
          temps <- loc.anom[d:(d+364)] + loc.clim[d:(d+364)]
          temp_list <- append(temp_list, temps, after = length(temp_list))
          d = d + 365
        }
        else { 
          temps <- loc.anom[d:(d+365)] + loc.clim[d:(d+365)]
          temp_list <- append(temp_list, temps, after = length(temp_list))
          d = d + 366
        }
      }
      else {
        temps <- loc.anom[d:(d+364)] + loc.clim[d:(d+364)]
        temp_list <- append(temp_list, temps, after = length(temp_list))
        d = d + 365
      }
      rep = rep + 1
    }
    
    ## make dataframe of date and corresponding temp values
    loc <- data.frame(date[], temp_list[])
    colnames(loc) <- c("date", "temp_value")
    
    ## add loc to loc_cumulative so one data frame contains data from all 10 year datasets
    loc_cumulative <- rbind(loc_cumulative, loc)
  }
  
  pop_id <- as.character(paste(unique_pairs$genus_species[num_unique], unique_pairs$latitude[num_unique], unique_pairs$longitude[num_unique], sep = "_"))
  
  ## add column for population to temperature data 
  if (num_unique == 1){
    temperature_data <- cbind(temperature_data, loc_cumulative)
    temperature_data <- temperature_data[,-1]
    colnames(temperature_data)[num_unique+1] <- pop_id
  }
  else {
    temperature_data <- cbind(temperature_data, loc_cumulative[,2])
    colnames(temperature_data)[num_unique+1] <- pop_id
  }
  
  num_unique = num_unique + 1;
  
}

temperature_data <- temperature_data[-1,]
temperature_data <- readRDS("./data-processed/precious_temps_original.rds")

## figure out which are NA and why:
isNA <- select(temperature_data, as.vector(which(colSums(is.na(temperature_data)) == nrow(temperature_data))))
indexNA <- as.vector(which(colSums(is.na(temperature_data)) == nrow(temperature_data)))
uniqueNA <- unique_pairs[indexNA-1,]

## get rid of marine and freshwater missing columns as they are unlikely to be available or accurate:
uniqueNA <- uniqueNA %>%
  filter(realm_general2 == "Terrestrial") %>%
  mutate(population_id = as.character(paste(genus_species, latitude, longitude, sep = "_")))

isNA <- isNA[,which(colnames(isNA) %in% uniqueNA$population_id)]

missing_long <- c()
missing_lat <- c()

## see what lat lon was and google maps it to investigate 
x=1
while (x < length(uniqueNA$genus_species) + 1) {
  missing_long <- append(missing_long, long[which.min(abs(long - uniqueNA$longitude[x]))])
  missing_lat <- append(missing_lat, lat[which.min(abs(lat - uniqueNA$latitude[x]))])
  x=x+1
}

#### find new grid squares by looking on google maps:
## move closer to Boston land mass:
missing_long[7] <- -70.5
missing_lat[7] <- 41.5

## five are on Martinique
missing_long[1:5] <- -60.5
missing_lat[1:5] <- 14.5

## no accurate temp data can be obtained since no data for island of St. Croix
missing_long[8] <- NA
missing_lat[8] <- NA

## Buchan, Victoria:
missing_long[9] <- 148.5
missing_lat[9] <- -37.5

## Murrindal, Victoria
missing_long[10] <- 148.5
missing_lat[10]<- -36.5

## Godsford, NSW:
missing_long[11] <- 150.5
missing_lat[11]<- -33.5
missing_long[13] <- 150.5
missing_lat[13]<- -33.5

## Wilson's Promontory National Park, Victoria
missing_long[12] <- 146.5
missing_lat[12] <- -38.5

## no accurate temp data can be obtained since island 
missing_long[6] <- NA
missing_lat[6] <- NA
missing_long[14] <- NA
missing_lat[14] <- NA

## Coconut grove miami:
missing_long[15:16] <- -80.5
missing_lat[14] <- -80.5


##assign new lats and longs to use when getting temp data 
uniqueNA$new_lat <- missing_lat
uniqueNA$new_long <- missing_long


## get temp data for populations that were NA
temperature_data_NA <- data.frame(matrix(nrow = 51043))
num_unique <- 1
while (num_unique < length(uniqueNA$population_id) + 1) {
  
  loc_cumulative <- data.frame(matrix(ncol=2))
  colnames(loc_cumulative) <- c("date", "temp_value")
  rep = 1880
  
  while (rep < 2020) {
    print(paste("On population number ", num_unique, ", getting temp data from ", rep, sep = ""))
    ## read in gridded data in nc file for the file_index from berkeley earth and store data in R workspace 
    filename <- paste("Gridded_daily_TMAX/Complete_TMAX_Daily_LatLong1_", rep, ".nc", sep = "")
    ncfile <- nc_open(filename)
    
    ## create variables for things needed to use data
    date <- ncvar_get(ncfile, "date_number")
    arr.anom <-ncvar_get(ncfile, "temperature")
    arr.clim <- ncvar_get(ncfile, "climatology")
    
    nc_close(ncfile)
    
    ## get clim and anom data for collection location of species 
    ## NaN here if location does not have data 
    loc_long_index <- which.min(abs(long - uniqueNA$new_long[num_unique]))
    loc_lat_index <- which.min(abs(lat - uniqueNA$new_lat[num_unique]))
    loc.anom <- arr.anom[loc_long_index,loc_lat_index,]
    loc.clim.365d <- arr.clim[loc_long_index,loc_lat_index,]
    
    ## account for leap year - duplicate day index added on feb 28 in clim array (this seems to be how they dealt with it when calculating anom)
    index_59 <- loc.clim.365d[59]
    loc.clim.366d <- append(loc.clim.365d, index_59, after = 59)
    
    ## repeat day list loc.clim.365d on normal years + loc.clim.366d on leap years 
    last_year = (rep + 10)
    loc.clim <- c()
    while (rep < last_year) {
      if (rep == 2019) {
        loc.clim <- append(loc.clim, loc.clim.365d[1:273], after = length(loc.clim))
      }
      else if (rep %% 4 == 0){
        if (rep == 1900) { 
          loc.clim <- append(loc.clim, loc.clim.365d, after = length(loc.clim))
        }
        else {
          loc.clim <- append(loc.clim, loc.clim.366d, after = length(loc.clim))
        }
      }
      else {
        loc.clim <- append(loc.clim, loc.clim.365d, after = length(loc.clim))
      }
      rep = rep + 1
    }
    
    ## create dataframe of actual temp values at location by adding anomaly and climatology values over the 10 years 
    temp_list <- c()
    max <- rep
    rep <- rep - 10
    d <- 1
    
    while (rep < max) {
      if (rep == 2019) {
        temps <- loc.anom[d:(d+272)] + loc.clim[d:(d+272)]
        temp_list <- append(temp_list, temps, after = length(temp_list))
        d = d + 273
      }
      else if (rep %% 4 == 0) {
        if (rep == 1900) {
          temps <- loc.anom[d:(d+364)] + loc.clim[d:(d+364)]
          temp_list <- append(temp_list, temps, after = length(temp_list))
          d = d + 365
        }
        else { 
          temps <- loc.anom[d:(d+365)] + loc.clim[d:(d+365)]
          temp_list <- append(temp_list, temps, after = length(temp_list))
          d = d + 366
        }
      }
      else {
        temps <- loc.anom[d:(d+364)] + loc.clim[d:(d+364)]
        temp_list <- append(temp_list, temps, after = length(temp_list))
        d = d + 365
      }
      rep = rep + 1
    }
    
    ## make dataframe of date and corresponding temp values
    loc <- data.frame(date[], temp_list[])
    colnames(loc) <- c("date", "temp_value")
    
    ## add loc to loc_cumulative so one data frame contains data from all 10 year datasets
    loc_cumulative <- rbind(loc_cumulative, loc)
  }
  
  pop_id <- as.character(paste(uniqueNA$genus_species[num_unique], uniqueNA$latitude[num_unique], uniqueNA$longitude[num_unique], sep = "_"))
  
  ## add column for population to temperature data 
  if (num_unique == 1){
    temperature_data_NA <- cbind(temperature_data_NA, loc_cumulative)
    temperature_data_NA <- temperature_data_NA[,-1]
    colnames(temperature_data_NA)[num_unique+1] <- pop_id
  }
  else {
    temperature_data_NA <- cbind(temperature_data_NA, loc_cumulative[,2])
    colnames(temperature_data_NA)[num_unique+1] <- pop_id
  }
  
  num_unique = num_unique + 1;
  
}


temperature_data_NA <- temperature_data_NA[-1,]
temperature_data_NA <- readRDS("data-processed/arr_precious_NA_temps.rds") %>%
  select(-date)

##combine:
temperature_data <- temperature_data[,-as.vector(which(colnames(temperature_data) %in% colnames(temperature_data_NA)))] ## remove NA columns
temperature_data <- cbind(temperature_data, temperature_data_NA) ## add columns that were NA back


## add rows for populations with the same lat and long but different species 
## will have same temperature data 
populations <- cadillac[!duplicated(cadillac[,c("genus_species", "latitude", "longitude")]),]

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



## write to file
write.csv(temperature_data, "data-processed/arr_temp-data.csv", row.names = FALSE)






## -------------------------------------------------------------------------------------
## get Tavg data 

## Tavg files for 1900 + 1910 + 1920 are corrupt - start from 1930 instead 
## read in gridded data in nc file for the file_index from berkeley earth and store data in R workspace 
filename <- paste("./Berkeley_Tavg/Complete_TAVG_Daily_LatLong1_1930.nc", sep = "")
ncfile <- nc_open(filename)

## create variables for things needed to use data
lat <- ncvar_get(ncfile, "latitude")
long <- ncvar_get(ncfile, "longitude")

## close the file
nc_close(ncfile)

## bring in data
cadillac <- read.csv("./data-raw/intratherm-may-2020-squeaky-clean.csv")

## filter out rows of data we cannot use  
cadillac <- cadillac %>%
  filter(!is.na(latitude)) %>%
  filter(!is.na(longitude)) %>%
  filter(!is.na(acclim_temp)) %>%
  filter(!is.na(lifespan_days)) %>%
  filter(lifespan_days != "unk" || "kunk") %>%
  filter(parameter_tmax_or_tmin == "tmax") 

## get rid of duplicate population rows since all will have the same temp data
unique_pairs <- cadillac[!duplicated(cadillac[,c("latitude", "longitude")]),]

## temps:
temperature_data <- data.frame(matrix(nrow = 32294))

## for each population:
num_unique <- 1
while (num_unique < nrow(unique_pairs) + 1) {
  
  loc_cumulative <- data.frame(matrix(ncol=2))
  colnames(loc_cumulative) <- c("date", "temp_value")
  rep = 1930
  
  while (rep < 2013) {
    print(paste("On population number ", num_unique, ", getting temp data from ", rep, sep = ""))
    ## read in gridded data in nc file for the file_index from berkeley earth and store data in R workspace 
    filename <- paste("./Berkeley_Tavg/Complete_TAVG_Daily_LatLong1_", rep, ".nc", sep = "")
    ncfile <- nc_open(filename)
    
    ## create variables for things needed to use data
    date <- ncvar_get(ncfile, "date_number")
    arr.anom <-ncvar_get(ncfile, "temperature")
    arr.clim <- ncvar_get(ncfile, "climatology")
    
    nc_close(ncfile)
    
    ## get clim and anom data for collection location of species 
    ## NaN here if location does not have data 
    loc_long_index <- which.min(abs(long - unique_pairs$longitude[num_unique]))
    loc_lat_index <- which.min(abs(lat - unique_pairs$latitude[num_unique]))
    loc.anom <- arr.anom[loc_long_index,loc_lat_index,]
    loc.clim.365d <- arr.clim[loc_long_index,loc_lat_index,]
    
    ## account for leap year - duplicate day index added on feb 28 in clim array (this seems to be how they dealt with it when calculating anom)
    index_59 <- loc.clim.365d[59]
    loc.clim.366d <- append(loc.clim.365d, index_59, after = 59)
    
    ## repeat day list loc.clim.365d on normal years + loc.clim.366d on leap years 
    last_year = (rep + 10)
    loc.clim <- c()
    
    if(last_year == 2020) {
      last_year = 2019
    }
    
    while (rep < last_year) {
      if (rep == 2018) {
        loc.clim <- append(loc.clim, loc.clim.365d[1:151], after = length(loc.clim))
      }
      else if (rep %% 4 == 0){
        if (rep == 1900) { 
          loc.clim <- append(loc.clim, loc.clim.365d, after = length(loc.clim))
        }
        else {
          loc.clim <- append(loc.clim, loc.clim.366d, after = length(loc.clim))
        }
      }
      else {
        loc.clim <- append(loc.clim, loc.clim.365d, after = length(loc.clim))
      }
      rep = rep + 1
    }
    
    ## create dataframe of actual temp values at location by adding anomaly and climatology values over the 10 years 
    temp_list <- c()
    max <- rep
    rep <- rep - 10
    d <- 1
    
    if (rep == 2009) {
      rep <- 2010
    }
    
    while (rep < max) {
      if (rep == 2018) {
        temps <- loc.anom[d:(d+150)] + loc.clim[d:(d+150)]
        temp_list <- append(temp_list, temps, after = length(temp_list))
        d = d + 150
      }
      else if (rep %% 4 == 0) {
        if (rep == 1900) {
          temps <- loc.anom[d:(d+364)] + loc.clim[d:(d+364)]
          temp_list <- append(temp_list, temps, after = length(temp_list))
          d = d + 365
        }
        else { 
          temps <- loc.anom[d:(d+365)] + loc.clim[d:(d+365)]
          temp_list <- append(temp_list, temps, after = length(temp_list))
          d = d + 366
        }
      }
      else {
        temps <- loc.anom[d:(d+364)] + loc.clim[d:(d+364)]
        temp_list <- append(temp_list, temps, after = length(temp_list))
        d = d + 365
      }
      rep = rep + 1
    }
    
    ## make dataframe of date and corresponding temp values
    loc <- data.frame(date[], temp_list[])
    colnames(loc) <- c("date", "temp_value")
    
    ## add loc to loc_cumulative so one data frame contains data from all 10 year datasets
    loc_cumulative <- rbind(loc_cumulative, loc)
  }
  
  pop_id <- as.character(paste(unique_pairs$genus_species[num_unique], unique_pairs$latitude[num_unique], unique_pairs$longitude[num_unique], sep = "_"))
  
  ## add column for population to temperature data 
  if (num_unique == 1){
    temperature_data <- cbind(temperature_data, loc_cumulative)
    temperature_data <- temperature_data[,-1]
    colnames(temperature_data)[num_unique+1] <- pop_id
  }
  else {
    temperature_data <- cbind(temperature_data, loc_cumulative[,2])
    colnames(temperature_data)[num_unique+1] <- pop_id
  }
  
  num_unique = num_unique + 1;
  
}


temperature_data <- temperature_data[-1,]
precious_temps_original_tavg = temperature_data

temperature_data <- readRDS("./data-processed/precious_temps_original_tavg.rds")

## add rows for populations with the same lat and long but different species 
## will have same temperature data 
populations <- cadillac[!duplicated(cadillac[,c("genus_species", "latitude", "longitude")]),]

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



## write to file
write.csv(temperature_data, "./data-processed/arr_temp-data_tavg.csv", row.names = FALSE)

