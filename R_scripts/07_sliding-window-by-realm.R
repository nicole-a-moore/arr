## sliding window by realm 
library(gdata)
library(evobiR)
library(tidyverse)

intratherm <- read.csv("./data-processed/intratherm_sliding-window-ready.csv") 

## subset: 
cadillac <- subset(intratherm, select = c(genus_species, population_id, acclim_temp, latitude, longitude, 
                                          genus, species, parameter_value, parameter_tmax_or_tmin, 
                                          realm_general2, lifespan_days,
                                          season_when_away_100km_start, season_when_away_100km_stop, 
                                          season_inactive_start,  season_inactive_stop,
                                          season_inactive_start2,  season_inactive_stop2,
                                          maximum_body_size_svl_hbl_cm, elevation_of_collection))


####      TERRESTRIAL       ####
################################
temp_data <- read_csv("./data-processed/arr_temp-data_tavg.csv")

## get rid of duplicate population rows and any marine data 
terrestrial <- cadillac %>%
  filter(realm_general2 == "Terrestrial")

unique_pairs <- subset(terrestrial, !duplicated(terrestrial$population_id)) 
  

sd_cumulative <- c()
experienced_var_mean <- c()
experienced_var_max <- c()

## for each terrestrial species:
num_unique <- 1
while (num_unique < nrow(unique_pairs)+1) {
  
  print(paste("Marking months when away for population ", num_unique, sep = ""))
  ## mark months that species is away as true 
  months <- c(jan <- FALSE, feb <- FALSE, mar <- FALSE, apr <- FALSE, may <- FALSE, jun <- FALSE, 
              jul <- FALSE, aug <- FALSE, sep <- FALSE, oct <- FALSE, nov <- FALSE, dec <- FALSE)
  
  if (!is.na(unique_pairs$season_when_away_100km_start[num_unique])) {
    away_start <- unique_pairs$season_when_away_100km_start[num_unique]
    away_stop <- unique_pairs$season_when_away_100km_stop[num_unique]
    
    months <- initialize_months(away_start, away_stop, months)
  }
  
  if (!is.na(unique_pairs$season_inactive_start[num_unique])) {
    away_start <- unique_pairs$season_inactive_start[num_unique]
    away_stop <- unique_pairs$season_inactive_stop[num_unique]
    
    months <- initialize_months(away_start, away_stop, months)
    if (!is.na(unique_pairs$season_inactive_start2[num_unique])) {
      away_start <- unique_pairs$season_inactive_start2[num_unique]
      away_stop <- unique_pairs$season_inactive_stop2[num_unique]
      
      months <- initialize_months(away_start, away_stop, months)
    }
  }
  
  print("Creating temp data frame for population...")
  ## make dataframe of date and temp values for collection location of that species 
  pop_id <- unique_pairs$population_id[num_unique]
  
  loc <- data.frame(matrix(nrow = 32293))
  loc[,1] <- temp_data$date
  loc[,2] <-  temp_data[,which(colnames(temp_data) == pop_id)]
  colnames(loc) <- c("date", "temp_value")
  
  print("Setting temps when away or inactive to NA...")
  ## set temps to NA:
  loc <- set_temps_to_NA(column = loc, months = months, realm = "Terrestrial")
  
  print("Sliding the window...")
  ## perform sliding window calculations:
  lifespan <- as.integer(unique_pairs$lifespan_days[num_unique])
  
  if(lifespan > nrow(loc)) {
    lifespan <- as.integer(nrow(loc)/365)
  }
  
  sd_vector <- SlidingWindow(FUN = sd_special, loc$temp_value, lifespan, 1)
  sd = sd_special(loc$temp_value)
  print(paste("Performed sliding window for ", num_unique, sep = ""))
  print(paste("Number of days iterating over: ", length(loc$temp_value),sep = ""))
  
  ## calculate mean sd and add to experienced variation vectors 
  pop_experienced_var_mean<- mean(sd_vector, na.rm = TRUE)
  pop_experienced_var_max<- max(sd_vector, na.rm = TRUE)
  ##LSV_max <- max(sd_vector, na.rm = TRUE)
  experienced_var_mean <- append(experienced_var_mean, pop_experienced_var_mean, 
                                 after = length(experienced_var_mean))
  experienced_var_max<- append(experienced_var_max, pop_experienced_var_max, 
                               after = length(experienced_var_max))
  sd_cumulative <- append(sd_cumulative, sd, after = length(sd_cumulative))
  
  print("Done! Moving on to next population :-)")
  ##move on to next species
  num_unique <- num_unique + 1
}

unique_pairs$experienced_var_mean <- experienced_var_mean
unique_pairs$experienced_var_max <- experienced_var_max

terrestrial_unique <- unique_pairs
saveRDS(terrestrial_unique, "./data-processed/terrestrial_expvar.rds")



####      FRESHWATER       ####
###############################
temp_data <- read_csv("./data-processed/arr_freshwater-temp-data.csv")

## get rid of duplicate population rows and any marine data 
freshwater <- cadillac %>%
  filter(realm_general2 == "Freshwater")

unique_pairs <- subset(freshwater, !duplicated(freshwater$population_id)) 


sd_cumulative <- c()
experienced_var_mean <- c()
experienced_var_max <- c()

## for each terrestrial species:
num_unique <- 1
while (num_unique < nrow(unique_pairs)+1) {
  
  print(paste("Marking months when away for population ", num_unique, sep = ""))
  ## mark months that species is away as true 
  months <- c(jan <- FALSE, feb <- FALSE, mar <- FALSE, apr <- FALSE, may <- FALSE, jun <- FALSE, 
              jul <- FALSE, aug <- FALSE, sep <- FALSE, oct <- FALSE, nov <- FALSE, dec <- FALSE)
  
  if (!is.na(unique_pairs$season_when_away_100km_start[num_unique])) {
    away_start <- unique_pairs$season_when_away_100km_start[num_unique]
    away_stop <- unique_pairs$season_when_away_100km_stop[num_unique]
    
    months <- initialize_months(away_start, away_stop, months)
  }
  
  if (!is.na(unique_pairs$season_inactive_start[num_unique])) {
    away_start <- unique_pairs$season_inactive_start[num_unique]
    away_stop <- unique_pairs$season_inactive_stop[num_unique]
    
    months <- initialize_months(away_start, away_stop, months)
    if (!is.na(unique_pairs$season_inactive_start2[num_unique])) {
      away_start <- unique_pairs$season_inactive_start2[num_unique]
      away_stop <- unique_pairs$season_inactive_stop2[num_unique]
      
      months <- initialize_months(away_start, away_stop, months)
    }
  }
  
  print("Creating temp data frame for population...")
  ## make dataframe of date and temp values for collection location of that species 
  pop_id <- unique_pairs$population_id[num_unique]
  
  loc <- data.frame(matrix(nrow = 16071))
  loc[,1] <- temp_data$date
  loc[,2] <-  temp_data[,which(colnames(temp_data) == pop_id)]
  colnames(loc) <- c("date", "temp_value")
  
  print("Setting temps when away or inactive to NA...")
  ## set temps to NA:
  loc <- set_temps_to_NA(column = loc, months = months, realm = "Freshwater")
  
  print("Sliding the window...")
  
  ## perform sliding window calculations:
  lifespan <- as.integer(unique_pairs$lifespan_days[num_unique])
  
  if(lifespan > nrow(loc)) {
    lifespan <- as.integer(nrow(loc)/365)
  }
  
  sd_vector <- SlidingWindow(FUN = sd_special, loc$temp_value, lifespan, 1)
  sd = sd_special(loc$temp_value)
  print(paste("Performed sliding window for ", num_unique, sep = ""))
  print(paste("Number of days iterating over: ", length(loc$temp_value),sep = ""))
  
  ## calculate mean sd and add to experienced variation vectors 
  pop_experienced_var_mean<- mean(sd_vector, na.rm = TRUE)
  pop_experienced_var_max<- max(sd_vector, na.rm = TRUE)
  ##LSV_max <- max(sd_vector, na.rm = TRUE)
  experienced_var_mean <- append(experienced_var_mean, pop_experienced_var_mean, 
                                 after = length(experienced_var_mean))
  experienced_var_max<- append(experienced_var_max, pop_experienced_var_max, 
                               after = length(experienced_var_max))
  sd_cumulative <- append(sd_cumulative, sd, after = length(sd_cumulative))
  
  print("Done! Moving on to next population :-)")
  ##move on to next species
  num_unique <- num_unique + 1
}

unique_pairs$experienced_var_mean <- experienced_var_mean
unique_pairs$experienced_var_max <- experienced_var_max

freshwater_unique <- unique_pairs
saveRDS(freshwater_unique, "./data-processed/freshwater_expvar.rds")





####         MARINE         ####
################################
temp_data <- read_csv("./data-processed/arr_marine-temp-data.csv")

## get rid of duplicate population rows and any marine data 
marine <- cadillac %>%
  filter(realm_general2 == "Marine")

unique_pairs <- subset(marine, !duplicated(marine$population_id)) 


sd_cumulative <- c()
experienced_var_mean <- c()
experienced_var_max <- c()

## for each terrestrial species:
num_unique <- 1
while (num_unique < nrow(unique_pairs)+1) {
  
  print(paste("Marking months when away for population ", num_unique, sep = ""))
  ## mark months that species is away as true 
  months <- c(jan <- FALSE, feb <- FALSE, mar <- FALSE, apr <- FALSE, may <- FALSE, jun <- FALSE, 
              jul <- FALSE, aug <- FALSE, sep <- FALSE, oct <- FALSE, nov <- FALSE, dec <- FALSE)
  
  if (!is.na(unique_pairs$season_when_away_100km_start[num_unique])) {
    away_start <- unique_pairs$season_when_away_100km_start[num_unique]
    away_stop <- unique_pairs$season_when_away_100km_stop[num_unique]
    
    months <- initialize_months(away_start, away_stop, months)
  }
  
  if (!is.na(unique_pairs$season_inactive_start[num_unique])) {
    away_start <- unique_pairs$season_inactive_start[num_unique]
    away_stop <- unique_pairs$season_inactive_stop[num_unique]
    
    months <- initialize_months(away_start, away_stop, months)
    if (!is.na(unique_pairs$season_inactive_start2[num_unique])) {
      away_start <- unique_pairs$season_inactive_start2[num_unique]
      away_stop <- unique_pairs$season_inactive_stop2[num_unique]
      
      months <- initialize_months(away_start, away_stop, months)
    }
  }
  
  print("Creating temp data frame for population...")
  ## make dataframe of date and temp values for collection location of that species 
  pop_id <- unique_pairs$population_id[num_unique]
  
  loc <- data.frame(matrix(nrow = 14096))
  loc[,1] <- temp_data$date
  loc[,2] <-  temp_data[,which(colnames(temp_data) == pop_id)]
  colnames(loc) <- c("date", "temp_value")
  
  print("Setting temps when away or inactive to NA...")
  ## set temps to NA:
  loc <- set_temps_to_NA(column = loc, months = months, realm = "Marine")
  
  print("Sliding the window...")
  ## perform sliding window calculations:
  sd_vector <- SlidingWindow(FUN = sd_special, loc$temp_value, as.integer(unique_pairs$lifespan_days[num_unique]), 1)
  sd = sd_special(loc$temp_value)
  print(paste("Performed sliding window for ", num_unique, sep = ""))
  print(paste("Number of days iterating over: ", length(loc$temp_value),sep = ""))
  
  ## calculate mean sd and add to experienced variation vectors 
  pop_experienced_var_mean<- mean(sd_vector, na.rm = TRUE)
  pop_experienced_var_max<- max(sd_vector, na.rm = TRUE)
  ##LSV_max <- max(sd_vector, na.rm = TRUE)
  experienced_var_mean <- append(experienced_var_mean, pop_experienced_var_mean, 
                                 after = length(experienced_var_mean))
  experienced_var_max<- append(experienced_var_max, pop_experienced_var_max, 
                               after = length(experienced_var_max))
  sd_cumulative <- append(sd_cumulative, sd, after = length(sd_cumulative))
  
  print("Done! Moving on to next population :-)")
  ##move on to next species
  num_unique <- num_unique + 1
}

unique_pairs$experienced_var_mean <- experienced_var_mean
unique_pairs$experienced_var_max <- experienced_var_max

marine_unique <- unique_pairs
saveRDS(marine_unique, "./data-processed/marine_expvar.rds")




## combine all of the output into one and then merge with intratherm ------------
combined_unique <- rbind(marine_unique, terrestrial_unique, freshwater_unique)
combined_unique <- combined_unique %>%
  select(population_id, experienced_var_mean, experienced_var_max)

merged <- left_join(intratherm, combined_unique)


## write to new file: 
write.csv(merged, "./data-processed/arr_sliding-window-by-realm-output.csv", row.names = FALSE)


#
###
#######
##############
##########################
################################################
# STANDARD DEVIATION FUNCTION TO DEAL WITH NA ##
################################################
sd_special <- function(x) {
  return(sd(x, na.rm = TRUE))
}
#
###
#######
##############
##########################
################################################
# QUANTILE FUNCTION TO DEAL WITH NA ############
################################################
IquantileR_special <- function(x) {
  return(quantile(x, probs=0.975, na.rm = TRUE) - quantile(x, probs=0.025, na.rm = TRUE))
}

#
###
#######
##############
##########################
############################################
###########################################################################
# FUNCTION FOR SETTING TEMP DATA TO NA DURING SEASON WHEN AWAY/INACTIVE ###
###########################################################################
## takes a temp data frame with columns for each population's temp data and a realm argument specifying which dataset the temp data came from 
## column names must contain the population's genus_species name
## returns a new version of the temp_data with temps during seasons away and seasons inactive set to NA
set_temps_to_NA <- function(column, realm, months) {

    ## set temps in column to NA 
    ## different conditions for each realm since all have different start and stop dates for temp data
    if (realm == "Terrestrial") {
      max <- 2019
      rep <- 1930
      day_index <- 1 
    }
    else if (realm == "Marine") {
      ## 1981-09-01 to 2020-04-04
      max = 2021
      rep = 1982
      day_index <- 1
    }
    else if (realm == "Freshwater") {
      max = 2002
      rep = 1958
      day_index <- 1
    }
    
    while (rep < max) {
      if (rep == 2019) {
        if (isTRUE(months[1])) {
          column[day_index:(day_index + 30),] <- NA
        }
        if (isTRUE(months[2])) {
          column[(day_index + 30):(day_index + 58),] <- NA
        }
        if (isTRUE(months[3])) {
          column[(day_index + 58):(day_index + 89),] <- NA
        }
        if (isTRUE(months[4])) {
          column[(day_index + 89):(day_index + 119),] <- NA
        }
        if (isTRUE(months[5])) {
          column[(day_index + 119):(day_index + 150),] <- NA
        }
        if (isTRUE(months[6])) {
          column[(day_index + 150):(day_index + 180),] <- NA
        }
        if (isTRUE(months[7])) {
          column[(day_index + 180):(day_index + 211),] <- NA
        }
        if (isTRUE(months[8])) {
          column[(day_index + 211):(day_index + 242),] <- NA
        }
        if (isTRUE(months[9])) {
          column[(day_index + 242):(day_index + 272),] <- NA
        }
        print(paste("Initialized ", rep, " months to NA - last year!", sep = ""))
        rep <- rep + 1
        day_index <- day_index + 365
      }
      else if (rep == 1900) {
        if (isTRUE(months[1])) {
          column[day_index:(day_index + 30),] <- NA
        }
        if (isTRUE(months[2])) {
          column[(day_index + 30):(day_index + 58),] <- NA
        }
        if (isTRUE(months[3])) {
          column[(day_index + 58):(day_index + 89),] <- NA
        }
        if (isTRUE(months[4])) {
          column[(day_index + 89):(day_index + 119),] <- NA
        }
        if (isTRUE(months[5])) {
          column[(day_index + 119):(day_index + 150),] <- NA
        }
        if (isTRUE(months[6])) {
          column[(day_index + 150):(day_index + 180),] <- NA
        }
        if (isTRUE(months[7])) {
          column[(day_index + 180):(day_index + 211),] <- NA
        }
        if (isTRUE(months[8])) {
          column[(day_index + 211):(day_index + 242),] <- NA
        }
        if (isTRUE(months[9])) {
          column[(day_index + 242):(day_index + 272),] <- NA
        }
        if (isTRUE(months[10])) {
          column[(day_index + 272):(day_index + 303),] <- NA
        }
        if (isTRUE(months[11])) {
          column[(day_index + 303):(day_index + 333),] <- NA
        }
        if (isTRUE(months[12])) {
          column[(day_index + 333):(day_index + 364),] <- NA
        }
        print(paste("Initialized ", rep, " months to NA (not a leap year)", sep = ""))
        rep <- rep + 1
        day_index <- day_index + 365
      }
      else if (rep %% 4 == 0) {
        if (isTRUE(months[1])) {
          column[day_index:(day_index + 30),] <- NA
        }
        if (isTRUE(months[2])) {
          column[(day_index + 30):(day_index + 59),] <- NA
        }
        if (isTRUE(months[3])) {
          column[(day_index + 59):(day_index + 90),] <- NA
        }
        if (isTRUE(months[4])) {
          column[(day_index + 90):(day_index + 120),] <- NA
        }
        if (isTRUE(months[5])) {
          column[(day_index + 120):(day_index + 151),] <- NA
        }
        if (isTRUE(months[6])) {
          column[(day_index + 151):(day_index + 181),] <- NA
        }
        if (isTRUE(months[7])) {
          column[(day_index + 181):(day_index + 212),] <- NA
        }
        if (isTRUE(months[8])) {
          column[(day_index + 212):(day_index + 243),] <- NA
        }
        if (isTRUE(months[9])) {
          column[(day_index + 243):(day_index + 273),] <- NA
        }
        if (isTRUE(months[10])) {
          column[(day_index + 273):(day_index + 304),] <- NA
        }
        if (isTRUE(months[11])) {
          column[(day_index + 304):(day_index + 334),] <- NA
        }
        if (isTRUE(months[12])) {
          column[(day_index + 334):(day_index + 365),] <- NA
        }
        print(paste("Initialized ", rep, " months to NA (leap year)", sep = ""))
        rep <- rep + 1
        day_index <- day_index + 366
      }
      else {
        if (isTRUE(months[1])) {
          column[day_index:(day_index + 30),] <- NA
        }
        if (isTRUE(months[2])) {
          column[(day_index + 30):(day_index + 58),] <- NA
        }
        if (isTRUE(months[3])) {
          column[(day_index + 58):(day_index + 89),] <- NA
        }
        if (isTRUE(months[4])) {
          column[(day_index + 89):(day_index + 119),] <- NA
        }
        if (isTRUE(months[5])) {
          column[(day_index + 119):(day_index + 150),] <- NA
        }
        if (isTRUE(months[6])) {
          column[(day_index + 150):(day_index + 180),] <- NA
        }
        if (isTRUE(months[7])) {
          column[(day_index + 180):(day_index + 211),] <- NA
        }
        if (isTRUE(months[8])) {
          column[(day_index + 211):(day_index + 242),] <- NA
        }
        if (isTRUE(months[9])) {
          column[(day_index + 242):(day_index + 272),] <- NA
        }
        if (isTRUE(months[10])) {
          column[(day_index + 272):(day_index + 303),] <- NA
        }
        if (isTRUE(months[11])) {
          column[(day_index + 303):(day_index + 333),] <- NA
        }
        if (isTRUE(months[12])) {
          column[(day_index + 333):(day_index + 364),] <- NA
        }
        print(paste("Initialized ", rep, " months to NA (not a leap year)", sep = ""))
        rep <- rep + 1
        day_index <- day_index + 365
      }
    }
  return(column)
}


#
###
#######
##############
##########################
############################################
#######################################################################
# FUNCTION FOR SETTING MONTHS LIST BASED ON START/STOP NUMERIC DATA ###
#######################################################################

initialize_months <- function(away_start, away_stop, months) {
  
  ##jan start
  if(away_start == 0.042) {
    if (away_stop == 0.125){
      months[1:2] = TRUE
    }
    else if (away_stop == 0.208){
      months[1:3] = TRUE
    }
    else if (away_stop == 0.292){
      months[1:4] = TRUE
    }
    else if (away_stop == 0.375){
      months[1:5] = TRUE
    }
    else if (away_stop == 0.458){
      months[1:6] = TRUE
    }
    else if (away_stop == 0.542){
      months[1:7] = TRUE
    }
    else if (away_stop == 0.625){
      months[1:8] = TRUE
    }
    else if (away_stop == 0.708){
      months[1:9] = TRUE
    }
    else if (away_stop == 0.792){
      months[1:10] = TRUE
    }
    else if (away_stop == 0.875){
      months[1:11] = TRUE
    }
    else if (away_stop == 0.958){
      months[1:12] = TRUE
    }
  }
  ##feb start
  else if(away_start == 0.125) {
    if (away_stop == 0.208){
      months[2:3] = TRUE
    }
    else if (away_stop == 0.292){
      months[2:4] = TRUE
    }
    else if (away_stop == 0.375){
      months[2:5] = TRUE
    }
    else if (away_stop == 0.458){
      months[2:6] = TRUE
    }
    else if (away_stop == 0.542){
      months[2:7] = TRUE
    }
    else if (away_stop == 0.625){
      months[2:8] = TRUE
    }
    else if (away_stop == 0.708){
      months[2:9] = TRUE
    }
    else if (away_stop == 0.792){
      months[2:10] = TRUE
    }
    else if (away_stop == 0.875){
      months[2:11] = TRUE
    }
    else if (away_stop == 0.958){
      months[2:12] = TRUE
    }
    else if (away_stop == 0.042){
      months[2:12] = TRUE
      months[1] = TRUE
    }
  }
  ## march start
  else if(away_start == 0.208) {
    if (away_stop == 0.292){
      months[3:4] = TRUE
    }
    else if (away_stop == 0.375){
      months[3:5] = TRUE
    }
    else if (away_stop == 0.458){
      months[3:6] = TRUE
    }
    else if (away_stop == 0.542){
      months[3:7] = TRUE
    }
    else if (away_stop == 0.625){
      months[3:8] = TRUE
    }
    else if (away_stop == 0.708){
      months[3:9] = TRUE
    }
    else if (away_stop == 0.792){
      months[3:10] = TRUE
    }
    else if (away_stop == 0.875){
      months[3:11] = TRUE
    }
    else if (away_stop == 0.958){
      months[3:12] = TRUE
    }
    else if (away_stop == 0.042){
      months[3:12] = TRUE
      months[1] = TRUE
    }
    else if (away_stop == 0.125){
      months[3:12] = TRUE
      months[1:2] = TRUE
    }
  }
  ## apr start
  else if(away_start == 0.292) {
    if (away_stop == 0.375){
      months[4:5] = TRUE
    }
    else if (away_stop == 0.458){
      months[4:6] = TRUE
    }
    else if (away_stop == 0.542){
      months[4:7] = TRUE
    }
    else if (away_stop == 0.625){
      months[4:8] = TRUE
    }
    else if (away_stop == 0.708){
      months[4:9] = TRUE
    }
    else if (away_stop == 0.792){
      months[4:10] = TRUE
    }
    else if (away_stop == 0.875){
      months[4:11] = TRUE
    }
    else if (away_stop == 0.958){
      months[4:12] = TRUE
    }
    else if (away_stop == 0.042){
      months[4:12] = TRUE
      months[1] = TRUE
    }
    else if (away_stop == 0.125){
      months[4:12] = TRUE
      months[1:2] = TRUE
    }
    else if (away_stop == 0.208){
      months[4:12] = TRUE
      months[1:3] = TRUE
    }
  }
  ## may start
  else if(away_start == 0.375) {
    if (away_stop == 0.458){
      months[5:6] = TRUE
    }
    else if (away_stop == 0.542){
      months[5:7] = TRUE
    }
    else if (away_stop == 0.625){
      months[5:8] = TRUE
    }
    else if (away_stop == 0.708){
      months[5:9] = TRUE
    }
    else if (away_stop == 0.792){
      months[5:10] = TRUE
    }
    else if (away_stop == 0.875){
      months[5:11] = TRUE
    }
    else if (away_stop == 0.958){
      months[5:12] = TRUE
    }
    else if (away_stop == 0.042){
      months[5:12] = TRUE
      months[1] = TRUE
    }
    else if (away_stop == 0.125){
      months[5:12] = TRUE
      months[1:2] = TRUE
    }
    else if (away_stop == 0.208){
      months[5:12] = TRUE
      months[1:3] = TRUE
    }
    else if (away_stop == 0.292){
      months[5:12] = TRUE
      months[1:4] = TRUE
    }
  }
  ## june start
  else if(away_start == 0.458) {
    if (away_stop == 0.542){
      months[6:7] = TRUE
    }
    else if (away_stop == 0.625){
      months[6:8] = TRUE
    }
    else if (away_stop == 0.708){
      months[6:9] = TRUE
    }
    else if (away_stop == 0.792){
      months[6:10] = TRUE
    }
    else if (away_stop == 0.875){
      months[6:11] = TRUE
    }
    else if (away_stop == 0.958){
      months[6:12] = TRUE
    }
    else if (away_stop == 0.042){
      months[6:12] = TRUE
      months[1] = TRUE
    }
    else if (away_stop == 0.125){
      months[6:12] = TRUE
      months[1:2] = TRUE
    }
    else if (away_stop == 0.208){
      months[6:12] = TRUE
      months[1:3] = TRUE
    }
    else if (away_stop == 0.292){
      months[6:12] = TRUE
      months[1:4] = TRUE
    }
    else if (away_stop == 0.375){
      months[6:12] = TRUE
      months[1:5] = TRUE
    }
  }
  ## july start
  else if(away_start == 0.542) {
    if (away_stop == 0.625){
      months[7:8] = TRUE
    }
    else if (away_stop == 0.708){
      months[7:9] = TRUE
    }
    else if (away_stop == 0.792){
      months[7:10] = TRUE
    }
    else if (away_stop == 0.875){
      months[7:11] = TRUE
    }
    else if (away_stop == 0.958){
      months[7:12] = TRUE
    }
    else if (away_stop == 0.042){
      months[7:12] = TRUE
      months[1] = TRUE
    }
    else if (away_stop == 0.125){
      months[7:12] = TRUE
      months[1:2] = TRUE
    }
    else if (away_stop == 0.208){
      months[7:12] = TRUE
      months[1:3] = TRUE
    }
    else if (away_stop == 0.292){
      months[7:12] = TRUE
      months[1:4] = TRUE
    }
    else if (away_stop == 0.375){
      months[7:12] = TRUE
      months[1:5] = TRUE
    }
    else if (away_stop == 0.458){
      months[7:12] = TRUE
      months[1:6] = TRUE
    }
  }
  ## aug start
  else if(away_start == 0.625) {
    if (away_stop == 0.708){
      months[8:9] = TRUE
    }
    else if (away_stop == 0.792){
      months[8:10] = TRUE
    }
    else if (away_stop == 0.875){
      months[8:11] = TRUE
    }
    else if (away_stop == 0.958){
      months[8:12] = TRUE
    }
    else if (away_stop == 0.042){
      months[8:12] = TRUE
      months[1] = TRUE
    }
    else if (away_stop == 0.125){
      months[8:12] = TRUE
      months[1:2] = TRUE
    }
    else if (away_stop == 0.208){
      months[8:12] = TRUE
      months[1:3] = TRUE
    }
    else if (away_stop == 0.292){
      months[8:12] = TRUE
      months[1:4] = TRUE
    }
    else if (away_stop == 0.375){
      months[8:12] = TRUE
      months[1:5] = TRUE
    }
    else if (away_stop == 0.458){
      months[8:12] = TRUE
      months[1:6] = TRUE
    }
    else if (away_stop == 0.542){
      months[8:12] = TRUE
      months[1:7] = TRUE
    }
  }
  ## sept start
  else if(away_start == 0.708) {
    if (away_stop == 0.792){
      months[9:10] = TRUE
    }
    else if (away_stop == 0.875){
      months[9:11] = TRUE
    }
    else if (away_stop == 0.958){
      months[9:12] = TRUE
    }
    else if (away_stop == 0.042){
      months[9:12] = TRUE
      months[1] = TRUE
    }
    else if (away_stop == 0.125){
      months[9:12] = TRUE
      months[1:2] = TRUE
    }
    else if (away_stop == 0.208){
      months[9:12] = TRUE
      months[1:3] = TRUE
    }
    else if (away_stop == 0.292){
      months[9:12] = TRUE
      months[1:4] = TRUE
    }
    else if (away_stop == 0.375){
      months[9:12] = TRUE
      months[1:5] = TRUE
    }
    else if (away_stop == 0.458){
      months[9:12] = TRUE
      months[1:6] = TRUE
    }
    else if (away_stop == 0.542){
      months[9:12] = TRUE
      months[1:7] = TRUE
    }
    else if (away_stop == 0.625){
      months[9:12] = TRUE
      months[1:8] = TRUE
    }
  }
  ## oct start
  else if(away_start == 0.792) {
    if (away_stop == 0.875){
      months[10:11] = TRUE
    }
    else if (away_stop == 0.958){
      months[10:12] = TRUE
    }
    else if (away_stop == 0.042){
      months[10:12] = TRUE
      months[1] = TRUE
    }
    else if (away_stop == 0.125){
      months[10:12] = TRUE
      months[1:2] = TRUE
    }
    else if (away_stop == 0.208){
      months[10:12] = TRUE
      months[1:3] = TRUE
    }
    else if (away_stop == 0.292){
      months[10:12] = TRUE
      months[1:4] = TRUE
    }
    else if (away_stop == 0.375){
      months[10:12] = TRUE
      months[1:5] = TRUE
    }
    else if (away_stop == 0.458){
      months[10:12] = TRUE
      months[1:6] = TRUE
    }
    else if (away_stop == 0.542){
      months[10:12] = TRUE
      months[1:7] = TRUE
    }
    else if (away_stop == 0.625){
      months[10:12] = TRUE
      months[1:8] = TRUE
    }
    else if (away_stop == 0.708){
      months[10:12] = TRUE
      months[1:9] = TRUE
    }
  }
  ## nov start
  else if(away_start == 0.875) {
    if (away_stop == 0.958){
      months[11:12] = TRUE
    }
    else if (away_stop == 0.042){
      months[11:12] = TRUE
      months[1] = TRUE
    }
    else if (away_stop == 0.125){
      months[11:12] = TRUE
      months[1:2] = TRUE
    }
    else if (away_stop == 0.208){
      months[11:12] = TRUE
      months[1:3] = TRUE
    }
    else if (away_stop == 0.292){
      months[11:12] = TRUE
      months[1:4] = TRUE
    }
    else if (away_stop == 0.375){
      months[11:12] = TRUE
      months[1:5] = TRUE
    }
    else if (away_stop == 0.458){
      months[11:12] = TRUE
      months[1:6] = TRUE
    }
    else if (away_stop == 0.542){
      months[11:12] = TRUE
      months[1:7] = TRUE
    }
    else if (away_stop == 0.625){
      months[11:12] = TRUE
      months[1:8] = TRUE
    }
    else if (away_stop == 0.708){
      months[11:12] = TRUE
      months[1:9] = TRUE
    }
    else if (away_stop == 0.792){
      months[11:12] = TRUE
      months[1:10] = TRUE
    }
  }
  ## dec start
  else if(away_start == 0.958) {
    if (away_stop == 0.042){
      months[12] = TRUE
      months[1] = TRUE
    }
    else if (away_stop == 0.125){
      months[12] = TRUE
      months[1:2] = TRUE
    }
    else if (away_stop == 0.208){
      months[12] = TRUE
      months[1:3] = TRUE
    }
    else if (away_stop == 0.292){
      months[12] = TRUE
      months[1:4] = TRUE
    }
    else if (away_stop == 0.375){
      months[12] = TRUE
      months[1:5] = TRUE
    }
    else if (away_stop == 0.458){
      months[12] = TRUE
      months[1:6] = TRUE
    }
    else if (away_stop == 0.542){
      months[12] = TRUE
      months[1:7] = TRUE
    }
    else if (away_stop == 0.625){
      months[12] = TRUE
      months[1:8] = TRUE
    }
    else if (away_stop == 0.708){
      months[12] = TRUE
      months[1:9] = TRUE
    }
    else if (away_stop == 0.792){
      months[12] = TRUE
      months[1:10] = TRUE
    }
    else if (away_stop == 0.875){
      months[12] = TRUE
      months[1:11] = TRUE
    }
  }
  return(months)
}

