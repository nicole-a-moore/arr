## script used to run sliding window analysis 
library(gdata)
library(evobiR)
library(tidyverse)


cadillac <- read.csv("./data-raw/intratherm-may-2020-squeaky-clean.csv")
temp_data <- read.csv("./data-processed/arr_terrestrial_temp_data.csv")

## filter out rows of data we cannot use for sliding window/ARR analysis: no location data, no acclimation, not tmax, no lifespan 

cadillac <- cadillac %>%
  subset(subset = !is.na(latitude)) %>%
  subset(subset = !is.na(longitude)) %>%
  subset(subset = !is.na(acclim_temp)) %>%
  subset(subset = !is.na(lifespan_days)) %>%
  subset(subset = lifespan_days != "unk") %>%
  subset(subset = lifespan_days != "kunk") %>%
  subset(subset = parameter_tmax_or_tmin == "tmax") 

cadillac <- drop.levels(cadillac)

## convert season_inactive and season_when_away_100km to start and end date numbers
##FIRST CLEAN 
## season_when_away_100km._start, season_when_away_100km._stop, season_inactive_start, season_inactive_stop



## make new dataframe with only each unique species + lat_long pair, remove columns unnecessary for later merging and analysis 
unique_pairs_unsubsetted <- cadillac[!duplicated(cadillac$population_id),]
unique_pairs <- unique_pairs_unsubsetted %>%
  select(c(genus_species, population_id, lifespan_days, season_when_away_100km._start, 
           season_when_away_100km._stop, season_inactive_start, 
           season_inactive_stop, season_inactive_start2,
           season_inactive_stop2, latitude, longitude, realm_general2,  maximum_body_size_svl_hbl_cm))

## merge so unique_pairs contains column matrix_row_index so the location can be linked to the temperature data file 
unique_pairs <- left_join(unique_pairs, study_locs, by = "lat_long")


experienced_var_mean_col <- c()
experienced_var_max_col <- c()

## for each species, initialize months when away or inactive to TRUE
num_unique <- 1
while (num_unique < 631) {
  ## mark months that species is away as true 
  months <- c(jan <- FALSE, feb <- FALSE, mar <- FALSE, apr <- FALSE, may <- FALSE, jun <- FALSE, jul <- FALSE, 
              aug <- FALSE, sep <- FALSE, oct <- FALSE, nov <- FALSE, dec <- FALSE)
  
  if (!is.na(unique_pairs$season_when_away_100km._start[num_unique])) {
    away_start <- unique_pairs$season_when_away_100km._start[num_unique]
    away_stop <- unique_pairs$season_when_away_100km._stop[num_unique]
    
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
  
  loc <- data.frame(matrix(nrow = 51042))
  loc[,1] <- temp_data[,1]
  loc[,2] <- temp_data[,num_unique + 1]
  colnames(loc) <- c("date", "temp_value")
  
  
  ## iterate through dates in loc by each month, setting temperature values to NA if month = TRUE (if species is away/inactive)
  ## for each year from 1775-2020 - so 245 iterations of 12 
  year_index = 1
  month_index = 1
  
  if (sum(is.na(loc$temp_value)) > 3000) {
    LSV_column_mean <- append(LSV_column_mean, NA, after = length(LSV_column_mean))
    LSV_column_max <- append(LSV_column_max, NA, after = length(LSV_column_max))
  }
  else {
    while (year_index < 246) {
      if (isTRUE(months[1])) {
        loc[month_index, 3] <- NA
      }
      if (isTRUE(months[2])) {
        loc[month_index + 1, 3] <- NA
      }
      if (isTRUE(months[3])) {
        loc[month_index + 2, 3] <- NA
      }
      if (isTRUE(months[4])) {
        loc[month_index + 3, 3] <- NA
      }
      if (isTRUE(months[5])) {
        loc[month_index + 4, 3] <- NA
      }
      if (isTRUE(months[6])) {
        loc[month_index + 5, 3] <- NA
      }
      if (isTRUE(months[7])) {
        loc[month_index + 6, 3] <- NA
      }
      if (isTRUE(months[8])) {
        loc[month_index + 7, 3] <- NA
      }
      if (isTRUE(months[9])) {
        loc[month_index + 8, 3] <- NA
      }
      if (isTRUE(months[10])) {
        loc[month_index + 9, 3] <- NA
      }
      if (isTRUE(months[11])) {
        loc[month_index + 10, 3] <- NA
      }
      if (isTRUE(months[12])) {
        loc[month_index + 11, 3] <- NA
      }
      year_index <- year_index + 1
      month_index <- month_index + 12
    }
    
    ## perform sliding window for location and store LSV in vector column
    ## convert lifespan to months from days (/30) so iteration works properly over month temp data 
    sd_vector <- SlidingWindow(FUN = sd_special, loc$temp_value, as.integer(unique_pairs$lifespan_days[num_unique]/30), 1)
    experienced_var_mean <- mean(sd_vector, na.rm = TRUE)
    experienced_var_max <- max(sd_vector, na.rm = TRUE)
    experienced_var_mean_col <- append(experienced_var_mean_col, experienced_var_mean, after = length(experienced_var_mean_col))
    experienced_var_max_col <- append(experienced_var_max_col, experienced_var_max, after = length(experienced_var_max_col))
  }
  
  num_unique <- num_unique + 1
  
}

unique_pairs$experienced_var_mean <- experienced_var_mean_col
unique_pairs$experienced_var_max <- experienced_var_max_col

## merge back 

## write to file
write.csv(write.csv(unique_pairs, "./data-processed/experienced_var.csv", row.names = FALSE))


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

