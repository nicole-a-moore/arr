## script used to run sliding window analysis
library(gdata)
library(evobiR)
library(tidyverse)

intratherm <- read.csv("./data-raw/intratherm-may-2020-squeaky-clean.txt")

## filter out rows of data we cannot use for tmax, sliding window/ARR analysis: no location data, no acclimation, tmin rows, no lifespan 
intratherm <- intratherm %>%
  filter(!is.na(latitude)) %>%
  filter(!is.na(longitude)) %>%
  filter(!is.na(acclim_temp)) %>%
  filter(!is.na(lifespan_days)) %>%
  filter(lifespan_days != "unk"|| "kunk") %>%
  filter(parameter_tmax_or_tmin == "tmax") 
intratherm <- drop.levels(intratherm)

## convert lifespan to numeric
intratherm$lifespan_days <- as.numeric(as.character(intratherm$lifespan_days))
  
## convert seasons when away and inactive 
intratherm <- convert_seasons_to_numeric(intratherm)

## create pop_id in the same form as temp data 
intratherm <- intratherm %>%
  mutate(population_id = as.character(paste(genus_species, latitude, longitude, sep = "_")))

## write to file:
write.csv(intratherm, "./data-processed/intratherm_sliding-window-ready.csv", row.names = FALSE)

## begin sliding window 
temp_data <- read_csv("./data-processed/arr_temp-data.csv")

cadillac <- intratherm 

## subset:
cadillac <- subset(cadillac, select = c(genus_species, population_id, acclim_temp, latitude, longitude, 
                                    genus, species, parameter_value, parameter_tmax_or_tmin, 
                                    realm_general2, lifespan_days,
                                    season_when_away_100km_start, season_when_away_100km_stop, 
                                    season_inactive_start,  season_inactive_stop, season_inactive_start2,
                                    season_inactive_stop2, maximum_body_size_svl_hbl_cm, elevation_of_collection))

## get rid of duplicate population rows and any marine data 
unique_pairs <- subset(cadillac, !duplicated(cadillac$population_id)) %>%
  filter(realm_general2 != "Marine")

sd_cumulative <- c()
experienced_var_mean <- c()
experienced_var_max <- c()

## for each species:
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
  
  loc <- data.frame(matrix(nrow = 51042))
  loc[,1] <- temp_data$date
  loc[,2] <-  temp_data[,which(colnames(temp_data) == pop_id)]
  colnames(loc) <- c("date", "temp_value")
    
  print("Setting temps when away or inactive to NA...")
  ## set temps to NA:
  loc <- set_temps_to_NA(max = 2020, rep = 1880, day_index = 1, loc = loc, months = months)
  
  print("Sliding the window...")
  ## perform sliding window calculations:
  lifespan <- as.integer(unique_pairs$lifespan_days[num_unique])
  
  if(!is.na(lifespan)) {
    if(lifespan > nrow(loc)) {
      lifespan <- as.integer(nrow(loc)/365)
    }
    sd_vector <- SlidingWindow(FUN = sd_special, loc$temp_value, lifespan, 1)
    print(paste("Performed sliding window for ", num_unique, sep = ""))
    print(paste("Number of days iterating over: ", length(loc$temp_value),sep = ""))
    ## calculate mean sd and add to experienced variation vectors 
    pop_experienced_var_mean<- mean(sd_vector, na.rm = TRUE)
    pop_experienced_var_max<- max(sd_vector, na.rm = TRUE)
  }
  else {
    pop_experienced_var_mean = NA
    pop_experienced_var_max = NA
  }
  sd = sd_special(loc$temp_value)
  
  experienced_var_mean <- append(experienced_var_mean, pop_experienced_var_mean, 
                                 after = length(experienced_var_mean))
  experienced_var_max<- append(experienced_var_max, pop_experienced_var_max, 
                               after = length(experienced_var_max))
  sd_cumulative <- append(sd_cumulative, sd, after = length(sd_cumulative))
  
  print("Done! Moving on to next population :-)")
  ##move on to next species
  num_unique <- num_unique + 1
}


experienced_var_mean <- ifelse((is.infinite(experienced_var_mean) | is.nan(experienced_var_mean)), "NA",
                               experienced_var_mean)
experienced_var_max <- ifelse((is.infinite(experienced_var_max) | is.nan(experienced_var_max)), "NA",
                               experienced_var_max)
sd_cumulative <- ifelse((is.infinite(sd_cumulative) | is.nan(sd_cumulative)), "NA",
                        sd_cumulative)

unique_pairs$experienced_var_mean <- experienced_var_mean
unique_pairs$experienced_var_max <- experienced_var_max
unique_pairs$sd_cumulative <- sd_cumulative

precious_sw_output <- unique_pairs

unique_pairs <- readRDS("./data-processed/precious_sw_output.rds")

## merge back to full intratherm database:

unique_pairs <- unique_pairs %>%
  select(-parameter_value, -parameter_tmax_or_tmin,-acclim_temp)

merged <- left_join(intratherm, unique_pairs)



## write to file
write.csv(merged, "./data-processed/arr_sliding-window-output.csv", row.names = FALSE)


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

#
###
#######
##############
##########################
############################################
#######################################################################
# FUNCTION FOR CONVERTING SEASONS TO START/STOP NUMERIC DATES #########
#######################################################################
convert_seasons_to_numeric <- function (data) {
  
  data <- data %>%
    filter(!is.na(latitude)) %>%
    filter(!is.na(longitude))
  
  ## if in southern hemisphere, mark:
  data <- data %>%
    mutate(is_in_south = ifelse(latitude < 0, "Y", "N"))
  
  
  ## convert season_inactive and season_when_away_100km to start and end date numbers
  swa <- data$season_when_away_100km
  
  if (length(swa) < 1) {
    print("It seems like something is wrong with the column 'season_when_away_100km'")
    return()
  }
  
  split <- str_split_fixed(swa, pattern = "-", n = 2)
  split <- as.data.frame(split)
  split$V1 <- as.character(split$V1)
  split$V2 <- as.character(split$V2)
  split$is_in_south <- data$is_in_south
  
  
  ## if is in south is true, flip the seasons:
  i = 1
  while (i < length(split$V1) + 1) {
    if(split$is_in_south[i] == 'Y') {
      if(split$V1[i] == "Summer") {
        split$V1[i] = "Winter"
      }
      else if(split$V1[i] == "Spring") {
        split$V1[i] = "Fall"
      }
      else if(split$V1[i] == "Winter") {
        split$V1[i] = "Summer"
      }
      else if(split$V1[i] == "Fall") {
        split$V1[i] = "Spring"
      }
      if (split$V2[i] == "Summer") {
        split$V2[i] = "Winter"
      }
      else if (split$V2[i] == "Spring") {
        split$V2[i] = "Fall"
      }
      else if (split$V2[i] == "Winter") {
        split$V2[i] = "Summer"
      }
      else if (split$V2[i] == "Fall") {
        split$V2[i] = "Spring"
      }
    }
    i = i + 1
  }
  
  split <- split %>% ## if only one month/season, put value in both columns
    mutate(V2 = if_else(V2=="", V1 , V2))
  
  season_when_away_100km_start <- split$V1 %>%
    str_replace(pattern = "Winter", replacement = "0.875") %>%
    str_replace(pattern = "Spring", replacement = "0.208") %>%
    str_replace(pattern = "Summer", replacement = "0.458") %>%
    str_replace(pattern = "Fall", replacement = "0.708") %>%
    str_replace(pattern = "Jan", replacement = "0.042") %>%
    str_replace(pattern = "Feb", replacement = "0.125") %>%
    str_replace(pattern = "Mar", replacement = "0.208") %>%
    str_replace(pattern = "Apr", replacement = "0.292") %>%
    str_replace(pattern = "May", replacement = "0.375") %>%
    str_replace(pattern = "Jun", replacement = "0.458") %>%
    str_replace(pattern = "Jul", replacement = "0.542") %>%
    str_replace(pattern = "Aug", replacement = "0.625") %>%
    str_replace(pattern = "Sep", replacement = "0.708") %>%
    str_replace(pattern = "Oct", replacement = "0.792") %>%
    str_replace(pattern = "Nov", replacement = "0.875") %>%
    str_replace(pattern = "Dec", replacement = "0.958") %>%
    str_replace(pattern = "kunk", replacement = "") %>%
    str_replace(pattern = "unk", replacement = "")
  
  season_when_away_100km_stop <- split$V2 %>%
    str_replace(pattern = "Winter", replacement = "0.125") %>%
    str_replace(pattern = "Spring", replacement = "0.375") %>%
    str_replace(pattern = "Summer", replacement = "0.625") %>%
    str_replace(pattern = "Fall", replacement = "0.792") %>%
    str_replace(pattern = "Jan", replacement = "0.042") %>%
    str_replace(pattern = "Feb", replacement = "0.125") %>%
    str_replace(pattern = "Mar", replacement = "0.208") %>%
    str_replace(pattern = "Apr", replacement = "0.292") %>%
    str_replace(pattern = "May", replacement = "0.375") %>%
    str_replace(pattern = "Jun", replacement = "0.458") %>%
    str_replace(pattern = "Jul", replacement = "0.542") %>%
    str_replace(pattern = "Aug", replacement = "0.625") %>%
    str_replace(pattern = "Sep", replacement = "0.708") %>%
    str_replace(pattern = "Oct", replacement = "0.792") %>%
    str_replace(pattern = "Nov", replacement = "0.875") %>%
    str_replace(pattern = "Dec", replacement = "0.958") %>%
    str_replace(pattern = "kunk", replacement = "") %>%
    str_replace(pattern = "unk", replacement = "")
  
  sia <- data$season_inactive
  
  if (length(sia) < 1) {
    print("It seems like something is wrong with the column 'season_inactive'")
    return()
  }
  
  split <- str_split_fixed(sia, pattern = "-", n = 2)
  split <- as.data.frame(split)
  split$V1 <- as.character(split$V1)
  split$V2 <- as.character(split$V2)
  split$is_in_south <- data$is_in_south
  
  ## set ones with multiple sia to first 
  split$V1[which(sia == "Summer and Winter")] = "Summer"
  split$V1[which(sia == "Spring and Winter")] = "Spring"
  split$V1[which(sia == "Fall and Winter")] = "Fall"
  
  ## if is in south is true, flip the seasons:
  i = 1
  while (i < length(split$V1) + 1) {
    if(split$is_in_south[i] == 'Y') {
      if(split$V1[i] == "Summer") {
        split$V1[i] = "Winter"
      }
      else if(split$V1[i] == "Spring") {
        split$V1[i] = "Fall"
      }
      else if(split$V1[i] == "Winter") {
        split$V1[i] = "Summer"
      }
      else if(split$V1[i] == "Fall") {
        split$V1[i] = "Spring"
      }
      if (split$V2[i] == "Summer") {
        split$V2[i] = "Winter"
      }
      else if (split$V2[i] == "Spring") {
        split$V2[i] = "Fall"
      }
      else if (split$V2[i] == "Winter") {
        split$V2[i] = "Summer"
      }
      else if (split$V2[i] == "Fall") {
        split$V2[i] = "Spring"
      }
    }
    i = i + 1
  }
  
  split <- split %>% ## if only one month/season, put value in both columns
    mutate(V2 = if_else(V2=="", V1 , V2))
  
  season_inactive_start <- split$V1 %>%
    str_replace(pattern = "Winter", replacement = "0.875") %>%
    str_replace(pattern = "Spring", replacement = "0.208") %>%
    str_replace(pattern = "Summer", replacement = "0.458") %>%
    str_replace(pattern = "Fall", replacement = "0.708") %>%
    str_replace(pattern = "Jan", replacement = "0.042") %>%
    str_replace(pattern = "Feb", replacement = "0.125") %>%
    str_replace(pattern = "Mar", replacement = "0.208") %>%
    str_replace(pattern = "Apr", replacement = "0.292") %>%
    str_replace(pattern = "May", replacement = "0.375") %>%
    str_replace(pattern = "Jun", replacement = "0.458") %>%
    str_replace(pattern = "Jul", replacement = "0.542") %>%
    str_replace(pattern = "Aug", replacement = "0.625") %>%
    str_replace(pattern = "Sep", replacement = "0.708") %>%
    str_replace(pattern = "Oct", replacement = "0.792") %>%
    str_replace(pattern = "Nov", replacement = "0.875") %>%
    str_replace(pattern = "Dec", replacement = "0.958") %>%
    str_replace(pattern = "kunk", replacement = "") %>%
    str_replace(pattern = "unk", replacement = "") %>%
    str_replace(pattern = "none", replacement = "")
  
  season_inactive_stop <- split$V2 %>%
    str_replace(pattern = "Winter", replacement = "0.125") %>%
    str_replace(pattern = "Spring", replacement = "0.375") %>%
    str_replace(pattern = "Summer", replacement = "0.625") %>%
    str_replace(pattern = "Fall", replacement = "0.792") %>%
    str_replace(pattern = "Jan", replacement = "0.042") %>%
    str_replace(pattern = "Feb", replacement = "0.125") %>%
    str_replace(pattern = "Mar", replacement = "0.208") %>%
    str_replace(pattern = "Apr", replacement = "0.292") %>%
    str_replace(pattern = "May", replacement = "0.375") %>%
    str_replace(pattern = "Jun", replacement = "0.458") %>%
    str_replace(pattern = "Jul", replacement = "0.542") %>%
    str_replace(pattern = "Aug", replacement = "0.625") %>%
    str_replace(pattern = "Sep", replacement = "0.708") %>%
    str_replace(pattern = "Oct", replacement = "0.792") %>%
    str_replace(pattern = "Nov", replacement = "0.875") %>%
    str_replace(pattern = "Dec", replacement = "0.958") %>%
    str_replace(pattern = "kunk", replacement = "") %>%
    str_replace(pattern = "unk", replacement = "") %>%
    str_replace(pattern = "none", replacement = "")
  
  ## set second season when inactive:
  split$V1 <- ""
  split$V2 <- ""
  split$V1[which(sia == "Summer and Winter")] = "Winter"
  split$V1[which(sia == "Spring and Winter")] = "Winter"
  split$V1[which(sia == "Fall and Winter")] = "Winter"
  
  ## if is in south is true, flip the seasons:
  i = 1
  while (i < length(split$V1) + 1) {
    if(split$is_in_south[i] == 'Y') {
      if(split$V1[i] == "Summer") {
        split$V1[i] = "Winter"
      }
      else if(split$V1[i] == "Spring") {
        split$V1[i] = "Fall"
      }
      else if(split$V1[i] == "Winter") {
        split$V1[i] = "Summer"
      }
      else if(split$V1[i] == "Fall") {
        split$V1[i] = "Spring"
      }
      if (split$V2[i] == "Summer") {
        split$V2[i] = "Winter"
      }
      else if (split$V2[i] == "Spring") {
        split$V2[i] = "Fall"
      }
      else if (split$V2[i] == "Winter") {
        split$V2[i] = "Summer"
      }
      else if (split$V2[i] == "Fall") {
        split$V2[i] = "Spring"
      }
    }
    i = i + 1
  }
  
  split <- split %>% ## if only one month/season, put value in both columns
    mutate(V2 = if_else(V2=="", V1 , V2))
  
  season_inactive_start2 <- split$V1 %>%
    str_replace(pattern = "Winter", replacement = "0.875")
  
  season_inactive_stop2 <- split$V2 %>%
    str_replace(pattern = "Winter", replacement = "0.125")
  
  data$season_when_away_100km_start <- season_when_away_100km_start
  data$season_when_away_100km_stop <- season_when_away_100km_stop
  data$season_inactive_start <- season_inactive_start
  data$season_inactive_stop <- season_inactive_stop
  data$season_inactive_start2 <- season_inactive_start2
  data$season_inactive_stop2 <- season_inactive_stop2
  
  return(data)
} 

#
###
#######
##############
##########################
############################################
################################################################
# FUNCTION TO SET TEMPS WHEN AWAY TO NA AT A LOCATION  #########
################################################################
set_temps_to_NA <- function (rep, max, day_index, loc, months) {
  ## max = last year in temp data + 1
  ## rep = first year
  ## day_index = 1

  ## iterate through dates by 365 days except on leap year, set temp_val to NA if months = true 
  ## for each year from 1880-2020, meaning 140 iterations of 365/366
   while (rep < max) {
    if (rep == 2019) {
      if (isTRUE(months[1])) {
        loc$temp_value[day_index:(day_index + 30)] <- NA
      }
      if (isTRUE(months[2])) {
        loc$temp_value[(day_index + 30):(day_index + 58)] <- NA
      }
      if (isTRUE(months[3])) {
        loc$temp_value[(day_index + 58):(day_index + 89)] <- NA
      }
      if (isTRUE(months[4])) {
        loc$temp_value[(day_index + 89):(day_index + 119)] <- NA
      }
      if (isTRUE(months[5])) {
        loc$temp_value[(day_index + 119):(day_index + 150)] <- NA
      }
      if (isTRUE(months[6])) {
        loc$temp_value[(day_index + 150):(day_index + 180)] <- NA
      }
      if (isTRUE(months[7])) {
        loc$temp_value[(day_index + 180):(day_index + 211)] <- NA
      }
      if (isTRUE(months[8])) {
        loc$temp_value[(day_index + 211):(day_index + 242)] <- NA
      }
      if (isTRUE(months[9])) {
        loc$temp_value[(day_index + 242):(day_index + 272)] <- NA
      }
      rep <- rep + 1
      day_index <- day_index + 365
    }
    else if (rep == 1900) {
      if (isTRUE(months[1])) {
        loc$temp_value[day_index:(day_index + 30)] <- NA
      }
      if (isTRUE(months[2])) {
        loc$temp_value[(day_index + 30):(day_index + 58)] <- NA
      }
      if (isTRUE(months[3])) {
        loc$temp_value[(day_index + 58):(day_index + 89)] <- NA
      }
      if (isTRUE(months[4])) {
        loc$temp_value[(day_index + 89):(day_index + 119)] <- NA
      }
      if (isTRUE(months[5])) {
        loc$temp_value[(day_index + 119):(day_index + 150)] <- NA
      }
      if (isTRUE(months[6])) {
        loc$temp_value[(day_index + 150):(day_index + 180)] <- NA
      }
      if (isTRUE(months[7])) {
        loc$temp_value[(day_index + 180):(day_index + 211)] <- NA
      }
      if (isTRUE(months[8])) {
        loc$temp_value[(day_index + 211):(day_index + 242)] <- NA
      }
      if (isTRUE(months[9])) {
        loc$temp_value[(day_index + 242):(day_index + 272)] <- NA
      }
      if (isTRUE(months[10])) {
        loc$temp_value[(day_index + 272):(day_index + 303)] <- NA
      }
      if (isTRUE(months[11])) {
        loc$temp_value[(day_index + 303):(day_index + 333)] <- NA
      }
      if (isTRUE(months[12])) {
        loc$temp_value[(day_index + 333):(day_index + 364)] <- NA
      }
      rep <- rep + 1
      day_index <- day_index + 365
    }
    else if (rep %% 4 == 0) {
      if (isTRUE(months[1])) {
        loc$temp_value[day_index:(day_index + 30)] <- NA
      }
      if (isTRUE(months[2])) {
        loc$temp_value[(day_index + 30):(day_index + 59)] <- NA
      }
      if (isTRUE(months[3])) {
        loc$temp_value[(day_index + 59):(day_index + 90)] <- NA
      }
      if (isTRUE(months[4])) {
        loc$temp_value[(day_index + 90):(day_index + 120)] <- NA
      }
      if (isTRUE(months[5])) {
        loc$temp_value[(day_index + 120):(day_index + 151)] <- NA
      }
      if (isTRUE(months[6])) {
        loc$temp_value[(day_index + 151):(day_index + 181)] <- NA
      }
      if (isTRUE(months[7])) {
        loc$temp_value[(day_index + 181):(day_index + 212)] <- NA
      }
      if (isTRUE(months[8])) {
        loc$temp_value[(day_index + 212):(day_index + 243)] <- NA
      }
      if (isTRUE(months[9])) {
        loc$temp_value[(day_index + 243):(day_index + 273)] <- NA
      }
      if (isTRUE(months[10])) {
        loc$temp_value[(day_index + 273):(day_index + 304)] <- NA
      }
      if (isTRUE(months[11])) {
        loc$temp_value[(day_index + 304):(day_index + 334)] <- NA
      }
      if (isTRUE(months[12])) {
        loc$temp_value[(day_index + 334):(day_index + 365)] <- NA
      }
      rep <- rep + 1
      day_index <- day_index + 366
    }
    else {
      if (isTRUE(months[1])) {
        loc$temp_value[day_index:(day_index + 30)] <- NA
      }
      if (isTRUE(months[2])) {
        loc$temp_value[(day_index + 30):(day_index + 58)] <- NA
      }
      if (isTRUE(months[3])) {
        loc$temp_value[(day_index + 58):(day_index + 89)] <- NA
      }
      if (isTRUE(months[4])) {
        loc$temp_value[(day_index + 89):(day_index + 119)] <- NA
      }
      if (isTRUE(months[5])) {
        loc$temp_value[(day_index + 119):(day_index + 150)] <- NA
      }
      if (isTRUE(months[6])) {
        loc$temp_value[(day_index + 150):(day_index + 180)] <- NA
      }
      if (isTRUE(months[7])) {
        loc$temp_value[(day_index + 180):(day_index + 211)] <- NA
      }
      if (isTRUE(months[8])) {
        loc$temp_value[(day_index + 211):(day_index + 242)] <- NA
      }
      if (isTRUE(months[9])) {
        loc$temp_value[(day_index + 242):(day_index + 272)] <- NA
      }
      if (isTRUE(months[10])) {
        loc$temp_value[(day_index + 272):(day_index + 303)] <- NA
      }
      if (isTRUE(months[11])) {
        loc$temp_value[(day_index + 303):(day_index + 333)] <- NA
      }
      if (isTRUE(months[12])) {
        loc$temp_value[(day_index + 333):(day_index + 364)] <- NA
      }
      rep <- rep + 1
      day_index <- day_index + 365
    }
  }
  
  return(loc)
}


