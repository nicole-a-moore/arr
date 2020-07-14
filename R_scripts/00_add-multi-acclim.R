## adding data from Globtherm that was excluded since only 1 population, multi acclim temps
## can be used for this analysis
library(tidyverse)
library(taxize)

intratherm <- read.csv("data-raw/intratherm-may-2020-squeaky-clean.csv", stringsAsFactors = FALSE) %>%
  mutate(population_id = paste(genus_species, latitude, elevation_of_collection, longitude, sep = "_")) 
multiacclim <- read.csv("data-raw/intratherm-multi-acclim.csv", stringsAsFactors = FALSE) 


## taxize old version of data base so species match 
###################################################
taxa <- data.frame(taxa = multiacclim$genus_species) ## create dataframe of names to check

syns <- unique(taxa)
tsn_search <- get_tsn(as.character(syns$taxa), accepted = FALSE) ## find tsn for each unique taxa
tsn_search <- readRDS("data-processed/tsn-search_multiacclim.rds")
tsns <- data.frame(tsn_search)
tsns$taxa <- syns$taxa
syns <- tsns

found <- syns %>%
  subset(match == "found") 

report <- lapply(found$ids, itis_acceptname)
report_df <- data.frame(matrix(unlist(report), nrow=247, byrow=T),stringsAsFactors=FALSE)
report_df <- readRDS("data-processed/report-df_multiacclim.rds")

found <- found %>%
  mutate(genus_species_corrected = report_df$X2)

## merge short unique list to long list of all taxa
merged_unique <- left_join(syns, found)
merged <- left_join(taxa, merged_unique) %>%
  mutate(taxa = as.character(taxa)) %>%
  mutate(genus_species_corrected = as.character(genus_species_corrected))

## if names found are not accepted names, then change to accepted name
i = 1
while (i < length(merged$taxa)) {
  if (!is.na(merged$genus_species_corrected[i])) {
    merged$taxa[i] <- as.character(merged$genus_species_corrected[i])
  }
  i = i+1
}

## create new genus and species columns, correct original dataset
split <- str_split_fixed(merged$taxa, pattern = " ", n = 2)
merged <- merged %>%
  mutate(genus = split[,1]) %>%
  mutate(species = split[,2])
  
multiacclim <- multiacclim %>%
  ungroup()%>%
  mutate(genus = merged$genus) %>%
  mutate(species = merged$species) %>%
  mutate(genus_species = merged$taxa) 

## change ones missing extra i:
## "Limnodynastes peroni"
## "Pseudophryne bibroni"
## change Takydromus hsuehshaner to Takydromus hsuehshanensis
## change Eleutherodactylus	portoricensus to Eleutherodactylus portoricensis
multiacclim$genus_species[which(multiacclim$genus_species == 
                                  "Limnodynastes peroni")] <- "Limnodynastes peronii"
multiacclim$genus[which(multiacclim$genus_species == "Limnodynastes peroni")] <- "Limnodynastes"
multiacclim$species[which(multiacclim$genus_species == "Limnodynastes peroni")] <- "peronii"

multiacclim$genus_species[which(multiacclim$genus_species == "Pseudophryne bibroni")] <- "Pseudophryne bibronii"
multiacclim$genus[which(multiacclim$genus_species == "Pseudophryne bibroni")] <- "Pseudophryne"
multiacclim$species[which(multiacclim$genus_species == "Pseudophryne bibroni")] <- "bibronii"

multiacclim$genus_species[which(multiacclim$genus_species == "Takydromus hsuehshaner")] <- "Takydromus hsuehshanensis"
multiacclim$genus[which(multiacclim$genus_species == "Takydromus hsuehshaner")] <- "Takydromus"
multiacclim$species[which(multiacclim$genus_species == "Takydromus hsuehshaner")] <- "hsuehshanensis"

multiacclim$genus_species[which(multiacclim$genus_species == "Eleutherodactylus portoricensus")] <- "Eleutherodactylus portoricensis"
multiacclim$genus[which(multiacclim$genus_species == "Eleutherodactylus portoricensus")] <- "Eleutherodactylus"
multiacclim$species[which(multiacclim$genus_species == "Eleutherodactylus portoricensus")] <- "portoricensis"

## get rid of duplicated rows, tmins
## remove rows of populations that are in intratherm
## get rid of non-existent Hyla alpina
## get rid of Cyprinodon sp.1
multiacclim <- multiacclim %>%
  unique() %>%
  filter(parameter_tmax_or_tmin == "tmax") %>%
  mutate(population_id = paste(genus_species, latitude, elevation_of_collection, longitude, sep = "_")) %>%
  filter(!genus_species %in% intratherm$genus_species) %>%
  filter(!genus_species == "Hyla alpina") %>%
  filter(!genus_species == "Cyprinodon sp.1") 


## see which species we do not have traits for 
traits <- read.csv("data-raw/intratherm-traits-clean-citations.csv")
species <- unique(multiacclim$genus_species)

havetraits <- species[which(species %in% paste(traits$genus, traits$species, sep = " "))]
newspecies <- species[which(!species %in% paste(traits$genus, traits$species, sep = " "))]

to_add <- multiacclim %>%
  filter(genus_species %in% newspecies) %>%
  filter(!duplicated(genus_species))

traits_to_add <- data.frame(matrix(ncol = 30, nrow = 40)) 
colnames(traits_to_add) <- colnames(traits)

traits_to_add <- traits_to_add %>%
  mutate(genus = to_add$genus, species = to_add$species)
  

## write rows in traits database for new species and write to file to fill in 
traits <- rbind(traits, traits_to_add)

write.csv(traits, "data-processed/intratherm-traits_with-new-species.csv", row.names = FALSE)




## filled in only necessary traits for new species, if lifespan unavailable skipped others (lifespan, season when away, etc.)

traits <- read.csv("data-raw/intratherm-traits_with-new-species_complete.csv")



## check higher taxa:
check <- multiacclim %>%
  select(genus_species, family, class, order, phylum) %>% View

## 	Chelon subviridis not accepted name according to world register of marine species
## change:
chelon <- which(multiacclim$genus_species == "Chelon subviridis")
multiacclim$genus_species[chelon] <- "Planiliza subviridis"
multiacclim$order[chelon] <- "Mugiliformes"
multiacclim$class[chelon] <- "Actinopterygii"
multiacclim$phylum[chelon] <- "Chordata"
multiacclim$genus[chelon] <- "Planiliza"
multiacclim$species[chelon] <- "subviridis"


## merge 
multiacclim <- multiacclim %>%
  rename(life_stage.x = life_stage) %>%
  mutate(elevation_of_collection = as.character(elevation_of_collection)) %>%
  mutate(acclim_temp = as.character(acclim_temp)) %>%
  mutate(ramping_rate = as.character(ramping_rate)) 
  
intratherm <- intratherm %>%
  mutate(elevation_of_collection = as.character(elevation_of_collection)) %>%
  mutate(acclim_temp = as.character(acclim_temp)) %>%
  mutate(ramping_rate = as.character(ramping_rate)) 
    

intratherm <- full_join(intratherm, multiacclim)


## get rid of unnecessary columns: 
intratherm <- intratherm %>%
  select(-c(
    "genus_from_study",
    "species_from_study",                                                                    
    "record_number",                                                                         
    "record",                                                                                
    "acclim_time_original",                                                                  
    "acclim_time_units",                                                                     
    "dispersal_type_walking_may_need_to_be_reworded_to_something_that_encaptures_slithering",
    "logic_source_for_dispersal_distance2",                                                  
    "age_maturity_days_female",                                                              
    "age_maturity_days_female_reference",                                                    
    "age_maturity_days_male",                                                                
    "age_maturity_days_male_reference",                                                      
    "age_maturity.days_unknown_sex",                                                        
    "dispersal_distance_category",                                                   
    "logic_source_for_dispersal_distance",                                                  
    "maximum_body_size_svl_hbl_cm",                                                          
    "source_for_maximum_body_size",                                                          
    "average_body_size_female",                                                             
    "average_body_size_female_reference",                                                    
    "average_body_size_male",                                                              
    "average_body_size_male_reference",                                                      
    "average_body_size_unknown_sex",                                                         
    "season_inactive",                                                                       
    "logic_source_for_season_inactive",                                                      
    "migratory",                                                                        
    "source_logic_for_migration_info",                                                       
    "season_when_away_10km",                                                           
    "season_when_away_100km",                                                           
    "logic_source_season_when_away",                                                     
    "Home.range.size_Km2",                                                                   
    "logic.source_for_Home.range.size",                                                      
    "data_gatherer",                                                           
    "lifespan_days",                                                                        
    "lifespan_days_reference",                                                              
    "is.nocturnal",                                                                          
    "is.nocturnal_source",                                                               
    "realm_general3",                                                                        
    "realm_specific",                                                                        
    "realm_general",                                                                         
    "notes",                                                                                 
    "n_cat",                                                                                 
    "family1" ,                                                                              
    "family2",                                                                             
    "genus2",                                                                                
    "species2",                                                                             
    "family3",                                                                               
    "genus3",                                                                               
    "species3" ,                                                                             
    "n_comment",                                                                            
    "raw_ctm2" ,                                                                            
    "acclim" ,                                                                             
    "acclim_time_life_120" ,                                                                 
    "log_acclim_time"  ,                                                                     
    "safety" ,                                                                       
    "photoperiod"  ,                                                                       
    "ep1"   ,                                                                                
    "ep2" ,                                                                                  
    "elevation",                                                                             
    "log_elevation" ,                                                                        
    "geog_cat" ,                                                                             
    "abs_lat",                                                                             
    "log_abs_lat" ,                                                                          
    "source"      ,                                                                          
    "comments" ,                                                                             
    "year"  ,                                                                             
    "citation" ,                                                                             
    "range",                                                                                 
    "log_range",                                                                             
    "introduced"  ,                                                                       
    "elev_min",                                                                             
    "elev_max",                                                                             
    "svl" ,                                                                             
    "log_svl",                                                                            
    "tmxslope"  ,                                                                           
    "bioclim5" ,                                                                             
    "iucn_2014",                                                                             
    "iucn_2104risk",                                                                         
    "sampling_habitat",                                                                      
    "acclimation",                                                                      
    "rate_acclimation_c_day"  ,                                                           
    "methodology" ,                                                                     
    "length_experiment_min",                                                               
    "body_length_mm",                                                                        
    "weight_g",                                                                              
    "phylogeny",                                                                             
    "reference" )
    )

intratherm <- left_join(intratherm, traits)


## add intratherm id for new ones: 
intratherm$intratherm_id[2791:3037] <- 3030:3276


## clean realm_general2 to get rid of aquatic
aquatic <- intratherm[which(intratherm$realm_general2 == "Aquatic"),]

aquatic$realm_general2[which(aquatic$intratherm_id %in% 3115:3119)] <- "Marine"
aquatic$realm_general2[which(aquatic$intratherm_id %in% 3175:3176)] <- "Terrestrial"
aquatic$realm_general2[which(aquatic$intratherm_id %in% 3182:3183)] <- "Terrestrial"

intratherm <- intratherm %>%
  filter(!intratherm$realm_general2 == "Aquatic") %>%
  rbind(., aquatic)



## clean traits:
## cleaning up season_when_away100km+ and season_when_inactive
## change names of columns to simplify:
colnames(intratherm)[which(names(intratherm) == "season_when_away_100km...migratory.only.")] <- "season_when_away_100km"
colnames(intratherm)[which(names(intratherm) == "season_when_away_10km...migratory.only.")] <- "season_when_away_10km"

## clean up 100km:
swa <- intratherm$season_when_away_100km
unique(swa)

swa <- str_replace(swa, pattern = "October/Novermber/December/January/February", replacement = "Oct-Feb") %>%
  str_replace(pattern = "spring-fall", replacement = "Spring-Fall") %>%
  str_replace(pattern = "fall/winter", replacement = "Fall-Winter") %>%
  str_replace(pattern = "March-April", replacement = "Mar-Apr") %>%
  str_replace(pattern = "August/September/October", replacement = "Aug-Oct") %>%
  str_replace(pattern = "summer", replacement = "Summer") %>%
  str_replace(pattern = "fall", replacement = "Fall") %>%
  str_replace(pattern = "autumn/winter", replacement = "Fall-Winter") %>%
  str_replace(pattern = "April-Aug/Sept", replacement = "Apr-Sep") %>%
  str_replace(pattern = "spring", replacement = "Spring") %>%
  str_replace(pattern = "late Spring/early Summer", replacement = "Spring-Summer") %>%
  str_replace(pattern = "3 months around April/May", replacement = "Apr-Jun") %>%
  str_replace(pattern = "late winter-early Spring", replacement = "Winter-Spring") %>%
  str_replace(pattern = "Fall \\- Winter", replacement = "Fall-Winter") %>%
  str_replace(pattern = "Jan \\- Jun", replacement = "Jan-Jun") 

swa[which(is.na(swa))] <- ""
unique(swa)

intratherm$season_when_away_100km <- swa



## clean up 10km:
swa <- intratherm$season_when_away_10km
unique(swa)

swa <- str_replace(swa, pattern = "spring-fall", replacement = "Spring-Fall") %>%
  str_replace(pattern = "winter(including Nov and Dec)", replacement = "Winter") %>%
  str_replace(pattern = "August/September/October", replacement = "Aug-Oct") %>%
  str_replace(pattern = "March-April", replacement = "Mar-Apr") %>%
  str_replace(pattern = "fall/winter", replacement = "Fall-Winter") %>%
  str_replace(pattern = "late spring/early summer", replacement = "Spring-Summer") %>%
  str_replace(pattern = "October/Novermber/December/January/February", replacement = "Oct-Feb") %>%
  str_replace(pattern = "spring", replacement = "Spring") %>%
  str_replace(pattern = "May-October", replacement = "May-Oct") %>%
  str_replace(pattern = "fall", replacement = "Fall") %>%
  str_replace(pattern = "winter", replacement = "Winter") %>%
  str_replace(pattern = "summer", replacement = "Summer") %>%
  str_replace(pattern = "autumn/winter", replacement = "Fall-Winter") %>%
  str_replace(pattern = "autumn/Winter", replacement = "Fall-Winter") %>%
  str_replace(pattern = "February/March/April", replacement = "Feb-Apr") %>%
  str_replace(pattern = "late winter-early spring", replacement = "Winter-Spring") %>%
  str_replace(pattern = "April-Aug/Sept", replacement = "Apr-Sep") %>%
  str_replace(pattern = "July-September", replacement = "Jul-Sep") %>%
  str_replace(pattern = "August/September/October/November/December/January/February/March/April", replacement = "Aug-Apr") %>%
  str_replace(pattern = "April/May/June", replacement = "Apr-Jun") %>%
  str_replace(pattern = "late Winter-early Spring", replacement = "Winter-Spring") %>%
  str_replace(pattern = "\\([^()]{0,}\\)", replacement = "") %>%
  str_replace(pattern = "3 months around April/May", replacement = "Apr-Jun") %>%
  str_replace(pattern = "Aug-Oct/November/December/January/Feb-Apr", replacement = "Aug-Apr") %>%
  str_replace(pattern = "Spring\\+ two months", replacement = "Spring-Summer") %>%
  str_replace(pattern = "Fall \\- Winter", replacement = "Fall-Winter") %>%
  str_replace(pattern = "Jan \\- Jun", replacement = "Jan-Jun") 

swa[which(is.na(swa))] <- ""
unique(swa)

intratherm$season_when_away_10km <- swa

## clean up season inactive:
## let dry == summer, hot == summer
sia <- intratherm$season_inactive
unique(sia)

sia <- str_replace(sia, pattern = "hot dry", replacement = "Summer") %>%
  str_replace(pattern = "Oct - Mar", replacement = "Oct-Mar") %>%
  str_replace(pattern = "summer/dry season", replacement = "Summer") %>%
  str_replace(pattern = "winter", replacement = "Winter") %>%
  str_replace(pattern = "Oct - Apr", replacement = "Oct-Apr") %>%
  str_replace(pattern = "June - Oct", replacement = "Jun-Oct") %>%
  str_replace(pattern = "fall-winter", replacement = "Fall-Winter") %>%
  str_replace(pattern = "spring\\+ Winter", replacement = "Spring and Winter") %>%
  str_replace(pattern = "hot", replacement = "Summer") %>%
  str_replace(pattern = "dry", replacement = "Summer") %>%
  str_replace(pattern = "summer", replacement = "Summer") %>%
  str_replace(pattern = "winter\\+hot", replacement = "Winter and Summer") %>%
  str_replace(pattern = "winter\\+dry", replacement = "Winter and Summer") %>%
  str_replace(pattern = "fall\\-Winter", replacement = "Fall and Winter") %>%
  str_replace(pattern = "Winter\\+Summer", replacement = "Summer and Winter") %>%
  str_replace(pattern = "none ", replacement = "") %>%
  str_replace(pattern = "none", replacement = "") %>%
  str_replace(pattern = "Winter\\?", replacement = "Winter")%>%
  str_replace(pattern = "Winter ", replacement = "Winter") 

sia[which(is.na(sia))] <- ""
unique(sia)
intratherm$season_inactive <- sia



## clean lifestage.x
unique(intratherm$life_stage.x)
intratherm <- intratherm %>%
  mutate(life_stage.x = ifelse(life_stage.x == "adult", "Adult", life_stage.x)) %>%
  mutate(life_stage.x = ifelse(life_stage.x == "adults", "Adult", life_stage.x)) %>%
  mutate(life_stage.x = ifelse(life_stage.x == "young", "Juvenile", life_stage.x)) %>%
  mutate(life_stage.x = ifelse(life_stage.x == "males_adults", "Adult", life_stage.x)) %>%
  mutate(life_stage.x = ifelse(life_stage.x == "juvenile", "Juvenile", life_stage.x)) %>%
  mutate(life_stage.x = ifelse(life_stage.x == "thalli", "Juvenile", life_stage.x))
  

## get rid of columns that will not be used in arr analysis 
intratherm <- intratherm %>%
  select(-c("Home.range.size_Km2",
            "logic.source_for_Home.range.size",
            "average_body_size_female",
            "average_body_size_male",
            "average_body_size_unknown_sex",
            "average_body_size_female_reference",
            "average_body_size_male_reference",
            "dispersal_distance_category",
            "dispersal_distance2_category",
            "logic_source_for_dispersal_distance2",
            "logic_source_for_dispersal_distance", 
            "age_maturity_days_female",
            "age_maturity_days_female_reference",
            "age_maturity_days_male",
            "age_maturity_days_male_reference",
            "age_maturity.days_unknown_sex")
  )


write.csv(intratherm, "data-processed/intratherm_version-for-arr.csv", row.names = FALSE)

## get temperature data and rerun analysis 

