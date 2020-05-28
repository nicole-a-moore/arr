## calculates ARR for each population
library(tidyverse)
library(broom)
library(lmodel2)
library(janitor)
library(RColorBrewer)
library(viridis)

## read in files
intratherm <- read.csv( "./data-processed/arr_sliding-window-output.csv") ## version of the database with only usable data and experienced variation column 

#make groups of each species for every lat and long
arr <- intratherm %>%
  mutate(acclim_temp = as.numeric(as.character(acclim_temp))) %>%
  group_by(population_id) %>%
  do(tidy(lm(parameter_value ~ acclim_temp, data = .), conf.int = TRUE)) %>% #fit lm to each group, the slope is ARR
  filter(term == "acclim_temp") %>% #extract just the slopes - get rid of intercept 
  dplyr::rename(ARR = estimate) #rename these ARR

## plot the ARR linear models to visualize:
ggplot(data = intratherm, aes(x = acclim_temp, y = parameter_value, col = population_id)) + 
  geom_point(size = 1, alpha = 0.1) +
  geom_errorbar(aes(ymin=parameter_value-error_estimate, ymax=parameter_value+error_estimate), 
                alpha = 0.1, width=0) +
  geom_smooth(method = "lm", se = FALSE, size = 0.4, alpha = 0.1) +
  theme_minimal() + 
  theme(legend.position = "none", panel.border = element_rect(colour = "dimgrey", fill = NA), 
        panel.grid.minor = element_line(colour = "white"), 
        panel.grid.major = element_line(colour = "white"), 
        plot.title = element_text(hjust = 0.5)) + 
  labs(x = "Acclimation temperature (°C)", y = "Upper critical thermal limit (°C)") + 
  scale_color_manual(values = magma(625))

## merge ARR results for each popultion back to intratherm
arr <- arr %>%
  select(-term)

intratherm <- left_join(intratherm, arr, by = "population_id")

## write to csv
write.csv(intratherm, "./data-processed/arr_sliding-window-output-with-arr.csv", row.names = FALSE)









##get rid of duplicated rows since acclimation temp no longer matters, get rid of negative ARRs
LSV_ARR_data <- LSV_ARR_data[!duplicated(LSV_ARR_data[,2]),]
LSV_ARR_data <- subset(LSV_ARR_data, subset = as.numeric(as.character(ARR)) > -0.15)

write.csv(LSV_ARR_data, "~/Documents/SUNDAY LAB/Intratherm/Data sheets/phylo_input.csv")

ggplot(data = cadillac, aes(x = acclim_temp, y = parameter_value, group = population_id, col = LSV)) + 
  geom_point(size = 1, alpha = 0.5) + geom_errorbar(aes(ymin=parameter_value-error_estimate,
                                                        ymax=parameter_value+error_estimate), width=0, alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, size = 0.4, alpha = 0.1) +
  theme_minimal() + 
  theme(panel.border = 
          element_rect(colour = "dimgrey", fill = NA), panel.grid.minor = 
          element_line(colour = "white"), panel.grid.major = 
          element_line(colour = "white"), plot.title = element_text(hjust = 0.5)) + 
  labs(x = "Acclimation temperature (°C)", y = "Upper critical thermal limit (°C)", colour = paste(" Mean", "\n", "experienced", "\n", "variation")) + scale_colour_viridis_c(option = "A")




## linear model of LSV weighted by 1/variance of ARR stope estimate
lm_output_all_mean <- LSV_ARR_data %>%
  drop_na(ARR,LSV_meanSD) %>%
  do(tidy(lm(ARR ~ LSV_meanSD, data = .), conf.int = TRUE))

lm_output_all_top<- LSV_ARR_data %>%
  drop_na(ARR,LSV_meanSD) %>%
  do(tidy(lme(ARR ~ LSV_meanSD, random = ~1|order/family/genus/species, 
              data = data, na.action= NULL), conf.int = TRUE))

ggplot(LSV_ARR_data, aes(x = LSV_meanSD, y = ARR, col = abs(latitude))) + geom_point() + geom_smooth(method = "lm", color = "black") + geom_errorbar(aes(ymin=ARR-std.error, ymax=ARR+std.error), width=0) + labs(x = "Mean lifespan specific variation (°C)", y = "Acclimation response ratio", col = "Absolute latitude (°N/S)" ) + scale_colour_viridis_c(option = "A") + geom_line(y = predict(y = ))




## write to file to make map + taxonomic breakdown from
write.csv(LSV_ARR_data, "~/Documents/SUNDAY LAB/Intratherm/Data sheets/LSV_ARR_output_critical.csv", row.names = FALSE)

## lm with max sd
lm_output_all_max<- LSV_ARR_data %>%
  drop_na(ARR,LSV_max) %>%
  do(tidy(lm(ARR ~ LSV_max, data = .), conf.int = TRUE))

ggplot(LSV_ARR_data, aes(x = LSV_max, y = ARR)) + geom_point() + geom_smooth(method = "lm") 


## lm by class mean
lm_output_class_mean <- LSV_ARR_data %>%
  drop_na(ARR,LSV_mean) %>%
  group_by(class.y) %>%
  do(tidy(lm(ARR ~ LSV_mean, data = .), conf.int = TRUE))

ggplot(LSV_ARR_data, aes(x = LSV_mean, y = ARR, col = class.y)) + geom_point() + geom_smooth(method = "lm") + labs(col = "Class")


## SMA regression:
sma_output_mean <- LSV_ARR_data %>%
  drop_na(ARR,LSV_mean) %>%
  lmodel2(ARR ~ LSV_mean, data = .) %>%
  plot("SMA") 


###### test whether sd and not sliding window works 
sd_vals <- c()

i <- 1
while (i < 409) {
  name <- paste("~/Documents/SUNDAY LAB/collection_location_tmax_series/loc", i, ".csv", sep = "")
  loc <- read.csv(name)
  
  standard_dev <- sd(loc$temp_value, na.rm = TRUE)
  sd_vals <- append(sd_vals, standard_dev, after = length(sd_vals))
  i = i+1
}

study_locs <- read.csv("~/Documents/SUNDAY LAB/Intratherm/Data sheets/locations_unique.csv")
study_locs$sd <- sd_vals




sd_data <- left_join(ARR_data, study_locs, by = c("latitude", "longitude"))  
sd_data <- sd_data[!duplicated(sd_data[,2]),]
sd_data <- subset(sd_data, subset = !is.na(ARR))

## does sd increase with increasing latitude as expected?
ggplot(sd_data, aes(x = abs(latitude), y = sd, col = (latitude >0))) + geom_point() 

## does sd at location over all years predict ARR? aka without sliding window
ggplot(sd_data, aes(x = sd, y = ARR, col = abs(latitude))) + geom_point()+ geom_smooth(method = "lm", se = TRUE, col = "black") + labs(x = "Standard deviation of monthly average TMAX", y = "Acclimation response ratio (ARR)") + scale_colour_viridis_c(option = "A") + geom_errorbar(aes(ymin=ARR-std.error, ymax=ARR+std.error), width=0) 

lm_output_all_sd <- sd_data %>%
  drop_na(ARR,sd) %>%
  do(tidy(lm(ARR ~ sd, data = .), conf.int = TRUE))


## does ARR vary predictably with latitude 
ggplot(sd_data, aes(x = abs(latitude), y = ARR)) + geom_point() + geom_smooth(method = "lm")


## does ARR vary predictably with maximum tmax value for a species? 
max_param <- merged %>%
  group_by(population_id) %>%
  arrange( desc(parameter_value) ) %>% 
  slice(1) 

ggplot(max_param, aes(x = parameter_value, y = ARR)) + geom_point() + geom_smooth(method = "lm")



### 
ggplot(LSV_ARR_data, aes(x = lifespan_days, y = LSV_mean)) + geom_point() + geom_smooth(method = "lm")

### LSV vs latitude 
ggplot(LSV_ARR_data, aes(x = abs(latitude), y = LSV_meanIQR)) + geom_point() + geom_smooth(method = "lm", col = "black") + labs(x = "abs(latitude)", y = "Mean range between 0.025 quantile and 0.975 quantile") 


### body size
lm_output_bs <- LSV_ARR_data %>%
  drop_na(ARR,maximum_body_size_svl_hbl_cm) %>%
  do(tidy(lm(ARR ~ maximum_body_size_svl_hbl_cm, data = .), conf.int = TRUE))

ggplot(LSV_ARR_data, aes(x = maximum_body_size_svl_hbl_cm, y = ARR)) + geom_point() + geom_smooth(method = "lm")


ggplot(LSV_ARR_data, aes(x = abs(latitude), y = lifespan_days)) + geom_point() + geom_smooth(method = "lm")

















