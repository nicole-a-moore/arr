## calculates ARR for each population
library(tidyverse)
library(broom)
library(lmodel2)
library(janitor)
library(RColorBrewer)
library(viridis)
library(nlme)

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
g = ggplot(data = intratherm, aes(x = acclim_temp, y = parameter_value, col = population_id)) + 
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

ggsave("./Figures/ARR-model-fits-all.png", g)

## merge ARR results for each popultion back to intratherm
arr <- arr %>%
  select(-term)

intratherm <- left_join(intratherm, arr, by = "population_id")

## write to csv
write.csv(intratherm, "./data-processed/arr_sliding-window-output-with-arr.csv", row.names = FALSE)






##get rid of duplicated rows since acclimation temp no longer matters, get rid of negative ARRs
intratherm2 <- intratherm %>%
  filter(!duplicated(intratherm$population_id)) %>%
  filter(as.numeric(as.character(ARR)) > -0.15)

write.csv(intratherm2, "~/Documents/SUNDAY LAB/Intratherm/Data sheets/phylo_input.csv")

## colour plot by experienced variation metric 
g = ggplot(data = intratherm, aes(x = acclim_temp, y = parameter_value, group = population_id, 
                              col = experienced_var_mean)) + 
  geom_point(size = 1, alpha = 0.5) + 
  geom_errorbar(aes(ymin=parameter_value-error_estimate,ymax=parameter_value+error_estimate),
                width=0, alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, size = 0.4, alpha = 0.1) +
  theme_minimal() + 
  theme(panel.border = element_rect(colour = "dimgrey", fill = NA), 
        panel.grid.minor = element_line(colour = "white"), panel.grid.major = 
        element_line(colour = "white"), plot.title = element_text(hjust = 0.5)) + 
  labs(x = "Acclimation temperature (°C)",y = "Upper critical thermal limit (°C)", 
       colour = paste(" Mean", "\n", "experienced", "\n", "variation")) +
  scale_colour_viridis_c(option = "A")

ggsave("./Figures/ARR-model-fits-col-is-exp-var.png", g)

## linear model of experienced_var_mean weighted by 1/variance of ARR stope estimate
## just kidding, scrap this because most ARR don't have an error estimate
lm_output_all_mean <- intratherm2 %>%
  drop_na(ARR, experienced_var_mean) %>%
  do(tidy(lm(ARR ~ experienced_var_mean, data = .), conf.int = TRUE))

lm_output_all_top<- intratherm2 %>%
  drop_na(ARR, experienced_var_mean) %>%
  do(tidy(lme(ARR ~ experienced_var_mean, random = ~1|order/family/genus/species,
              data = ., na.action= NULL), conf.int = TRUE))

ggplot(intratherm2, aes(x = experienced_var_mean, y = ARR, col = abs(latitude))) + 
  geom_point() + 
  geom_smooth(method = "lm", color = "black") + 
  geom_errorbar(aes(ymin=ARR-std.error, ymax=ARR+std.error), width=0) + 
  labs(x = "Mean experienced variation (°C)", y = "Acclimation response ratio", 
       col = "Absolute latitude (°N/S)" ) + 
  scale_colour_viridis_c(option = "A")



## write to file to make map + taxonomic breakdown from
intratherm2 <- intratherm2 %>%
  filter(realm_general2 != "Marine") 

write.csv(intratherm2, "./data-processed/arr_map-and-model-input.csv", row.names = FALSE)

## lm with max sd
lm_output_all_max<- intratherm2 %>%
  drop_na(ARR,experienced_var_max) %>%
  do(tidy(lm(ARR ~ experienced_var_max, data = .), conf.int = TRUE))

ggplot(intratherm2, aes(x = experienced_var_max, y = ARR)) + geom_point() + geom_smooth(method = "lm") 


## lm by class mean
lm_output_class_mean <- intratherm2 %>%
  drop_na(ARR,experienced_var_mean) %>%
  group_by(class) %>%
  mutate(experienced_var_mean = mean(experienced_var_mean)) %>%
  do(tidy(lm(ARR ~ experienced_var_mean, data = .), conf.int = TRUE))

ggplot(intratherm2, aes(x = experienced_var_mean, y = ARR, col = class)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  labs(col = "Class")


## SMA regression:
sma_output_mean <- intratherm2 %>%
  drop_na(ARR,experienced_var_mean) %>%
  lmodel2(ARR ~ experienced_var_mean, data = .) %>%
  plot("SMA") 


###### EXPLORE!! ------------------------------------------

## does sd increase with increasing latitude as expected?
ggplot(intratherm2, aes(x = abs(latitude), y = sd_cumulative, col = (latitude >0))) +
  geom_point() 

## does sd at location over all years predict ARR? aka without sliding window
ggplot(intratherm2, aes(x = sd_cumulative, y = ARR, col = abs(latitude))) + 
  geom_point() + 
  geom_smooth(method = "lm", se = TRUE, col = "black") +
  labs(x = "Standard deviation of monthly average TMAX", y = "Acclimation response ratio (ARR)") + 
  scale_colour_viridis_c(option = "A") + 
  geom_errorbar(aes(ymin=ARR-std.error, ymax=ARR+std.error), width=0) 

lm_output_all_sd <- intratherm2 %>%
  drop_na(ARR, sd_cumulative) %>%
  do(tidy(lm(ARR ~ sd_cumulative, data = .), conf.int = TRUE))


## does ARR vary predictably with latitude 
ggplot(intratherm2, aes(x = abs(latitude), y = ARR)) + 
  geom_point() + 
  geom_smooth(method = "lm")


## does ARR vary predictably with maximum tmax value for a species? 
max_param <- intratherm %>%
  group_by(population_id) %>%
  arrange(desc(parameter_value)) %>% 
  slice(1) 

ggplot(max_param, aes(x = parameter_value, y = ARR)) + 
  geom_point() + 
  geom_smooth(method = "lm")



### 
ggplot(intratherm2, aes(x = lifespan_days, y = experienced_var_mean)) + 
  geom_point() + 
  geom_smooth(method = "lm")

### experienced_var vs latitude 
ggplot(intratherm2, aes(x = abs(latitude), y = experienced_var_mean)) + 
  geom_point() + 
  geom_smooth(method = "lm", col = "black") + 
  labs(x = "abs(latitude)", y = "experienced_var_mean") 


### body size
lm_output_bs <- intratherm2 %>%
  drop_na(ARR, maximum_body_size_svl_hbl_cm) %>%
  do(tidy(lm(ARR ~ maximum_body_size_svl_hbl_cm, data = .), conf.int = TRUE))

ggplot(intratherm2, aes(x = maximum_body_size_svl_hbl_cm, y = ARR)) + 
  geom_point() + 
  geom_smooth(method = "lm")


ggplot(intratherm2, aes(x = abs(latitude), y = lifespan_days)) + 
  geom_point() + 
  geom_smooth(method = "lm")