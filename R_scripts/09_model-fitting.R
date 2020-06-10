## fitting models 
library(nlme)
library(MuMIn)
library(lme4)
library(tidyverse)
library(broom)

data <- read.csv("./data-processed/arr_sliding-window-output-with-arr.csv")

## get rid of multiple rows per population and ARRs less than 0.15
model_input <- data %>%
  filter(!duplicated(population_id)) %>%
  filter(as.numeric(as.character(ARR)) > -0.15) %>%
  filter(!is.na(experienced_var_mean)) %>%
  filter(!is.na(ARR))

## figure making for paper:
data$ARR[which.min(data$experienced_var_mean)]

## add line representing other model fit:
xs = c(1.797841, 15.37186)
y1 = 0.005559816*1.797841 + 0.224658873 
y2 = 0.005559816*15.37186 + 0.224658873 
beta = c(y1, y2)

geom_abline(a = 0.224658873, b = 0.005559816, from = 0, to = 12) 

g = ggplot(model_input, aes(x = experienced_var_mean, y = ARR, col = abs(latitude))) + 
  geom_point() + 
  geom_smooth(method = "lm", color = "black") + 
  geom_errorbar(aes(ymin=ARR-std.error, ymax=ARR+std.error), width=0) + 
  labs(x = "Mean experienced variation (°C)", y = "Acclimation response ratio", col = paste(" Absolute ", "\n", "latitude (°N/S)" )) +
  scale_colour_viridis_c(option = "A") + 
  geom_segment(aes(x = xs[1], xend = xs[2], y = beta[1], yend = beta[2]), colour = "black",lty = "dashed")

ggsave("./Figures/experienced-var-vs-arr.png", g)
  



### basic linear model:
# experienced_var_mean
fit1 <- lm(ARR ~ experienced_var_mean, data = data)


### linear mixed effect models accounting for relatedness using a nested grouping structure: 
# experienced_var_mean
fit5 <- lme(ARR ~ experienced_var_mean, random = ~1|order/family/genus/species, 
            data = model_input, na.action= NULL)
# experienced_var_mean + body size
fit6 <- lme(ARR ~ experienced_var_mean*maximum_body_size_svl_hbl_cm, random = ~1|order/family/genus/species, 
            data = model_input, na.action= NULL)
# experienced_var_mean + realm
fit7 <- lme(ARR ~ experienced_var_mean*realm_general2, random = ~1|order/family/genus/species, 
            data = model_input, na.action= NULL)
# experienced_var_mean + body size + realm
fit8 <- lme(ARR ~ experienced_var_mean*maximum_body_size_svl_hbl_cm+realm_general2, 
            random = ~1|order/family/genus/species, data = model_input, na.action= NULL)

model.sel(fit1,fit5, fit6, fit7, fit8)

plot(residuals(fit5), col = data$species)

plot(resid(fit5))
qqnorm(resid(fit5))
qqline(resid(fit5))


### comparing models with and without error propagation:
with_err <- subset(data, std.error != "NaN")
with_err <- subset(with_err, std.error != "0")

## all data
fit9 <- lm(ARR ~ experienced_var_mean, data = data) 
## data with std error estimates
fit10 <- lm(ARR ~ experienced_var_mean, data =  with_err) 
## propagating error for data with std error estimates

fit11 <- lm(ARR ~ experienced_var_mean, data =  with_err, weights = (std.error^2)) 
## propagating error for all data, NaNs = lowest weight in dataset 
weights = (data$std.error^2)
lowest = weights[which.min(weights)]
weights = replace(weights, which(is.na(weights)), lowest)
weights = replace(weights, which(is.infinite(weights)), lowest)
fit12 <- lm(ARR ~ experienced_var_mean, data =  data, weights = weights) 

model.sel(fit10, fit11)
model.sel(fit9, fit12)




means <- data %>%
  group_by(phylum) %>%   #make groups of each species for every lat and long
  summarize(meanARR = mean(ARR, na.rm=TRUE), sdARR = sd(ARR, na.rm=TRUE), meanexpvar = mean(experienced_var_mean, na.rm=TRUE) , sdexpvar = sd(experienced_var_mean, na.rm=TRUE)) 

ggplot(means, aes(x = meanexpvar, y = meanARR, col = phylum)) + 
  geom_point()  + 
  geom_smooth(method = "lm")

means <- data %>%
  group_by(class) %>%   #make groups of each species for every lat and long
  summarize(meanARR = mean(ARR, na.rm=TRUE), sdARR = sd(ARR, na.rm=TRUE), meanexpvar = mean(experienced_var_mean, na.rm=TRUE) , sdexpvar = sd(experienced_var_mean, na.rm=TRUE)) 

ggplot(means, aes(x = meanexpvar, y = meanARR, col = class)) + 
  geom_point()  + 
  geom_smooth(method = "lm")

means <- data %>%
  group_by(order) %>%   #make groups of each species for every lat and long
  summarize(meanARR = mean(ARR, na.rm=TRUE), sdARR = sd(ARR, na.rm=TRUE), meanexpvar = mean(experienced_var_mean, na.rm=TRUE) , sdexpvar = sd(experienced_var_mean, na.rm=TRUE)) 

ggplot(means, aes(x = meanexpvar, y = meanARR, col = order)) + 
  geom_point()  + 
  geom_smooth(method = "lm")

means <- data %>%
  group_by(family) %>%   #make groups of each species for every lat and long
  summarize(meanARR = mean(ARR, na.rm=TRUE), sdARR = sd(ARR, na.rm=TRUE), meanexpvar = mean(experienced_var_mean, na.rm=TRUE) , sdexpvar = sd(experienced_var_mean, na.rm=TRUE)) 

ggplot(means, aes(x = meanexpvar, y = meanARR, col = family)) + 
  geom_point()  + 
  geom_smooth(method = "lm")

means <- data %>%
  group_by(genus) %>%   #make groups of each species for every lat and long
  summarize(meanARR = mean(ARR, na.rm=TRUE), sdARR = sd(ARR, na.rm=TRUE), meanexpvar = mean(experienced_var_mean, na.rm=TRUE) , sdexpvar = sd(experienced_var_mean, na.rm=TRUE)) 

ggplot(means, aes(x = meanexpvar, y = meanARR, col = genus)) + 
  geom_point()  + 
  geom_smooth(method = "lm")
  
means <- data %>%
  group_by(species) %>%   #make groups of each species for every lat and long
  summarize(meanARR = mean(ARR, na.rm=TRUE), sdARR = sd(ARR, na.rm=TRUE), meanexpvar = mean(experienced_var_mean, na.rm=TRUE) , sdexpvar = sd(experienced_var_mean, na.rm=TRUE)) 

ggplot(means, aes(x = meanexpvar, y = meanARR, col = species)) + 
  geom_point()  + 
  geom_smooth(method = "lm")