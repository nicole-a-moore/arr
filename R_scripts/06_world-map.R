## map making
## creates a plot of all collection locations 
library(ggplot2)
library(sf)
theme_set(theme_bw())
library(rnaturalearth)
library(rnaturalearthdata)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(ggplotify)
library(base2grob)


world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)
#> [1] "sf" "data.frame"

##plot world data using geometry stored in sf function
##all data must be in sf format


data <- read.csv("./data-processed/arr_sliding-window-output-with-arr.csv")
## get rid of multiple rows per population and ARRs less than 0.15
model_input <- data %>%
  filter(!duplicated(population_id)) %>%
  filter(as.numeric(as.character(ARR)) > -0.15) %>%
  filter(!is.na(experienced_var_mean)) %>%
  filter(!is.na(ARR)) 

model_input <-droplevels(model_input)

loc <- model_input %>%
  select(latitude, longitude, class, genus_species, realm_general2)

colnames(loc) <- c("Latitude", "Longitude", "Class", "Genus_species", "Realm")

bar <- as.data.frame(table(loc$Class))
colnames(bar) <- c("Class", "Frequency")
loc <- left_join(loc, bar)

## reorder classes by frequency
loc$Class <- reorder(loc$Class, loc$Frequency)


## make colour palette
myPalette <- brewer.pal(n = length(unique(loc$Class)), name = "Dark2") 


map = ggplot(data=world) + 
  geom_sf(color = "dimgrey", fill = "white", size = 0.1) + 
  theme(panel.grid.major = element_line(colour = "light grey", size = 0.05),
        panel.border = element_rect(colour = "transparent"), 
        panel.background = element_rect(fill = "grey96"), legend.position = "none") + 
  coord_sf(expand = FALSE) +
  scale_x_continuous(breaks = c(-150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150)) + 
  scale_y_continuous(breaks = c(-90, -60, -30, 0, 30, 60, 90))  +
  geom_point(data = loc, aes(y = Latitude, x = Longitude, 
                             col = Class), size = 0.5) +
  scale_colour_manual(values = myPalette)


## horizontal bar chart
bchart <- ggplot(data = loc, aes(x = Class, y = Frequency, fill = Class )) + 
  geom_bar(stat = "identity", position = 'dodge', width = 0.5) + 
  coord_flip() + 
  scale_fill_manual(values = myPalette) + 
  theme(legend.position = "none", panel.grid.major = element_line(colour = "transparent"), 
        panel.grid.minor = element_line(colour = "transparent"),
        panel.border = element_rect(colour = "transparent"), 
        panel.background = element_rect(fill = "transparent"), axis.title.x = element_blank(), 
        axis.title.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        axis.ticks.y = element_blank(), plot.title = element_text(hjust = -0.7, size = 12)) + 
  labs(title = "Class") 

blank <- ggplot() + geom_blank() + theme_minimal()

g1 <- as.grob(map)
g2 <- as.grob(bchart)
g3 <- as.grob(blank)

grobs <- list(g1, g2, g3)

lay <- rbind(c(1,1,1,1,3),
             c(1,1,1,1,2),
             c(1,1,1,1,2),
             c(1,1,1,1,2),
             c(1,1,1,1,3))


beautiful_map <- grid.arrange(grobs = grobs, layout_matrix = lay)

ggsave(beautiful_map, height = 5, width = 5*2.350109, path = "./Figures/", filename = "map-of-collection-locs.png")


