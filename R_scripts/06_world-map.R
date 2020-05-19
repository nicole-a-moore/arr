## map making
## creates a plot of all collection locations 

library(ggplot2)
library(sf)
theme_set(theme_bw())
#> Linking to GEOS 3.7.2, GDAL 2.4.2, PROJ 5.2.0
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

loc <- read.csv("")
loc <- data.frame(loc$latitude, loc$longitude, loc$class, loc$genus_species, loc$realm_general2)
colnames(loc) <- c("Latitude", "Longitude", "Class", "Genus_species", "Realm")

## make colour palette
myPalette <- brewer.pal(n = 4, name = "Dark2")

map = ggplot(data=world) + 
  geom_sf(color = "dimgrey", fill = "white", size = 0.1) + 
  theme(panel.grid.major = element_line(colour = "light grey", size = 0.05),
        panel.border = element_rect(colour = "transparent"), 
        panel.background = element_rect(fill = "grey96"), legend.position = "none") + coord_sf(expand = FALSE) +
  scale_x_continuous(breaks = c(-150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150)) + 
  scale_y_continuous(breaks = c(-90, -60, -30, 0, 30, 60, 90))  +
  geom_point(data = loc, aes(y = Latitude, x = Longitude, col = Class), size = 0.3) +
  scale_color_manual(values = myPalette)

## horizontal bar chart
bar <- as.data.frame(table(loc$Class))
colnames(bar) <- c("Class", "Frequency")
bchart <- ggplot(data = bar, aes(x = Class, y = Frequency, fill = Class )) + 
  geom_bar(stat = "identity", position = 'dodge', width = 0.5) + 
  coord_flip() + 
  scale_color_manual(values = myPalette) + 
  theme(legend.position = "none", panel.grid.major = element_line(colour = "transparent"), 
        panel.grid.minor = element_line(colour = "transparent"),
        panel.border = element_rect(colour = "transparent"), 
        panel.background = element_rect(fill = "transparent"), axis.title.x = element_blank(), 
        axis.title.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        axis.ticks.y = element_blank(), plot.title = element_text(hjust = -0.7, size = 12)) + labs(title = "Class") 

blank <- ggplot() + geom_blank() + theme_minimal()

g1 <- as.grob(map)
g2 <- as.grob(bchart)
g3 <- as.grob(blank)

grobs <- list(g1, g2, g3)

lay <- rbind(c(1,1,1,1,3),
             c(1,1,1,1,3),
             c(1,1,1,1,2),
             c(1,1,1,1,3),
             c(1,1,1,1,3))

grid.arrange(grobs = grobs, layout_matrix = lay)


## figure out data proportions by class and species 
table(loc$Class)
length(unique(loc$Genus_species))

## by realm 
table(loc$Realm)

## by classes
actino <- subset(loc, subset = Class == "Actinopterygii")
length(actino$Genus_species)
length(unique(actino$Genus_species))

amph <- subset(loc, subset = Class == "Amphibia")
length(amph$Genus_species)
length(unique(amph$Genus_species))

rept <- subset(loc, subset = Class == "Reptilia")
length(rept$Genus_species)
length(unique(rept$Genus_species))

insect <- subset(loc, subset = Class == "Insecta")
length(insect$Genus_species)
length(unique(insect$Genus_species))

