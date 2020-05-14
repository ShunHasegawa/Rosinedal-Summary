library(maps)
library(mapdata)
library(maptools)  #for shapefiles
library(scales)  #for transparency
pcontorta <- readShapePoly("pinucont.shp")   #layer of data for species range
samps <- read.csv("FieldSamples.csv")   #my data for sampling sites, contains a column of "lat" and a column of "lon" with GPS points in decimal degrees
map("worldHires","Sweden", xlim=c(-140,-110),ylim=c(48,64), col="gray90", fill=TRUE)  #plot the region of Canada I want
map("worldHires","Sweden", col="gray90", fill=TRUE)  
# map.cities("Sweden", "Västerbotten")
map.cities(world.cities, country = "Sweden", "Umea")

range(transect_cor$lon)
range(transect_cor$lat)

data(world.cities)
world.cities %>% 
  filter(country.etc == "Sweden")
map("worldHires","usa", xlim=c(-140,-110),ylim=c(48,64), col="gray95", fill=TRUE, add=TRUE)  #add the adjacent parts of the US; can't forget my homeland
plot(pcontorta, add=TRUE, xlim=c(-140,-110),ylim=c(48,64), col=alpha("darkgreen", 0.6), border=FALSE)  #plot the species range
points(samps$lon, samps$lat, pch=19, col="red", cex=0.5)  #plot my sample sites
