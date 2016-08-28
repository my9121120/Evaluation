library(ggmap)
library(ggplot2)
library(sp)
load("DQData")
properties = names(DQpropertyall)
plot(DQpropertyall$X, DQpropertyall$Y)
plot(DQpropertyall$OM.g.kg.1., DQpropertyall$pH)

DQCluster = DQpropertyall[DQpropertyall$X < 517600,]
DQFurther = DQpropertyall[DQpropertyall$X > 517600,]

coordinates(DQCluster)<-c("X", "Y")
proj4string(DQCluster)<-CRS("+proj=tmerc +lat_0=0 +lon_0=120
                                +k=1 +x_0=500000 +y_0=0 +ellps=krass
                                +units=m +no_defs")
longlat <- CRS("+proj=longlat +ellps=WGS84")
DQCluster_GPS <- spTransform(DQCluster, longlat)

DQCluster_GPS <- as.data.frame(DQCluster_GPS)
lon = mean(DQCluster_GPS[,"X"])
lat = mean(DQCluster_GPS[,"Y"])
map <- get_googlemap('Shanxi', zoom = 4, markers = data.frame(longitude=lon, latitude=lat)
                     , language = "en-EN")

ggmap(map) + xlab('Longitude') + ylab('Latitude')
