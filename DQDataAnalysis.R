library("gstat")
library(sp)
library(rgdal)
rm(list = ls())
a = read.csv("DQpropertyall.csv", numerals = "no.loss")
load("DQData")
properties = names(DQpropertyall)
plot(DQpropertyall$X, DQpropertyall$Y)
plot(DQpropertyall$OM.g.kg.1., DQpropertyall$pH)

DQCluster = DQpropertyall[DQpropertyall$X < 517600,]
DQCluster = DQCluster[DQCluster$Y > 3383000,]
DQFurther = DQpropertyall[DQpropertyall$X > 517600,]
DQFurther = DQFurther[DQCluster$Y < 3383000,]

coordinates(DQCluster)<-c("X", "Y")
proj4string(DQCluster)<-CRS("+proj=tmerc +lat_0=0 +lon_0=120
                                +k=1 +x_0=500000 +y_0=0 +ellps=krass
                                +units=m +no_defs")
hscat(DQCluster$pH~1, DQCluster, breaks = (0:9)*10)
vgm = variogram(DQCluster$pH~1, DQCluster)
fittingvgm = fit.variogram(vgm, vgm("Sph"))
plot(vgm, fittingvgm)
plot(vgm, fit.variogram(vgm, vgm("Exp")))

plot(variogram(DQCluster$pH~1, DQCluster))
plot(variogram(DQCluster$pH~1, DQCluster, 
               alpha = c(0, 45, 90, 135)))

loc = data.frame(DQpropertyall)[,c("X", "Y")]
loccluster = loc[loc$X<517600 & loc$Y>3383000,]
plot(loccluster$X, loccluster$Y)
distance = dist(loc)
maxdist = max(distance)
mindist = min(distance)
distance = dist(loccluster)
maxdist = max(distance)
mindist = min(distance)




coordinates(DQpropertyall)<-c("X", "Y")
spplot(DQpropertyall,"pH", cex = 1)


#The coordinate system might be EPSG:2412 
#Beijing 1954 / 3-degree Gauss-Kruger zone 36
proj4string(DQpropertyall)<-CRS("+proj=tmerc +lat_0=0 +lon_0=120
                                +k=1 +x_0=500000 +y_0=0 +ellps=krass
                                +units=m +no_defs")

longlat = CRS("+proj=longlat +ellps=WGS84")
DQpropertyall = spTransform(DQpropertyall, longlat)
DQpropertyall = as.data.frame(DQpropertyall)
loc = DQpropertyall[,c("X", "Y")]
plot(loc$X, loc$Y, cex = 0.5)


hist(DQpropertyall$pH)
hist(DQpropertyall$OM.g.kg.1.)
hist(DQpropertyall$silt.)
hist(DQpropertyall$clay.)
hist(DQpropertyall$sand.)