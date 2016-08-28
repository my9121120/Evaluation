library(ggmap)
library(ggplot2)
library(sp)
library(gstat)
load("DQData")

DQCluster = DQpropertyall[DQpropertyall$X < 517600,]
DQFurther = DQpropertyall[DQpropertyall$X > 517600,]

xmin <- min(DQCluster$X)
xmax <- max(DQCluster$X)

ymin <- min(DQCluster$Y)
ymax <- max(DQCluster$Y)

x <- seq(xmin, xmax, by = 2)
y <- seq(ymin, ymax, by = 2)

DQGrid <- expand.grid(X = x, Y = y)
#DQGrid[,"mark"] = 0

coordinates(DQCluster)<-c("X", "Y")
proj4string(DQCluster)<-CRS("+proj=tmerc +lat_0=0 +lon_0=120
                                +k=1 +x_0=500000 +y_0=0 +ellps=krass
                                +units=m +no_defs")


#coordinates(DQGrid)<-c("X", "Y")
gridded(DQGrid) <- ~X+Y
proj4string(DQGrid) <- CRS("+proj=tmerc +lat_0=0 +lon_0=120
                                +k=1 +x_0=500000 +y_0=0 +ellps=krass
                                +units=m +no_defs")

x <- gstat(formula=DQCluster$pH ~ 1, data = DQCluster, degree = 1)
kt <- predict(x, newdata = DQGrid)
spplot(kt[1])
