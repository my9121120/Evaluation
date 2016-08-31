library(ggmap)
library(ggplot2)
library(sp)
library(gstat)
library(graphics)
load("DQData")

DQCluster = DQpropertyall[DQpropertyall$X < 517600,]
DQFurther = DQpropertyall[DQpropertyall$X > 517600,]

xmin <- min(DQCluster$X)
xmax <- max(DQCluster$X)

ymin <- min(DQCluster$Y)
ymax <- max(DQCluster$Y)

length <- xmax - xmin
width <- ymax - ymin
area <- length * width

x <- seq(xmin, xmax, by = 2)
y <- seq(ymin, ymax, by = 2)

DQGrid <- expand.grid(X = x, Y = y)
#DQGrid[,"mark"] = 0
#DQCluster$pH = 5
coordinates(DQCluster)<-c("X", "Y")
proj4string(DQCluster)<-CRS("+proj=tmerc +lat_0=0 +lon_0=120
                                +k=1 +x_0=500000 +y_0=0 +ellps=krass
                                +units=m +no_defs")


#coordinates(DQGrid)<-c("X", "Y")
gridded(DQGrid) <- ~X+Y
proj4string(DQGrid) <- CRS("+proj=tmerc +lat_0=0 +lon_0=120
                                +k=1 +x_0=500000 +y_0=0 +ellps=krass
                                +units=m +no_defs")

scale = list("SpatialPolygonsRescale", layout.scale.bar(),
             offset = c(517205,3383090), scale = 20, fill=c("transparent","black"))
text1 = list("sp.text", c(517205,3383095), "0", cex = 0.7)
text2 = list("sp.text", c(517225,3383095), "20 m", cex = 0.7)
arrow = list("SpatialPolygonsRescale", layout.north.arrow(),
             offset = c(517210,3383140), scale = 20)
pts = list("sp.points", DQCluster, cex= 0.5, pch = 4, col = "black")


s <- spplot(DQCluster, 'pH', scales=list(draw=T), cex = 0.1, cuts = 2,
       sp.layout = list(text1, text2, scale, arrow, pts),
       key.space = "right")
s$legend <- NULL
s





# x <- gstat(formula=DQCluster$pH ~ 1, data = DQCluster, degree = 1)
# kt <- predict(x, newdata = DQGrid)
# spplot(kt[1])
