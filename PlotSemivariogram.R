library("gstat")
library(sp)
library(e1071)
library(stats)
library(graphics)
library(rgdal)
rm(list = ls())
#a = read.csv("DQpropertyall.csv", numerals = "no.loss")
source("multicrossvalidation.R")
load("DQData")
properties = names(DQpropertyall)


DQCluster = DQpropertyall[DQpropertyall$X < 517600,]

DQCluster[, "lnOM"] = log(DQCluster$OM.g.kg.1.)
crs <- CRS("+proj=tmerc +lat_0=0 +lon_0=120
           +k=1 +x_0=500000 +y_0=0 +ellps=krass
           +units=m +no_defs")
loc <- c("X", "Y")
target <- "pH"

Cir_vgm <- generate_vgm(target, DQCluster, "Cir", crs, loc, 0.5)
Sph_vgm <- generate_vgm(target, DQCluster, "Sph", crs, loc, 0.5)
Exp_vgm <- generate_vgm(target, DQCluster, "Exp", crs, loc, 0.5)

target <- "OM.g.kg.1."

Cir_vgm_OM <- generate_vgm(target, DQCluster, "Cir", crs, loc, 0.5)
Sph_vgm_OM <- generate_vgm(target, DQCluster, "Sph", crs, loc, 0.5)
Exp_vgm_OM <- generate_vgm(target, DQCluster, "Exp", crs, loc, 0.5)

target <- "lnOM"

Cir_vgm_lnOM <- generate_vgm(target, DQCluster, "Cir", crs, loc, 0.5)
Sph_vgm_lnOM <- generate_vgm(target, DQCluster, "Sph", crs, loc, 0.5)
Exp_vgm_lnOM <- generate_vgm(target, DQCluster, "Exp", crs, loc, 0.5)

par(mfrow = c(3,1))

plot(gamma~dist, Cir_vgm$vgm,ylim = c(0, 1.05*max(Cir_vgm$vgm$gamma)),col="blue", 
     ylab = 'semivariance', xlab= 'distance(m)', main = "(a)",font.main = 1)
maxdist <- round(max(Cir_vgm$vgm$dist),0)
lines(variogramLine(Cir_vgm$fitted_vgm,maxdist),lty = 1, col ="black")
lines(variogramLine(Sph_vgm$fitted_vgm,maxdist),lty = 3, col ="green")
lines(variogramLine(Exp_vgm$fitted_vgm,maxdist),lty = 5, col = "red")
legend(x = 80,y =0.1, legend=c("Circular model","Spherical model","Exponential mode"),
       lty = c(1, 3, 5), col = c("black","green","red"), bg = "transparent")

plot(gamma~dist, Cir_vgm_OM$vgm,ylim = c(0, 1.05*max(Cir_vgm_OM$vgm$gamma)),col="blue", 
     ylab = 'semivariance', xlab= 'distance(m)', main = "(b)",font.main = 1)
maxdist <- round(max(Cir_vgm_OM$vgm$dist),0)
lines(variogramLine(Cir_vgm_OM$fitted_vgm,maxdist),lty = 1, col ="black")
lines(variogramLine(Sph_vgm_OM$fitted_vgm,maxdist),lty = 3, col ="green")
lines(variogramLine(Exp_vgm_OM$fitted_vgm,maxdist),lty = 5, col = "red")
legend(x = 80,y =20, legend=c("Circular model","Spherical model","Exponential mode"),
       lty = c(1, 3, 5), col = c("black","green","red"), bg = "transparent")

plot(gamma~dist, Cir_vgm_lnOM$vgm,ylim = c(0, 1.05*max(Cir_vgm_lnOM$vgm$gamma)),col="blue", 
     ylab = 'semivariance', xlab= 'distance(m)', main = "(c)",font.main = 1)
maxdist <- round(max(Cir_vgm_lnOM$vgm$dist),0)
lines(variogramLine(Cir_vgm_lnOM$fitted_vgm,maxdist),lty = 1, col ="black")
lines(variogramLine(Sph_vgm_lnOM$fitted_vgm,maxdist),lty = 3, col ="green")
lines(variogramLine(Exp_vgm_lnOM$fitted_vgm,maxdist),lty = 5, col = "red")
legend(x = 80,y =0.06, legend=c("Circular model","Spherical model","Exponential mode"),
       lty = c(1, 3, 5), col = c("black","green","red"), bg = "transparent")



# maxdist <- round(max(Cir_vgm$vgm$dist),0)
# mypanel = function(x,y,...) {                                                 
#   vgm.panel.xyplot(x,y,...)
#   panel.lines(variogramLine(Sph_vgm$fitted_vgm,maxdist),lty = 2, col ="green")
#   panel.lines(variogramLine(Exp_vgm$fitted_vgm,maxdist),lty = 3, col = "red")
# }
# plot(Cir_vgm$vgm, model=Cir_vgm$fitted_vgm, panel=mypanel)
# 
# maxdist <- round(max(Cir_vgm_OM$vgm$dist),0)
# mypanel = function(x,y,...) {                                                 
#   vgm.panel.xyplot(x,y,...)
#   panel.lines(variogramLine(Sph_vgm_OM$fitted_vgm,maxdist),lty = 2, col ="green")
#   panel.lines(variogramLine(Exp_vgm_OM$fitted_vgm,maxdist),lty = 3, col = "red")
# }
# print(plot(Cir_vgm_OM$vgm, model=Cir_vgm_OM$fitted_vgm,
#            col = "black", panel=mypanel))
# 
# maxdist <- round(max(Cir_vgm_lnOM$vgm$dist),0)
# mypanel = function(x,y,...) {                                                 
#   vgm.panel.xyplot(x,y,...)
#   panel.lines(variogramLine(Sph_vgm_lnOM$fitted_vgm,maxdist),lty = 2, col ="green")
#   panel.lines(variogramLine(Exp_vgm_lnOM$fitted_vgm,maxdist),lty = 3, col = "red")
# }
# print(plot(Cir_vgm_lnOM$vgm, model=Cir_vgm_lnOM$fitted_vgm,
#            col = "black", panel=mypanel))