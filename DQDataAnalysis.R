library("gstat")
library(sp)
library(rgdal)
rm(list = ls())
#a = read.csv("DQpropertyall.csv", numerals = "no.loss")
source("multicrossvalidation.R")
load("DQData")
properties = names(DQpropertyall)
plot(DQpropertyall$X, DQpropertyall$Y)
plot(DQpropertyall$OM.g.kg.1., DQpropertyall$pH)

DQCluster = DQpropertyall[DQpropertyall$X < 517600,]
DQFurther = DQpropertyall[DQpropertyall$X > 517600,]

ratio = 0.8
set.seed(1)
DQCluster["mark"] = 0
row_counts <- dim(DQCluster)[1]
row_selected <- sample(1:row_counts, floor(row_counts * ratio))
DQCluster[row_selected, "mark"] = 1
train_data <- DQCluster[DQCluster["mark"] == 1,]
test_data <- DQCluster[DQCluster["mark"] == 0,]

crs <- CRS("+proj=tmerc +lat_0=0 +lon_0=120
                                +k=1 +x_0=500000 +y_0=0 +ellps=krass
                                +units=m +no_defs")
loc <- c("X", "Y")

train_data_sp <- train_data
test_data_sp <-test_data
coordinates(train_data_sp) <- loc
proj4string(train_data_sp) <-crs
coordinates(test_data_sp) <- loc
proj4string(test_data_sp) <-crs
row_names <- c("RMSEBasedPara", "MAEBasedPara", "MEBasedPara")

k <- 5
target <- "pH"
neighbors <- 4:20
powers <- 1:4

idw_para <- idw.para(target, train_data, neighbors, powers,
                     crs, loc, k)

prediction_RMSE <- idw(pH~1, train_data_sp, test_data_sp, idp = idw_para$RMSE$power,
                      nmax =  idw_para$RMSE$neighbor)
measure_RMSE <- measures(prediction_RMSE$var1.pred, test_data[,target])

prediction_MAE <- idw(pH~1, train_data_sp, test_data_sp, idp = idw_para$MAE$power,
                     nmax =  idw_para$MAE$neighbor)
measure_MAE <- measures(prediction_MAE$var1.pred, test_data[,target])

prediction_ME <- idw(pH~1, train_data_sp, test_data_sp, idp = idw_para$ME$power,
                    nmax =  idw_para$ME$neighbor)
measure_ME <- measures(prediction_ME$var1.pred, test_data[,target])

idw_prediction_dataframe <- data.frame("RMSEBasedPara" = prediction_RMSE$var1.pred,
                                       "MAEBasedPara" = prediction_MAE$var1.pred,
                                       "MEBasedPara" = prediction_ME$var1.pred,
                                       "Observation" = test_data[,target])

idw_measure_dataframe <- rbind(as.data.frame(measure_RMSE),
                               as.data.frame(measure_MAE),
                               as.data.frame(measure_ME))
rownames(idw_measure_dataframe) <- row_names


degrees <- 1:3
surf_para <- surf.para(target, train_data, degrees, crs, loc, k)

prediction_RMSE <- krige(pH~1, train_data_sp, test_data_sp, 
                         degree = surf_para$RMSE$degree)
measure_RMSE <- measures(prediction_RMSE$var1.pred, test_data[,target])

prediction_MAE <- krige(pH~1, train_data_sp, test_data_sp, 
                        degree = surf_para$MAE$degree)
measure_MAE <- measures(prediction_MAE$var1.pred, test_data[,target])

prediction_ME <- krige(pH~1, train_data_sp, test_data_sp, 
                       degree = surf_para$ME$degree)
measure_ME <- measures(prediction_ME$var1.pred, test_data[,target])

surf_prediction_dataframe <- data.frame("RMSEBasedPara" = prediction_RMSE$var1.pred,
                                       "MAEBasedPara" = prediction_MAE$var1.pred,
                                       "MEBasedPara" = prediction_ME$var1.pred,
                                       "Observation" = test_data[,target])

surf_measure_dataframe <- rbind(as.data.frame(measure_RMSE),
                               as.data.frame(measure_MAE),
                               as.data.frame(measure_ME))
rownames(surf_measure_dataframe) <- row_names


Sph_vgm <- generate_vgm(target, train_data, "Sph", crs, loc)
plot(Sph_vgm$vgm, Sph_vgm$fitted_vgm)

Sph_para <- krige.para(target, train_data, neighbors, Sph_vgm$fitted_vgm, crs, loc, k)

prediction_RMSE <- krige(pH~1, train_data_sp, test_data_sp, 
                         model = Sph_vgm$fitted_vgm,
                         nmax = Sph_para$RMSE$neighbor)
measure_RMSE <- measures(prediction_RMSE$var1.pred, test_data[,target])

prediction_MAE <- krige(pH~1, train_data_sp, test_data_sp, 
                        model = Sph_vgm$fitted_vgm,
                        nmax = Sph_para$MAE$neighbor)
measure_MAE <- measures(prediction_MAE$var1.pred, test_data[,target])

prediction_ME <- krige(pH~1, train_data_sp, test_data_sp, 
                       model = Sph_vgm$fitted_vgm,
                       nmax = Sph_para$ME$neighbor)
measure_ME <- measures(prediction_ME$var1.pred, test_data[,target])

Sph_prediction_dataframe <- data.frame("RMSEBasedPara" = prediction_RMSE$var1.pred,
                                        "MAEBasedPara" = prediction_MAE$var1.pred,
                                        "MEBasedPara" = prediction_ME$var1.pred,
                                        "Observation" = test_data[,target])

Sph_measure_dataframe <- rbind(as.data.frame(measure_RMSE),
                                as.data.frame(measure_MAE),
                                as.data.frame(measure_ME))
rownames(Sph_measure_dataframe) <- row_names


Exp_vgm <- generate_vgm(target, train_data, "Exp", crs, loc)
plot(Exp_vgm$vgm, Exp_vgm$fitted_vgm)

Exp_para <- krige.para(target, train_data, neighbors, Exp_vgm$fitted_vgm, crs, loc, k)

prediction_RMSE <- krige(pH~1, train_data_sp, test_data_sp, 
                         model = Exp_vgm$fitted_vgm,
                         nmax = Exp_para$RMSE$neighbor)
measure_RMSE <- measures(prediction_RMSE$var1.pred, test_data[,target])

prediction_MAE <- krige(pH~1, train_data_sp, test_data_sp, 
                        model = Exp_vgm$fitted_vgm,
                        nmax = Exp_para$MAE$neighbor)
measure_MAE <- measures(prediction_MAE$var1.pred, test_data[,target])

prediction_ME <- krige(pH~1, train_data_sp, test_data_sp, 
                       model = Exp_vgm$fitted_vgm,
                       nmax = Exp_para$ME$neighbor)
measure_ME <- measures(prediction_ME$var1.pred, test_data[,target])

Exp_prediction_dataframe <- data.frame("RMSEBasedPara" = prediction_RMSE$var1.pred,
                                       "MAEBasedPara" = prediction_MAE$var1.pred,
                                       "MEBasedPara" = prediction_ME$var1.pred,
                                       "Observation" = test_data[,target])

Exp_measure_dataframe <- rbind(as.data.frame(measure_RMSE),
                               as.data.frame(measure_MAE),
                               as.data.frame(measure_ME))
rownames(Exp_measure_dataframe) <- row_names



Cir_vgm <- generate_vgm(target, train_data, "Cir", crs, loc)
plot(Cir_vgm$vgm, Cir_vgm$fitted_vgm)

Cir_para <- krige.para(target, train_data, neighbors, Cir_vgm$fitted_vgm, crs, loc, k)

prediction_RMSE <- krige(pH~1, train_data_sp, test_data_sp, 
                         model = Cir_vgm$fitted_vgm,
                         nmax = Cir_para$RMSE$neighbor)
measure_RMSE <- measures(prediction_RMSE$var1.pred, test_data[,target])

prediction_MAE <- krige(pH~1, train_data_sp, test_data_sp, 
                        model = Cir_vgm$fitted_vgm,
                        nmax = Cir_para$MAE$neighbor)
measure_MAE <- measures(prediction_MAE$var1.pred, test_data[,target])

prediction_ME <- krige(pH~1, train_data_sp, test_data_sp, 
                       model = Cir_vgm$fitted_vgm,
                       nmax = Cir_para$ME$neighbor)
measure_ME <- measures(prediction_ME$var1.pred, test_data[,target])

Cir_prediction_dataframe <- data.frame("RMSEBasedPara" = prediction_RMSE$var1.pred,
                                       "MAEBasedPara" = prediction_MAE$var1.pred,
                                       "MEBasedPara" = prediction_ME$var1.pred,
                                       "Observation" = test_data[,target])

Cir_measure_dataframe <- rbind(as.data.frame(measure_RMSE),
                               as.data.frame(measure_MAE),
                               as.data.frame(measure_ME))
rownames(Cir_measure_dataframe) <- row_names


idw_para_dataframe <- para_dataframe(idw_para)
surf_para_dataframe <- para_dataframe(surf_para)
Sph_para_dataframe <- para_dataframe(Sph_para)
Exp_para_dataframe <- para_dataframe(Exp_para)
Cir_para_dataframe <- para_dataframe(Cir_para)

measures_frame = dataframe_list_convertion(list(idw_measure_dataframe,
                               surf_measure_dataframe,
                               Sph_measure_dataframe,
                               Exp_measure_dataframe,
                               Cir_measure_dataframe),
                               c("idw", "surf", "Sph",
                                 "Exp", "Cir"))
titles <- c("RMSE", "MAE", "ME")
titles <- paste(titles, "minmized based parameters selection")
measure_barplot(measures_frame, titles)


###################################################################################
#Predicting OM.g.kg.1. with various interpolation methods
target <- "OM.g.kg.1."
idw_para <- idw.para(target, train_data, neighbors, powers,
                     crs, loc, k)

prediction_RMSE <- idw(OM.g.kg.1.~1, train_data_sp, test_data_sp, 
                       idp = idw_para$RMSE$power,
                       nmax =  idw_para$RMSE$neighbor)
measure_RMSE <- measures(prediction_RMSE$var1.pred, test_data[,target])

prediction_MAE <- idw(OM.g.kg.1.~1, train_data_sp, test_data_sp, 
                      idp = idw_para$MAE$power,
                      nmax =  idw_para$MAE$neighbor)
measure_MAE <- measures(prediction_MAE$var1.pred, test_data[,target])

prediction_ME <- idw(OM.g.kg.1.~1, train_data_sp, test_data_sp, 
                     idp = idw_para$ME$power,
                     nmax =  idw_para$ME$neighbor)
measure_ME <- measures(prediction_ME$var1.pred, test_data[,target])

idw_prediction_dataframe <- data.frame("RMSEBasedPara" = prediction_RMSE$var1.pred,
                                       "MAEBasedPara" = prediction_MAE$var1.pred,
                                       "MEBasedPara" = prediction_ME$var1.pred,
                                       "Observation" = test_data[,target])

idw_measure_dataframe <- rbind(as.data.frame(measure_RMSE),
                               as.data.frame(measure_MAE),
                               as.data.frame(measure_ME))
rownames(idw_measure_dataframe) <- row_names


degrees <- 1:3
surf_para <- surf.para(target, train_data, degrees, crs, loc, k)

prediction_RMSE <- krige(OM.g.kg.1.~1, train_data_sp, test_data_sp, 
                         degree = surf_para$RMSE$degree)
measure_RMSE <- measures(prediction_RMSE$var1.pred, test_data[,target])

prediction_MAE <- krige(OM.g.kg.1.~1, train_data_sp, test_data_sp, 
                        degree = surf_para$MAE$degree)
measure_MAE <- measures(prediction_MAE$var1.pred, test_data[,target])

prediction_ME <- krige(OM.g.kg.1.~1, train_data_sp, test_data_sp, 
                       degree = surf_para$ME$degree)
measure_ME <- measures(prediction_ME$var1.pred, test_data[,target])

surf_prediction_dataframe <- data.frame("RMSEBasedPara" = prediction_RMSE$var1.pred,
                                        "MAEBasedPara" = prediction_MAE$var1.pred,
                                        "MEBasedPara" = prediction_ME$var1.pred,
                                        "Observation" = test_data[,target])

surf_measure_dataframe <- rbind(as.data.frame(measure_RMSE),
                                as.data.frame(measure_MAE),
                                as.data.frame(measure_ME))



Sph_vgm <- generate_vgm(target, train_data, "Sph", crs, loc)
plot(Sph_vgm$vgm, Sph_vgm$fitted_vgm)

Sph_para <- krige.para(target, train_data, neighbors, Sph_vgm$fitted_vgm, crs, loc, k)

prediction_RMSE <- krige(OM.g.kg.1.~1, train_data_sp, test_data_sp, 
                         model = Sph_vgm$fitted_vgm,
                         nmax = Sph_para$RMSE$neighbor)
measure_RMSE <- measures(prediction_RMSE$var1.pred, test_data[,target])

prediction_MAE <- krige(OM.g.kg.1.~1, train_data_sp, test_data_sp, 
                        model = Sph_vgm$fitted_vgm,
                        nmax = Sph_para$MAE$neighbor)
measure_MAE <- measures(prediction_MAE$var1.pred, test_data[,target])

prediction_ME <- krige(OM.g.kg.1.~1, train_data_sp, test_data_sp, 
                       model = Sph_vgm$fitted_vgm,
                       nmax = Sph_para$ME$neighbor)
measure_ME <- measures(prediction_ME$var1.pred, test_data[,target])

Sph_prediction_dataframe <- data.frame("RMSEBasedPara" = prediction_RMSE$var1.pred,
                                       "MAEBasedPara" = prediction_MAE$var1.pred,
                                       "MEBasedPara" = prediction_ME$var1.pred,
                                       "Observation" = test_data[,target])

Sph_measure_dataframe <- rbind(as.data.frame(measure_RMSE),
                               as.data.frame(measure_MAE),
                               as.data.frame(measure_ME))
rownames(Sph_measure_dataframe) <- row_names


Exp_vgm <- generate_vgm(target, train_data, "Exp", crs, loc)
plot(Exp_vgm$vgm, Exp_vgm$fitted_vgm)

Exp_para <- krige.para(target, train_data, neighbors, Exp_vgm$fitted_vgm, crs, loc, k)

prediction_RMSE <- krige(OM.g.kg.1.~1, train_data_sp, test_data_sp, 
                         model = Exp_vgm$fitted_vgm,
                         nmax = Exp_para$RMSE$neighbor)
measure_RMSE <- measures(prediction_RMSE$var1.pred, test_data[,target])

prediction_MAE <- krige(OM.g.kg.1.~1, train_data_sp, test_data_sp, 
                        model = Exp_vgm$fitted_vgm,
                        nmax = Exp_para$MAE$neighbor)
measure_MAE <- measures(prediction_MAE$var1.pred, test_data[,target])

prediction_ME <- krige(OM.g.kg.1.~1, train_data_sp, test_data_sp, 
                       model = Exp_vgm$fitted_vgm,
                       nmax = Exp_para$ME$neighbor)
measure_ME <- measures(prediction_ME$var1.pred, test_data[,target])

Exp_prediction_dataframe <- data.frame("RMSEBasedPara" = prediction_RMSE$var1.pred,
                                       "MAEBasedPara" = prediction_MAE$var1.pred,
                                       "MEBasedPara" = prediction_ME$var1.pred,
                                       "Observation" = test_data[,target])

Exp_measure_dataframe <- rbind(as.data.frame(measure_RMSE),
                               as.data.frame(measure_MAE),
                               as.data.frame(measure_ME))
rownames(Exp_measure_dataframe) <- row_names



Cir_vgm <- generate_vgm(target, train_data, "Cir", crs, loc)
plot(Cir_vgm$vgm, Cir_vgm$fitted_vgm)

Cir_para <- krige.para(target, train_data, neighbors, Cir_vgm$fitted_vgm, crs, loc, k)

prediction_RMSE <- krige(OM.g.kg.1.~1, train_data_sp, test_data_sp, 
                         model = Cir_vgm$fitted_vgm,
                         nmax = Cir_para$RMSE$neighbor)
measure_RMSE <- measures(prediction_RMSE$var1.pred, test_data[,target])

prediction_MAE <- krige(OM.g.kg.1.~1, train_data_sp, test_data_sp, 
                        model = Cir_vgm$fitted_vgm,
                        nmax = Cir_para$MAE$neighbor)
measure_MAE <- measures(prediction_MAE$var1.pred, test_data[,target])

prediction_ME <- krige(OM.g.kg.1.~1, train_data_sp, test_data_sp, 
                       model = Cir_vgm$fitted_vgm,
                       nmax = Cir_para$ME$neighbor)
measure_ME <- measures(prediction_ME$var1.pred, test_data[,target])

Cir_prediction_dataframe <- data.frame("RMSEBasedPara" = prediction_RMSE$var1.pred,
                                       "MAEBasedPara" = prediction_MAE$var1.pred,
                                       "MEBasedPara" = prediction_ME$var1.pred,
                                       "Observation" = test_data[,target])

Cir_measure_dataframe <- rbind(as.data.frame(measure_RMSE),
                               as.data.frame(measure_MAE),
                               as.data.frame(measure_ME))
rownames(Cir_measure_dataframe) <- row_names


idw_para_dataframe <- para_dataframe(idw_para)
surf_para_dataframe <- para_dataframe(surf_para)
Sph_para_dataframe <- para_dataframe(Sph_para)
Exp_para_dataframe <- para_dataframe(Exp_para)
Cir_para_dataframe <- para_dataframe(Cir_para)

measures_frame = dataframe_list_convertion(list(idw_measure_dataframe,
                                                surf_measure_dataframe,
                                                Sph_measure_dataframe,
                                                Exp_measure_dataframe,
                                                Cir_measure_dataframe),
                                           c("idw", "surf", "Sph",
                                             "Exp", "Cir"))
titles <- c("RMSE", "MAE", "ME")
titles <- paste(titles, "minmized based parameters selection")
measure_barplot(measures_frame, titles)

  


  
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