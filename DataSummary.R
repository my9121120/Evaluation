# library("gstat")
# library(sp)
# library(rgdal)
# library(ggplot2)
# library(ggmap)
# library(maps)
library(e1071)
library(stats)
library(graphics)
library(lattice)
#rm(list = ls())
#a = read.csv("DQpropertyall.csv", numerals = "no.loss")
#source("multicrossvalidation.R")
load("DQData")
DQCluster = DQpropertyall[DQpropertyall$X < 517600,]
DQCluster[, "lnOM"] = log(DQCluster$OM.g.kg.1.)
statistics <- summary(DQCluster[, c("pH", "OM.g.kg.1.", "lnOM")])

sd_pH <- sd(DQCluster$pH)
skew_pH <- skewness(DQCluster$pH, type = 1)
kurt_pH <- kurtosis(DQCluster$pH, type = 1)

sd_SOM <- sd(DQCluster$OM.g.kg.1.)
skew_SOM <- skewness(DQCluster$OM.g.kg.1., type = 1)
kurt_SOM <- kurtosis(DQCluster$OM.g.kg.1., type = 1)

sd_lnSOM <- sd(DQCluster$lnOM)
skew_lnSOM <- skewness(DQCluster$lnOM, type = 1)
kurt_lnSOM <- kurtosis(DQCluster$lnOM, type = 1)

stat_data <- data.frame("pH"=round(c(sd_pH,skew_pH,kurt_pH),digits = 2),
                        "SOM"=round(c(sd_SOM,skew_SOM,kurt_SOM),digits = 2),
                        "lnOM"=round(c(sd_lnSOM,skew_lnSOM,kurt_lnSOM),digits = 2))
rownames(stat_data) <- c("sd", "skew", "kurt")

plot(DQCluster$pH, ylim=c(4,10))
plot(DQCluster$OM.g.kg.1., ylim=c(9,40))

par(mfrow = c(1,3))
pH <- DQCluster$pH
mean_pH <- mean(pH)
hist(pH, xlab = "pH", ylab = "Density", freq = FALSE,
     main = "(a)", font.main = 1)
curve(dnorm(x, mean = mean_pH, sd = sd_pH), 
      add = TRUE, col = "darkblue")


OM <- DQCluster$OM.g.kg.1.
hist(OM, xlab = "OM", ylab = "Density",
     freq = FALSE, main = "(b)", font.main = 1)
curve(dnorm(x, mean = mean(OM), sd = sd_SOM),
      add = TRUE, col = "darkblue")

lnOM <- DQCluster$lnOM
hist(lnOM, xlab = "ln(OM)", ylab = "Density", 
     freq = FALSE, main = "(c)", font.main = 1)
curve(dnorm(x, mean = mean(lnOM), sd = sd_lnSOM), 
      add = TRUE, col = "darkblue")
