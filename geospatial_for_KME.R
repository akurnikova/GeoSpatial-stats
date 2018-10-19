## Created by
## A. Kurnikova
## (2016)

## First read in the file
dat <- read.csv(file="~/Downloads/Site5_Conduct_1500.csv",head=TRUE, sep = ",")

## Import (need to install first) necessary libraries
library(ggplot2) # For plotting
library(reshape) # For plotting

## by-hand gridding
sY <- seq(min(dat$Y),max(dat$Y), by=0.00015)
sX <- seq(min(dat$X),max(dat$X), by=0.00015)
Z <- dat$Conductivity.500TIME
X <- dat$X
Y <- dat$Y

grid.mean <- matrix(data = NA, nrow = (NROW(sX)-1), ncol = (NROW(sY)-1))
grid.sd <- matrix(data = NA, nrow = (NROW(sX)-1), ncol = (NROW(sY)-1))

for(i in 1:(NROW(sY)-1)){
  iy1 = which(Y<sY[i+1])
  iy2 = which(Y>=sY[i])
  iy = intersect(iy1,iy2)
  for(j in 1:(NROW(sX)-1)){
    ix1 <- which(X<sX[j+1])
    ix2 <- which(X>=sX[j])
    ix <- intersect(ix1,ix2)
    grid.mean[j,i] <- mean(Z[intersect(ix,iy)])
    grid.sd[j,i] <- sd(Z[intersect(ix,iy)])
  }
}

dg <- as.data.frame(grid.mean)

colnames(dg) <- sY[-1]
dg["xvals"]  <- sX[-1]
dg.melt <- melt(dg,id.vars = "xvals",variable.name="z")

colnames(dg.melt) <- c("X","Y","Z")
## Moran's I - by grid
library(ape)
moransI.dists <- as.matrix(dist(cbind(dg.melt$X, dg.melt$Y)))
moransI.dists.inv <- 1/moransI.dists
diag(moransI.dists.inv) <- 0
moransI.dists.inv[is.infinite(moransI.dists.inv)] <- 0
Moran.I(dg.melt$Z, moransI.dists.inv,na.rm = TRUE)

## Plot!!!!
ggplot(dg.melt ,aes(x = xvals,y = variable,z= value)) + geom_tile(aes(fill=dg.melt$value))


## Try a variogram
library(geoR)
dists <- dist(dat[,1:2])
summary(dists)
breaks = seq(0, 0.0025, l = 20)
v1 <- variog(coords = dat[,1:2], data = dat[,15], breaks = breaks)
plot(v1, type = "b") 


## Moran's I - by raw points
library(ape)
moransI.dists <- as.matrix(dist(cbind(dat$X, dat$Y)))
moransI.dists.inv <- 1/moransI.dists
diag(moransI.dists.inv) <- 0
moransI.dists.inv[is.infinite(moransI.dists.inv)] <- 0
Moran.I(dat$Conductivity.500TIME, moransI.dists.inv,na.rm = TRUE)


SourceX <- -101.199722
SourceY <- 48.808333
Source2X <- -101.199444
Source2Y <- 48.811111
dat_conv <-geoXY(dat$Y, dat$X)
Source1_conv <- geoXY(SourceY, SourceX, min(dat$Y), min(dat$X))
Source2_conv <- geoXY(Source2Y, Source2X, min(dat$Y), min(dat$X))

dat$StationDist <- sqrt((dat_conv[,1]-Source1_conv[,1])^2+(dat_conv[,2]-Source1_conv[,2])^2)
dat$StationDist2 <- sqrt((dat_conv[,1]-Source2_conv[,1])^2+(dat_conv[,2]-Source2_conv[,2])^2)
dat$CombinedDist <- 1/(1/(dat$StationDist)+(1/(dat$StationDist2)))

  plot(dat$StationDist,  dat$Conductivity.500TIME)
points(dat$StationDist2,  dat$Conductivity.500TIME,col = 2)
points(dat$CombinedDist,  dat$Conductivity.500TIME,col = 3)


##
library("PASWR")
Distance_Threshold <- 0.001
ECA_elevation_threshold_test <- 60
SIGN.test(dat$Conductivity.500TIME[dat$CombinedDist>Distance_Threshold], md = ECA_elevation_threshold_test, alternative = "less", conf.level = 0.95)
SIGN.test(dat$Conductivity.500TIME[dat$CombinedDist<Distance_Threshold], md = ECA_elevation_threshold_test, alternative = "less", conf.level = 0.95)



Distance_Threshold <- 0.001
## Trying to exclude the values around the station
sY <- seq(min(dat$Y),max(dat$Y), by=0.00015)
sX <- seq(min(dat$X),max(dat$X), by=0.00015)
Z <- dat$Conductivity.500TIME[dat$CombinedDist>Distance_Threshold]
X <- dat$X[dat$CombinedDist>Distance_Threshold]
Y <- dat$Y[dat$CombinedDist>Distance_Threshold]

grid.mean <- matrix(data = NA, nrow = (NROW(sX)-1), ncol = (NROW(sY)-1))
grid.sd <- matrix(data = NA, nrow = (NROW(sX)-1), ncol = (NROW(sY)-1))

for(i in 1:(NROW(sY)-1)){
  iy1 = which(Y<sY[i+1])
  iy2 = which(Y>=sY[i])
  iy = intersect(iy1,iy2)
  for(j in 1:(NROW(sX)-1)){
    ix1 <- which(X<sX[j+1])
    ix2 <- which(X>=sX[j])
    ix <- intersect(ix1,ix2)
    grid.mean[j,i] <- mean(Z[intersect(ix,iy)])
    grid.sd[j,i] <- sd(Z[intersect(ix,iy)])
  }
}

dg <- as.data.frame(grid.mean)

colnames(dg) <- sY[-1]
dg["xvals"]  <- sX[-1]
dg.melt <- melt(dg,id.vars = "xvals",variable.name="z")

colnames(dg.melt) <- c("X","Y","Z")
## Moran's I - by grid
library(ape)
moransI.dists <- as.matrix(dist(cbind(dg.melt$X, dg.melt$Y)))
moransI.dists.inv <- 1/moransI.dists
diag(moransI.dists.inv) <- 0
moransI.dists.inv[is.infinite(moransI.dists.inv)] <- 0
Moran.I(dg.melt$Z, moransI.dists.inv,na.rm = TRUE)
