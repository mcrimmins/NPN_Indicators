# nClimGrid seasonal averages
# MAC 12/20/21

library(raster)
#library(dplyr)
#library(ggplot2)

# set rasteroptions
rasterOptions(progress = 'text')

# functions
perc.rank <- function(x) trunc(rank(x))/length(x)

# load data into stack
temp <-stack("/scratch/crimmins/climgrid/netcdf/nclimgrid_tavg.nc")

# find last complete year
dates<-as.data.frame(names(temp))

dates<-seq(as.Date(strsplit(names(temp[[1]]), "X")[[1]][2],"%Y.%m.%d"),
           as.Date(strsplit(names(temp[[nlayers(temp)]]), "X")[[1]][2],"%Y.%m.%d"), by="month")

# subset data to shorter period like NN period 2009-2021
temp<-subset(temp, which(dates>="2009-01-01"))
dates<-dates[which(dates>="2009-01-01")]

# 3-month moving average
#x <- calc(st1, function(x) movingFun(x, 3, mean))
win <- 2
fMov <- function(x) movingFun(x, win, mean, type = 'to')
beginCluster(7)
  temp3mo <- clusterR(temp, calc, args=list(fun=fMov), export='win')
endCluster()

names(temp3mo)<-dates

# get climo of 3mo averages
dates<-as.data.frame(dates)
dates$month<-as.numeric(format(dates$dates, "%m"))
temp3mo_avg<- stackApply(temp3mo, dates$month, fun = mean)

# get anom
fun <- function(x, y) {
  x - y
}

temp3mo_anom <- overlay(x = temp3mo, y = temp3mo_avg, fun = fun)
names(temp3mo_anom)<-dates$dates

writeRaster(temp3mo_anom,filename="/scratch/crimmins/climgrid/processed/conus/CONUSmonthly_nClimGrid_2mo_tmean_anom_011895_092021.grd", overwrite=TRUE )
