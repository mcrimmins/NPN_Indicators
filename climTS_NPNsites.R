# nClimGrid point extract for NPN sites
# MAC 10/26/21

library(raster)
library(dplyr)
library(ggplot2)

# set rasteroptions
rasterOptions(progress = 'text')

# functions
perc.rank <- function(x) trunc(rank(x))/length(x)

# load data into stack
temp1 <-stack("/scratch/crimmins/climgrid/netcdf/nclimgrid_tmin.nc")
temp2 <-stack("/scratch/crimmins/climgrid/netcdf/nclimgrid_tmax.nc")

#library(readxl)
npnData <- read.csv("spp-pp_7ormore_yrs.csv")

# only unique sites
npnDataSites <- npnData[!duplicated(npnData[ , c("site_id")]), ]

# get sites into spdf
xy <- npnDataSites[,c(9,8)]
spdf <- SpatialPointsDataFrame(coords = xy, data = npnDataSites,
                               proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

# raster getValues approach https://newbedev.com/increasing-speed-of-crop-mask-extract-raster-by-many-polygons-in-r
library(foreach)
library(doParallel)

registerDoParallel(6)

climDataMin<-list()
climDataMax<-list()

for(i in 1:nrow(npnDataSites)){
  cell<-cellFromXY(temp1,spdf[i,])
  # check for empty locations
  if(is.na(cell)==TRUE){
    climDataMin[[i]]<-NA
    climDataMax[[i]]<-NA
      }else{
        r <- rasterFromCells(temp1, cell,values=F)
        climDataMin[[i]] <- foreach(j = 1:dim(temp1)[3],.packages='raster',.combine=rbind,.inorder=T) %dopar% {
          #get value and store
          getValues(crop(temp1[[j]],r))
        }  
        r <- rasterFromCells(temp2, cell,values=F)
        climDataMax[[i]] <- foreach(j = 1:dim(temp2)[3],.packages='raster',.combine=rbind,.inorder=T) %dopar% {
            #get value and store
            getValues(crop(temp2[[j]],r))
        }
      }
  print(i)
  # add error catch with lat/lon outside of nClimGrid
}
#proc.time() - ptm
endCluster()

# calculate monthly mean
climDataMean<-list()
climDataMean<-lapply(seq_along(climDataMin), function(i) rowMeans(data.frame((climDataMin[[i]]),(climDataMax[[i]])), na.rm = TRUE))

save(climDataMin,climDataMax,climDataMean,npnDataSites, file = "climDataNPN_MinMax.RData")

