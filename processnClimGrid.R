# process nClimGrid data
# MAC 08/30/21
# update netcdf using downloadnClimGrid_netcdf.R

library(raster)

# set rasteroptions
rasterOptions(progress = 'text')

# load data into stack
temp <-stack("/scratch/crimmins/climgrid/netcdf/nclimgrid_tavg.nc")

# western US only
#e <- extent(-115.5, -102, 30.5, 39)
#prec <- crop(prec, e)	

# find last complete year
dates<-as.data.frame(names(temp))

dates<-seq(as.Date(strsplit(names(temp[[1]]), "X")[[1]][2],"%Y.%m.%d"),
           as.Date(strsplit(names(temp[[nlayers(temp)]]), "X")[[1]][2],"%Y.%m.%d"), by="month")

# subset to complete years
# temp<-subset(prec, which(dates<="2020-12-01"))

writeRaster(temp,filename="/scratch/crimmins/climgrid/processed/conus/CONUSmonthly_nClimGrid_prec_011895_092021.grd", overwrite=TRUE )