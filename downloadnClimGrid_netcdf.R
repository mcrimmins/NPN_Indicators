# download nClimGrid data
# MAC 08/30/21
# get data from https://www.ncei.noaa.gov/data/nclimgrid-monthly/access/
#
# nclimgrid_prcp.nc	2021-10-07 14:52	1.3G	 
# nclimgrid_tavg.nc	2021-10-07 14:52	972M	 
# nclimgrid_tmax.nc	2021-10-07 14:52	1.0G	 
# nclimgrid_tmin.nc	2021-10-07 14:52	1.0G	
# https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.ncdc:C00332

library(RCurl)
library(raster)

# get precip grid ~ 1 GB
URL <- "https://www.ncei.noaa.gov/data/nclimgrid-monthly/access/nclimgrid_prcp.nc"
download.file(URL, destfile = "/scratch/crimmins/climgrid/netcdf/nclimgrid_prcp.nc", method="curl")

# get tmin grid ~ 1 GB
URL <- "https://www.ncei.noaa.gov/data/nclimgrid-monthly/access/nclimgrid_tmin.nc"
download.file(URL, destfile = "/scratch/crimmins/climgrid/netcdf/nclimgrid_tmin.nc", method="curl")

# get tmax grid ~ 1 GB
URL <- "https://www.ncei.noaa.gov/data/nclimgrid-monthly/access/nclimgrid_tmax.nc"
download.file(URL, destfile = "/scratch/crimmins/climgrid/netcdf/nclimgrid_tmax.nc", method="curl")

# get tavg grid ~ 1 GB
URL <- "https://www.ncei.noaa.gov/data/nclimgrid-monthly/access/nclimgrid_tavg.nc"
download.file(URL, destfile = "/scratch/crimmins/climgrid/netcdf/nclimgrid_tavg.nc", method="curl")
