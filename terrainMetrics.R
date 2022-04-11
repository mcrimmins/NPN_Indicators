# develop terrain metrics from topo data for NPN analysis
# MAC 01/19/22

library(raster)

# get elevation data 30arc seconds ~ 1km
elevation <- raster::getData('alt', country='USA')
  elevation<-elevation[[1]]
TPI <- raster::terrain(elevation, opt='TPI')
TRI <- raster::terrain(elevation, opt='TRI')
slope<- raster::terrain(elevation, opt='slope')
aspect<- raster::terrain(elevation, opt='aspect')

# larger TPI neighborhood
# TPI for different neighborhood size:
tpiw <- function(x, w=5) {
  m <- matrix(1/(w^2-1), nc=w, nr=w)
  m[ceiling(0.5 * length(m))] <- 0
  f <- focal(x, m)
  x - f
}
TPI5 <- tpiw(elevation, w=5)
# }

anthrome<-raster("anthromes2015AD.asc")

save.image(file="terrainData.RData")

# 11	11: Urban 
# 12	12: Mixed settlements 
# 21	21: Rice villages 
# 22	22: Irrigated villages 
# 23	23: Rainfed villages 
# 24	24: Pastoral villages 
# 31	31: Residential irrigated croplands 
# 32	32: Residential rainfed croplands 
# 33	33: Populated croplands 
# 34	34: Remote croplands 
# 41	41: Residential rangelands 
# 42	42: Populated rangelands 
# 43	43: Remote rangelands 
# 51	51: Residential woodlands 
# 52	52: Populated woodlands 
# 53	53: Remote woodlands 
# 54	54: Inhabited treeless & barren lands 
# 61	61: Wild woodlands
# 62	62: Wild treeless and barren lands
# 63	63: Ice, uninhabited 
# 70	NODATA
