# Analyze NPN sites with nClimGrid data 
# MAC 10/27/21

library(raster)
library(dplyr)
library(ggplot2)

# set rasteroptions
rasterOptions(progress = 'text')

# functions
perc.rank <- function(x) trunc(rank(x))/length(x)

# load gridded data into stack
temp <-stack("/scratch/crimmins/climgrid/netcdf/nclimgrid_tavg.nc")
# load climate data from climTS_NPNsites.R
load("~/RProjects/NPNClimate/climDataNPN_MinMax.RData")
# set climate data type
  climData<-climDataMean

# load original data file
#npnData <- read.csv("spp-pp_7ormore_yrs.csv")
npnData <- read.csv("npnData_climate_analysis_120921-filtered.csv")
# trim to previous file format
npnData<-npnData[,2:45]

# only unique sites
#npnDataSites <- npnData[!duplicated(npnData[ , c("site_id")]), ]

# get sites into spdf
xy <- npnDataSites[,c(9,8)]
spdf <- SpatialPointsDataFrame(coords = xy, data = npnDataSites,
                               proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

# look for first_yes_doy
obsCols<-which(grepl('first_yes_doy', colnames(npnData))==TRUE)
firstYes<-strsplit(colnames(npnData)[which(grepl('first_yes_doy', colnames(npnData))==TRUE)],"\\.")
firstYes<- do.call(rbind.data.frame,firstYes)
colnames(firstYes)<-c("name","year")
firstYes$year<-as.numeric(as.character(firstYes$year))

# get observing period
yrsPOR<-apply(npnData[,obsCols],1,na.contiguous)
#yrsPOR<-do.call(rbind.data.frame,lapply(seq_along(yrs), function(i) attr(yrs[[i]], "tsp")[1:2]+2008))
yrsPOR<-lapply(seq_along(yrsPOR), function(i) attr(yrsPOR[[i]], "tsp")[1:2]+firstYes$year[1]-1)
# find empty
empty<-which(!lengths(yrsPOR))
  yrsPOR[empty]<-NA
# convert df and bind to sites
yrsPOR <- do.call(rbind.data.frame,yrsPOR)
colnames(yrsPOR) <- c("firstYr","lastYr")
  yrsPOR$firstYr[is.na(yrsPOR$firstYr)]<-firstYes$year[1]
  yrsPOR$lastYr[is.na(yrsPOR$lastYr)]<-firstYes$year[nrow(firstYes)]
npnData<-cbind.data.frame(npnData,yrsPOR)

npnData[npnData==-9999] <- NA

# filter obs based on prior no
noCols<-which(grepl('numdays_since_prior_no', colnames(npnData))==TRUE)
for(k in 1:length(obsCols)){
  npnData[,obsCols[k]]<-ifelse(npnData[,noCols[k]]>15 | is.na(npnData[,noCols[k]]), NA, npnData[,obsCols[k]])
}
# filter out plants with all missing based on prior nos
f1 <- function(x) length(which(!is.na(x)))
npnData$validYrs<-apply(npnData[,obsCols], MARGIN=1, FUN=f1)
npnData<-subset(npnData, validYrs>=6)


# extract climate data for NPN obs

# dates
dates<-names(temp)
dates<-as.Date(dates,format = "X%Y.%m.%d")

npnData$climMean<-NA
npnData$climSD<-NA
npnData$percSD<-NA
npnData$percMean<-NA
npnData$doyMean<-NA
npnData$doyMedian<-NA
npnData$doySD<-NA  
npnData$climCorr<-NA
npnData$doy_slope<-NA
npnData$doy_pval<-NA
npnData$clim_slope<-NA
npnData$clim_pval<-NA
npnData$meanMo<-NA

#Kendall::MannKendall(as.numeric(npnData[1,19:31]))[1]

for(i in 1:nrow(npnData)){
  # get data for each site from climate data list
  climTS<-climData[[which(npnDataSites$site_id==npnData$site_id[i])]]
  
  # get DOY stats
  npnData$doyMean[i]<-mean(as.numeric(npnData[i,obsCols]), na.rm=TRUE)
  npnData$doyMedian[i]<-median(as.numeric(npnData[i,obsCols]), na.rm=TRUE)
  npnData$doySD[i]<-sd(as.numeric(npnData[i,obsCols]), na.rm=TRUE)
  meanMo<-as.numeric(format(as.Date(paste0(round(npnData$doyMedian[i],0),"-2000"), format = "%j-%Y"),"%m"))
  
  npnData$meanMo[i]<-meanMo
  
  #  por<-13
  # get NPN obPeriod
  # por<-testData$yrs[i]
  por<-(npnData$lastYr[i]-npnData$firstYr[i])+1
  
  # create dataframe
  climTS<-cbind.data.frame(dates,climTS)
  climTS$month<-as.numeric(format(climTS$dates,"%m"))
  climTS$year<-as.numeric(format(climTS$dates, "%Y"))
  # create seasons
  #climTS$seas <- cut(climTS$month,c(0,3,6,9,12))
  #levels(climTS$seas) = c("JFM","AMJ","JAS","OND")
  
  # seasonal stats
  # seasStats<-climTS %>% group_by(seas, year) %>%
  #   summarize(climVar=mean(climTS, na.rm=TRUE)) %>%
  #   filter(month==meanMo)
  
  # roll mean for seasonal temp
  climTS$climTS<-zoo::rollapply(climTS$climTS, 3, mean, na.rm = TRUE, align="right", fill=NA)
  
  # seasonal stats
  seasStats<-subset(climTS, month==meanMo)
  
  # calc rolling stats
  #seasStats$rollSD<-(zoo::rollapply(seasStats$climVar, por, sd, na.rm = TRUE, align="right", fill=NA))
  #seasStats$rollMean<-(zoo::rollapply(seasStats$climVar, por, mean, na.rm = TRUE, align="right", fill=NA))
  # 
  seasStats$rollSD<-(zoo::rollapply(seasStats$climTS, por, sd, na.rm = TRUE, align="right", fill=NA))
  seasStats$rollMean<-(zoo::rollapply(seasStats$climTS, por, mean, na.rm = TRUE, align="right", fill=NA))
  # trim to non NA
  seasStats<-na.omit(seasStats)
  
  # perc ranks of rolling sd and means
  seasStats$percSD<-perc.rank(seasStats$rollSD)
  seasStats$percMean<-perc.rank(seasStats$rollMean)
  # look for empty ts
  
  length(seasStats$percSD[which(seasStats$year==npnData$lastYr[i])])
  
  if(length(seasStats$percSD[nrow(seasStats)])==0 | length(seasStats$percSD[which(seasStats$year==npnData$lastYr[i])])==0 ){
    npnData$climSD[i]<-NA
    npnData$climMean[i]<-NA
    npnData$percSD[i]<-NA
    npnData$percMean[i]<-NA
    npnData$climCorr[i]<-NA
    npnData$doy_slope[i]<-NA
    npnData$doy_pval[i]<-NA
    npnData$clim_slope[i]<-NA
    npnData$clim_pval[i]<-NA
  }else{
    #npnData$percSD[i]<-seasStats$percSD[nrow(seasStats)]
    #npnData$percMean[i]<-seasStats$percMean[nrow(seasStats)]
    npnData$percSD[i]<-seasStats$percSD[which(seasStats$year==npnData$lastYr[i])]
    npnData$percMean[i]<-seasStats$percMean[which(seasStats$year==npnData$lastYr[i])]
    npnData$climSD[i]<-seasStats$rollSD[which(seasStats$year==npnData$lastYr[i])]
    npnData$climMean[i]<-seasStats$rollMean[which(seasStats$year==npnData$lastYr[i])]
    npnData$climCorr[i]<-cor(as.numeric(npnData[i,obsCols]),
                             as.numeric(seasStats$climTS[(nrow(seasStats)-12):nrow(seasStats)]),
                             use="na.or.complete", method="spearman")
    # mk test on doy
    npnData$doy_slope[i]<-rkt::rkt(seq(min(npnData$firstYr),max(npnData$lastYr),1), as.numeric(npnData[i,obsCols]))$B
    npnData$doy_pval[i] <-rkt::rkt(seq(min(npnData$firstYr),max(npnData$lastYr),1), as.numeric(npnData[i,obsCols]))$sl
    # mk test on clim
    npnData$clim_slope[i]<-rkt::rkt(seq(min(npnData$firstYr),max(npnData$lastYr),1),
                                 as.numeric(seasStats$climTS[(nrow(seasStats)-12):nrow(seasStats)]))$B
    npnData$clim_pval[i]<-rkt::rkt(seq(min(npnData$firstYr),max(npnData$lastYr),1),
                                as.numeric(seasStats$climTS[(nrow(seasStats)-12):nrow(seasStats)]))$sl
  }
  
}

npnData$doyCOV<-(npnData$doySD/npnData$doyMean)*100
npnData$climCOV<-(npnData$climSD/abs(npnData$climMean))*100


#####
# Analyses and plots

# map data
library(ggplot2)
# load spatial data
states <- map_data("state")

# subset 
#subData<-subset(npnData,phenophase_description=="Breaking leaf buds")
# or all data
subData<-npnData
subData<-subset(subData, state!="AK")

# climate metrics at sites
#subDataSites<-subData[!duplicated(subData[ , c("site_id")]), ]

subDataSites<- subData %>% group_by(common_name,phenophase_description) %>%
                    mutate(n=n()) %>%
                    filter(n>=6)

subDataSites$phase_name<-paste(subDataSites$phenophase_description,subDataSites$common_name)
subDataSites$sig_doy<-ifelse(subDataSites$doy_pval<=0.05, "sig","non-sig")
subDataSites$sig_clim<-ifelse(subDataSites$clim_pval<=0.05, "sig","non-sig")

# facet phenophase, plot correlation
ggplot() +
  geom_polygon(data = states, aes(x = long, y = lat, group = group), fill=NA, color="black", size=0.1)  +
  #coord_fixed(xlim=c(out$meta$ll[1]-zoomLev, out$meta$ll[1]+zoomLev), ylim=c(out$meta$ll[2]-zoomLev, out$meta$ll[2]+zoomLev), ratio = 1) +
  #coord_fixed(xlim=c(-118, -102), ylim=c(28.75, 37.5), ratio = 1) +
  geom_point(data = subDataSites, aes(x = longitude, y = latitude, color=doy_slope))+
  scale_colour_gradient2(low = "blue",mid = "white", high = "red")+
  facet_wrap(.~phase_name)

ggplot(subDataSites, aes(y=climCorr))+
  geom_boxplot(varwidth = TRUE)+
  facet_wrap(.~phase_name)

ggplot() +
  geom_polygon(data = states, aes(x = long, y = lat, group = group), fill=NA, color="black", size=0.1)  +
  #coord_fixed(xlim=c(out$meta$ll[1]-zoomLev, out$meta$ll[1]+zoomLev), ylim=c(out$meta$ll[2]-zoomLev, out$meta$ll[2]+zoomLev), ratio = 1) +
  #coord_fixed(xlim=c(-118, -102), ylim=c(28.75, 37.5), ratio = 1) +
  geom_point(data = subDataSites, aes(x = longitude, y = latitude, color=clim_slope, shape=sig_clim), size=4)+
  scale_colour_gradient2(low = "blue",mid = "white", high = "red",limits=c(-0.1, 0.1), oob=scales::squish)+
  scale_shape_manual(values=c(4,20))+
  facet_wrap(.~phase_name)

ggplot() +
  geom_polygon(data = states, aes(x = long, y = lat, group = group), fill=NA, color="black", size=0.1)  +
  #coord_fixed(xlim=c(out$meta$ll[1]-zoomLev, out$meta$ll[1]+zoomLev), ylim=c(out$meta$ll[2]-zoomLev, out$meta$ll[2]+zoomLev), ratio = 1) +
  #coord_fixed(xlim=c(-118, -102), ylim=c(28.75, 37.5), ratio = 1) +
  geom_point(data = subDataSites, aes(x = longitude, y = latitude, color=percMean, size=percSD, alpha=0.5))+
  #scale_colour_gradient(low = "yellow", high = "red")+
  scale_colour_gradient2(low = "blue",mid = "white", high = "red", midpoint = 0.5)+
  facet_wrap(.~phase_name)

ggplot() +
  geom_polygon(data = states, aes(x = long, y = lat, group = group), fill=NA, color="black", size=0.1)  +
  #coord_fixed(xlim=c(out$meta$ll[1]-zoomLev, out$meta$ll[1]+zoomLev), ylim=c(out$meta$ll[2]-zoomLev, out$meta$ll[2]+zoomLev), ratio = 1) +
  #coord_fixed(xlim=c(-118, -102), ylim=c(28.75, 37.5), ratio = 1) +
  geom_point(data = subDataSites, aes(x = longitude, y = latitude, color=meanMo))+
  #scale_colour_gradient(low = "yellow", high = "red")+
  scale_colour_gradient2(low = "blue",mid = "green", high = "red")+
  facet_wrap(.~phase_name)

# relationships between metrics
plot(subData$percSD,subData$doyCOV)
cor(subData$percSD,subData$doyCOV, use = "na.or.complete")

library("PerformanceAnalytics")
my_data <- subData[,47:55]
chart.Correlation(my_data, histogram=TRUE, pch=19, method = "spearman")

ggplot(subData, aes(percSD,doyCOV, color=climCorr))+
  geom_point()+
  scale_colour_gradient2(low = "blue",mid = "grey", high = "red", midpoint = 0)

ggplot(npnData, aes(climCorr))+
  geom_histogram()+
  facet_wrap(.~phenophase_description)

ggplot() +
  geom_polygon(data = states, aes(x = long, y = lat, group = group), fill=NA, color="black", size=0.1)  +
  #coord_fixed(xlim=c(out$meta$ll[1]-zoomLev, out$meta$ll[1]+zoomLev), ylim=c(out$meta$ll[2]-zoomLev, out$meta$ll[2]+zoomLev), ratio = 1) +
  #coord_fixed(xlim=c(-118, -102), ylim=c(28.75, 37.5), ratio = 1) +
  geom_point(data = subData, aes(x = longitude, y = latitude, color=climCorr, size=doyCOV))+
  scale_colour_gradient(low = "yellow", high = "red")

ggplot() +
  geom_polygon(data = states, aes(x = long, y = lat, group = group), fill=NA, color="black", size=0.1)  +
  #coord_fixed(xlim=c(out$meta$ll[1]-zoomLev, out$meta$ll[1]+zoomLev), ylim=c(out$meta$ll[2]-zoomLev, out$meta$ll[2]+zoomLev), ratio = 1) +
  #coord_fixed(xlim=c(-118, -102), ylim=c(28.75, 37.5), ratio = 1) +
  geom_point(data = subData, aes(x = longitude, y = latitude, color=climCorr, size=doySD))+
  scale_colour_gradient2(low = "blue",mid = "grey", high = "red", midpoint = 0)



# write out data file
write.csv(npnData, file="npnData_climate_analysis_120921.csv")

write.csv(subDataSites, file="subDataSites_climate_analysis_121721.csv")
