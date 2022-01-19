# Site-level lat model
# adapted from lat_model.R
# MAC 01/12/22
# Pull data from NN
#
library(rnpn)
library(dplyr)
library(plyr)


NNObs <- npn_download_site_phenometrics(request_source = "TCrimmins", c(2009:2021),
                                        species_id = c(12), phenophase_id = c(501),
                                        additional_fields = c("species_functional_type", "species_category", "Observed_Status_Conflict_Flag",
                                                              "site_name", "Partner_Group", "Dataset_ID"))
# species_id = c(12,79), phenophase_id = c(483,501),

# remove >1st instances of "yes" on the same individual in the same year
leaf_data_1yes <- NNObs %>%
  group_by(species, site_id, mean_first_yes_year) %>%
  filter(mean_first_yes_doy == min(mean_first_yes_doy)) %>%
  slice(1) %>% # takes the first occurrence if there is a tie
  ungroup()
# remove any observations with no "no" observations within 30 days prior to a "yes" observation
leaf_data_1yes_prevno <- leaf_data_1yes[which(leaf_data_1yes$mean_numdays_since_prior_no >=0),]
leaf_data_1yes_prevno <- leaf_data_1yes_prevno[which(leaf_data_1yes_prevno$mean_numdays_since_prior_no <=15),]

NNObs2<-leaf_data_1yes_prevno[which(leaf_data_1yes_prevno$mean_first_yes_doy <=180),]

NNObs2<-subset(NNObs2,longitude<0)
NNObs2<-subset(NNObs2,longitude>-100)
NNObs2<-subset(NNObs2,latitude<=50)

# compare beech and flowering dogwood
    # NNObs2$name_pheno<-paste0(NNObs2$common_name,"-",NNObs2$phenophase_description)
    # temp<-subset(NNObs2, name_pheno %in% c("American beech-Leaves","flowering dogwood-Open flowers"))
    # 
    # ggplot(NNObs2, aes(latitude,mean_first_yes_doy, color=common_name))+
    #   geom_point(size=1, alpha=0.4)+
    #   facet_wrap(.~mean_first_yes_year)+
    #   stat_smooth(method = lm)+
    #   theme_bw()

# lat/elev model
latModel <- lm(mean_first_yes_doy ~ latitude+elevation_in_meters, data = NNObs2)
summary(latModel)

# doy~lat plot
library(ggplot2)
ggplot(NNObs2)+
  geom_point(aes(latitude, mean_first_yes_doy, color=as.factor(mean_first_yes_year)))+
  geom_smooth(data=NNObs2,aes(latitude, mean_first_yes_doy), method=lm, se=TRUE)+
  ggtitle(paste0(NNObs2$common_name[1],"-",NNObs2$phenophase_description[1]))+
  ggpubr::stat_cor(aes(latitude, mean_first_yes_doy),method = "spearman")+
  theme_bw()


NNObs2$latResids<-latModel$residuals
NNObs2$year<-NNObs2$mean_first_yes_year

ggplot(NNObs2)+
  geom_point(aes(latitude, latResids, color=as.factor(mean_first_yes_year)))+
  geom_smooth(data=NNObs2,aes(latitude, latResids), method=lm, se=TRUE)+
  ggtitle(paste0(NNObs2$common_name[1],"-",NNObs2$phenophase_description[1]))+
  theme_bw()


# red maple
#latModel <- lm(mean_first_yes_doy ~ latitude, data = NNObs2)
#summary(latModel)

# lm vs quantreg residual analysis
latModel_lm <- lm(mean_first_yes_doy ~ latitude+elevation_in_meters, data = NNObs2)
summary(latModel_lm)
latModel_rq <- quantreg::rq(mean_first_yes_doy ~ latitude+elevation_in_meters, data = NNObs2)
summary(latModel_rq)
  temp1<-as.data.frame(latModel_lm$residuals)
    colnames(temp1)<-"resid"
    temp1$model<-"lm"
  temp2<-as.data.frame(latModel_rq$residuals)
    colnames(temp2)<-"resid"
    temp2$model<-"rq"
  temp1<-rbind.data.frame(temp1,temp2)  
ggplot(temp1, aes(model, resid))+
  geom_boxplot()
  
  
#latModel <- lm(first_yes_doy ~ latitude, data = NNObs2)
# median regression
#latModel <- quantreg::rq(first_yes_doy ~ latitude+elevation_in_meters, data = NNObs2)
#summary(latModel)
# use model to find first day of spring/last between US lat range
predDate<-predict(latModel, newdata=data.frame(latitude = c(30,50), elevation_in_meters = c(50,50)))
as.Date(predDate, origin = "2016-01-01")


NNObs2$latResids<-latModel$residuals
NNObs2$year<-NNObs2$mean_first_yes_year

#screen residuals
hist(NNObs2$latResids)
#NNObs2<-subset(NNObs2, latResids<=30 & latResids>=-30)

# add in month
NNObs2$month<-as.numeric(format(as.Date(paste0(NNObs2$mean_first_yes_doy,"-",NNObs2$mean_first_yes_year),format="%j-%Y"),"%m"))

#####
library(raster)
#library(rasterVis)
library(ggplot2)
# ALL YEARS TOGETHER
tmean_anom_full<-stack("/scratch/crimmins/climgrid/processed/conus/CONUSmonthly_nClimGrid_2mo_tmean_anom_011895_092021.grd")
# make 2009-2021 anom grid
tmean_anom_npn<-stack("/scratch/crimmins/climgrid/processed/conus/CONUSmonthly_nClimGrid_2mo_tmean_anom_012009_092021.grd")

# find last complete year -- FULL
dates_full<-as.data.frame(names(tmean_anom_full))
dates_full<-as.data.frame(seq(as.Date(strsplit(names(tmean_anom_full[[1]]), "X")[[1]][2],"%Y.%m.%d"),
                              as.Date(strsplit(names(tmean_anom_full[[nlayers(tmean_anom_full)]]), "X")[[1]][2],"%Y.%m.%d"), by="month"))
colnames(dates_full)<-"date"
dates_full$month<-as.numeric(format(dates_full$date, "%m"))
dates_full$year<-as.numeric(format(dates_full$date, "%Y"))

# find last complete year -- NPN PoR
dates_npn<-as.data.frame(names(tmean_anom_npn))
dates_npn<-as.data.frame(seq(as.Date(strsplit(names(tmean_anom_npn[[1]]), "X")[[1]][2],"%Y.%m.%d"),
                             as.Date(strsplit(names(tmean_anom_npn[[nlayers(tmean_anom_npn)]]), "X")[[1]][2],"%Y.%m.%d"), by="month"))
colnames(dates_npn)<-"date"
dates_npn$month<-as.numeric(format(dates_npn$date, "%m"))
dates_npn$year<-as.numeric(format(dates_npn$date, "%Y"))

#####
# # static month period
# NNObsClim <- list()
# mo=5 # end month spring 3-mo period
# 
# for(i in min(NNObs2$mean_first_yes_year):max(NNObs2$mean_first_yes_year)){
#   
#   # subset raster for each -- full 
#   subAnom<-tmean_anom_full[[which(dates_full$month==mo & dates_full$year==i)]]
#   tempObs<-subset(NNObs2, year==i)
#   temp_full<-extract(subAnom[[1]], SpatialPoints(tempObs[,c("longitude","latitude")]), sp = T)
#   colnames(temp_full@data)<-"climAnom_full"
#   
#   # subset raster for each -- npn por 
#   subAnom<-tmean_anom_npn[[which(dates_npn$month==mo & dates_npn$year==i)]]
#   tempObs<-subset(NNObs2, year==i)
#   temp_npn<-extract(subAnom[[1]], SpatialPoints(tempObs[,c("longitude","latitude")]), sp = T)
#   colnames(temp_npn@data)<-"climAnom_npn"
#   
#   tempObs<-cbind.data.frame(tempObs, temp_npn@data, temp_full@data)
#   
#   NNObsClim[[i]] <- tempObs # add it to your list
#   
# }
# ######

######
# self adjust clim window

NNObsClim <- list()
# end month spring 3-mo period

for(i in min(NNObs2$mean_first_yes_year):max(NNObs2$mean_first_yes_year)){
  
  tempObs<-subset(NNObs2, year==i)
  # full extract
  subAnom<-tmean_anom_full[[which(dates_full$year==i)]]
  temp_full<-extract(subAnom[[tempObs$month]], SpatialPoints(tempObs[,c("longitude","latitude")]), sp = T)
  temp_full<-as.data.frame(diag(as.matrix(temp_full@data)))
  colnames(temp_full)<-"climAnom_full"
            
  # npn extract -- npn por 
  subAnom<-tmean_anom_npn[[which(dates_npn$year==i)]]
  temp_npn<-extract(subAnom[[tempObs$month]], SpatialPoints(tempObs[,c("longitude","latitude")]), sp = T)
  temp_npn<-as.data.frame(diag(as.matrix(temp_npn@data)))
  colnames(temp_npn)<-"climAnom_npn"
  
  tempObs<-cbind.data.frame(tempObs, temp_npn, temp_full)
  
  NNObsClim[[i]] <- tempObs # add it to your list
  print(i)
}
#####



NNObsClim <- do.call(rbind, NNObsClim)
NNObsClim$anomType_full<-ifelse(NNObsClim$climAnom_full<0,"neg","pos")
NNObsClim$anomType_npn<-ifelse(NNObsClim$climAnom_npn<0,"neg","pos")

cor(NNObsClim$latResids, NNObsClim$climAnom_npn, "na.or.complete", method="spearman")
cor(NNObsClim$latResids, NNObsClim$climAnom_full, "na.or.complete", method="spearman")

states <- map_data("state")

ggplot(NNObsClim) + 
  geom_polygon(data = states, aes(x = long, y = lat, group = group), fill=NA, color="black", size=0.1)  +
  geom_jitter(data = NNObsClim, aes(x = longitude, y = latitude, color=latResids),size = 1, width=1, height=1, alpha=0.8)+ # shape=anomType_npn
  #geom_point(data = NNObsClim, aes(x = longitude, y = latitude, color=latResids),size = 1, alpha=0.8)+ # shape=anomType_npn
  scale_color_gradient2(low = "red",mid = "grey", high = "blue", limits=c(-10, 10), oob=scales::squish)+
  #scale_shape_manual(values=c(16, 3))+
  facet_wrap(~ year) +
  coord_equal( xlim = c(-100, -60), ylim = c(25, 50))+
  theme_bw()+
  #ggtitle(paste0(NNObs2$common_name[1],"-",NNObs2$phenophase_description[1],": model resids vs local 2-mo clim anom (2009-2021 base)"))
  ggtitle(paste0(NNObs2$common_name[1],"-",NNObs2$phenophase_description[1],": model resids"))

# subset a month
#test<-subset(NNObsClim, mean_first_yes_year==2017)
# ggplot(subset(NNObsClim, mean_first_yes_year==2017))+
#   geom_polygon(data = states, aes(x = long, y = lat, group = group), fill=NA, color="black", size=0.1)  +
#   #geom_jitter(data = NNObsClim, aes(x = longitude, y = latitude, color=latResids),size = 1, width=1, height=1, alpha=0.8)+ # shape=anomType_npn
#   geom_point(data = subset(NNObsClim, mean_first_yes_year==2017), aes(x = longitude, y = latitude, color=latResids),size = 1, alpha=0.8)+ # shape=anomType_npn
#   scale_color_gradient2(low = "red",mid = "grey", high = "blue", limits=c(-20, 20), oob=scales::squish)+
#   coord_equal( xlim = c(-100, -60), ylim = c(25, 50))+
#   theme_bw()
  

ggplot(NNObsClim)+
  geom_point(aes(latResids, climAnom_npn, color=as.factor(mean_first_yes_year)))+
  geom_smooth(data=NNObsClim,aes(latResids,climAnom_npn),method=lm,se=FALSE)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept=0)+
  ggtitle(paste0(NNObs2$common_name[1],"-",NNObs2$phenophase_description[1]))+
  ggpubr::stat_cor(aes(latResids,climAnom_npn),method = "spearman")+
  theme_bw()

ggplot(NNObsClim)+
  geom_point(aes(climAnom_npn, climAnom_full, color=as.factor(mean_first_yes_year)))+
  geom_smooth(data=NNObsClim,aes(climAnom_npn,climAnom_full),method=lm,se=TRUE)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept=0)+
  ggtitle(paste0(NNObs2$common_name[1],"-",NNObs2$phenophase_description[1],": climate period comparison"))+
  ggpubr::stat_cor(aes(climAnom_npn,climAnom_full),method = "spearman")+
  theme_bw()+
  geom_abline(color = "grey")

ggplot(NNObsClim)+
  geom_point(aes(latResids, climAnom_npn, color=climAnom_full))+
  geom_smooth(data=NNObsClim,aes(latResids,climAnom_npn),method=lm,se=FALSE)+
  scale_color_gradient2(low = "blue",mid = "grey", high = "red", midpoint=0)+ # limits=c(-20, 20), oob=scales::squish
  geom_hline(yintercept = 0)+
  geom_vline(xintercept=0)+
  ggtitle(paste0(NNObs2$common_name[1],"-",NNObs2$phenophase_description[1]))+
  ggpubr::stat_cor(aes(latResids,climAnom_npn),method = "spearman")+
  theme_bw()


hist(NNObsClim$mean_first_yes_doy)

ggplot(NNObsClim)+
  geom_point(aes(latResids, climAnom_npn, color=climAnom_full))+
  geom_smooth(data=NNObsClim,aes(latResids,climAnom_npn),method=lm,se=FALSE)+
  scale_color_gradient2(low = "blue",mid = "grey", high = "red", midpoint = 0)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept=0)+
  facet_wrap(.~mean_first_yes_year)+
  ggtitle(paste0(NNObs2$common_name[1],"-",NNObs2$phenophase_description[1]))+
  ggpubr::stat_cor(aes(latResids,climAnom_npn),method = "spearman")+
  theme_bw()

ggplot(NNObsClim)+
  geom_point(aes(latResids, latitude, color=climAnom_npn))+
  #geom_smooth(data=NNObsClim,aes(latResids,climAnom_npn),method=lm,se=FALSE)+
  scale_color_gradient2(low = "blue",mid = "grey", high = "red", midpoint = 0)+
  #geom_hline(yintercept = 0)+
  geom_vline(xintercept=0)+
  facet_wrap(.~mean_first_yes_year)+
  ggtitle(paste0(NNObs2$common_name[1],"-",NNObs2$phenophase_description[1]))+
  #ggpubr::stat_cor(aes(latResids,climAnom_npn),method = "spearman")+
  theme_bw()

hist(NNObsClim$latResids)

ggplot(NNObsClim)+
  scale_color_gradient2(low = "blue",mid = "grey", high = "red", midpoint = 0)+
  geom_point(aes(latResids,mean_first_yes_doy, color=climAnom_npn))+
  facet_wrap(.~mean_first_yes_year)

# annual stats on residuals

annClim<-NNObsClim %>% group_by(mean_first_yes_year) %>% 
                        dplyr::summarize(count=n(),
                                  meanResid=mean(latResids),
                                  sdResid=sd(latResids),
                                  meanclimAnom_npn=mean(climAnom_npn, na.rm=TRUE),
                                  meanclimAnom_full=mean(climAnom_full, na.rm=TRUE))

ggplot(annClim, aes(mean_first_yes_year, meanResid))+
    geom_line()+
    geom_ribbon(data=annClim, mapping=aes(x=mean_first_yes_year, ymin=meanResid-sdResid, ymax=meanResid+sdResid),alpha=0.4)+
  geom_hline(yintercept = 0)+
  theme_bw()

temp<-tidyr::gather(annClim, anomType, value, -mean_first_yes_year, -count, -meanResid,-sdResid)
ggplot(temp, aes(mean_first_yes_year, value,color=anomType))+
  geom_line()+
  #geom_ribbon(data=annClim, mapping=aes(x=mean_first_yes_year, ymin=meanResid-sdResid, ymax=meanResid+sdResid),alpha=0.4)+
  geom_hline(yintercept = 0)+
  theme_bw()

plot(annClim$meanResid, annClim$meanclimAnom_npn)
cor(annClim$meanResid, annClim$meanclimAnom_npn, method = "spearman")

plot(annClim$meanResid, annClim$meanclimAnom_full)
cor(annClim$meanResid, annClim$meanclimAnom_full)

ggplot(annClim, aes(mean_first_yes_year,count))+
  geom_bar(stat='identity')


# modeling lat, elev and local clim
climModel <- lm(mean_first_yes_doy ~ latitude+elevation_in_meters, data = NNObsClim)
summary(climModel)
climModel <- lm(mean_first_yes_doy ~ latitude+elevation_in_meters+climAnom_npn, data = NNObsClim)
summary(climModel)

jtools::summ(climModel, scale=TRUE)
jtools::summ(climModel, confint = TRUE, digits = 3)
jtools::effect_plot(climModel, pred = climAnom_npn, interval = TRUE, plot.points = TRUE)
jtools::plot_summs(climModel, scale = TRUE)

#NNObsClim$climResids<-climModel$residuals

# slope time series

temp<-subset(NNObs2, latitude>=35 & latitude<=45)

pred<-temp %>%
  dplyr::group_by(mean_first_yes_year,common_name) %>%
  dplyr::summarize(latSlope = summary(lm(mean_first_yes_doy ~ latitude+elevation_in_meters))$coef['latitude','Estimate'],
                   confLow  = confint(lm(mean_first_yes_doy ~ latitude+elevation_in_meters), 'latitude', level=0.95)[1,1],
                   confHigh = confint(lm(mean_first_yes_doy ~ latitude+elevation_in_meters), 'latitude', level=0.95)[1,2])


ggplot(pred, aes(mean_first_yes_year,latSlope,color=common_name))+
  geom_line()+
  ggtitle("Latitude Slope")+
  geom_hline(yintercept = summary(lm(mean_first_yes_doy ~ latitude+elevation_in_meters,data=NNObs2))$coef['latitude','Estimate'])+
  geom_ribbon(data=pred, aes(ymin = confLow, ymax = confHigh, x = mean_first_yes_year, group=common_name, fill=common_name), alpha = 0.4) 

ggplot(temp, aes(latitude,mean_first_yes_doy,color=common_name))+
  geom_point(size=1, alpha=0.4)+
  facet_wrap(.~mean_first_yes_year)+
  stat_smooth(method = lm)+
  theme_bw()

hist(NNObs2$latitude)

