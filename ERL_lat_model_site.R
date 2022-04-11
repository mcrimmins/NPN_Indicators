# GENERATE FIGURES 1-3 for ERL paper
# Site-level lat model
# adapted from lat_model.R and lat_model_site.R
# MAC 01/31/22

library(rnpn)
library(dplyr)
library(plyr)
library(raster)
#library(rasterVis)
library(ggplot2)

# load terrain data
load("terrainData.RData")

# DOWNLOAD DATA 
##### change to spp id 79, pheno id 483 for supp fig
# set spp/phenophase
NNObs <- npn_download_site_phenometrics(request_source = "TCrimmins", c(2009:2021),
                                        species_id = c(12), phenophase_id = c(501),
                                        additional_fields = c("species_functional_type", "species_category", "Observed_Status_Conflict_Flag",
                                                              "site_name", "Partner_Group", "Dataset_ID"))
  ## CLEAN up data
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
  ##

# DEVELOP LAT/ELEV MODEL
latModel <- lm(mean_first_yes_doy ~ latitude+elevation_in_meters, data = NNObs2)
summary(latModel)
  # add residuals to dataframe
  NNObs2$latResids<-latModel$residuals
  NNObs2$year<-NNObs2$mean_first_yes_year
  # add in month
  NNObs2$month<-as.numeric(format(as.Date(paste0(NNObs2$mean_first_yes_doy,"-",NNObs2$mean_first_yes_year),format="%j-%Y"),"%m"))

# EXTRACT CLIMATE ANOMS (local 2-month) add to new dataframe
#####
  # CLIMATE GRIDS
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

# local clim window
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
  
  # terrain data
  temp_elev<-extract(elevation, SpatialPoints(tempObs[,c("longitude","latitude")]), sp = T)@data
    colnames(temp_elev)<-"elev"
  temp_slope<-extract(slope, SpatialPoints(tempObs[,c("longitude","latitude")]), sp = T)@data
    colnames(temp_slope)<-"slope"  
  temp_asp<-extract(aspect, SpatialPoints(tempObs[,c("longitude","latitude")]), sp = T)@data
    colnames(temp_asp)<-"aspect"
  temp_tpi<-extract(TPI, SpatialPoints(tempObs[,c("longitude","latitude")]), sp = T)@data
    colnames(temp_tpi)<-"TPI"
  temp_tri<-extract(TRI, SpatialPoints(tempObs[,c("longitude","latitude")]), sp = T)@data
    colnames(temp_tri)<-"TRI"
  temp_anth<-extract(anthrome, SpatialPoints(tempObs[,c("longitude","latitude")]), sp = T)@data
    colnames(temp_anth)<-"anthrome"  
  temp_tpi5<-extract(TPI5, SpatialPoints(tempObs[,c("longitude","latitude")]), sp = T)@data
    colnames(temp_tpi5)<-"TPI5"   
  
  tempObs<-cbind.data.frame(tempObs, temp_npn, temp_full, temp_elev, temp_slope,temp_asp, temp_tpi, temp_tpi5, temp_tri, temp_anth)
  
  NNObsClim[[i]] <- tempObs # add it to your list
  print(i)
}
#####

# add climate anoms back in to dataframe
NNObsClim <- do.call(rbind, NNObsClim)
NNObsClim$anomType_full<-ifelse(NNObsClim$climAnom_full<0,"neg","pos")
NNObsClim$anomType_npn<-ifelse(NNObsClim$climAnom_npn<0,"neg","pos")

#cor(NNObsClim$latResids, NNObsClim$climAnom_npn, "na.or.complete", method="pearson")
#cor(NNObsClim$latResids, NNObsClim$climAnom_full, "na.or.complete", method="spearman")

# annual stats on residuals

annClim<-NNObsClim %>% group_by(mean_first_yes_year) %>% 
                        dplyr::summarize(count=n(),
                                  meanResid=mean(latResids),
                                  sdResid=sd(latResids),
                                  meanclimAnom_npn=mean(climAnom_npn, na.rm=TRUE),
                                  meanclimAnom_full=mean(climAnom_full, na.rm=TRUE))
#annClim$label<-paste0(round(annClim$meanResid,1)," (",round(annClim$sdResid,1),")")
# annClim$label<-paste0(round(annClim$meanResid,1)," (",round(annClim$sdResid,1),")\n",
#                       round(annClim$meanclimAnom_npn,1)," (",round(annClim$meanclimAnom_full,1),")")
annClim$label<-paste0("n=",annClim$count)

#####
# lat-elev vs doy scatterplot - FIGURE 1

p<-ggplot(NNObsClim)+
  geom_point(aes(latitude, mean_first_yes_doy, color=elevation_in_meters))+
  geom_smooth(aes(latitude, mean_first_yes_doy), method=lm, se=FALSE)+
  scale_colour_gradientn(colours = terrain.colors(10), name="Elevation (m)")+
  #ggtitle(paste0(NNObs2$common_name[1],"-",NNObs2$phenophase_description[1]))+
  #ggpubr::stat_cor(aes(latitude, mean_first_yes_doy),method = "pearson")+
  theme_bw()+
  xlab("Latitude (deg N)")+
  ylab("Day of Year")

tiff("/home/crimmins/RProjects/NPNClimate/figs/fig1.tif",
     width = 7.5, height = 3.5, units = "in", res = 300L)
print(p)
dev.off()

# flipped
p<-ggplot(NNObsClim)+
  geom_point(aes(mean_first_yes_doy,latitude, color=elevation_in_meters))+
  geom_smooth(aes(mean_first_yes_doy,latitude), method=lm, se=FALSE)+
  scale_colour_gradientn(colours = terrain.colors(10), name="Elevation (m)")+
  #ggtitle(paste0(NNObs2$common_name[1],"-",NNObs2$phenophase_description[1]))+
  #ggpubr::stat_cor(aes(latitude, mean_first_yes_doy),method = "pearson")+
  theme_bw()+
  ylab("Latitude (deg N)")+
  xlab("Day of Year")

tiff("/home/crimmins/RProjects/NPNClimate/figs/fig1v2.tif",
     width = 7.5, height = 3.5, units = "in", res = 300L)
print(p)
dev.off()

plot(NNObsClim$latitude,NNObsClim$mean_first_yes_doy)
abline(lm(NNObsClim$mean_first_yes_doy ~ NNObsClim$latitude))

#####


#####
# residual map - FIGURE 2
states <- map_data("state")
countries<-map_data("world")
countries<-subset(countries, region %in% c("Canada","Mexico"))

library(cowplot)
library(ggplotify)
library(grid)

pMap<-ggplot(NNObsClim) + 
          geom_polygon(data = countries, aes(x = long, y = lat, group = group), fill="grey", color="black", size=0.1)  +
          geom_polygon(data = states, aes(x = long, y = lat, group = group), fill="white", color="grey", size=0.1)  +
          #geom_jitter(data = NNObsClim, aes(x = longitude, y = latitude, color=latResids),size = 1, width=1, height=1, alpha=0.8)+ # shape=anomType_npn
          geom_point(data = NNObsClim, aes(x = longitude, y = latitude, fill=latResids),colour="grey", shape=21, size = 1.5, position="jitter", alpha=0.8)+
          #geom_point(data = NNObsClim, aes(x = longitude, y = latitude, color=latResids),size = 1, alpha=0.8)+ # shape=anomType_npn
          scale_fill_gradient2(low = "#b2182b",mid = "#f7f7f7", high = "#2166ac", limits=c(-10, 10), oob=scales::squish,
                               name="Phenophase\nanom (days)",
                               labels = c("<-10", "-5", "0", "5",">10"),
                               breaks = c(-10, -5, 0, 5, 10))+
          #scale_fill_brewer(type="div",palette="RdBu")+
          #scale_shape_manual(values=c(16, 3))+
          facet_wrap(~ mean_first_yes_year) +
          coord_equal( xlim = c(-98.5, -66.5), ylim = c(25, 50))+
          theme_bw()+
          theme(panel.background = element_rect(fill = "powderblue"),
                axis.text.x = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks = element_blank(),
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                line = element_blank())+
          #ggtitle(paste0(NNObs2$common_name[1],"-",NNObs2$phenophase_description[1],": model resids vs local 2-mo clim anom (2009-2021 base)"))
          #ggtitle(paste0(NNObs2$common_name[1],"-",NNObs2$phenophase_description[1],": model resids"))+
           geom_text(data = annClim, aes(label=label), size=3,
                        x = Inf, y = -Inf, hjust=1.1, vjust=-0.5,
                        inherit.aes = FALSE)

pScat<-ggplot(NNObsClim)+
          #geom_point(aes(latResids, climAnom_npn, color=as.factor(mean_first_yes_year)))+
          geom_point(aes(latResids, climAnom_npn), size=1)+
          geom_smooth(data=NNObsClim,aes(latResids,climAnom_npn),method=lm,se=FALSE)+
          geom_hline(yintercept = 0)+
          geom_vline(xintercept=0)+
          #ggtitle(paste0(NNObs2$common_name[1],"-",NNObs2$phenophase_description[1]))+
          #ggpubr::stat_cor(aes(latResids,climAnom_npn),method = "spearman")+
          scale_colour_discrete("Year")+
          theme_bw()+
          theme(legend.position="bottom")+
          ylab("Temp anom (deg C)")+
          xlab("Phenophase anom (days)")

cor.test(NNObsClim$latResids, NNObsClim$climAnom_npn, method = "pearson")
#cor(NNObsClim$latResids, NNObsClim$climAnom_npn, method = "pearson", use = "pairwise.complete.obs")

# as inset outside plot
pScat<-as.grob(pScat)
p<-ggdraw(pMap) +
  draw_plot(pScat, 0.05, -0.37, scale=0.25)+
  draw_plot_label(
    c("(a)", "(b)"),
    c(0.01, 0.37),
    c(0.955, 0.25),
    size = 14
  )
#draw_plot(pMap, 0.3, -0.41, scale=0.35)
#p<-p+annotation_custom(grobWest)+annotation_custom(grobCent)+annotation_custom(grobEast)
# plot to png
tiff("/home/crimmins/RProjects/NPNClimate/figs/fig2.tif", width = 7.5, height = 6.5, units = "in", res = 300L)
#png("/home/crimmins/RProjects/NPNClimate/figs/fig2.png", width = 7.5, height = 6.5, units = "in", res = 300L)
#grid.newpage()
print(p, newpage = FALSE)
dev.off()

#####

#####
# yearly indicators time series FIGURE 3

p<-ggplot(annClim, aes(x=mean_first_yes_year)) +
  geom_line( aes(y=meanResid), color="blue") + 
  geom_line( aes(y=meanclimAnom_npn*5), color="red") + # Divide by 10 to get the same range than the temperature
  geom_line( aes(y=meanclimAnom_full*5), color="orange", linetype = "dashed") +
  geom_ribbon( mapping=aes(ymin=meanResid-sdResid, ymax=meanResid+sdResid),alpha=0.2, fill="blue")+
  scale_y_continuous(
    # Features of the first axis
    name = "Phenophase Anom (days)",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*1/5, name="Temp Anom (deg C)"))+
  geom_hline(yintercept = 0)+
  xlab("Year")+
  theme_bw()

scale_color_manual(name="Temp and DOY anomalies",values=c("DOY"="blue","2009-2021"="red","1895-2021"))
ggplot()+
  geom_line(data=Summary,aes(y=Y1,x= X,colour="Y1"),size=1 )+
  geom_line(data=Summary,aes(y=Y2,x= X,colour="Y2"),size=1) +
  scale_color_manual(name = "Y series", values = c("Y1" = "darkblue", "Y2" = "red"))

cor.test(annClim$meanResid, annClim$meanclimAnom_npn, method = "pearson")

tiff("/home/crimmins/RProjects/NPNClimate/figs/fig3.tif", width = 7.5, height = 3.5, units = "in", res = 300L)
#png("/home/crimmins/RProjects/NPNClimate/figs/fig2.png", width = 7.5, height = 6.5, units = "in", res = 300L)
#grid.newpage()
print(p, newpage = FALSE)
dev.off()


#####



