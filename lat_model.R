# pulling NN obs to show how to use data in a particular year as an indicator
# T Crimmins, M Crimmins, theresa@usanpn.org
# 12-18-21
#
# Pull data from NN
#
library(rnpn)
library(dplyr)
library(plyr)
# species IDs, phenophaseIDs
# red maple = 3; flowers or flower buds = 500
# choke cherry = XX
# balsam fir
# use this to figure out what phenophase IDs to pull - put in speciesID and year
#test <- npn_phenophases_by_species(3, 2010)
# download "individual phenometrics" data; red maple spp=3, open flowers=501

# models to consider 61/371 (sugar maple), 

NNObs <- npn_download_individual_phenometrics(request_source = "TCrimmins", c(2009:2021),
                                              species_id = 28, phenophase_id = 371,
                                              additional_fields = c("species_functional_type", "species_category", "Observed_Status_Conflict_Flag",
                                                                    "site_name", "Partner_Group", "Dataset_ID"))
# remove >1st instances of "yes" on the same individual in the same year
leaf_data_1yes <- NNObs %>%
  group_by(species, individual_id, first_yes_year) %>%
  filter(first_yes_doy == min(first_yes_doy)) %>%
  slice(1) %>% # takes the first occurrence if there is a tie
  ungroup()
# remove any observations with no "no" observations within 30 days prior to a "yes" observation
leaf_data_1yes_prevno <- leaf_data_1yes[which(leaf_data_1yes$numdays_since_prior_no >=0),]
leaf_data_1yes_prevno <- leaf_data_1yes_prevno[which(leaf_data_1yes_prevno$numdays_since_prior_no <=15),]

# ##### NEEDED?
# # remove individual observations that have been flagged with the Observer_Status_Conflict_Flag
# leaf_data_1yes_flag <- leaf_data_1yes[which(leaf_data_1yes$observed_status_conflict_flag==-9999),]
# # remove sites for which greater than 5 percent of observations have a status conflict record
# # to do this - calculate the percentage of conflicts per site
# # get total records per site
# sitetotals <- leaf_data_1yes_flag %>%
#   group_by(site_id) %>%
#   dplyr::summarise(n=n())
# sitetotals <- as.data.frame(sitetotals)
# # filter out observations with conflict flags
# conflict <- leaf_data_1yes_flag[which(leaf_data_1yes_flag$observed_status_conflict_flag== "MultiObserver-StatusConflict"),]
# # use that to get total conflicts per site
# conflict_totals <- conflict %>%
#   group_by(site_id) %>%
#   dplyr::summarise(nflag=n())
# conflict_totals <- as.data.frame(conflict_totals)
# # merge the total records and total conflicts (replace NAs with zeroes for sites
# # with no conflicts)
# hi_conflict <- merge(conflict_totals, sitetotals, by="site_id", all=TRUE)
# hi_conflict[is.na(hi_conflict)] <- 0
# # calculate percentage of conflict records
# hi_conflict$percentage <- hi_conflict$nflag/hi_conflict$n*100
# # identify sites for which conflict records make up more than 5 percent of total records
# lowconflict_sites <- hi_conflict[which(hi_conflict$percentage<=5),]
# lowconflict_sites <- as.data.frame(lowconflict_sites[,1])
# names(lowconflict_sites) <- "site_id"
# # remove the high-conflict sites from the previously filtered data (with single conflict flags removed) to see the difference
# leaf_data_1yes_flagsite <- merge(leaf_data_1yes_flag, lowconflict_sites, by="site_id")
# # Remove spring observations after DOY 180
# NNObs2 <- leaf_data_1yes_flagsite[which(leaf_data_1yes_flagsite$first_yes_doy <=180),]
# # count up #years by individuals*phenophases
# sitesperYr <- NNObs %>%
#   group_by(first_yes_year) %>%
#   dplyr::summarize(count = n())
# ##################################################################################3

NNObs2<-leaf_data_1yes_prevno[which(leaf_data_1yes_prevno$first_yes_doy <=180),]

NNObs2<-subset(NNObs2,longitude<0)

#### look for lat effect

plot(NNObs2$latitude,NNObs2$first_yes_doy)

library(ggplot2)
ggplot(NNObs2, aes(latitude, first_yes_doy, color=(first_yes_year)))+
  geom_point()+
  #facet_wrap(.~first_yes_year)+
  stat_smooth(method = lm)
  #stat_smooth()

latModel <- lm(first_yes_doy ~ latitude+elevation_in_meters, data = NNObs2)
#latModel <- lm(first_yes_doy ~ latitude, data = NNObs2)

summary(latModel)

NNObs2$latResids<-latModel$residuals
NNObs2$year<-NNObs2$first_yes_year

ggplot(NNObs2, aes(latitude, latResids, color=first_yes_year))+
  geom_point()+
  stat_smooth(method = lm)+
  stat_smooth()

states <- map_data("state")


ggplot() +
  geom_polygon(data = states, aes(x = long, y = lat, group = group), fill=NA, color="black", size=0.1)  +
  #coord_fixed(xlim=c(out$meta$ll[1]-zoomLev, out$meta$ll[1]+zoomLev), ylim=c(out$meta$ll[2]-zoomLev, out$meta$ll[2]+zoomLev), ratio = 1) +
  #coord_fixed(xlim=c(-118, -102), ylim=c(28.75, 37.5), ratio = 1) +
  geom_point(data = NNObs2, aes(x = longitude, y = latitude, fill=latResids),colour="grey", shape=21, size = 2, position="jitter")+
  scale_fill_gradient2(low = "red",mid = "white", high = "blue", limits=c(-20, 20), oob=scales::squish)+
  lims(x = c(min(NNObs2$longitude), max(NNObs2$longitude)), y = c(min(NNObs2$latitude), max(NNObs2$latitude)))+
  facet_wrap(.~first_yes_year)+
  theme_bw()+
  ggtitle("Black Cherry - 371")

ggplot(NNObs2, aes(latResids)) +
  geom_density()

ggplot() +
  geom_polygon(data = states, aes(x = long, y = lat, group = group), fill=NA, color="black", size=0.1)  +
  #coord_fixed(xlim=c(out$meta$ll[1]-zoomLev, out$meta$ll[1]+zoomLev), ylim=c(out$meta$ll[2]-zoomLev, out$meta$ll[2]+zoomLev), ratio = 1) +
  #coord_fixed(xlim=c(-118, -102), ylim=c(28.75, 37.5), ratio = 1) +
  geom_point(data = NNObs2, aes(x = longitude, y = latitude, color=first_yes_doy))+
  #scale_colour_gradient2(low = "red",mid = "white", high = "blue", limits=c(-20, 20), oob=scales::squish)+
  scale_colour_gradient2(low = "red",mid = "green", high = "blue")+
  facet_wrap(.~first_yes_year)+
  ggtitle("Red Maple - 371")

ggplot() +
  geom_polygon(data = states, aes(x = long, y = lat, group = group), fill=NA, color="black", size=0.1)  +
  #coord_fixed(xlim=c(out$meta$ll[1]-zoomLev, out$meta$ll[1]+zoomLev), ylim=c(out$meta$ll[2]-zoomLev, out$meta$ll[2]+zoomLev), ratio = 1) +
  #coord_fixed(xlim=c(-118, -102), ylim=c(28.75, 37.5), ratio = 1) +
  geom_point(data = NNObs2, aes(x = longitude, y = latitude, color=first_yes_month))+
  #scale_colour_gradient2(low = "red",mid = "white", high = "blue", limits=c(-20, 20), oob=scales::squish)+
  scale_colour_gradient2(low = "red",mid = "green", high = "blue")+
  facet_wrap(.~first_yes_year)

ggplot()+
  geom_histogram(data=NNObs2, aes(first_yes_month))+
  facet_wrap(.~year)

# add in rasters, trying https://gis.stackexchange.com/questions/377444/plotting-a-raster-stack-with-ggplot2
library(raster)
#library(rasterVis)
library(ggplot2)

#tmean_anom<-stack("/scratch/crimmins/climgrid/processed/conus/CONUSmonthly_nClimGrid_tmean_anom_011895_092021.grd")
# make 2009-2021 anom grid
tmean_anom<-stack("/scratch/crimmins/climgrid/processed/conus/CONUSmonthly_nClimGrid_tmean_anom_012009_092021.grd")



# find last complete year
dates<-as.data.frame(names(tmean_anom))
dates<-as.data.frame(seq(as.Date(strsplit(names(tmean_anom[[1]]), "X")[[1]][2],"%Y.%m.%d"),
           as.Date(strsplit(names(tmean_anom[[nlayers(tmean_anom)]]), "X")[[1]][2],"%Y.%m.%d"), by="month"))
colnames(dates)<-"date"
dates$month<-as.numeric(format(dates$date, "%m"))
dates$year<-as.numeric(format(dates$date, "%Y"))

# subset anom stack 
subAnom<-tmean_anom[[which(dates$month==4 & dates$year>=min(NNObs2$first_yes_year))]]
#names(subAnom)<-seq(2009,2021,1)

coords <- xyFromCell(subAnom, seq_len(ncell(subAnom)))
subAnom <- stack(as.data.frame(getValues(subAnom)))
names(subAnom) <- c('value', 'variable')
subAnom$year<-as.numeric(format(as.Date(subAnom$variable,format="X%Y.%m.%d"),"%Y"))

subAnom <- cbind(coords, subAnom)
ggplot(subAnom) + 
  geom_tile(aes(x, y, fill = value)) +
  geom_polygon(data = states, aes(x = long, y = lat, group = group), fill=NA, color="black", size=0.1)  +
  geom_jitter(data = NNObs2, aes(x = longitude, y = latitude, color=latResids),size = 2, width=1, height=1)+
  scale_color_gradient2(low = "red",mid = "white", high = "blue", limits=c(-20, 20), oob=scales::squish)+
  scale_fill_gradient2(low="blue",mid = "white",high = "red", midpoint = 0) +
  facet_wrap(~ year) +
  coord_equal()

ggplot(NNObs2, aes(latResids))+
  geom_histogram()+
  facet_wrap(~ year)
  

yearStats<-NNObs2 %>% 
      dplyr::group_by(first_yes_year) %>%
        dplyr::summarise(meanResid=mean(latResids, na.rm=TRUE),
                         medResid=median(latResids, na.rm=TRUE),
                         sdResid=sd(latResids, na.rm=TRUE),
                         n = n())

ggplot(yearStats, aes(first_yes_year, meanResid, size=sdResid))+
  geom_point()+
  geom_hline(yintercept = 0)

# plot resid against climate

# subset anom stack 
subAnom<-tmean_anom[[which(dates$month==5 & dates$year>=min(NNObs2$first_yes_year))]]

# loop through years
tempObs<-subset(NNObs2, year==2021)
temp<-extract(subAnom[[1]], SpatialPoints(tempObs[,c("longitude","latitude")]), sp = T)
colnames(temp@data)<-"climAnom"
tempObs<-cbind.data.frame(tempObs, temp@data)

cor(tempObs$latResids, tempObs$climAnom, "na.or.complete", method="spearman")
plot(tempObs$latResids, tempObs$climAnom)

ggplot(tempObs, aes(latResids, climAnom))+
  geom_point()+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept=0)

ggplot(NNObs2, aes(first_yes_year))+
  geom_bar()

ggplot(NNObs2, aes(y=latResids))+
  geom_boxplot(varwidth = TRUE)+
  facet_wrap(.~year)

# bootstrap the regression coefficients

sample_coef_intercept <- NULL
sample_coef_lat <- NULL
sample_coef_elev <- NULL

for (i in 1:1000) {
  sample_d = NNObs2[sample(1:nrow(NNObs2), nrow(NNObs2), replace = TRUE), ]
  
  model_bootstrap <-lm(first_yes_doy ~ latitude+elevation_in_meters, data = sample_d)
  
  sample_coef_intercept <-
    c(sample_coef_intercept, model_bootstrap$coefficients[1])
  
  sample_coef_lat <-
    c(sample_coef_lat, model_bootstrap$coefficients[2])
  
  sample_coef_elev <-
    c(sample_coef_elev, model_bootstrap$coefficients[3])
}

coefs <- t(rbind(sample_coef_intercept, sample_coef_lat,sample_coef_elev))

# standard model
lm(first_yes_doy ~ latitude+elevation_in_meters, data = NNObs2)

summary(coefs)

# visualizing regression results
fit1<-lm(first_yes_doy ~ latitude+elevation_in_meters, data = NNObs2)

equation1=function(x){coef(fit1)[2]*x+coef(fit1)[1]}
equation2=function(x){coef(fit1)[2]*x+coef(fit1)[1]+coef(fit1)[3]}

ggplot(NNObs2,aes(y=first_yes_doy,x=latitude,color=elevation_in_meters))+geom_point()+
  stat_function(fun=equation1,geom="line",color=scales::hue_pal()(2)[1])+
  stat_function(fun=equation2,geom="line",color=scales::hue_pal()(2)[2])

#https://cran.r-project.org/web/packages/jtools/vignettes/summ.html
jtools::summ(fit1, scale=TRUE)
jtools::summ(fit1, confint = TRUE, digits = 3)

jtools::effect_plot(fit1, pred = elevation_in_meters, interval = TRUE, plot.points = TRUE)

jtools::plot_summs(fit1, scale = TRUE)
