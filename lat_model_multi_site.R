# Multi-station lat model
# adapted from lat_model.R
# MAC 01/10/22

#
library(rnpn)
library(dplyr)
library(plyr)
library(ggplot2)

NNObs <- npn_download_site_phenometrics(request_source = "TCrimmins", c(2009:2021),
                                        species_id = c(12,79), phenophase_id = c(483,501),
                                        additional_fields = c("species_functional_type", "species_category", "Observed_Status_Conflict_Flag",
                                                              "site_name", "Partner_Group", "Dataset_ID"))
# species_id = c(12,79), phenophase_id = c(483,501),

# remove >1st instances of "yes" on the same individual in the same year
leaf_data_1yes <- NNObs %>%
  group_by(species, phenophase_id, site_id, mean_first_yes_year) %>%
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

# NNObs
NNObs2$name_pheno<-paste0(NNObs2$common_name,"-",NNObs2$phenophase_description)
temp<-subset(NNObs2, name_pheno %in% c("American beech-Leaves","flowering dogwood-Open flowers"))

  ggplot(temp, aes(latitude,mean_first_yes_doy,color=name_pheno))+
    geom_point(size=1, alpha=0.4)+
    facet_wrap(.~mean_first_yes_year)+
    stat_smooth(method = lm)+
    theme_bw()
  
  # doy~lat by phenophase
  ggplot(NNObs2, aes(latitude,mean_first_yes_doy, color=phenophase_description))+
    geom_point()+
    stat_smooth(method = lm)

# subset to specific spp/pheno 
NNObs2<-subset(NNObs2, name_pheno %in% c("American beech-Leaves","flowering dogwood-Open flowers"))  
  

# slope coefficient 
test<-lm(mean_first_yes_doy ~ latitude+elevation_in_meters, data=NNObs2)
summary(test)
summary(lm(mean_first_yes_doy ~ latitude+elevation_in_meters,data=NNObs2))$coef['latitude','Estimate']

pred<-NNObs2 %>%
  dplyr::group_by(name_pheno) %>%
  dplyr::summarize(latSlope = summary(lm(mean_first_yes_doy ~ latitude+elevation_in_meters))$coef['latitude','Estimate'])

pred<-NNObs2 %>%
  dplyr::group_by(mean_first_yes_year,name_pheno) %>%
  dplyr::summarize(latSlope = summary(lm(mean_first_yes_doy ~ latitude+elevation_in_meters))$coef['latitude','Estimate'],
                   confLow  = confint(lm(mean_first_yes_doy ~ latitude+elevation_in_meters), 'latitude', level=0.95)[1,1],
                   confHigh = confint(lm(mean_first_yes_doy ~ latitude+elevation_in_meters), 'latitude', level=0.95)[1,2])

pred<-subset(pred, latSlope<10 & latSlope>=-5)                   

ggplot(pred, aes(mean_first_yes_year,latSlope,color=name_pheno))+
  geom_line()+
  ggtitle("Latitude Slope")+
  geom_hline(yintercept = summary(lm(mean_first_yes_doy ~ latitude+elevation_in_meters,data=NNObs2))$coef['latitude','Estimate'])+
  geom_ribbon(data=pred, aes(ymin = confLow, ymax = confHigh, x = mean_first_yes_year, group=name_pheno, fill=name_pheno), alpha = 0.4) 

ggplot(NNObs2, aes(latitude,mean_first_yes_doy,color=name_pheno))+
  geom_point(size=1, alpha=0.4)+
  facet_wrap(.~mean_first_yes_year)+
  stat_smooth(method = lm)+
  theme_bw()+
  ggtitle("spp-phenophase")
  

predLow = predict(lm(mean_first_yes_doy ~ latitude+elevation_in_meters, data = NNObs2),
                  newdata=data.frame(latitude = c(35), elevation_in_meters = c(350)), interval="prediction")[1,2]


# multiple models
pred1<-NNObs2 %>%
      dplyr::group_by(mean_first_yes_year, name_pheno) %>%
      dplyr::summarize(pred = predict(lm(mean_first_yes_doy ~ latitude+elevation_in_meters),
                               newdata=data.frame(latitude = c(35), elevation_in_meters = c(350))),
                       predLow = predict(lm(mean_first_yes_doy ~ latitude+elevation_in_meters),
                                      newdata=data.frame(latitude = c(35), elevation_in_meters = c(350)),interval="prediction")[2],
                       predHigh = predict(lm(mean_first_yes_doy ~ latitude+elevation_in_meters),
                                         newdata=data.frame(latitude = c(35), elevation_in_meters = c(350)),interval="prediction")[3]) 

pred1$lat<-"lat35"
temp<-tidyr::spread(pred1[,1:3], name_pheno, pred)
cor(temp[,c(2:3)], method = "spearman")

pred2<-NNObs2 %>%
  dplyr::group_by(mean_first_yes_year, name_pheno) %>%
  dplyr::summarize(pred = predict(lm(mean_first_yes_doy ~ latitude+elevation_in_meters),
                                  newdata=data.frame(latitude = c(45), elevation_in_meters = c(350))),
                   predLow = predict(lm(mean_first_yes_doy ~ latitude+elevation_in_meters),
                                     newdata=data.frame(latitude = c(45), elevation_in_meters = c(350)),interval="prediction")[2],
                   predHigh = predict(lm(mean_first_yes_doy ~ latitude+elevation_in_meters),
                                      newdata=data.frame(latitude = c(45), elevation_in_meters = c(350)),interval="prediction")[3]) 
pred2$lat<-"lat45"
temp<-tidyr::spread(pred2[,1:3], name_pheno, pred)
cor(temp[,c(2:3)], method = "spearman")


# pred3<-NNObs2 %>%
#   dplyr::group_by(mean_first_yes_year, name_pheno) %>%
#   dplyr::summarize(pred = predict(lm(mean_first_yes_doy ~ latitude+elevation_in_meters),
#                                   newdata=data.frame(latitude = c(45), elevation_in_meters = c(350))),
#                    predLow = predict(lm(mean_first_yes_doy ~ latitude+elevation_in_meters),
#                                      newdata=data.frame(latitude = c(45), elevation_in_meters = c(350)),interval="prediction")[2],
#                    predHigh = predict(lm(mean_first_yes_doy ~ latitude+elevation_in_meters),
#                                       newdata=data.frame(latitude = c(45), elevation_in_meters = c(350)),interval="prediction")[3]) 
# pred3$lat<-"lat50"
# temp<-tidyr::spread(pred3, name_pheno, pred)
# cor(temp[,c(3:4)], method = "spearman")


pred<-rbind.data.frame(pred1,pred2)
pred<-subset(pred, pred<180 & pred>0)

library(ggplot2)

ggplot(pred, aes(mean_first_yes_year,pred,color=lat))+
  geom_line()+
  facet_wrap(.~name_pheno)

ggplot(pred, aes(mean_first_yes_year,pred,color=name_pheno))+
  geom_line()+
  facet_wrap(.~lat)

ggplot(pred, aes(mean_first_yes_year,pred,color=lat))+
  geom_line()+
  #ggtitle(paste0(NNObs2$common_name[1],"-",NNObs2$phenophase_description[1], ":Start of Spring"))+
  #geom_hline(yintercept = meanDate[1])+
  #geom_hline(yintercept = meanDate[2])+
  geom_ribbon(data=pred, aes(ymin = predLow, ymax = predHigh, x = mean_first_yes_year, group=lat,fill=lat), alpha = 0.4)+
  facet_wrap(.~name_pheno)+
  ggtitle("DOY predictions")

ggplot(pred, aes(mean_first_yes_year,pred,color=name_pheno))+
  geom_line()+
  #ggtitle(paste0(NNObs2$common_name[1],"-",NNObs2$phenophase_description[1], ":Start of Spring"))+
  #geom_hline(yintercept = meanDate[1])+
  #geom_hline(yintercept = meanDate[2])+
  geom_ribbon(data=pred, aes(ymin = predLow, ymax = predHigh, x = mean_first_yes_year, group=name_pheno,fill=name_pheno), alpha = 0.4)+
  facet_wrap(.~lat)+
  ggtitle("DOY predictions")

# stem plot - spring length

temp<-tidyr::spread(pred[,c(1:3,6)], lat, pred)
temp$length<-temp$lat45-temp$lat35

ggplot(temp, aes(mean_first_yes_year,pred,color=name_pheno))+
  #geom_line()+
  #ggtitle(paste0(NNObs2$common_name[1],"-",NNObs2$phenophase_description[1], ":Start of Spring"))+
  #geom_hline(yintercept = meanDate[1])+
  #geom_hline(yintercept = meanDate[2])+
  #geom_ribbon(data=pred, aes(ymin = predLow, ymax = predHigh, x = mean_first_yes_year, group=lat,fill=lat), alpha = 0.4)+
  geom_segment(aes(x=mean_first_yes_year,y=lat35,xend=mean_first_yes_year,yend=length+lat35))+
  facet_wrap(.~name_pheno)+
  ggtitle("Season length (prediction at 35 and 45 deg lat)")





