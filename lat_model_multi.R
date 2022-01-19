# Multi-station lat model
# adapted from lat_model.R
# MAC 01/10/22

#
library(rnpn)
library(dplyr)
library(plyr)
library(ggplot2)

# models to consider 61/371 (sugar maple), 
NNObs <- npn_download_individual_phenometrics(request_source = "TCrimmins", c(2009:2021),
                                              species_id = c(12,79), phenophase_id = c(483,501),
                                              additional_fields = c("species_functional_type", "species_category", "Observed_Status_Conflict_Flag",
                                                                    "site_name", "Partner_Group", "Dataset_ID"))
# SCREEN DATA
# remove >1st instances of "yes" on the same individual in the same year
leaf_data_1yes <- NNObs %>%
  group_by(species, individual_id, first_yes_year) %>%
  filter(first_yes_doy == min(first_yes_doy)) %>%
  slice(1) %>% # takes the first occurrence if there is a tie
  ungroup()
# remove any observations with no "no" observations within 30 days prior to a "yes" observation
leaf_data_1yes_prevno <- leaf_data_1yes[which(leaf_data_1yes$numdays_since_prior_no >=0),]
leaf_data_1yes_prevno <- leaf_data_1yes_prevno[which(leaf_data_1yes_prevno$numdays_since_prior_no <=15),]
# more screening based on season length
NNObs2<-leaf_data_1yes_prevno[which(leaf_data_1yes_prevno$first_yes_doy <=180),]
# geography
NNObs2<-subset(NNObs2,longitude<0)
NNObs2<-subset(NNObs2,latitude<=50)

# NNObs
NNObs2$name_pheno<-paste0(NNObs2$common_name,"-",NNObs2$phenophase_description)
temp<-subset(NNObs2, name_pheno %in% c("American beech-Leaves","flowering dogwood-Open flowers"))

ggplot(temp, aes(latitude,first_yes_doy,color=name_pheno))+
  geom_point(size=1, alpha=0.4)+
  facet_wrap(.~first_yes_year)+
  stat_smooth(method = lm)+
  theme_bw()

# doy~lat by phenophase
ggplot(NNObs2, aes(latitude,first_yes_doy, color=phenophase_description))+
  geom_point()+
  stat_smooth(method = lm)


# slope coefficient 
test<-lm(first_yes_doy ~ latitude+elevation_in_meters, data=NNObs2)
summary(test)
summary(lm(first_yes_doy ~ latitude+elevation_in_meters,data=NNObs2))$coef['latitude','Estimate']

pred<-NNObs2 %>%
  dplyr::group_by(common_name) %>%
  dplyr::summarize(latSlope = summary(lm(first_yes_doy ~ latitude+elevation_in_meters))$coef['latitude','Estimate'])

pred<-NNObs2 %>%
  dplyr::group_by(first_yes_year,common_name) %>%
  dplyr::summarize(latSlope = summary(lm(first_yes_doy ~ latitude+elevation_in_meters))$coef['latitude','Estimate'],
                   confLow  = confint(lm(first_yes_doy ~ latitude+elevation_in_meters), 'latitude', level=0.95)[1,1],
                   confHigh = confint(lm(first_yes_doy ~ latitude+elevation_in_meters), 'latitude', level=0.95)[1,2])
                   

ggplot(pred, aes(first_yes_year,latSlope,color=common_name))+
  geom_line()+
  ggtitle("Latitude Slope - leaves")+
  geom_hline(yintercept = summary(lm(first_yes_doy ~ latitude+elevation_in_meters,data=NNObs2))$coef['latitude','Estimate'])+
  geom_ribbon(data=pred, aes(ymin = confLow, ymax = confHigh, x = first_yes_year, group=common_name, fill=common_name), alpha = 0.4) 

ggplot(NNObs2, aes(latitude,first_yes_doy,color=common_name))+
  geom_point(size=1, alpha=0.4)+
  facet_wrap(.~first_yes_year)+
  stat_smooth(method = lm)+
  theme_bw()+
  ggtitle("leaves phenophase")
  


# multiple models
pred1<-NNObs2 %>%
      dplyr::group_by(first_yes_year, phenophase_description) %>%
      dplyr::summarize(pred = predict(lm(first_yes_doy ~ latitude+elevation_in_meters),
                               newdata=data.frame(latitude = c(35), elevation_in_meters = c(350)))) 
pred1$lat<-"lat35"
temp<-tidyr::spread(pred1, phenophase_description, pred)
cor(temp[,c(3:4)], method = "spearman")

pred2<-NNObs2 %>%
  dplyr::group_by(first_yes_year, phenophase_description) %>%
  dplyr::summarize(pred = predict(lm(first_yes_doy ~ latitude+elevation_in_meters),
                                  newdata=data.frame(latitude = c(45), elevation_in_meters = c(350)))) 
pred2$lat<-"lat45"
temp<-tidyr::spread(pred2, phenophase_description, pred)
cor(temp[,c(3:4)], method = "spearman")



pred3<-NNObs2 %>%
  dplyr::group_by(first_yes_year, common_name) %>%
  dplyr::summarize(pred = predict(lm(first_yes_doy ~ latitude+elevation_in_meters),
                                  newdata=data.frame(latitude = c(50), elevation_in_meters = c(350)))) 
pred3$lat<-"lat50"
temp<-tidyr::spread(pred3, common_name, pred)
cor(temp[,c(3:5)], method = "spearman")


pred<-rbind.data.frame(pred1,pred2)
pred<-subset(pred, pred<180 & pred>0)

library(ggplot2)

ggplot(pred, aes(first_yes_year,pred,color=lat))+
  geom_line()+
  facet_wrap(.~phenophase_description)

ggplot(pred, aes(first_yes_year,pred,color=phenophase_description))+
  geom_line()+
  facet_wrap(.~lat)


