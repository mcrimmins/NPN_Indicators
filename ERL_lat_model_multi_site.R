# ERL Figure 4
# Multi-station lat model
# adapted from lat_model.R, lat_model_multi_site.R
# MAC 01/31/22

library(rnpn)
library(dplyr)
library(plyr)
library(ggplot2)

# DOWNLOAD DATA 
#####
# set spp/phenophases
NNObs <- npn_download_site_phenometrics(request_source = "TCrimmins", c(2009:2021),
                                        species_id = c(12,79), phenophase_id = c(483,501),
                                        additional_fields = c("species_functional_type", "species_category", "Observed_Status_Conflict_Flag",
                                                              "site_name", "Partner_Group", "Dataset_ID"))

## CLEAN up data
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

# create spp~phenophase combo labels
NNObs2$name_pheno<-paste0(NNObs2$common_name,"-",NNObs2$phenophase_description)
# subset to specific spp/pheno 
NNObs2<-subset(NNObs2, name_pheno %in% c("American beech-Leaves","flowering dogwood-Open flowers"))  
  
# spp~phenophase models
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

# combine 35/45 results into common df
pred<-rbind.data.frame(pred1,pred2)
pred<-subset(pred, pred<180 & pred>0)

# FIGURE 4
# spring progression at 35 to 45 lat
library(ggplot2)
# stem plot - spring length

temp<-tidyr::spread(pred[,c(1:3,6)], lat, pred)
temp$length<-temp$lat45-temp$lat35
temp$name_pheno<-factor(temp$name_pheno, levels=c("flowering dogwood-Open flowers","American beech-Leaves"))
# remove 2009
temp<-subset(temp, mean_first_yes_year>2009)

p<-ggplot(temp, aes(ymin = lat35, ymax = lat45, x = mean_first_yes_year)) + 
  geom_linerange(aes(color = name_pheno),position = position_dodge(width = 0.4), size = 2)+
  scale_x_continuous(breaks=seq(2009,2021,1))+
  scale_colour_manual(values = c( "#8da0cb","#66c2a5"), name="")+
  theme_bw()+
  ylab("Phenophase estimate (Day of Year)")+
  xlab("Year")+
  theme(legend.position="bottom")+
  theme(legend.text=element_text(size=11))

#ggsave(p, file="./figs/fig4.eps", device="eps")

tiff("/home/crimmins/RProjects/NPNClimate/figs/fig4.tif",
     width = 7.5, height = 3.5, units = "in", res = 300L)
print(p)
dev.off()
