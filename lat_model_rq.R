# NPN indicators - quantreg method
# adapted from lat_model.R
# MAC 01/10/22

#
library(rnpn)
library(dplyr)
library(plyr)
library(ggplot2)

# models to consider 61/371 (sugar maple), 
NNObs <- npn_download_individual_phenometrics(request_source = "TCrimmins", c(2009:2021),
                                              species_id = 3, phenophase_id = 371,
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


# lat+elevation model
latModel <- quantreg::rq(first_yes_doy ~ latitude+elevation_in_meters, data = NNObs2)
  summary(latModel)
  summary(latModel, se='boot')
predDate<-predict(latModel, newdata=data.frame(latitude = c(30,50), elevation_in_meters = c(50,50)), interval="confidence")
# compare to lm 
# latModel <- lm(first_yes_doy ~ latitude+elevation_in_meters, data = NNObs2)
#   summary(latModel)
#   predict(latModel, newdata=data.frame(latitude = c(30,50), elevation_in_meters = c(50,50)), interval="prediction")



# goodness of fit https://stats.stackexchange.com/questions/129200/r-squared-in-quantile-regression
latNull <- quantreg::rq(first_yes_doy ~ 1, data = NNObs2)
  rho <- function(u,tau=.5)u*(tau - (u < 0))
  V <- sum(rho(latModel$resid, latModel$tau))
  V0 <- sum(rho(latNull$resid, latNull$tau))
R1 <- 1-V/V0

# use model to find first day of spring/last between US lat range
predDate<-predict(latModel, newdata=data.frame(latitude = c(30,50), elevation_in_meters = c(50,50)))
meanDate<-as.Date(predDate, origin = "2016-01-01")

# range of tau
latModel <- quantreg::rq(first_yes_doy ~ latitude+elevation_in_meters, data = NNObs2, tau=c(0.5,0.025,0.95))
summary(latModel)
predDate<-predict(latModel, newdata=data.frame(latitude = c(30,50), elevation_in_meters = c(50,50)))

# quantreg models by year SoS
#qrFit <- NNObs2 %>% group_by(first_yes_year) %>%
#  do(fitLat = quantreg::rq(first_yes_doy ~ latitude+elevation_in_meters, data = .))

sos<-(lapply(split(NNObs2, NNObs2$first_yes_year), function(x) {
  m<-quantreg::rq(first_yes_doy ~ latitude+elevation_in_meters, data=x)
  predict(m, newdata=data.frame(latitude = c(30,50), elevation_in_meters = c(350,350)))
}))

# pull out confidence interval
sos <- data.frame(matrix(unlist(sos), nrow=length(sos), byrow=TRUE))
colnames(sos)<-c("begin","end")
sos$length<-sos$end-sos$begin
sos$begin<-as.Date(sos$begin, origin = "2016-01-01")
sos$end<-as.Date(sos$end, origin = "2016-01-01")
#sos$begin<-format(as.Date(sos$begin, origin = "2016-01-01"),"%b-%d")
#sos$end<-format(as.Date(sos$end, origin = "2016-01-01"),"%b-%d")
sos$year<-seq(2009,2021,1)

sosG<-tidyr::gather(sos[,c("begin","end","year")], "point", "date", -year)

ggplot(sosG, aes(year,date,color=point))+
  geom_line()+
  ggtitle("Sugar maple 371 - Start of Spring")+
  geom_hline(yintercept = meanDate[1])+
  geom_hline(yintercept = meanDate[2])

#######
# plot with conf intervals
sos<-(lapply(split(NNObs2, NNObs2$first_yes_year), function(x) {
  m<-quantreg::rq(first_yes_doy ~ latitude+elevation_in_meters, data=x)
  predict(m, newdata=data.frame(latitude = c(30,50), elevation_in_meters = c(50,50)), interval="confidence")
}))

# pull out confidence interval
sos <- data.frame(matrix(unlist(sos), nrow=length(sos), byrow=TRUE))
colnames(sos)<-c("fitBegin","fitEnd","lwrBegin","lwrEnd","hgrBegin","hgrEnd")
  cols<-1:6
  sos[,cols] = apply(sos[,cols], 2, function(x) as.Date(x, origin = "2016-01-01",format="%Y-%m-%d"))
  myfun <- function(x) as.Date(x, format="%Y-%m-%d", origin="1970-01-01")
  sos <- data.frame(lapply(sos, myfun))
  # seas length
  sos$length<-sos$fitEnd-sos$fitBegin
  # add years
  sos$year<-seq(2009,2021,1)

sosG<-tidyr::gather(sos[,c("fitBegin","fitEnd","lwrBegin","lwrEnd","hgrBegin","hgrEnd","year")], "point", "date", -year)

ggplot(sosG, aes(year,date,color=point))+
  geom_line()+
  ggtitle("Sugar maple 371 - Start of Spring")+
  geom_hline(yintercept = meanDate[1])+
  geom_hline(yintercept = meanDate[2])

# grouped conf intervals
sosG1<-tidyr::gather(sos[,c("fitBegin","lwrBegin","hgrBegin","year")], "point", "date", -year,-lwrBegin,-hgrBegin)
  colnames(sosG1)[1:2]<-c("lwr","hgr")
sosG2<-tidyr::gather(sos[,c("fitEnd","lwrEnd","hgrEnd","year")], "point", "date", -year,-lwrEnd,-hgrEnd)
  colnames(sosG2)[1:2]<-c("lwr","hgr")
sosG1<-rbind.data.frame(sosG1,sosG2)  
  
ggplot(sosG1, aes(year,date,color=point))+
  geom_line()+
  ggtitle("Sugar maple 371 - Start of Spring")+
  geom_hline(yintercept = meanDate[1])+
  geom_hline(yintercept = meanDate[2])+
  geom_ribbon(data=sosG1, aes(ymin = lwr, ymax = hgr, x = year, group=point,fill=point), alpha = 0.4) 

## pred intervals
#####
# plot with prediction intervals base on tau
# plot with conf intervals
sos<-(lapply(split(NNObs2, NNObs2$first_yes_year), function(x) {
  m<-quantreg::rq(first_yes_doy ~ latitude+elevation_in_meters, data=x, tau = c(0.5,0.025,0.95))
  predict(m, newdata=data.frame(latitude = c(30,50), elevation_in_meters = c(50,50)))
}))

# pull out pred interval
sos <- data.frame(matrix(unlist(sos), nrow=length(sos), byrow=TRUE))
  colnames(sos)<-c("fitBegin","fitEnd","lwrBegin","lwrEnd","hgrBegin","hgrEnd")
  cols<-1:6
  sos[,cols] = apply(sos[,cols], 2, function(x) as.Date(x, origin = "2016-01-01",format="%Y-%m-%d"))
  myfun <- function(x) as.Date(x, format="%Y-%m-%d", origin="1970-01-01")
  sos <- data.frame(lapply(sos, myfun))
  # seas length
  sos$length<-sos$fitEnd-sos$fitBegin
  # add years
  sos$year<-seq(2009,2021,1)


# grouped conf intervals
sosG1<-tidyr::gather(sos[,c("fitBegin","lwrBegin","hgrBegin","year")], "point", "date", -year,-lwrBegin,-hgrBegin)
colnames(sosG1)[1:2]<-c("lwr","hgr")
sosG2<-tidyr::gather(sos[,c("fitEnd","lwrEnd","hgrEnd","year")], "point", "date", -year,-lwrEnd,-hgrEnd)
colnames(sosG2)[1:2]<-c("lwr","hgr")
sosG1<-rbind.data.frame(sosG1,sosG2) 

ggplot(sosG1, aes(year,date,color=point))+
  geom_line()+
  ggtitle("Sugar maple 371 - Start of Spring")+
  geom_hline(yintercept = meanDate[1])+
  geom_hline(yintercept = meanDate[2])+
  geom_ribbon(data=sosG1, aes(ymin = lwr, ymax = hgr, x = year, group=point,fill=point), alpha = 0.4) 

ggplot(NNObs2, aes(latitude, first_yes_doy, color=(elevation_in_meters)))+
  geom_point()+
  facet_wrap(.~first_yes_year)+
  stat_smooth(method = lm)+
  geom_quantile(quantiles = c(0.5,0.025,0.975))
  #stat_smooth()+
  #ggtitle("N Red Oak-371")

ggplot(NNObs2, aes(latitude, first_yes_doy, color=(elevation_in_meters)))+
  geom_point()+
  facet_wrap(.~first_yes_year)+
  stat_smooth(method = lm)+
  geom_quantile(quantiles = c(0.5,0.025,0.975))



#####

## pred intervals with lm
#####
# http://www.sthda.com/english/articles/40-regression-analysis/166-predict-in-r-model-predictions-and-confidence-intervals/
# use model to find first day of spring/last between US lat range
latModel<-lm(first_yes_doy ~ latitude+elevation_in_meters, data=NNObs2)
predDate<-predict(latModel, newdata=data.frame(latitude = c(30,50), elevation_in_meters = c(350,350)))
meanDate<-as.Date(predDate, origin = "2016-01-01")

sos<-(lapply(split(NNObs2, NNObs2$first_yes_year), function(x) {
  m<-lm(first_yes_doy ~ latitude+elevation_in_meters, data=x)
  predict(m, newdata=data.frame(latitude = c(30,50), elevation_in_meters = c(350,350)),interval="prediction")
}))

# pull out pred interval
sos <- data.frame(matrix(unlist(sos), nrow=length(sos), byrow=TRUE))
colnames(sos)<-c("fitBegin","fitEnd","lwrBegin","lwrEnd","hgrBegin","hgrEnd")
cols<-1:6
sos[,cols] = apply(sos[,cols], 2, function(x) as.Date(x, origin = "2016-01-01",format="%Y-%m-%d"))
myfun <- function(x) as.Date(x, format="%Y-%m-%d", origin="1970-01-01")
sos <- data.frame(lapply(sos, myfun))
# seas length
sos$length<-sos$fitEnd-sos$fitBegin
# add years
sos$year<-seq(2009,2021,1)


# grouped conf intervals
sosG1<-tidyr::gather(sos[,c("fitBegin","lwrBegin","hgrBegin","year")], "point", "date", -year,-lwrBegin,-hgrBegin)
colnames(sosG1)[1:2]<-c("lwr","hgr")
sosG2<-tidyr::gather(sos[,c("fitEnd","lwrEnd","hgrEnd","year")], "point", "date", -year,-lwrEnd,-hgrEnd)
colnames(sosG2)[1:2]<-c("lwr","hgr")
sosG1<-rbind.data.frame(sosG1,sosG2) 

ggplot(sosG1, aes(year,date,color=point))+
  geom_line()+
  ggtitle(paste0(NNObs2$common_name[1],"-",NNObs2$phenophase_description[1], ":Start of Spring"))+
  geom_hline(yintercept = meanDate[1])+
  geom_hline(yintercept = meanDate[2])+
  geom_ribbon(data=sosG1, aes(ymin = lwr, ymax = hgr, x = year, group=point,fill=point), alpha = 0.4) 

# Display conf and pred intervals 
model <- lm(first_yes_doy ~ latitude+elevation_in_meters, data=NNObs2)
# 1. Add predictions 
pred.int <- predict(model, interval = "prediction")
mydata <- cbind(NNObs2, pred.int)
# 2. Regression line + confidence intervals
library("ggplot2")
p <- ggplot(mydata, aes(latitude, first_yes_doy)) +
  geom_point() +
  stat_smooth(method = lm)
# 3. Add prediction intervals
p + geom_line(aes(y = lwr), color = "red", linetype = "dashed")+
  geom_line(aes(y = upr), color = "red", linetype = "dashed")

#####




