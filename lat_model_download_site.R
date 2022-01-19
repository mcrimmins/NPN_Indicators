# download data for site level
# analyze NN data along lat gradient
# adapted from lat_model2.R
# MAC 01/13/22
#
# Pull data from NN
#
library(rnpn)
library(dplyr)
library(plyr)

# spp id list
sspID<-c(3,61,82,102,79,36,81,702,708,12,705,976,27,52,100,28,795,704,117,35,131,60,915,90)
phenoID<-c(371,483,501)
#j=1

tempList<-list()
NPNlist <- list()

# download "individual phenometrics" data; red maple spp=3, open flowers=501
for(j in 1:length(phenoID)){  # length(phenoID)
  for(i in 1:length(sspID)){

  tempList[[i]] <- npn_download_site_phenometrics(request_source = "TCrimmins", c(2009:2021),
                                          species_id = sspID[i], phenophase_id = phenoID[j],
                                          additional_fields = c("species_functional_type", "species_category", "Observed_Status_Conflict_Flag",
                                                                "site_name", "Partner_Group", "Dataset_ID"))
  }
  NPNlist[[j]] <- tempList
}

save(NPNlist, file="NPNData_site_for_lat_model.RData")

#######

load("NPNData_site_for_lat_model.RData")

# loop through list to get stats

j=3 # phenophase

NNstats<-data.frame(name=character(),
                    species_id=integer(),
                    phenophase_id=integer(),
                    n=integer(),
                    minLat=integer(),
                    maxLat=integer(),
                    medLat=integer(),
                    medDOY=integer(),
                    formula=character(),
                    R2=integer())
# NNstats$n<-NA
# NNstats$minLat<-NA
# NNstats$maxLat<-NA
# NNstats$medLat<-NA
# NNstats$formula<-NA
# NNstats$R2<-NA

for(i in 1:length(NPNlist[[j]])){
  
  NNObs<-NPNlist[[j]][[i]]
  
  if(length(NNObs)!=0){
  # remove >1st instances of "yes" on the same individual in the same year
  leaf_data_1yes <- NNObs %>%
    group_by(species, site_id, mean_first_yes_year) %>%
    filter(mean_first_yes_doy == min(mean_first_yes_doy)) %>%
    slice(1) %>% # takes the first occurrence if there is a tie
    ungroup()
  # remove any observations with no "no" observations within 30 days prior to a "yes" observation
  leaf_data_1yes_prevno <- leaf_data_1yes[which(leaf_data_1yes$mean_numdays_since_prior_no >=0),]
  leaf_data_1yes_prevno <- leaf_data_1yes_prevno[which(leaf_data_1yes_prevno$mean_numdays_since_prior_no <=15),]
  
  NNObs<-leaf_data_1yes_prevno[which(leaf_data_1yes_prevno$mean_first_yes_doy <=180),]
  
  NNObs<-subset(NNObs,longitude<0)
  NNObs<-subset(NNObs,latitude<=50)
  
  ##
  
  if(nrow(NNObs)!=0){
    # develop models
    intercept_only <- lm(mean_first_yes_doy ~ 1, data=NNObs)
    all <- lm(mean_first_yes_doy ~ latitude+elevation_in_meters, data=NNObs)
    backward <- step(all, direction='backward', scope=formula(all), trace=0)
    
    # summarize in table
    tempdf<-cbind.data.frame(NNObs$common_name[1],NNObs$species_id[1],NNObs$phenophase_id[1],nrow(NNObs),min(NNObs$latitude),max(NNObs$latitude),median(NNObs$latitude),
                             median(NNObs$mean_first_yes_doy),as.character(backward$call)[2],summary(backward)$r.squared)
    
    NNstats<-rbind.data.frame(NNstats,tempdf)
    }
      else{
    }
  #as.character(backward$call)[2], 
  
  }
  else{
    
  }
}

colnames(NNstats)<-c("name","spp_id","pheno_id","n","minLat","maxLat","medLat","medDOY","formula","R2")

NNstats$range<-NNstats$maxLat-NNstats$minLat







