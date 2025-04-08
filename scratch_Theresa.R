# download & analysis of Caryn's leaf-out and flowering sequence/constancy
# Updated July 2022
# T Crimmins, M Crimmins
library(rnpn)
library(readr)
library(dplyr)
library(ggplot2)
library(data.table)
# download individual phenometrics from NE US
# LEAF-OUT: download BLB data
npnData <- npn_download_individual_phenometrics(
  request_source="TCrimmins",
  years=c(2009:2022),
  states = c("ME", "NY", "NH", "VT", "MA", "CT", "RI", "MN", "OH", "MI", "WI", "PA", "NJ"),
  phenophase_ids = c(371, 373), #only breaking leaf buds
  additional_fields = c("Observed_Status_Conflict_Flag", "Partner_Group"))
####################################
# BLOOM: download open flower data
# npnData <- npn_download_individual_phenometrics(
#   request_source="TCrimmins",
#   years=c(2009:2022),
#   states = c("ME", "NY", "NH", "VT", "MA", "CT", "RI", "MN", "OH", "MI", "WI", "PA", "NJ"),
#   phenophase_ids = c(205, 210, 494, 501), #just open flowers
#   additional_fields = c("Observed_Status_Conflict_Flag", "Partner_Group"))
####################################
# remove >1st instances of "yes" on the same individual in the same year
leaf_data_1yes <- npnData %>% 
  group_by(site_id, common_name, phenophase_id, first_yes_year) %>% 
  filter(first_yes_doy == min(first_yes_doy)) %>% 
  slice(1) %>% # takes the first occurrence if there is a tie
  ungroup()
# remove any observations with no "no" observations within 7 days prior to a "yes" observation
leaf_data_1yes_prevno <- leaf_data_1yes[which(leaf_data_1yes$numdays_since_prior_no >=0),]
leaf_data_1yes_prevno <- leaf_data_1yes_prevno[which(leaf_data_1yes_prevno$numdays_since_prior_no <=7),]
# remove individual observations that have been flagged with the Observer_Status_Conflict_Flag
leaf_data_1yes_flag <- leaf_data_1yes_prevno[which(leaf_data_1yes_prevno$observed_status_conflict_flag==-9999),]
# remove sites for which greater than 5 percent of observations have a status conflict record
# to do this - calculate the percentage of conflicts per site
# get total records per site
sitetotals <- leaf_data_1yes_flag %>% 
  group_by(site_id) %>%
  dplyr::summarise(n=n())
sitetotals <- as.data.frame(sitetotals)
# filter out observations with conflict flags
conflict <- leaf_data_1yes_flag[which(leaf_data_1yes_flag$observed_status_conflict_flag== "MultiObserver-StatusConflict"),]
# use that to get total conflicts per site
conflict_totals <- conflict %>% 
  group_by(site_id) %>%
  dplyr::summarise(nflag=n())
conflict_totals <- as.data.frame(conflict_totals)
# merge the total records and total conflicts (replace NAs with zeroes for sites
# with no conflicts)
hi_conflict <- merge(conflict_totals, sitetotals, by="site_id", all=TRUE)
hi_conflict[is.na(hi_conflict)] <- 0
# calculate percentage of conflict records
hi_conflict$percentage <- hi_conflict$nflag/hi_conflict$n*100
# identify sites for which conflict records make up more than 5 percent of total records
lowconflict_sites <- hi_conflict[which(hi_conflict$percentage<=5),]
lowconflict_sites <- as.data.frame(lowconflict_sites[,1])
names(lowconflict_sites) <- "site_id"
# remove the high-conflict sites from the previously filtered data (with single conflict flags removed) to see the difference
NNObs <- merge(leaf_data_1yes_flag, lowconflict_sites, by="site_id")
## Check to see if any sites have >1 individual per species in a year
spp_by_site_yr <- NNObs %>% 
  group_by(site_id, first_yes_year, common_name) %>% 
  dplyr::summarize(ct_inviduals_by_site = n())
# remove site-species with < 4 years of observations
# calculate number of years of obs per site-species combo
yrs_by_site_species <- NNObs %>% 
  group_by(site_id, common_name) %>% 
  dplyr::summarize(ct_years_by_site_species = n())
# join this table to NNObssub
join3 <- merge(NNObs, yrs_by_site_species, by.x = c('site_id', 'common_name'), by.y = c('site_id', 'common_name'), all.x = TRUE) 
# drop species with fewer than 5 years of observations at a site
NNObssub <- subset(join3, ct_years_by_site_species > 4)
# remove sites with < 2 spp*pp
# calculate number of spp by site
spp_by_site <- NNObssub %>% 
  group_by(site_id, first_yes_year) %>% 
  dplyr::summarize(ct_spp_by_site = n())
# join this table to NNObs
join2 <- merge(NNObssub, spp_by_site, by.x = c('site_id', 'first_yes_year'), by.y = c('site_id', 'first_yes_year'), all.x = TRUE) 
# drop years from sites with fewer than 2 spp
NNObssub2 <- subset(join2, ct_spp_by_site > 1)
## BEFORE FINALIZING:
# do one more check to ensure no sites with <2 species
# write out csvS
#write.csv(NNObssub2, "~/RProjects/CarynB/NNobs_leaf_NY7d_7-28-22.csv")
#write.csv(NNObssub2, "~/RProjects/CarynB/NNobs_flowering_NY7d_7-28-22.csv")
########################################################################################################
# Data prep:
#NNObssub2 <- read_csv("~/RProjects/CarynB/NNobs_leaf_NY7d_7-28-22.csv")
#NNObssub2 <- read_csv("~/RProjects/CarynB/NNobs_flowering_NY7d_7-28-22.csv")
temp<-NNObssub2[,c("site_id","first_yes_year","common_name","first_yes_doy")]
# spread out obs years to cols
temp<-tidyr::spread(temp, first_yes_year, first_yes_doy)
# get site ids
ids<-unique(temp$site_id)
# loop through ids
allData <- list()
idx <- 1
pairList<-list()
for(i in 1:length(ids)){
  temp1 <- subset(temp, site_id==ids[i])
  # force common name sort
  #temp1 <- temp1[order(temp1$common_name), ]
  # get combos
  combos <- t(combn(temp1$common_name,2))
  combos <- rbind(combos,combos[,c(2,1)]) # I think this is what generates all of the duplicates
  for(j in 1:nrow(combos)){
    temp2 <- na.omit(as.data.frame(t(subset(temp1, common_name %in% combos[j,]))))
    # year/site data
    temp3<-temp2
    temp3$site_id<-temp3$V1[1]
    temp3$spp1<-temp3$V1[2]
    temp3$spp2<-temp3$V2[2]
    temp3<-temp3[3:nrow(temp3),]
    temp3$year<-as.numeric(rownames(temp3))
    colnames(temp3)[1:2]<-c("spp1DOY","spp2DOY")
    pairList[[idx]]<-temp3
    ##    
    #temp2 <- temp2[!(temp2$`1`==temp2$`2`),] # deletes years with equal DOY
    
    #pairList[[idx]]<-temp2
    temp2$first <- max.col(temp2)
    allData[[idx]] <- cbind.data.frame(ids[i],as.data.frame(t(combos[j,])),length(which(temp2$first==1))/(nrow(temp2)-1),(nrow(temp2)-1))
    idx <- idx+1
  }
}
allData<-do.call(rbind.data.frame, allData)
colnames(allData) <- c("site_id","spp1","spp2","prop","n")
allPairs<-do.call(rbind.data.frame, pairList)
# normalize values in AllData: take everything <0.5 and calc inverse of it (1-x)
allData$prop2 <- ifelse(allData$prop < 0.5, (1-allData$prop), allData$prop)    
# drop unneeded columns
final <- subset(allData, select = -c(prop))
# add lat/long back in
sites <- NNObssub2[,c("site_id","latitude","longitude")]
sites <- dplyr::distinct(sites)
# rename "scale"
names(final)[names(final) == "prop2"] <- "consistency"
final2 <- merge(final, sites, by = "site_id")
# calculate mean, SD of first yes DOY by species - can't use orig input file (NNobs...) because ties were dropped - shouldn't retain these
# instead - count #yrs per row in df "temp"
# subset NNobs by that one, THEN proceed w/steps below
# count up number of years with values for each row
temp$countyrs <- rowSums(temp>0, na.rm = TRUE)
temp$countyrs <- temp$countyrs - 2
# subset df temp to n > 4
temp <- temp[temp$countyrs > 4, ]
# subset NNobs_leafout_rank_NY7d_v4 by temp$countyrs
subset_leafout <- merge(x = NNObssub2, y = temp, by.x = c("site_id", "common_name"), by.y = c("site_id", "common_name"), all.y = TRUE)
# calculate site-level mean, SD of first yes DOY by species
SiteDOY <- subset_leafout %>%
  group_by(common_name, site_id) %>%
  dplyr::summarize(n_obs = n(),
                   SiteMeanDOY = mean(first_yes_doy),
                   SiteSD_DOY = sd(first_yes_doy))
# join this to final2 - first for spp1, then spp2
final3 <- merge(x = final2, y = SiteDOY, by.x = c("site_id", "spp1"), by.y = c("site_id", "common_name"), all.x = TRUE)
names(final3)[names(final3) == 'SiteMeanDOY'] <- "sp1meanDOY"
names(final3)[names(final3) == "SiteSD_DOY"] <- "sp1SD"
final4 <- merge(x = final3, y = SiteDOY, by.x = c("site_id", "spp2"), by.y = c("site_id", "common_name"), all.x = TRUE)
names(final4)[names(final4) == 'SiteMeanDOY'] <- "sp2meanDOY"
names(final4)[names(final4) == "SiteSD_DOY"] <- "sp2SD"
# drop unwanted columns
final4 <- subset(final4, select = -c(n_obs.x, n_obs.y))
# calculate DOYDiff 
final4$DOYdiff <- abs((final4$sp1meanDOY)-(final4$sp2meanDOY))
### I think we need to filter this on the "n" field - I think that is the count of years for a species pair at a site!!
# filter final4 to n > 4
final4 <- subset(final4, n > 4)
###########################################################################
# write out csv - ### THIS is file used in analyses - is leaner file because it removes "ties" 
# These ARE the final files used in the analyses
readr::write_csv(final4, "~/RProjects/CarynB/leaf_consistency_7-28-22.csv")
readr::write_csv(final4, "~/RProjects/CarynB/leaf_consistency_7-28-22-FULL.csv")
readr::write_csv(final4, "~/RProjects/CarynB/flowerining_consistency_7-28-22.csv")
readr::write_csv(final4, "~/RProjects/CarynB/flowerining_consistency_7-28-22-FULL.csv")
###################################################
# trying to generate data for parallel line plots
##### THESE SUMMARIES AREN"T QUITE RIGHT - STILL CAN BE IMPROVED TO FINALIZE
##### PARALLEL LINE PLOTS
# go back to NNObssub2
NNObssub2 <- NNObssub2[,c("site_id","first_yes_year","common_name","first_yes_doy")]
# need to drop ties, sites w/only 1 species, and species*site combos w/<5 *shared* years
# drop ties
NNObssub3 <- unique(setDT(NNObssub2), by = c("site_id", "first_yes_year", "first_yes_doy"))
# remove sites with < 2 spp*pp
#calculate number of spp by site
species_by_site <- NNObssub3 %>% 
  group_by(site_id, first_yes_year) %>% 
  dplyr::summarize(ct_species_by_site = n())
# join this table to NNObssub3
NNObssub4 <- merge(NNObssub3, species_by_site, by.x = c('site_id', 'first_yes_year'), by.y = c('site_id', 'first_yes_year'), all.x = TRUE) 
# drop sites with < 2 species
NNObssub4 <- subset(NNObssub4, ct_species_by_site > 1)
# and then remove site*species with < 5y of observations remaining
#calculate number of years of observations per species per site
yrs_by_site2 <- NNObssub4 %>% 
  group_by(site_id, common_name) %>% 
  dplyr::summarize(ct_yrs_by_site_species = n())
# join this table to NNObssub4
NNObssub5 <- merge(NNObssub4, yrs_by_site2, by.x = c('site_id', 'common_name'), by.y = c('site_id', 'common_name'), all.x = TRUE) 
# drop rows with years < 4 
NNObssub5 <- subset(NNObssub5, ct_yrs_by_site_species > 4)
# write out this file as a .csv - THIS IS THE FILE FOR FIG 1 - PARALLEL LINE PLOTS
readr::write_csv(NNObssub5, "~/RProjects/CarynB/leafpair_ordering_8-3-22.csv")
readr::write_csv(NNObssub5, "~/RProjects/CarynB/flowerpair_ordering_8-3-22.csv")
###########################################################################
#### combine leaf & flower pair sequence files to evaluate Q5
# WE DROPPED ALL OF THIS BECAUSE IT WASN"T WORKING
# have to slim down AllPairs using same criteria as above before can combine & analyze
# subset "allPairs" based on "final4"
# FIRST, generate "final4" above using FULL (don't remove inverse species pairs) - do this separately for leaf & bloom
# then, will join final4 to allPairs
# then, drop all rows in allPairs without info from final4 - this will be input "leaf" or "flower" dataframe to then merge 
#leafPairs <- merge(allPairs, final4, by.x = c('site_id', 'spp1', 'spp2'), by.y = c('site_id', 'spp1', 'spp2'), all.x = TRUE) 
#bloomPairs <- merge(allPairs, final4, by.x = c('site_id', 'spp1', 'spp2'), by.y = c('site_id', 'spp1', 'spp2'), all.x = TRUE) 
#subset based on "n" field - drop all NA's - do separately for leaf & bloom
#leafPairs2 <- leafPairs[!is.na(leafPairs$n),]
#bloomPairs2 <- bloomPairs[!is.na(bloomPairs$n),]
#leafPairs2$phenophase <- "leaf"
#bloomPairs2$phenophase <- "flower"
#leafPairs2$unique <- paste(leafPairs2$site_id,"-",leafPairs2$year,"-",leafPairs2$spp1,"-",leafPairs2$spp2)
#bloomPairs2$unique <- paste(bloomPairs2$site_id,"-",bloomPairs2$year,"-",bloomPairs2$spp1,"-",bloomPairs2$spp2)
#inner join
#join = merge(x = leafPairs2, y = bloomPairs2, by = "unique", all = FALSE)
#write out file as .csv
#readr::write_csv(join, "~/RProjects/CarynB/joined_phenophase_sequence_7-28-22.csv")
############## I'm NOT 100% CONFIDENT THIS GETS EVERYONE, BECAUSE SPECIES PAIRS COULD BE REVERSED......
# SO, trying this a second way - more like previously - join the two giant tables, THEN worry about removing
# drop ties, drop duplicates, sites w/<5y of obs, sites w/<2 species
# regenerate allPairs above, don't trim at all
#allPairsleaf <- allPairs
#allPairsbloom <- allPairs
#allPairsleaf$phenophase <- "leaf"
#allPairsbloom$phenophase <- "flower"
#allPairsleaf$unique <- paste(allPairsleaf$site_id,"-",allPairsleaf$year,"-",allPairsleaf$spp1,"-",allPairsleaf$spp2)
#allPairsbloom$unique <- paste(allPairsbloom$site_id,"-",allPairsbloom$year,"-",allPairsbloom$spp1,"-",allPairsbloom$spp2)
#join = merge(x = allPairsleaf, y = allPairsbloom, by = "unique", all = TRUE)
#join = merge(x = allPairsleaf, y = allPairsbloom, by = "unique", all = FALSE)
#readr::write_csv(join, "~/RProjects/CarynB/joined_phenophase_sequence_7-28-22-BIG.csv")
############################################################################
# further analyses
# comparing consistency to DOY
# revisited 8-2-22
# 
# read file back in, if needed
leafdata <- read_csv("~/RProjects/CarynB/leaf_consistency-dropped_sites-for_consis_vs_DOYdiff_in_R.csv")
# calc averages where >1 site*spp pair instances exist (only forsythia-flowering dogwood for leaf-out?)
#leafdata2 <- leafdata %>%
#              group_by(spp2, spp1) %>%
#              dplyr::summarize(n_obs = n(),
#                               mean_consistency = mean(consistency),
#                               mean_sp1meanDOY = mean(sp1meanDOY),
#                               mean_DOYdiff = mean(DOYdiff),
#                               mean_latitude = mean(latitude))

### comparing consistency to diffDOY
# for consistency values, convert <0 values to 0
leafdata$binary <- ifelse(leafdata$consistency <1, 0, 1)
# logistic regression between consistency and DiffDOY - chose log reg because data so skewed
model <- glm(binary ~ DOYdiff, data = leafdata, family = binomial)
summary(model)$coef
confint(model)
# plot
leafdata %>%
  ggplot(aes(DOYdiff, binary)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  labs(
    title = "Logistic Regression Model", 
    x = "DOYdiff",
    y = "consistency"
  )
### Try log transforming consistency - doesn't really do the trick
leafdata$logconsis <- log(leafdata$consistency)
hist(leafdata$logconsis)
###############################
# consistency by latitude
model <- glm(binary ~ mean_latitude, data = leafdata2, family = binomial)
summary(model)$coef
confint(model)
# plot
leafdata2 %>%
  ggplot(aes(mean_latitude, binary)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  labs(
    title = "Logistic Regression Model", 
    x = "mean_DOYdiff",
    y = "consistency"
  )
# comparing consistency to DOY
# 
# read file back in, if needed
bloomdata <- read_csv("~/RProjects/CarynB/floweringconsistency_4-25-22.csv")
# calc averages where >1 site*spp pair instances exist 
bloomdata2 <- bloomdata %>%
  group_by(spp2, spp1) %>%
  dplyr::summarize(n_obs = n(),
                   mean_consistency = mean(consistency),
                   mean_sp1meanDOY = mean(sp1meanDOY),
                   mean_DOYdiff = mean(DOYdiff),
                   mean_latitude = mean(latitude))
# for consistency values, convert <0 values to 0
bloomdata2$binary <- ifelse(bloomdata2$mean_consistency <1, 0, 1)
### comparing consistency to diffDOY
# logistic regression between consistency and DOY (= site*spp-pair mean DOY)
model2 <- glm(binary ~ mean_DOYdiff, data = bloomdata2, family = binomial)
summary(model2)$coef
confint(model2)
## plot
bloomdata2 %>%
  ggplot(aes(mean_DOYdiff, binary)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  labs(
    title = "Logistic Regression Model", 
    x = "mean_DOYdiff",
    y = "consistency"
  )
###################################
# consistency by latitude
model <- glm(binary ~ mean_latitude, data = bloomdata2, family = binomial)
summary(model)$coef
confint(model)
# plot
bloomdata2 %>%
  ggplot(aes(mean_latitude, binary)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  labs(
    title = "Logistic Regression Model", 
    x = "mean_latitude",
    y = "consistency"
  )
