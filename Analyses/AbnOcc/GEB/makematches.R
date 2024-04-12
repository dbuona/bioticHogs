###environmental matching####
rm(list=ls())
options(stringsAsFactors = FALSE)
graphics.off()

library(dplyr)
library(tidyr)
library("MatchIt")

d<-read.csv("Analyses/AbnOcc/GEB/envplotdata.csv") ### read in data
d<-d %>% group_by(Plot) %>% filter(Year==max(Year)) %>% ungroup() 2###include only more recent surveys from repeat sampled plots

d$status<-ifelse(d$RelCov_I>0,1,0) ## assign native status

d.n<-filter(d,status==0) ### filter out natives

d.i<-filter(d,status==1) ### filter out non-natives
d.i<-filter(d.i,RelCov_I>5) ## include only plots with more than 5% non-native cover

d<-rbind(d.i,d.n) ## rejoin native and non-native plots

counters<-d %>% group_by(L4_KEY)%>%count() ## count the number of plots per L4 ecoregion
counters<-filter(counters,n>1) ### select only ecoregions with at least 2 plots

d<-filter(d,L4_KEY %in% counters$L4_KEY) ### filter out regions with less than 2 plots
d<- d[complete.cases(d), ] # select only compelte cases (needed for matching) 57,535

set.seed(111)
m2<-matchit(status~L4_KEY+
              NDMI_median_Apr1_Sept30_1985_2020_90m.tif+
              TMIN_Dec_Feb_1981_2018_mean.tif+
              SG_n_tot_M_sl2_100m.tif+
              VCF_percent_treeCover_2000_2016_mean_integer.tif+
              gHM.tif,data=d,method = "nearest",distance="glm",exact=~L4_KEY) ## match plots based on environmental variables

mdat2<-match.data(m2) ## covert to data frame #23,198 plot matches


#### make plot S1#########
plot(m2, type = "density", interactive = FALSE,
     which.xs = ~NDMI_median_Apr1_Sept30_1985_2020_90m.tif+ TMIN_Dec_Feb_1981_2018_mean.tif + SG_n_tot_M_sl2_100m.tif+VCF_percent_treeCover_2000_2016_mean_integer.tif+ gHM.tif)
############################




write.csv(mdat2,"Analyses/AbnOcc/GEB/matchedplots.csv",row.names = FALSE)
