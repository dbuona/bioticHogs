###start cleaning concept papers, match plots by characteristics

rm(list=ls())
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())
#rstan_options(auto_write = TRUE)

graphics.off()
library(dplyr)
library(tidyr)

library(ggplot2)

library(ape)
library(geodist)
library("MatchIt")

setwd("/Users/danielbuonaiuto/Documents/git/bioticHogs/")
### atempted plot matching
##data prep

##3read in plot summeries
d<-read.csv("Data/FULLDatabase_10272022_diversity_PctCov_100_10272022_INHABIT.csv")
d2<-dplyr::select(d,gHM.tif,Dataset,Plot,Year,RelCov_I,NDMI_median_Apr1_Sept30_1985_2020_90m.tif,SG_n_tot_M_sl2_100m.tif,TMIN_Dec_Feb_1981_2018_mean.tif,VCF_percent_treeCover_2000_2016_mean_integer.tif,gHM.tif,L4_KEY)
write.csv(d2,"Analyses/AbnOcc/GEB/envplotdata.csv",row.names = FALSE)
d<-d %>%
  group_by(Plot) %>%
  filter(Year==max(Year)) %>%
  ungroup() 

d$status<-ifelse(d$RelCov_I>0,1,0)

d.n<-filter(d,status==0)
d.i<-filter(d,status==1)
d.i<-filter(d.i,RelCov_I>5)
d<-rbind(d.i,d.n)


d$Richness_tot<-d$Richness_I+d$Richness_N+d$Richness_Unk


#qants<-d %>%group_by (L4_KEY) %>% summarise(Qs=round(quantile(RelCov_I,na.rm=TRUE,prob=c(0,.25)),1))
#qants$Q<-rep(c("Q25","Q75"),length(unique(d$L4_KEY)))
#qants<-tidyr::spread(qants,"Q","Qs")
#qants<-filter(qants,Q25!=Q75)

#d<-left_join(d,qants)
#d<-filter(d,!is.na(Q75))
#d$status<-NA
#d$status[which(d$RelCov_I<=d$Q25)]<- 0
#d$status[which(d$RelCov_I>=d$Q75)]<- 1

d2<-dplyr::select(d,Site,Lat,Long,gHM.tif,Dataset,Plot,Year,status,Richness_tot,NDMI_median_Apr1_Sept30_1985_2020_90m.tif,SG_n_tot_M_sl2_100m.tif,TMIN_Dec_Feb_1981_2018_mean.tif,VCF_percent_treeCover_2000_2016_mean_integer.tif,gHM.tif,L4_KEY)
#d2<-dplyr::select(d,Site,Lat,Long,Dataset,Plot,Year,status,Richness_tot,NDMI_median_Apr1_Sept30_1985_2020_90m.tif,TMIN_Dec_Feb_1981_2018_mean.tif,gHM.tif,L4_KEY)
d2<- d2[complete.cases(d2), ] ##55809

counters<-d2 %>% group_by(L4_KEY)%>%count()
counterss<-filter(counters,n>1)
d2<-filter(d2,L4_KEY %in% counterss$L4_KEY)

set.seed(111)
m2<-matchit(status~L4_KEY+
                  NDMI_median_Apr1_Sept30_1985_2020_90m.tif+
                  TMIN_Dec_Feb_1981_2018_mean.tif+
                  SG_n_tot_M_sl2_100m.tif+
                  VCF_percent_treeCover_2000_2016_mean_integer.tif+
                  gHM.tif,data=d2,method = "nearest",distance="glm",exact=~L4_KEY)

mdat2<-match.data(m2)###22,627 plots


?matchit()
pdf("Analyses/AbnOcc/propplot.pdf")
plot(m2, type = "density", interactive = FALSE,
     which.xs = ~NDMI_median_Apr1_Sept30_1985_2020_90m.tif+ TMIN_Dec_Feb_1981_2018_mean.tif + SG_n_tot_M_sl2_100m.tif+VCF_percent_treeCover_2000_2016_mean_integer.tif+ gHM.tif)
dev.off()

envpair<-dplyr::select(mdat2,Plot,L4_KEY,status,subclass)
write.csv(envpair,"Analyses/AbnOcc/Input/env_matched_plots.csv")

