###match making for knb


rm(list=ls())
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())
#rstan_options(auto_write = TRUE)

graphics.off()
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)

library(ape)
library(geodist)
library("MatchIt")

setwd("/Users/danielbuonaiuto/Documents/git/bioticHogs/")
d<-read.csv("Analyses/AbnOcc/GEB/envplotdata.csv") ### read in data

set.seed(111)
m2<-matchit(status~L4_KEY+Dataset+
              NDMI_median_Apr1_Sept30_1985_2020_90m.tif+
              TMIN_Dec_Feb_1981_2018_mean.tif+
              SG_n_tot_M_sl2_100m.tif+
              VCF_percent_treeCover_2000_2016_mean_integer.tif+
              gHM.tif,data=d,method = "nearest",distance="glm",exact=~L4_KEY+Dataset) ## match plots based on environmental variables

mdat2<-match.data(m2) ## covert to data frame #20,900 plot matches


#### make plot S1#########
plot(m2, type = "density", interactive = FALSE,
     which.xs = ~NDMI_median_Apr1_Sept30_1985_2020_90m.tif+ TMIN_Dec_Feb_1981_2018_mean.tif + SG_n_tot_M_sl2_100m.tif+VCF_percent_treeCover_2000_2016_mean_integer.tif+ gHM.tif)

write.csv(mdat2,"Analyses/AbnOcc/GEB/env_matchedplots.csv",row.names = FALSE) ##write matched data frame

rm(d)
#############################
d<-read.csv("Analyses/AbnOcc/GEB/spatplotdata.csv") ### read in data
d$unit<-paste(d$L4_KEY,d$Dataset,sep="_")
counters2<-d %>% group_by(L4_KEY,Dataset)%>%count()
counterss2<-filter(counters2,n>1)
counterss2$unit<-paste(counterss2$L4_KEY,counterss2$Dataset,sep="_")
d3<-filter(d,unit %in% counterss2$unit)

d.n<-filter(d,status==0)
d.i<-filter(d,status==1)

pairz.data<-data.frame()

ecoregionz<-sort(unique(counterss2$unit),decreasing = TRUE)
head(ecoregionz,10)
#ecoregionz<-ecoregionz[7:1199]


for (z in c(1:length(ecoregionz))){
  goo<-filter(d3,unit==ecoregionz[z]) 
  cordz<-data.frame(Lon=goo$Long,Lat=goo$Lat)
  disty1<-geodist(cordz,measure="haversine")
  rownames(disty1)<-goo$Plot
  colnames(disty1)<-goo$Plot
  diss<-melt(disty1)
  diss<-filter(diss,Var1!=Var2)
  diss<-filter(diss,Var2!=Var1)
  colnames(diss)<-c("plot1","plot2","distance")
  diss<-filter(diss,plot1 %in% unique(d.n$Plot))###make plot1 native
  diss<-filter(diss,plot2 %in% unique(d.i$Plot)) 
  diss<-diss %>% dplyr::group_by(plot2) %>% dplyr::slice(which.min(distance)) %>%ungroup() ###
  #diss<-filter(diss,plot1 %in% d.n$Plot) ### fitler this column to only native
  #diss<-filter(diss,plot2%in% d.i$Plot) ### filter this column to only invaded
  
  
  diss$Distance<-diss$distance/1000
 
  diss$unit<-ecoregionz[z]

  #pairz2<-filter(pairz2,isna==FALSE)
  
  
  pairz.data<-rbind(diss,pairz.data)
  print(paste("done",which(ecoregionz==ecoregionz[z]))) 
}

pairz2<-pairz.data %>% dplyr::filter(Distance<200) #se the min distance to every invaded
mean(pairz2$Distance)
median(pairz2$Distance)
sd(pairz2$Distance)


pairz2<-tibble::rownames_to_column(pairz2,"match2")
pairz2$match<-paste("match",pairz2$match2,sep="_")


pairzclose<-pairz2
median(pairzclose$Distance)

pairz.t<-dplyr::select(pairzclose,match,unit,plot1)
pairz.t$status<-0
pairz.t2<-dplyr::select(pairzclose,match,unit,plot2)
pairz.t2$status<-1

colnames(pairz.t)[3]<-"Plot"
colnames(pairz.t2)[3]<-"Plot"

distpair<-rbind(pairz.t,pairz.t2)
nrow(distpair)/2

write.csv(distpair,"Analyses/AbnOcc/GEB/spat_matchedplots.csv",row.names = FALSE) ##write matched data frame
