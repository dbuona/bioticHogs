#abundance vs. occurance

rm(list=ls())
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

graphics.off()
library(rstan)
library(dplyr)
library(vegan)
library("geodist")
library(ggplot2)
library(hillR)
library(betapart)
library(reshape2)


setwd("~/PowellCenter/")# for PC
load("46ers.Rda")
set.seed(3)

###read in data
d1<-read.csv("Data/FULLDatabase_05272022.csv")
d.env1<-read.csv("Data/FULLDatabase_05272022_plotsenv4Jul2022.csv")
d.ref<-read.csv("Data/FullDatabase_diversity_July4.csv")

counters<- d.env1%>% group_by(US_L3NAME) %>% count() %>% arrange(n)
hundo<-filter(counters,n>100)
hundo<-filter(hundo,n<1000)

d.env<-filter(d.env1, US_L3NAME %in% c(hundo$US_L3NAME))

d<-filter(d1, Plot %in% unique(d.env$Plot)) 

table(d$Resampled)

d<- d%>% group_by(Plot) %>% 
  filter(Year==max(Year)) %>%
  ungroup() ### this take resapled plots and chooses the most recent

d.ref2<-filter(d.ref,Plot %in% c(unique(d$Plot))) ### this filter the invasion levels data

d.ref2<-d.ref2 %>%
  group_by(Plot) %>%
  filter(Year==max(Year)) %>%
  ungroup() #  Same as 

matricize<-function(x){
  temp<-dplyr::select(x,Plot,AcceptedTaxonName,PctCov_100) ### select useful columns
  temp<-dplyr::filter(temp,!is.na(temp$AcceptedTaxonName))
  temp<- tidyr::spread(temp,AcceptedTaxonName,PctCov_100)
  temp[is.na(temp)] <- 0 
  temp<-as.data.frame(temp)
  rownames(temp)<-temp$Plot # convert plot names to row names
  temp<-dplyr::select(temp,-Plot)
  temp<-as.matrix.data.frame(temp)# remove plot column
}


mydat<-data.frame(site1=character(),site2=character(),TD_beta=numeric(),local_similarity=numeric(),Ecoregion=character())
ecoregionz<-c(unique(d.env$US_L3NAME))


for (z in c(1:length(ecoregionz))){
  
  tlands<-dplyr::filter(d.env,US_L3NAME==ecoregionz[z]) ##
  dat<-filter(d, Plot %in% c(unique(tlands$Plot)))
  comms<-matricize(dat)
  
  daty<-hill_taxa_parti_pairwise(comms,q = 0)
  daty<-dplyr::select(daty,site1,site2, TD_beta,local_similarity)
  daty$Ecoregion<-ecoregionz[z]
  mydat<-rbind(daty,mydat)
  print(paste("done",which(ecoregionz==ecoregionz[z]))) 
}
save.image("~/PowellCenter/46ers.Rda")

d.native<-filter(d,NativeStatus!="I")

mydat2<-data.frame(site1=character(),site2=character(),TD_beta=numeric(),local_similarity=numeric(),Ecoregion=character())
for (z in c(1:length(ecoregionz))){
  
  tlands<-dplyr::filter(d.env,US_L3NAME==ecoregionz[z]) ##
  dat<-filter(d.native, Plot %in% c(unique(tlands$Plot)))
  comms<-matricize(dat)
  
  daty<-hill_taxa_parti_pairwise(comms,q = 0)
  daty<-dplyr::select(daty,site1,site2,TD_beta,local_similarity)
  daty$Ecoregion<-ecoregionz[z]
  mydat2<-rbind(daty,mydat2)
  print(paste("done",which(ecoregionz==ecoregionz[z]))) 
}



mydat3<-data.frame(site1=character(),site2=character(),TD_beta=numeric(),local_similarity=numeric(),Ecoregion=character())
ecoregionz<-c(unique(d.env$US_L3NAME))


for (z in c(1:length(ecoregionz))){
  
  tlands<-dplyr::filter(d.env,US_L3NAME==ecoregionz[z]) ##
  dat<-filter(d, Plot %in% c(unique(tlands$Plot)))
  comms<-matricize(dat)
  
  daty<-hill_taxa_parti_pairwise(comms,q = 2)
  daty<-dplyr::select(daty,site1,site2, TD_beta,local_similarity)
  daty$Ecoregion<-ecoregionz[z]
  mydat3<-rbind(daty,mydat3)
  print(paste("done",which(ecoregionz==ecoregionz[z]))) 
}


mydat4<-data.frame(site1=character(),site2=character(),TD_beta=numeric(),local_similarity=numeric(),Ecoregion=character())
for (z in c(1:length(ecoregionz))){
  
  tlands<-dplyr::filter(d.env,US_L3NAME==ecoregionz[z]) ##
  dat<-filter(d.native, Plot %in% c(unique(tlands$Plot)))
  comms<-matricize(dat)
  
  daty<-hill_taxa_parti_pairwise(comms,q = 2)
  daty<-dplyr::select(daty,site1,site2,TD_beta,local_similarity)
  daty$Ecoregion<-ecoregionz[z]
  mydat4<-rbind(daty,mydat4)
  print(paste("done",which(ecoregionz==ecoregionz[z]))) 
}

colnames(mydat)[3]<-"Beta_all_species"
colnames(mydat2)[3]<-"Beta_native_only"

colnames(mydat)[4]<-"all_species"
colnames(mydat2)[4]<-"native_only"


mydatOcc<-left_join(mydat2,mydat)

#ggplot(mydatOcc,aes(native_only,all_species))+geom_point(size=.1)+geom_abline(intercept=0,slope=1,color="royalblue")+facet_wrap(~Ecoregion)



mydatOcc$H<-log((mydatOcc$all_species+.001)/(mydatOcc$native_only+.001))

tail(mydatOcc$H)





colnames(mydat3)[3]<-"Beta_all_species"
colnames(mydat4)[3]<-"Beta_native_only"

colnames(mydat3)[4]<-"all_species"
colnames(mydat4)[4]<-"native_only"
mydatAbn<-left_join(mydat4,mydat3)
mydatAbn$H<-log(((mydatAbn$all_species+.001)/(mydatAbn$native_only+.001)))

mydatAbn$metric<-"Abn"
mydatOcc$metric<-"Occ"

mydatAll<-rbind(mydatAbn,mydatOcc)

coords<-select(d.native,Plot,Lat,Long)
coords<-distinct(coords)
colnames(coords)<-c("site1", "Lat.1", "Lon.1")
coords2<-coords
colnames(coords2)<-c("site2", "Lat.2", "Lon.2")

mydatAll<-left_join(mydatAll,coords)
mydatAll<-left_join(mydatAll,coords2)

cord1<-data.frame(Lon=mydatAll$Lon.1,Lat=mydatAll$Lat.1)
cord2<-data.frame(Lon=mydatAll$Lon.2,Lat=mydatAll$Lat.2)

disy1<-geodist(x=cord1,y=cord2,measure="haversine",paired = TRUE)

mydatAll$Distance<-disy1/100

ggplot(mydatAll,aes(metric,H,color=Ecoregion))+stat_summary(fun.data="mean_cl_boot")+
  theme(legend.position="none")

library(brms)

mod1<-brm(H~Ecoregion*metric,data=mydatAll)
