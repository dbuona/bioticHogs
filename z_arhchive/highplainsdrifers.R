##High Plains Drifters data generation.

rm(list=ls())
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#graphics.off()

library(dplyr)
library(vegan)
library("geodist")
library(ggplot2)
library(hillR)
library(betapart)
library(reshape2)

setwd("~/Desktop/Powell/")# for mac
set.seed(3)

###read in data
d1<-read.csv("Data/FULLDatabase_05272022.csv")
d.env1<-read.csv("Data/FULLDatabase_05272022_plotsenv4Jul2022.csv")
d.ref<-read.csv("Data/FullDatabase_diversity_July4.csv")

unique(d.env1$US_L3NAME)
d.env<-filter(d.env1, US_L3NAME %in% c("High Plains","Driftless Area", "Erie Drift Plain","Southern Michigan/Northern Indiana Drift Plains"))
d<-filter(d1, Plot %in% unique(d.env$Plot)) ### filter the main dataset

#### deal with resampled plots
table(d$Resampled)

d<- d%>% group_by(Plot) %>% 
  filter(Year==max(Year)) %>%
  ungroup() ### this take resapled plots and chooses the most recent

d.ref2<-filter(d.ref,Plot %in% c(unique(d$Plot))) ### this filter the invasion levels data

d.ref2<-d.ref2 %>%
  group_by(Plot) %>%
  filter(Year==max(Year)) %>%
  ungroup() #  Same as above, only one on observation per plot
#################################

table(d.ref2$Dataset) ### here in theory you could select for a dataset though i dont

###NWCA are likly dulicatied

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

distancize<-function(x,y){
  temp2<-dplyr::select(x,Plot,Original.Long,Original.Lat) ## subset to lat longs
  temp2<-dplyr::distinct(temp2) ## remove dplicates casue each species was a row
  remove<-setdiff(temp2$Plot,rownames(y))
  temp2<-filter(temp2,!Plot %in% c(remove))
  setNames(data.frame(temp2),c("Plot","Long","Lat"))
  
}
datacize<-function(x,y,z){
  test<-as.matrix(x$beta.bray)
    test[upper.tri(test)] <- NA
   test<- setNames(melt(test,na.rm = TRUE),c("site1","site2", "beta_div"))

  test2<-as.matrix(x$beta.bray.bal)
  test2[upper.tri(test2)] <- NA
  test2<-setNames(melt(test2,na.rm=TRUE),c("site1","site2", "turnover"))
  
  test3<-as.matrix(x$beta.bray.gra)
  test3[upper.tri(test3)] <- NA
    
  test3<-setNames(melt(test3,na.rm=TRUE),c("site1","site2", "nestedness"))
  test<-left_join(test,test2)
  test<-left_join(test,test3)
  
  testD<-as.matrix(y)
  testD[upper.tri(testD)] <- NA
  testD<-setNames(melt(testD,na.rm=TRUE),c("site1","site2", "Distance"))
  test$Distance<-testD$Distance
  test$metric<-"abundance"
  
  testo<-as.matrix(z$beta.bray)
  testo[upper.tri(testo)] <- NA
  testo<- setNames(melt(testo,na.rm=TRUE),c("site1","site2", "beta_div"))
  
  testo2<-as.matrix(z$beta.bray.bal)
  testo2[upper.tri(testo2)] <- NA
  testo2<-setNames(melt(testo2,na.rm=TRUE),c("site1","site2", "turnover"))
  
  testo3<-as.matrix(z$beta.bray.gra)
  testo3[upper.tri(testo3)] <- NA
  
  testo3<-setNames(melt(testo3,na.rm=TRUE),c("site1","site2", "nestedness"))
  testo<-left_join(testo,testo2)
  testo<-left_join(testo,test3)
  
  testoD<-as.matrix(y)
  testoD[upper.tri(testoD)] <- NA
  testoD<-setNames(melt(testoD,na.rm=TRUE),c("site1","site2", "Distance"))
  testo$Distance<-testoD$Distance
  testo$metric<-"occurance"
  
 
  
  
  test<-rbind(test,testo)
}


mydat<-data.frame(site1=factor(), site2=factor(), beta_div=factor(), turnover=numeric(), nestedness=numeric(),
                  distance=numeric() ,metric=character(),status=character(),status=character,Ecoregion=character())

###loop 1
ecoregionz<-unique(d.env$NA_L3NAME) ##one has no numeric data
for (z in c(1:length(ecoregionz))){
  
  tlands<-dplyr::filter(d.env,US_L3NAME==ecoregionz[z]) ##
  
  dat<-filter(d, Plot %in% c(unique(tlands$Plot)))
  comms<-matricize(dat)
  occs<-ifelse(comms==0,comms,1)
  comms.out.abun<-beta.pair.abund(comms)
  comms.out.occs<-beta.pair.abund(occs)
  disty<-distancize(dat,comms)
  b<-geodist(disty,paired = TRUE,measure = "haversine")
  daty<-datacize(comms.out.abun,b,comms.out.occs)
  daty$status<-"introduced species included"
  daty$Ecoregion=ecoregionz[z]
  mydat<-rbind(daty,mydat)
  print(paste("done",which(ecoregionz==ecoregionz[z])))    
}

###601724

d.native<-filter(d,NativeStatus!="I")

mydat2<-data.frame(site1=factor(), site2=factor(), beta_div=factor(), turnover=numeric(), nestedness=numeric(),
                   distance=numeric() ,metric=character(),status=character(),status=character,Ecoregion=character())

###loop 2

for (z in c(1:length(ecoregionz))){
  
  tlands<-dplyr::filter(d.env,US_L3NAME==ecoregionz[z]) ##
  
  dat<-filter(d.native, Plot %in% c(unique(tlands$Plot)))
  comms<-matricize(dat)
  occs<-ifelse(comms==0,comms,1)
  comms.out.abun<-beta.pair.abund(comms)
  comms.out.occs<-beta.pair.abund(occs)
  disty<-distancize(dat,comms)
  b<-geodist(disty,paired = TRUE,measure = "haversine")
  daty<-datacize(comms.out.abun,b,comms.out.occs)
  daty$status<-"native species only"
  daty$Ecoregion=ecoregionz[z]
  mydat2<-rbind(daty,mydat2)
  print(paste("done",which(ecoregionz==ecoregionz[z])))    
}


#### how important are the zeros


rm(mydat3)
mydatPlot<-rbind(mydat,mydat2)
mydatPlot$Distance<-mydatPlot$Distance/100
mydatPlot<-filter(mydatPlot,site1!=site2)

write.csv(mydatPlot,"highplainsdrifterdata.csv",row.names = FALSE)




#mydatPlot$DistanceR<-ifelse(mydatPlot$Distance==0,runif(length(mydatPlot$Distance==0),0,1),mydatPlot$Distance)

ggplot(mydatPlot,aes(Distance,beta_div))+geom_smooth(aes(color=status,fill=status,linetype=metric),method="glm",method.args = list(family = "quasibinomial"),size=1.2)+
  facet_wrap(~Ecoregion,scales="free")+ggthemes::theme_few(base_size = 11)+ylab("Bray-Curtis Dissimilarity")+xlab("Distance (km)")+
  scale_color_viridis_d(begin = .1,end=.7)+scale_fill_viridis_d(begin = .1,end=.7)

ggplot(mydatPlot,aes(DistanceR,beta_div))+geom_smooth(aes(color=status,fill=status,linetype=metric),method="glm",method.args = list(family = "quasibinomial"),size=1.2)+
  facet_wrap(~Ecoregion,scales="free")+ggthemes::theme_few(base_size = 11)+ylab("Bray-Curtis Dissimilarity")+xlab("Distance (km)")+
  scale_color_viridis_d(begin = .1,end=.7)+scale_fill_viridis_d(begin = .1,end=.7)

ggplot(mydatPlot2,aes(Distance,beta_div))+geom_smooth(aes(color=status,fill=status,linetype=metric),method="glm",method.args = list(family = "quasibinomial"),size=1.2)+
  facet_wrap(~Ecoregion,scales="free")+ggthemes::theme_few(base_size = 11)+ylab("Bray-Curtis Dissimilarity")+xlab("Distance (km)")+
  scale_color_viridis_d(begin = .1,end=.7)+scale_fill_viridis_d(begin = .1,end=.7)


head(mydat)
head(mydat2)
mydat3<-mydat2
mydat3<-select(mydat3,-status,-Distance)
colnames(mydat3)[3:5]<-paste(colnames(mydat3)[3:5],"native",sep = ".")

mydat4<-left_join(mydat3,mydat,by=c("site1","site2","metric","Ecoregion"))

ggplot(mydat4,aes(beta_div.native,beta_div))+geom_point()+
  geom_smooth(method="lm",color="royalblue")+
  facet_grid(metric~Ecoregion)+geom_abline(intercept = 0,slope=1,color="red")+ggthemes::theme_few()+
  ylab("Pairwise dissimilarity native and introdiced")+xlab("Pairwise dissimilarity native only")

ggplot(mydat4,aes(turnover.native,turnover))+geom_point()+
  geom_smooth(method="lm",color="royalblue")+
  facet_grid(metric~Ecoregion)+geom_abline(intercept = 0,slope=1,color="red")

ggplot(mydat4,aes(nestedness.native,nestedness))+geom_point()+
  geom_smooth(method="lm",color="royalblue")+
  facet_grid(metric~Ecoregion)+geom_abline(intercept = 0,slope=1,color="red")
