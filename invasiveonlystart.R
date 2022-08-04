rm(list=ls())
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

graphics.off()

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
d<-read.csv("Data/FULLDatabase_05272022.csv")

di<-filter(d,NativeStatus=="I")

table(di$Resampled)

di<- di%>% group_by(Plot) %>% 
  filter(Year==max(Year)) %>%
  ungroup() ### this take resapled plots and chooses the most recent


nrow(d)/nrow(di)

d.env<-read.csv("Data/FULLDatabase_05272022_plotsenv4Jul2022.csv")
d.ref<-read.csv("Data/FullDatabase_diversity_July4.csv")

di.env<-filter(d.env, Plot %in% unique(di$Plot))
#counters<- di.env %>% group_by(US_L3NAME) %>% count() %>% arrange(n)



ecoregionz<-unique(di.env$US_L3NAME)


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
  test<-setNames(melt(as.matrix(x$beta.bray)),c("site1","site2", "beta_div"))
  test2<-setNames(melt(as.matrix(x$beta.bray.bal)),c("site1","site2", "turnover"))
  test3<-setNames(melt(as.matrix(x$beta.bray.gra)),c("site1","site2", "nestedness"))
  test<-left_join(test,test2)
  test<-left_join(test,test3)
  testD<-setNames(melt(as.matrix(y)),c("site1","site2", "Distance"))
  test$Distance<-testD$Distance
  test$metric<-"abundance"
  
  testo<-setNames(melt(as.matrix(z$beta.bray)),c("site1","site2", "beta_div"))
  testo2<-setNames(melt(as.matrix(z$beta.bray.bal)),c("site1","site2", "turnover"))
  testo3<-setNames(melt(as.matrix(z$beta.bray.gra)),c("site1","site2", "nestedness"))
  testo<-left_join(testo,testo2)
  testo<-left_join(testo,testo3)
  testo$Distance<-test$Distance
  testo$metric<-"occurance"
  
  test<-rbind(test,testo)
}


mydat<-data.frame(site1=factor(), site2=factor(), beta_div=factor(), turnover=numeric(), nestedness=numeric(),
                  distance=numeric() ,metric=character(),status=character(),Ecoregion=character())

###loop 1
for (z in c(1:length(ecoregionz))){
  
  tlands<-dplyr::filter(di.env,US_L3NAME==ecoregionz[z]) ##
  
  dat<-filter(di, Plot %in% c(unique(tlands$Plot)))
  comms<-matricize(dat)
  occs<-ifelse(comms==0,comms,1)
  comms.out.abun<-beta.pair.abund(comms)
  comms.out.occs<-beta.pair.abund(occs)
  disty<-distancize(dat,comms)
  b<-geodist(disty,paired = TRUE,measure = "haversine")
  daty<-datacize(comms.out.abun,b,comms.out.occs)
  daty$Ecoregion=unique(tlands$US_L3NAME)
        mydat<-rbind(daty,mydat)
      print(paste("done",which(ecoregionz==ecoregionz[z])))    
      }

mydat$Distance<-mydat$Distance/100
ggplot(mydat,aes(Distance,beta_div))+geom_smooth(aes(linetype=metric),method="glm",method.args = list(family = "quasibinomial"),size=1.2)+
  facet_wrap(~Ecoregion,scales="free_x")+ggthemes::theme_few()


colnames(mydat)
colnames(daty)
length(as.vector(b))
length(as.vector(goober$beta.bray.gra))
ggplot(mydat,aes(Distance,beta_div))+geom_smooth(method="glm",method.args = list(family = "quasibinomial"),size=1.2)
  geom_smooth(aes(Distance,turnover),method="glm",method.args = list(family = "quasibinomial"),size=1.2,color="red")+
  geom_smooth(aes(Distance,nestedness),method="glm",method.args = list(family = "quasibinomial"),size=1.2,color="green")


