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
set.seed(Togetheree3) 
load("NOREASTERN.Rda")
###read in data
d1<-read.csv("Data/FULLDatabase_05272022.csv")
d.env1<-read.csv("Data/FULLDatabase_05272022_plotsenv4Jul2022.csv")
d.ref<-read.csv("Data/FullDatabase_diversity_July4.csv")

###FYI what datsets are resampled
d.env.re<-filter(d.env1,Resampled=="Y")
resampledd<-d.env.re%>% group_by(Plot,Dataset) %>% count()

table(resampledd$Dataset)
counters<- d.env1%>% group_by(NA_L1NAME,NA_L2NAME,US_L3NAME) %>% count() %>% arrange(n)

d.env<-filter(d.env1,NA_L2CODE %in% c(8.1,5.3,5.2,8.2))
#d.env<-filter(d.env,Dataset=="AIM") ### for now we wont deal with inter

counters2<- d.env%>% group_by(US_L4NAME) %>% count() %>% arrange(n)

d<-filter(d1, Plot %in% unique(d.env$Plot)) 
unique(d.env$US_L3NAME)
library(maps)

jpeg("~/git/bioticHogs/desertsolitaire/NEforestmap.jpeg",width = 9, height= 8, units="in",res=300 )
par(mar = c(.1, .1, .1, .1))

map("state", interior = FALSE)
points(d.env$Lon, d.env$Lat, pch=15, cex=.2,col=as.factor(d.env$NA_L1NAME))
par(xpd=TRUE)
legend("topright",title="Ecoregions: Level I", legend=unique(d.env$NA_L1NAME),inset=c(.5,0.1), ncol=1,pch=15,cex=0.7, col=unique(as.factor(d.env$NA_L1NAME)))
dev.off()
d<- d%>% group_by(Plot) %>% 
  filter(Year==max(Year)) %>%
  ungroup() ### this take resapled plots and chooses the most recent

d.ref2<-filter(d.ref,Plot %in% c(unique(d$Plot))) ### this filter the invasion levels data

d.ref2<-d.ref2 %>%
  group_by(Plot) %>%
  filter(Year==max(Year)) %>%
  ungroup() #  Same as above, only one on observation per plot
#################################

table(d.ref2$Dataset)
colnames(d)

matricize<-function(x){
  temp<-dplyr::select(x,Plot,AcceptedTaxonName,PctCov_100) ### select useful columns
  temp<-dplyr::filter(temp,!is.na(temp$AcceptedTaxonName))
  temp<- tidyr::spread(temp,AcceptedTaxonName,PctCov_100)
  temp[is.na(temp)] <- 0 
  temp<-as.data.frame(temp)
  rownames(temp)<-temp$Plot # convert plot names to row names
  temp<-dplyr::select(temp,-Plot)
  temp<-as.matrix.data.frame(temp)#

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


colnames(mydat)[3]<-"Beta_all_species"
colnames(mydat2)[3]<-"Beta_native_only"

colnames(mydat)[4]<-"all_species"
colnames(mydat2)[4]<-"native_only"


mydatOcc<-left_join(mydat2,mydat)

mydatOcc$H_sim<-log((mydatOcc$all_species+.001)/(mydatOcc$native_only+.001))
mydatOcc$H_disim<-log((mydatOcc$Beta_all_species)/(mydatOcc$Beta_native_only))

dev.off()
ggplot(mydatOcc,aes(1,H_sim))+stat_summary(aes(color=Ecoregion))+geom_hline(yintercept=0,linetype="dashed")+scale_color_viridis_d()


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
save.image("~/PowellCenter/NOREASTERN.Rda")


colnames(mydat3)[3]<-"Beta_all_species"
colnames(mydat4)[3]<-"Beta_native_only"

colnames(mydat3)[4]<-"all_species"
colnames(mydat4)[4]<-"native_only"


mydatAbn<-left_join(mydat4,mydat3)

mydatAbn$H_sim<-log((mydatAbn$all_species+.001)/(mydatAbn$native_only+.001))
mydatAbn$H_disim<-log((mydatAbn$Beta_all_species)/(mydatAbn$Beta_native_only))

ggpubr::ggarrange(ggplot(mydatOcc,aes(1,H_sim))+stat_summary(aes(color=Ecoregion))+geom_hline(yintercept=0,linetype="dashed")+scale_color_viridis_d(),
ggplot(mydatAbn,aes(1,H_sim))+stat_summary(aes(color=Ecoregion))+geom_hline(yintercept=0,linetype="dashed")+scale_color_viridis_d(),common.legend = TRUE)

mydatOcc$metric<-0
mydatAbn$metric<-1

mydatALL<-rbind(mydatOcc,mydatAbn)

mydatALL$region1<-ifelse(mydatALL$Ecoregion %in% c("North Central Appalachians","Northeastern Highlands","Northern Minnesota Wetlands","Northern Lakes and Forests"),"Northern Forests","Temperate Forests")


coords<-select(d.native,Plot,Lat,Long)
coords<-distinct(coords)
colnames(coords)<-c("site1", "Lat.1", "Lon.1")
coords2<-coords
colnames(coords2)<-c("site2", "Lat.2", "Lon.2")

mydatALL<-left_join(mydatALL,coords)
mydatALL<-left_join(mydatALL,coords2)

cord1<-data.frame(Lon=mydatALL$Lon.1,Lat=mydatALL$Lat.1)
cord2<-data.frame(Lon=mydatALL$Lon.2,Lat=mydatALL$Lat.2)

disy1<-geodist(x=cord1,y=cord2,measure="haversine",paired = TRUE)

mydatALL$Distance<-disy1/100

ggplot()+
  geom_smooth(data=mydatALL,aes(Distance,native_only,color=as.factor(metric)),method="glm",method.args = list(family = "quasibinomial"),size=.8,alpha=0.3)+
geom_smooth(data=mydatALL,aes(Distance,all_species,color=as.factor(metric)),method="glm",method.args = list(family = "quasibinomial"),size=.8,alpha=0.3,linetype="dashed")+facet_wrap(~region1)

  facet_grid(metric~Ecoregion,scales="free_x")+ylab("Taxonomic Similarity")+ggthemes::theme_few(base_size = 9)

  Apal<-filter(mydatALL,Ecoregion=="Northern Lakes and Forests")
  ggplot()+
    geom_smooth(data=Apal,aes(Distance,native_only,color=as.factor(metric)),method="glm",method.args = list(family = "quasibinomial"),size=.8,alpha=0.3,linetype="dashed")+
    geom_smooth(data=Apal,aes(Distance,all_species,color=as.factor(metric)),method="glm",method.args = list(family = "quasibinomial"),size=.8,alpha=0.3)
  

  
  jpeg("~/git/bioticHogs/desertsolitaire/Northeastforests1.jpeg",width = 8, height= 9, units="in",res=300 )
ggplot(mydatALL,aes(metric,H_sim))+stat_summary(aes(color=Ecoregion,shape=region1))+geom_smooth(method="lm",aes(color=Ecoregion),se = FALSE,size=0.1,linetype="dashed")+
  geom_hline(yintercept=0)+scale_color_viridis_d()+scale_fill_viridis_d()+ggthemes::theme_few(base_size = 11)+scale_x_continuous(breaks=c(0,1))+ylab("Homoginization Index")
dev.off()

mod.metric<-brms::brm(H_sim~metric,data=mydatALL)
mod.metric.eco<-brms::brm(H_sim~metric*region1,data=mydatALL)
mod.metric.eco.partypool<-brms::brm(H_sim~metric+(metric|Ecoregion),data=mydatALL)

brms::conditional_effects(mod.metric)
brms::conditional_effects(mod.metric.eco)

newdata<-data.frame(metric=c(0,1,0,1))
examp<-cbind(newdata,fitted(mod.metric,newdata=newdata,probs=(c(.025,.25,.75,.975))))

ggplot(examp,aes(metric,Estimate))+geom_point(size=2,shape=0)+geom_errorbar(aes(ymin=`Q2.5`,ymax=`Q97.5`),width=0)+
  scale_x_continuous(breaks = c(0,1),labels = c("Occ","Abn"))+ggthemes::theme_few()+geom_hline(yintercept=0,linetype="dashed")

jpeg("~/git/bioticHogs/desertsolitaire/Northeastforests1to1.jpeg",width = 8, height= 9, units="in",res=300 )
ggplot(mydatALL,aes(native_only,all_species))+geom_point(aes(color=Ecoregion),size=0.1)+facet_grid(Ecoregion~as.factor(metric))+
  geom_abline(intercept=0,slope=1)+ggthemes::theme_few()+scale_colour_viridis_d()
dev.off()
