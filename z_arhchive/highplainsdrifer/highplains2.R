#highplains drifers

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

setwd("~/PowellCenter/")# for PC
set.seed(3)


###read in data
d1<-read.csv("Data/FULLDatabase_05272022.csv")
d.env1<-read.csv("Data/FULLDatabase_05272022_plotsenv4Jul2022.csv")
d.ref<-read.csv("Data/FullDatabase_diversity_July4.csv")

unique(d.env1$US_L3NAME)

d.env<-filter(d.env1, US_L3NAME %in% c("High Plains","Driftless Area"))
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

mydat<-data.frame(site1=character(),site2=character(),local_similarity=numeric(),Ecoregion=character())
ecoregionz<-c(unique(d.env$US_L3NAME))


for (z in c(1:length(ecoregionz))){
  
  tlands<-dplyr::filter(d.env,US_L3NAME==ecoregionz[z]) ##
  dat<-filter(d, Plot %in% c(unique(tlands$Plot)))
  comms<-matricize(dat)
  
  daty<-hill_taxa_parti_pairwise(comms,q = 2)
  daty<-dplyr::select(daty,site1,site2,local_similarity)
  daty$Ecoregion<-ecoregionz[z]
  mydat<-rbind(daty,mydat)
  print(paste("done",which(ecoregionz==ecoregionz[z]))) 
}



d.native<-filter(d,NativeStatus!="I")

mydat2<-data.frame(site1=character(),site2=character(),local_similarity=numeric(),Ecoregion=character())
for (z in c(1:length(ecoregionz))){
  
  tlands<-dplyr::filter(d.env,US_L3NAME==ecoregionz[z]) ##
  dat<-filter(d.native, Plot %in% c(unique(tlands$Plot)))
  comms<-matricize(dat)
  
  daty<-hill_taxa_parti_pairwise(comms,q = 2)
  daty<-dplyr::select(daty,site1,site2,local_similarity)
  daty$Ecoregion<-ecoregionz[z]
  mydat2<-rbind(daty,mydat2)
  print(paste("done",which(ecoregionz==ecoregionz[z]))) 
}

colnames(mydat)[3]<-"all_species"
colnames(mydat2)[3]<-"native_only"

mydat3<-left_join(mydat2,mydat)

ggplot(mydat3,aes(native_only,all_species))+geom_point(size=.1)+geom_abline(intercept=0,slope=1,color="royalblue")+facet_wrap(~Ecoregion)


mydat3$H<-log((mydat3$all_species+.001)/(mydat3$native_only+.001))


ggplot(mydat3,aes(H))+geom_histogram(bins=60)+facet_wrap(~Ecoregion,scales="free")+geom_vline(xintercept=0)
library(brms)
#mod1<-brm(H~(1|Ecoregion),data=mydat3)


coords<-select(d.native,Plot,Lat,Long)
coords<-distinct(coords)
colnames(coords)<-c("site1", "Lat.1", "Lon.1")
coords2<-coords
colnames(coords2)<-c("site2", "Lat.2", "Lon.2")

mydat3<-left_join(mydat3,coords)
mydat3<-left_join(mydat3,coords2)

cord1<-data.frame(Lon=mydat3$Lon.1,Lat=mydat3$Lat.1)
cord2<-data.frame(Lon=mydat3$Lon.2,Lat=mydat3$Lat.2)

disty1<-geodist(x=cord1,y=cord2,measure="haversine",paired = TRUE)
mydat3$Distance<-disty1/100

jpeg("~/git/bioticHogs/highplainsdrifer/hilRdistance.jpeg")
ggplot()+
  geom_smooth(data=mydat3,aes(Distance,native_only),method="glm",method.args = list(family = "quasibinomial"),size=1,color="navyblue",fill="navyblue",alpha=0.2)+
  geom_smooth(data=mydat3,aes(Distance,all_species),method="glm",method.args = list(family = "quasibinomial"),size=1,color="firebrick",fill="firebrick",alpha=0.2)+
  facet_wrap(~Ecoregion,scales="free_x")+ylab("Taxonomic Similarity")+ggthemes::theme_few()
dev.off()


mydato<-data.frame(site1=character(),site2=character(),local_similarity=numeric(),Ecoregion=character())
ecoregionz<-c(unique(d.env$US_L3NAME))


for (z in c(1:length(ecoregionz))){
  
  tlands<-dplyr::filter(d.env,US_L3NAME==ecoregionz[z]) ##
  dat<-filter(d, Plot %in% c(unique(tlands$Plot)))
  comms<-matricize(dat)
  
  daty<-hill_taxa_parti_pairwise(comms,q = 0)
  daty<-dplyr::select(daty,site1,site2,local_similarity)
  daty$Ecoregion<-ecoregionz[z]
  mydato<-rbind(daty,mydato)
  print(paste("done",which(ecoregionz==ecoregionz[z]))) 
}

mydat2o<-data.frame(site1=character(),site2=character(),local_similarity=numeric(),Ecoregion=character())
for (z in c(1:length(ecoregionz))){
  
  tlands<-dplyr::filter(d.env,US_L3NAME==ecoregionz[z]) ##
  dat<-filter(d.native, Plot %in% c(unique(tlands$Plot)))
  comms<-matricize(dat)
  
  daty<-hill_taxa_parti_pairwise(comms,q = 0)
  daty<-dplyr::select(daty,site1,site2,local_similarity)
  daty$Ecoregion<-ecoregionz[z]
  mydat2o<-rbind(daty,mydat2o)
  print(paste("done",which(ecoregionz==ecoregionz[z]))) 
}

colnames(mydato)[3]<-"all_species"
colnames(mydat2o)[3]<-"native_only"

mydat3o<-left_join(mydat2o,mydato)

ggplot(mydat3o,aes(native_only,all_species))+geom_point(size=.1)+geom_abline(intercept=0,slope=1,color="royalblue")+facet_wrap(~Ecoregion)


mydat3o$H<-log((mydat3o$all_species+.001)/(mydat3o$native_only+.001))

mydat3o$Distance<-mydat3$Distance

ggplot()+
  geom_smooth(data=mydat3o,aes(Distance,native_only),method="glm",method.args = list(family = "quasibinomial"),size=1,color="navyblue",fill="navyblue",alpha=0.2)+
  geom_smooth(data=mydat3o,aes(Distance,all_species),method="glm",method.args = list(family = "quasibinomial"),size=1,color="firebrick",fill="firebrick",alpha=0.2)+
  facet_wrap(~Ecoregion,scales="free_x")+ylab("Taxonomic Similarity")+ggthemes::theme_few()+ggtitle("Occurance")

ggplot()+
  geom_smooth(data=mydat3,aes(Distance,native_only),method="glm",method.args = list(family = "quasibinomial"),size=1,color="navyblue",fill="navyblue",alpha=0.2)+
  geom_smooth(data=mydat3,aes(Distance,all_species),method="glm",method.args = list(family = "quasibinomial"),size=1,color="firebrick",fill="firebrick",alpha=0.2)+
  facet_wrap(~Ecoregion,scales="free_x")+ylab("Taxonomic Similarity")+ggthemes::theme_few()+ggtitle("Abundance")


mydat3$metric<-"abundance"
mydat3o$metric<-"occurance"

mydat3<-select(mydat3,-Lat.1,-Lon.1,-Lat.2,-Lon.2)
mydat4<-(rbind(mydat3,mydat3o))

mydat4a<-select(mydat4,-native_only)
mydat4b<-select(mydat4,-all_species)

colnames(mydat4a)[4]<-"local_similarity"
mydat4a$status<-"all species"

colnames(mydat4b)[3]<-"local_similarity"
mydat4b$status<-"native only"

mydat4b<-mydat4b[, c(1,2,4,3,5,6,7,8)]

mydat5<-rbind(mydat4a,mydat4b)


ggplot()+
  geom_smooth(data=mydat5,aes(Distance,local_similarity,linetype=status),method="glm",method.args = list(family = "quasibinomial"),size=.8, color="black",alpha=0.3)+
 
  facet_grid(metric~Ecoregion,scales="free_x")+ylab("Taxonomic Similarity")+ggthemes::theme_few(base_size = 9)
 


colnames(mydato)
save.image("~/PowellCenter/Data/highplainsdrift2.Rda")
#####Helen and )
checkers<-mydato %>% group_by(site1,Ecoregion)%>% summarise(mean_sim=mean(all_species))

d.ref2<-select(d.ref2,Plot,TotalPctCover_I,RelCov_I)
colnames(d.ref2)[1]<-"site1"
checkers<-left_join(checkers,d.ref2)

ggplot(checkers,aes(RelCov_I,mean_sim))+geom_point()+geom_smooth(method="lm")+facet_wrap(~Ecoregion,scales="free_x")+
  ggthemes::theme_few()


ggplot(checkers,aes(TotalPctCover_I,mean_sim))+geom_point()+geom_smooth(method="lm")+facet_wrap(~Ecoregion,scales="free_x")+
  ggthemes::theme_few()
colnames(checkers)

mydat5<-left_join(mydat5,d.ref2)
colnames(mydat5)[c(9,10)]<-c("pctcov1","relcov1")
d.ref3<-d.ref2
colnames(d.ref3)<-c("site2", "pctcov2","relcov2")

mydat5<-left_join(mydat5,d.ref3)

mydat5$I_index<-(mydat5$relcov1+mydat5$relcov2)/100

ggplot(mydat5,aes(I_index,local_similarity))+geom_point()+geom_smooth()+facet_wrap(~Ecoregion)

mydat5$I_bin<-NA
mydat5$I_bin<-ifelse(mydat5$I_index<1.5,"high","medium")
mydat5$I_bin<-ifelse(mydat5$I_index<1,"low",mydat5$I_bin)
mydat5$I_bin<-ifelse(mydat5$I_index<.5,"very low",mydat5$I_bin)

ggplot()+
  geom_smooth(data=mydat5,aes(Distance,local_similarity,linetype=status,color=I_bin),method="glm",method.args = list(family = "quasibinomial"),size=.8,alpha=0.3)+
  
  facet_grid(metric~Ecoregion,scales="free_x")+ylab("Taxonomic Similarity")+ggthemes::theme_few(base_size = 9)

###plots#####################################
#1 map

library(maps)
jpeg("~/git/bioticHogs/highplainsdrifer/map.jpeg",width=8, height=8, unit= "in", res=200)
map("state", interior = FALSE)
points(cord1$Lon, cord1$Lat, pch=.5, cex=1,col="purple")
dev.off()

#####site 1
un<-ggplot(checkers,aes(RelCov_I,mean_sim))+geom_point()+geom_smooth(method="lm")+facet_wrap(~Ecoregion,scales="free_x")+
  ggthemes::theme_few()

DU<-ggplot(checkers,aes(TotalPctCover_I,mean_sim))+geom_point()+geom_smooth(method="lm")+facet_wrap(~Ecoregion,scales="free_x")+
  ggthemes::theme_few()
colnames(checkers)

jpeg("~/git/bioticHogs/highplainsdrifer/mean_COVERS.jpeg",width=8, height=8, unit= "in", res=200)
ggpubr::ggarrange(un,DU,nrow=2)
dev.off()
#2 One to ones
onetoone_o<-ggplot(mydat3o,aes(native_only,all_species))+geom_point(size=.1)+geom_abline(intercept=0,slope=1,color="royalblue",size=1.5)+facet_wrap(~Ecoregion)+
  ggthemes::theme_few()
onetoone_a<-ggplot(mydat3,aes(native_only,all_species))+geom_point(size=.1)+geom_abline(intercept=0,slope=1,color="royalblue",size=1.5)+facet_wrap(~Ecoregion)+
  ggthemes::theme_few()

jpeg("~/git/bioticHogs/highplainsdrifer/onetoones.jpeg",width=8, height=8, unit= "in", res=200)
ggpubr::ggarrange(onetoone_o,onetoone_a,nrow=2,labels = c("Occ","Abn"))

dev.off()
### distance decay

jpeg("~/git/bioticHogs/highplainsdrifer/distance_decay.jpeg",width=8, height=8, unit= "in", res=200)
ggplot()+
  geom_smooth(data=mydat5,aes(Distance,local_similarity,linetype=status),method="glm",method.args = list(family = "quasibinomial"),size=.8, color="black",alpha=0.3)+
  
  facet_grid(metric~Ecoregion,scales="free_x")+ylab("Taxonomic Similarity")+ggthemes::theme_few(base_size = 9)
dev.off()


ggplot(mydat4,aes(metric,H))+geom_boxplot()+facet_wrap(~Ecoregion,scales="free")

library(brms)
mod1<-brm(H~metric*Ecoregion,warmup=3000,iter=4000,data=mydat4,control=list(adapt_delta = .9))
summary(mod1)
conditional_effects(mod1)

cony<-conditional_effects(mod1)

jpeg("~/git/bioticHogs/highplainsdrifer/H_mod.jpeg",width=8, height=8, unit= "in", res=200)
plot(cony, points = F,plot = F)[[3]]+ggthemes::theme_few()+
  scale_color_viridis_d(option = "B",begin = .7, end=.05)+ylab("Homoginization Index")+
  geom_hline(yintercept=0,linetype="dotted")
dev.off()
