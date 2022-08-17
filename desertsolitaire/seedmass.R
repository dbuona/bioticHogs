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

setwd("~/PowellCenter/")# for pc
set.seed(3)

###read in data
#d1<-read.csv("Data/FULLDatabase_05272022.csv")
d.env1<-read.csv("Data/FULLDatabase_05272022_plotsenv4Jul2022.csv")
#d.ref<-read.csv("Data/FullDatabase_diversity_July4.csv")

d<-read.csv("Data/spcisFullEco_withTRY_24May2022 (1).csv")
colnames(d)
unique(d.env1$NA_L2NAME)
d.env<-filter(d.env1, NA_L3NAME %in% c("High Plains"))
table(d.env$US_L3NAME)

d.env<-filter(d.env, NA_L3NAME %in% c("High Plains"))
d<-filter(d, Plot %in% unique(d.env$Plot))    

d<-filter(d,!is.na(SLA))




d<- d%>% group_by(Plot) %>% 
  filter(Year==max(Year)) %>%
  ungroup() ### this take resapled plots and chooses the most recent
d<-filter(d,!is.na(AcceptedTaxonName))

dtrait<-select(d,AcceptedTaxonName,SLA)
dtrait<-distinct(dtrait)
library(tibble)

dtrait<-arrange(dtrait,AcceptedTaxonName)

dtrait<-column_to_rownames(dtrait, var = "AcceptedTaxonName")


numbers<-d %>%group_by(Plot)%>% count() 
numbers<-filter(numbers,n>1)

d<-filter(d,Plot %in% c(numbers$Plot))
matricize<-function(x){
  temp<-dplyr::select(x,Plot,AcceptedTaxonName,PctCov) ### select useful columns
  temp<-dplyr::filter(temp,!is.na(temp$AcceptedTaxonName))
  temp<- tidyr::spread(temp,AcceptedTaxonName,PctCov)
  temp[is.na(temp)] <- 0 
  temp<-as.data.frame(temp)
  rownames(temp)<-temp$Plot # convert plot names to row names
  temp<-dplyr::select(temp,-Plot)
  temp<-as.matrix.data.frame(temp)# remove plot column
}
dat<-matricize(d)


M1<-hill_func_parti_pairwise(dat,dtrait,q=0)
M2<-hill_taxa_parti_pairwise(dat,q=0)


M1<-select(M1,site1,site2,local_similarity)
M2<-select(M2,site1,site2,local_similarity)
colnames(M1)[3]<-"Funct"
colnames(M2)[3]<-"Taxon"

M3<-left_join(M1,M2)
#M4<-left_join(M1,M2)


#ggplot(M4,aes(Taxon,Funct))+geom_point(size=0.1)+geom_abline(intercept=0,slope=1,color="royalblue",size=1.5)+
 # ggthemes::theme_few()+ylim(0,1)




coords<-select(d,Plot,Original.Lat,Original.Long)
coords<-distinct(coords)
colnames(coords)<-c("site1", "Lat.1", "Lon.1")
coords2<-coords
colnames(coords2)<-c("site2", "Lat.2", "Lon.2")

M3<-left_join(M3,coords)
M3<-left_join(M3,coords2)

cord1<-data.frame(Lon=M3$Lon.1,Lat=M3$Lat.1)
cord2<-data.frame(Lon=M3$Lon.2,Lat=M3$Lat.2)

disty1<-geodist(x=cord1,y=cord2,measure="haversine",paired = TRUE)
M3$Distance<-disty1/100
M3$Funct<-ifelse(M3$Funct<0,0,M3$Funct)

M3a<-select(M3,-Funct)
M3b<-select(M3,-Taxon)
colnames(M3a)[3]<-"Similarity"
colnames(M3b)[3]<-"Similarity"

M3a$class<-"taxomonic"
M3b$class<-"trait"

M5<-rbind(M3a,M3b)

ggplot()+
  geom_smooth(data=M5,aes(Distance,Similarity,linetype=class),method="glm",method.args = list(family = "quasibinomial"),size=1,alpha=0.2,color="black",fill="black")+
 
  ggthemes::theme_few()+scale_color_manual(values=c("navyblue","firebrick"))+scale_fill_manual(values=c("navyblue","firebrick"))+ theme(legend.position="top")

ggplot(M3,aes(Taxon,Funct))+geom_point(size=0.1)+geom_abline(intercept=0,slope=1,color="royalblue",size=1.5)+
  ggthemes::theme_few()

#jpeg("~/git/bioticHogs/desertsolitaire/sonoran_seedmass.jpeg")
#ggpubr::ggarrange(a,b)
#dev.off()
#save.image("~/git/bioticHogs/highplainsdrifer/SLA.Rda")
d.native<-d.native<-filter(d,NativeStatus!="I")

numbers.n<-d.native %>%group_by(Plot)%>% count() 
numbers.n<-filter(numbers.n,n>1)

d.native<-filter(d.native,Plot %in% c(numbers.n$Plot))


dat.native<-matricize(d.native)


M1.native<-hill_func_parti_pairwise(dat.native,dtrait,q=0)
M2.native<-hill_taxa_parti_pairwise(dat.native,q=0)

M1.native<-select(M1.native,site1,site2,local_similarity)
M2.native<-select(M2.native,site1,site2,local_similarity)
colnames(M1.native)[3]<-"Funct"
colnames(M2.native)[3]<-"Taxon"

M3.native<-left_join(M1.native,M2.native)

M3native.a<-select(M3.native,-Funct)
M3native.b<-select(M3.native,-Taxon)
colnames(M3native.a)[3]<-"Similarity"
colnames(M3native.b)[3]<-"Similarity"

M3native.a$class<-"taxomonic"
M3native.b$class<-"trait"

M5.native<-rbind(M3native.a,M3native.b)

colnames(M5.native)[3]<-"Sim_native_only"

M6<-left_join(M5.native,M5,by=c("site1","site2","class"))
ggplot()+
stat_smooth(data=M6,aes(Distance,Similarity,linetype=class),method="glm",method.args = list(family = "quasibinomial"),size=1,alpha=0.2,color="black",fill="black")

range(M6$Sim_native_only)
M6$Sim_native_only<-ifelse(M6$Sim_native_only<0,0,M6$Sim_native_only)

M6native.a<-select(M6,-Similarity)
M6native.b<-select(M6,-Sim_native_only)

colnames(M6native.a)[3]<-"Similarity"


M6native.a$comm<-"native only"
M6native.b$comm<-"all species"

M7<-rbind(M6native.a,M6native.b)  

jpeg("~/git/bioticHogs/highplainsdrifer/sla_occ_highplains.jpeg")
ggplot()+stat_smooth(data=M7,aes(Distance,Similarity,linetype=class,color=comm),method="glm",method.args = list(family = "quasibinomial"),size=1,alpha=0.2)+
  ggthemes::theme_few()+scale_color_manual(values=c("black","skyblue"))
dev.off()  
