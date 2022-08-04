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
d.env<-filter(d.env1, US_L3NAME %in% c("High Plains"))
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

abun<-matricize(d)

abun.out<-hill_taxa_parti_pairwise(abun,q = 2)

d.native<-filter(d,NativeStatus!="I")

abun.native<-matricize(d.native)

abun.out.native<-hill_taxa_parti_pairwise(abun.native,q = 2)


abun.out<-select(abun.out,site1,site2,local_similarity)
abun.out.native<-select(abun.out.native,site1,site2,local_similarity)
colnames(abun.out.native)[3]<-"native"
colnames(abun.out)[3]<-"invaders"

outy<-left_join(abun.out.native,abun.out)
outy<-filter(outy,site1!=site2)

ggplot(outy,aes(native,invaders))+
  geom_point()+geom_abline(slope = 1,intercept = 0,color="royalblue")

outy$H<-log((outy$invaders+.001)/(outy$native+.001))




abun.out.occ<-hill_taxa_parti_pairwise(abun,q = 0)
abun.out.native.occ<-hill_taxa_parti_pairwise(abun.native,q = 0)

abun.out.occ<-select(abun.out.occ,site1,site2,local_similarity)
abun.out.native.occ<-select(abun.out.native.occ,site1,site2,local_similarity)
colnames(abun.out.native.occ)[3]<-"native"
colnames(abun.out.occ)[3]<-"invaders"

outy.occ<-left_join(abun.out.native.occ,abun.out.occ)
outy.occ<-filter(outy.occ,site1!=site2)
outy.occ$H<-log((outy.occ$invaders+.001)/(outy.occ$native+.001))



ggplot(outy,aes(native,invaders))+
  geom_point()+geom_abline(slope = 1,intercept = 0,color="royalblue")
ggplot(outy.occ,aes(native,invaders))+
  geom_point()+geom_abline(slope = 1,intercept = 0,color="royalblue")

ggpubr::ggarrange(a,b)
outy<-filter(outy,H!=0)
outy.occ<-filter(outy.occ,H!=0)
ggplot(outy.occ,aes(H,fill=H.bin))+geom_histogram(bins=1000)
ggplot(outy,aes(H,fill=H.bin))+geom_histogram(bins=1000)




outy$H.bin<-ifelse(outy$H<=0,"D","H")
outy$H.bin<-ifelse(outy$H==0,"N",outy$H.bin)
ggplot(outy,aes(H.bin))+geom_bar()

outy.occ$H.bin<-ifelse(outy.occ$H<=0,"D","H")
outy.occ$H.bin<-ifelse(outy.occ$H==0,"N",outy.occ$H.bin)
ggplot(outy.occ,aes(H.bin))+geom_bar()



mod1<-brm(invaders~native,data=outy)
mod2<-brm(invaders~native,data=outy.occ)

pp_check(mod2)
pp_check(mod1)       
