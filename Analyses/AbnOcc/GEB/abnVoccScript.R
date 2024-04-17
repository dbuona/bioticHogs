###main analysis for abn vs. occ biotic homgenization analyses####

rm(list=ls())
options(stringsAsFactors = FALSE)
graphics.off()

##load librarys
library(dplyr)
library(tidyr)
library(hillR)
library(ggplot2)


setwd("~/Documents/git/bioticHogs/")# set working directory

d1<-read.csv("Data/SPCIS_plant_taxa.csv") ## public dataset from

use.list<-read.csv("Analyses/AbnOcc/GEB/matchedplots.csv") ## environmental paired plots from makematches.R
use.list<-dplyr::select(use.list,Plot,L4_KEY,status,subclass)## drop unneeded columns

colnames(use.list)[4]<-"pair" ## reassign column name to represent nonative-native plot pairings

d<-d1 %>% group_by(Plot) %>% filter(Year==max(Year)) %>% ungroup() ###include only more recent surveys from repeat sampled plots

d<-left_join(d,use.list)## join pair assignments to main data

d<-dplyr::filter(d,!is.na(d$AcceptedTaxonName)) ## remove species with no names

###generate function  for coverting data frame to wide matrix to use in beta diversity calculation
matricize<-function(x){
  temp<-dplyr::select(x,pair,AcceptedTaxonName,PctCov_100) ### select useful columns
  temp<-dplyr::filter(temp,!is.na(temp$AcceptedTaxonName))
  temp<- tidyr::spread(temp,AcceptedTaxonName,PctCov_100)
  temp[is.na(temp)] <- 0 
  temp<-as.data.frame(temp)
  rownames(temp)<-temp$pair# convert plot names to row names
  temp<-dplyr::select(temp,-pair)
  temp<-as.matrix.data.frame(temp)#
  
}


count.use<-use.list%>%group_by(L4_KEY,status) %>% count() ## count number of invaded and univaded plots in each ecoregion
count.use<-filter(count.use,n>1) ##include only regions with at least 2 native and non-native plots

##seperate native and non-native plots from pairing data
list.i<-filter(use.list, status==1)### invaded
list.p<-filter(use.list, status==0) ### native

##make 4 blank data frame to input beta diversity calcualtions into
natyA<-data.frame() # native plots with abudance based calculation
natyO<-data.frame() # native plots with occurence based calculation
inyA<-data.frame() # non-native plots with abudance based calculation
inyO<-data.frame() # non-native plots with occurence based calculation

ecoregionz<-sort(unique(count.use$L4_KEY)) ### vector of ecoregions to calculate beta diversity within

##calculate beta diversity within each ecoregion. this takes a long time, so alternatively, load output in line 94
if(FALSE){ #skip the time consuming loop and load put
for (z in c(1:length(ecoregionz))){
  tlands<-dplyr::filter(d,L4_KEY==ecoregionz[z])

  dat.i<-filter(tlands,Plot %in% list.i$Plot) ## non-native plots within ecoregion[z]
  dat.p<-filter(tlands,Plot %in% list.p$Plot) ## native plots within ecoregion[z]
  
  commsall<-matricize(dat.i) # reformat data
  commsnative<-matricize(dat.p) # reformat data
  
  datyall.occ<-hill_taxa_parti_pairwise(commsall,q = 0) # calculate occurence-based pairwise beta diversity
  datyall.occ$L4_KEY<-ecoregionz[z] #append ecoregion label

  datyall.abn<-hill_taxa_parti_pairwise(commsall,q = 1) # calculate abundance-based pairwise beta diversity
  datyall.abn$L4_KEY<-ecoregionz[z] #append ecoregion label
 
  daty.naty.occ<-hill_taxa_parti_pairwise(commsnative,q = 0) # calculate occurence-based pairwise beta diversity
  daty.naty.occ$L4_KEY<-ecoregionz[z] #append ecoregion label

  daty.naty.abn<-hill_taxa_parti_pairwise(commsnative,q = 1) # calculate abundance-based pairwise beta diversity
  daty.naty.abn$L4_KEY<-ecoregionz[z] #append ecoregion label
 
  natyA<-rbind(daty.naty.abn,natyA) #bind each ecoregion in the loop to main data
  natyO<-rbind(daty.naty.occ,natyO) #bind each ecoregion in the loop to main data
  
  inyA<-rbind(datyall.abn,inyA) #bind each ecoregion in the loop to main data
  inyO<-rbind(datyall.occ,inyO) #bind each ecoregion in the loop to main data
  
  print(paste("done",which(ecoregionz==ecoregionz[z]))) #finish loop   
}

## the column "local_similarity" represents the sorensen (occ) and horn similarity (abn) index
## this is what we use in our analyses

InyA<-dplyr::select(inyA,site1,site2,local_similarity,L4_KEY) ##select useful columns
NatyA<-dplyr::select(natyA,site1,site2,local_similarity,L4_KEY) ##select useful columns
colnames(InyA)[3]<-"invA" #rename local similariy
colnames(NatyA)[3]<-"natA" #rename local similariy
Abn.env<-left_join(InyA,NatyA) # join abundance native and non-native comparisions

InyO<-dplyr::select(inyO,site1,site2,local_similarity,L4_KEY) ##select useful columns
NatyO<-dplyr::select(natyO,site1,site2,local_similarity,L4_KEY) ##select useful columns
colnames(InyO)[3]<-"invO" #rename local similariy
colnames(NatyO)[3]<-"natO" #rename local similariy
Occ.env<-left_join(InyO,NatyO) # join occurence native and non-native comparisions
ENV<-left_join(Occ.env,Abn.env) # join abundance and occurence based calculations

ENV$H.A<-ENV$invA-ENV$natA ## calculate H index for abundance
ENV$H.O<-ENV$invO-ENV$natO## calculate H index for occurence

write.csv(ENV,"Analyses/AbnOcc/GEB/ABNvOCC.csv",row.names = FALSE)
}

env<-read.csv("Analyses/AbnOcc/GEB/ABNvOCC.csv")
source("Analyses/AbnOcc/GEB/addons.R") ## add whether or not plots are from same dataset and have same dominant invader

##How many paris fall into each quantdrent (Figure 3)
env$quad<-NA
env$quad[which(env$H.A>0 &env$H.O>0)]<-1
env$quad[which(env$H.A<=0&env$H.O>0)]<-4  
env$quad[which(env$H.A<0&env$H.O<0)]<-3 
env$quad[which(env$H.A>0&env$H.O<=0)]<-2
env<-filter(env,!is.na(quad))

## descriptive statistics
env$dif<-abs(env$H.O-env$H.A)
round(mean(env$dif,na.rm=TRUE),3)
round(sd(env$dif,na.rm=TRUE),3)
env$number<-ifelse(env$dif>=.5,1,0)
round(table(env$number)[2]/(table(env$number)[1]+table(env$number)[2]),2)

cor.test(env$H.A,env$H.O ,method="pearson") ## corelatin coefficient

countit<-env %>%group_by(quad) %>% count()    ### count pairs in each quadrent 
countit$perc<-round(countit$n/ sum(countit$n),3)### percentage of point in each quadrent

#make heatmap
heat.time<-env
heat.time$H.A<-round(heat.time$H.A,1)
heat.time$H.O<-round(heat.time$H.O,1)

heat.test<-heat.time %>% group_by(H.O,H.A) %>% count() #count how many are in each .1 x.1 grid cell
heat.test$frequency<-heat.test$n/sum(heat.test$n) ## what percentage of totdal data are in each cell?
heat.test$frequency<-round(heat.test$frequency,2)

#jpeg("Analyses/AbnOcc/heatmap.jpeg",width=8,height=8,units = 'in',res=300) Figure 3
ggplot()+
  geom_point(data=env,aes(x=H.O,y=H.A),size=0.01)+
  geom_tile(data=heat.test,aes(x=H.O,y=H.A,fill=frequency),alpha=.8)+
  tidybayes::theme_tidybayes()+
  scale_fill_distiller(palette = "PuRd",direction=1 )+
  geom_hline(yintercept=0,color="darkgray")+geom_vline(xintercept=0,color="darkgray")+
  ylab("abundance")+xlab("occurence")+
  #geom_text(data=heat.test,aes(x=H.O,y=H.A,label=frequency3),size=2,color="navy")+
  #geom_smooth(method="lm",color="yellow")
  geom_abline(slope = 1,intercept = 0,color="black",linetype="dotdash",size=1)
dev.off()

#do pairs include the same dominant invader, and are they from the same original dataset
env$same.dom<-ifelse(env$dom.invader.1==env$dom.invader.2,1,0) ## are the dominant invaders the same?
env$same.dat.i<-ifelse(env$dataset.i.plot.1==env$dataset.i.plot.2,1,0) #are the non-native plots from the same data sets?
env$same.dat.u<-ifelse(env$dataset.u.plot.1==env$dataset.u.plot.2,1,0)#are the native plots from the same data sets?
env$same.all<-ifelse(env$dataset.u.plot.1==env$dataset.u.plot.2& env$dataset.u.plot.1==env$dataset.i.plot.1&env$dataset.i.plot.1==env$dataset.i.plot.2,1,0) # are all 4 plots that make up the comparision from the same dataset


table(env$same.all)
803541/969488
##make plot 4
ex1<-ggplot(env,aes(as.factor(same.all),dif))+stat_summary(size=0.001,alpha=0.2,aes(group=L4_KEY),position="jitter",shape=1)+
  stat_summary(
    fun.data= "mean_sdl", fun.args = list(mult = 1),size=.75)+coord_cartesian(ylim=c(0,.5))+ggthemes::theme_few(base_size = 9)+ylab("absolute difference abn vs. occ")+
  scale_x_discrete(name="",labels=c("multiple datasets","same dataset"))

ex2<-ggplot(env,aes(as.factor(same.dom),dif))+stat_summary(size=0.001,alpha=0.2,aes(group=L4_KEY),position="jitter",shape=1)+
  stat_summary(
    fun.data= "mean_sdl", fun.args = list(mult = 1),size=.75)+coord_cartesian(ylim=c(0,.5))+ggthemes::theme_few(base_size = 9)+ylab("")+
  scale_x_discrete(name="",labels=c("different invaders","same invader"))+theme(axis.text.y = element_blank(),axis.ticks.y = element_blank())

ex4<-ggplot(env,aes(as.factor(quad)))+geom_bar(position="fill",aes(fill=as.factor(same.dom)),width = .4)+
  scale_fill_manual(name="",values=c("grey40","lightgray"),labels=c("different invaders","same invader"))+ggthemes::theme_few(base_size = 9)+
  scale_x_discrete(name="",labels=c("H abn \nH occ","H abn \nD occ","D abn \nD occ","D abn \nH occ"))+theme(legend.position = "top")+
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank())+ylab("")

ex3<-ggplot(env,aes(as.factor(quad)))+geom_bar(position="fill",aes(fill=as.factor(same.all)),width = .4,color="black")+
  scale_fill_manual(name="",values=c("black","white"),labels=c("multiple datasets","same dataset"))+ggthemes::theme_few(base_size = 9)+
  scale_x_discrete(name="",labels=c("H abn \nH occ","H abn \nD occ","D abn \nD occ","D abn \nH occ"))+theme(legend.position = "top")+ylab("proportion")

#jpeg("Analyses/AbnOcc/explaination1.jpeg",width=7,height=7,units = 'in',res=200)  Figure 4
ggpubr::ggarrange(ex1,ex2,ex3,ex4,widths=c(.35,.3),heights=c(.3,.35),labels=c("a)","b)","c)","d)"))
dev.off()  

