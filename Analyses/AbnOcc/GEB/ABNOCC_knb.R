#### anaylsis for knb

rm(list=ls())
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())
#rstan_options(auto_write = TRUE)
graphics.off()

library(dplyr)
library(tidyr)
library(hillR)
library(ggplot2)
setwd("~/Documents/git/bioticHogs/")# mac


d1<-read.csv("Data/SPCIS_plant_taxa.csv") ## public dataset from

use.list<-read.csv("Analyses/AbnOcc/GEB/env_matchedplots.csv") ## environmental paired plots from makematches.R
use.list<-dplyr::select(use.list,Plot,L4_KEY,status,subclass)## drop unneeded columns

colnames(use.list)[4]<-"pair" ## reassign column name to represent nonative-native plot pairings

d<-d1 %>% group_by(Plot) %>% filter(Year==max(Year)) %>% ungroup() ###include only more recent surveys from repeat sampled plots

d<-left_join(d,use.list)## join pair assignments to main data

d<-dplyr::filter(d,!is.na(d$AcceptedTaxonName)) ## remove species with no names



env<-read.csv("Analyses/AbnOcc/Input/abn_v_occ_environmental.csv")

invaded<-filter(use.list,status==1) ###subset to invaded plots
uninvaded<-filter(use.list,status!=1) ###subset to uninvaded plots

invaded<-dplyr::select(invaded,Plot,pair) ## select relevent columns
uninvaded<-dplyr::select(uninvaded,Plot,pair)

uninvaded1<-uninvaded ## dupicated data
uninvaded2<-uninvaded  ## dupicated data

colnames(uninvaded1)<-c("uninvaded.plot.1","site1")
colnames(uninvaded2)<-c("uninvaded.plot.2","site2") 

invaded1<-invaded
invaded2<-invaded  
colnames(invaded1)<-c("invaded.plot.1","site1")
colnames(invaded2)<-c("invaded.plot.2","site2") 

env<-left_join(env,uninvaded1)
env<-left_join(env,uninvaded2)

env<-left_join(env,invaded1)
env<-left_join(env,invaded2)

##now get dominate invader



d.inv<-filter(d,NativeStatus=="I")

d.inv<-dplyr::select(d.inv,Plot,AcceptedTaxonName,PctCov_100)  

dom<-d.inv %>% group_by(Plot) %>% slice(which.max(PctCov_100))

dom1<-dom
dom2<-dom

colnames(dom1)<-c("invaded.plot.1","dom.invader.1","cov.dom.1")
colnames(dom2)<-c("invaded.plot.2","dom.invader.2","cov.dom.2")

env<-left_join(env,dom1)
env<-left_join(env,dom2)




###bromus
env$same.dom<-ifelse(env$dom.invader.1==env$dom.invader.2,"same dominant invader","different dominant invaders")
d.refy<-read.csv("Data/FULLDatabase_10272022_diversity_PctCov_100_10272022_INHABIT.csv")
d.ref<-d.refy %>%
  group_by(Plot) %>%
  filter(Year==max(Year)) %>%
  ungroup() 
d.key<-dplyr::select(d.ref,L4_KEY,L1_KEY)
d.key<-distinct(d.key)
env<-left_join(env,d.key)


####output env data
colnames(env)
env<-dplyr::select(env,-X)
colnames(env)
df<-env[,c(1,2,10,11,12,13,14,15,16,17,18,19,4,3,5,6,7,8,9)]
write.csv(df,"Analyses/AbnOcc/GEB/KNB/ABNvOCC_data.csv",row.names = FALSE)





br1<-dplyr::filter(env,dom.invader.1=="Bromus tectorum")
br2<-dplyr::filter(env,dom.invader.2=="Bromus tectorum")
br<-rbind(br1,br2)
br<-distinct(br)


brD<-dplyr::filter(br,L1_KEY=="10  NORTH AMERICAN DESERTS")
#bralt<-dplyr::filter(br,L2_KEY=="9.3  WEST-CENTRAL SEMI-ARID PRAIRIES")
onsps<-dplyr::filter(d.ref,Richness_I==1)
twosps<-dplyr::filter(d.ref,Richness_I>=2)
nrow(onsps)/(nrow(onsps)+nrow(twosps))

brsp1<-filter(brD,invaded.plot.1 %in% onsps$Plot)
brsp1<-filter(brD,invaded.plot.2 %in% onsps$Plot)


br$n_invaders<-ifelse(br$site1 %in% brsp1$site1 & br$site2 %in%brsp1$site2,"one","multiple")

table(br$n_invaders)





env$quad<-NA
env$quad[which(env$H.A>0 &env$H.O>0)]<-1
env$quad[which(env$H.A<=0&env$H.O>0)]<-4  
env$quad[which(env$H.A<0&env$H.O<0)]<-3 
env$quad[which(env$H.A>0&env$H.O<=0)]<-2
env<-filter(env,!is.na(quad))






env$dif<-abs(env$H.O-env$H.A)
round(mean(env$dif,na.rm=TRUE),3)
round(sd(env$dif,na.rm=TRUE),3)
env$number<-ifelse(env$dif>=.5,1,0)
round(table(env$number)[2]/(table(env$number)[1]+table(env$number)[2]),2)


cor.test(env$H.A,env$H.O ,method="pearson")

countit<-env %>%group_by(quad) %>% count()    
countit$perc<-round(countit$n/ sum(countit$n),3)


heat.time<-env
heat.time$H.A<-round(heat.time$H.A,1)
heat.time$H.O<-round(heat.time$H.O,1)


heat.time$catcher<-abs(heat.time$H.A)+abs(heat.time$H.O)


heat.test<-heat.time %>% group_by(H.O,H.A) %>% count()
heat.test$frequency<-heat.test$n/sum(heat.test$n)
heat.test$frequency<-round(heat.test$frequency,2)

heat.test$frequency2<-ifelse(heat.test$frequency==0,NA,heat.test$frequency)
heat.test$frequency3<-ifelse(is.na(heat.test$frequency2),"",paste(heat.test$frequency2*100,"%",sep=""))

#FIGURE 2
ggplot()+
  geom_point(data=env,aes(x=H.O,y=H.A),size=0.01)+
  geom_tile(data=heat.test,aes(x=H.O,y=H.A,fill=frequency),alpha=.8)+
  ggthemes::theme_few(base_size = 14)+
  scale_fill_distiller(palette = "Blues",direction=1 )+
  geom_hline(yintercept=0,color="darkgray")+geom_vline(xintercept=0,color="darkgray")+
  ylab("abundance")+xlab("occurence")+
  #geom_text(data=heat.test,aes(x=H.O,y=H.A,label=frequency3),size=2,color="navy")+
  #geom_smooth(method="lm",color="yellow")
  geom_abline(slope = 1,intercept = 0,color="black",linetype="dotdash",size=1)
dev.off()

d.refy<-read.csv("Data/FULLDatabase_10272022_diversity_PctCov_100_10272022_INHABIT.csv")
d.ref<-d.refy %>%
  group_by(Plot) %>%
  filter(Year==max(Year)) %>%
  ungroup() 

br1<-dplyr::filter(df,dom.invader.1=="Bromus tectorum") #select plots where B. tectorum is dominant in plot 1
br2<-dplyr::filter(df,dom.invader.2=="Bromus tectorum") #select plots where B. tectorum is dominant in plot 2
br<-rbind(br1,br2) # this make a data frame where either B. tectorum is dominant in one or both of the invaded plots
br<-distinct(br)### remove any dumplicates plots

brD<-dplyr::filter(br,L1_KEY=="10  NORTH AMERICAN DESERTS")  

onsps<-dplyr::filter(d.ref,Richness_I==1)
twosps<-dplyr::filter(d.ref,Richness_I>=2)
nrow(onsps)/(nrow(onsps)+nrow(twosps))

brsp1<-filter(brD,invaded.plot.1 %in% onsps$Plot)
brsp1<-filter(brD,invaded.plot.2 %in% onsps$Plot)


br$n_invaders<-ifelse(br$site1 %in% brsp1$site1 & br$site2 %in%brsp1$site2,"one","multiple")


coords<-dplyr::select(d.ref,Plot,Lat,Long)
coodsp1<-coords
coordsp2<-coords

colnames(coodsp1)<-c("invaded.plot.1","invaded.Lat.1","invaded.Long.1")

colnames(coordsp2)<-c("invaded.plot.2","invaded.Lat.2","invaded.Long.2")
coodsp1<-filter(coodsp1,invaded.plot.1 %in% unique(df$invaded.plot.1))
coodsp1<-distinct(coodsp1)

coordsp2<-filter(coordsp2,invaded.plot.2 %in% unique(df$invaded.plot.2))
coordsp2<-distinct(coordsp2)

envgoo<-dplyr::left_join(brsp1,coodsp1)
envgoo<-dplyr::left_join(envgoo,coordsp2)

coorz1<-dplyr::select(envgoo,invaded.Lat.1,invaded.Long.1)
coorz2<-dplyr::select(envgoo,invaded.Lat.2,invaded.Long.2)
library(geodist)
colnames(coorz1)<-c("Lat","Lon")
colnames(coorz2)<-c("Lat","Lon")
goingthedistance<-geodist(x=coorz1,coorz2,method,measure="haversine",paired = TRUE)
envgoo$plot_distances<-goingthedistance
envgoo$plot_distances<-envgoo$plot_distances/1000
envgoo$logdist<-log(envgoo$plot_distances)
envgoo<-filter(envgoo,plot_distances>0)
table(is.na(envgoo$logdist))

envgoo$domers<-NA
envgoo$domers[which(envgoo$dom.invader.1=="Bromus tectorum" & envgoo$dom.invader.2=="Bromus tectorum")]<- "mono-specific"
envgoo$domers[which(is.na(envgoo$domers))]<- "hetero-specific"

envgoo$cover_diff<-NA
envgoo$cover_diff[which(envgoo$cov.dom.1>=15 & envgoo$cov.dom.2>=15)]<- "high"
envgoo$cover_diff[which(envgoo$cov.dom.1<15 & envgoo$cov.dom.2<15)]<- "low"
envgoo$cover_diff[which(is.na(envgoo$cover_diff))]<- "mixed"

write.csv(envgoo,"Analyses/AbnOcc/GEB/KNB/BROMUS_data.csv")
