rm(list=ls())
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())
#rstan_options(auto_write = TRUE)

graphics.off()
library(dplyr)
library(tidyr)
library(hillR)
library(ggplot2)
#require(V.PhyloMaker2)
library(ape)

setwd("~/Documents/git/bioticHogs/")
#(3) 

d1<-read.csv("Data/FULLDatabase_10272022.csv")


d.ref<-read.csv("Data/FULLDatabase_10272022_diversity_PctCov_100_10272022_INHABIT.csv")

d.ref.late<-d.ref %>%
  group_by(Plot) %>%
  filter(Year==max(Year)) %>%
  ungroup()


d.ref.early<-d.ref %>%
  group_by(Plot) %>%
  filter(Year==min(Year)) %>%
  ungroup()


d.key<-dplyr::select(d.ref,Plot, NA_L1NAME,NA_L2NAME,NA_L3NAME)
d.key<-distinct(d.key)
colnames(d1)[9]<-"NA_L1NAME"

d1<-left_join(d1,d.key) ### join
  
d.res<-filter(d1,Resampled=="Y")
  length(unique(d.res$Plot)) #1894 plots 

 
  
  d.early<-d.res %>%
    group_by(Plot) %>%
    filter(Year==min(Year)) %>%
    ungroup()

  d.late<-d.res %>%
    group_by(Plot) %>%
    filter(Year==max(Year)) %>%
    ungroup()

  
  
###are there plots that startout uninvaded and end invaded  
startI<-d.early %>% group_by(Plot,NativeStatus) %>% count()   
endI<-d.late %>% group_by(Plot,NativeStatus) %>% count()     

startI2<-filter(startI,NativeStatus=="I")  
endI2<-filter(endI,NativeStatus=="I")  

startpris<-filter(startI,!Plot %in% startI2$Plot)

prisy<-unique(startpris$Plot)
invaded<-filter(endI2,Plot %in%prisy)
####plots that start uninvaded and end with at least 5% invasive conver

i.cov<-filter(d.ref.late,Plot %in% unique(invaded$Plot))
i.cov<-filter(i.cov, RelCov_I>5)
#####
d.start<-filter(d.early, Plot %in% unique(i.cov$Plot))
d.end<-filter(d.late, Plot %in% unique(i.cov$Plot))

######## find plots sampled over 20 year
year1<-select(d.start,Plot,Year)
year1<-distinct(year1)


year2<-select(d.end,Plot,Year)
year2<-distinct(year2)
colnames(year1)[2]<-"start"
colnames(year2)[2]<-"end"
years<-left_join(year1,year2)
years$return<-years$end-years$start
ggplot(years,aes(return))+geom_histogram()


twents<-filter(years,return==20)

d.start.5<-filter(d.start,Plot %in%twents$Plot)
d.end.5<-filter(d.end,Plot %in%twents$Plot)





  matricize.time<-function(x){
    temp<-dplyr::select(x,Plot,AcceptedTaxonName,PctCov_100) ### select useful columns
    temp<-dplyr::filter(temp,!is.na(temp$AcceptedTaxonName))
    temp<- tidyr::spread(temp,AcceptedTaxonName,PctCov_100)
    temp[is.na(temp)] <- 0 
    temp<-as.data.frame(temp)
    rownames(temp)<-temp$Plot# convert plot names to row names
    temp<-dplyr::select(temp,-Plot)
    temp<-as.matrix.data.frame(temp)#
    
  } 

###do beta diversity within all L2 Ecoregions    
cownt<-d.end.5 %>% group_by(Plot,NA_L3NAME) %>% count() 
cownt2<-cownt %>%group_by(NA_L3NAME) %>% count()   
cownt2<-filter(cownt2,n>1)  
ecoregionz<-intersect(cownt2$NA_L3NAME,d.start$NA_L3NAME) #38 plots from 3 L3 regions
  
  lateA<-data.frame()
  lateO<-data.frame()
  earlyA<-data.frame()
  earlyO<-data.frame()
  
  
  for (i in c(1:length(ecoregionz))){
    dat.early<-dplyr::filter(d.start.5,NA_L3NAME==ecoregionz[i])
    dat.late<-dplyr::filter(d.end.5,NA_L3NAME==ecoregionz[i])  
    
    comms.early<-matricize.time(dat.early)
    comms.late<-matricize.time(dat.late)
    
    daty.early.O<-hill_taxa_parti_pairwise(comms.early,q = 0)
    daty.early.O$L3_KEY<-ecoregionz[i]
    
    daty.early.A<-hill_taxa_parti_pairwise(comms.early,q = 1)
    daty.early.A$L3_KEY<-ecoregionz[i]
    
    daty.late.O<-hill_taxa_parti_pairwise(comms.late,q = 0)
    daty.late.O$L3_KEY<-ecoregionz[i]
    
    daty.late.A<-hill_taxa_parti_pairwise(comms.late,q = 1)
    daty.late.A$L3_KEY<-ecoregionz[i]
    
    lateA<-rbind(daty.late.A,lateA)
    earlyA<-rbind(daty.early.A,earlyA)
    
    lateO<-rbind(daty.late.O,lateO)
    earlyO<-rbind(daty.early.O,earlyO)
    print(ecoregionz[i])
  }
  length(unique(lateA$L3_KEY))

  LateA<-dplyr::select(lateA,site1,site2,local_similarity)
  EarlyA<-dplyr::select(earlyA,site1,site2,local_similarity)
  colnames(LateA)[3]<-"lateA"
  colnames(EarlyA)[3]<-"earlyA"
  Abn<-left_join(LateA,EarlyA)
  
  LateO<-dplyr::select(lateO,site1,site2,local_similarity)
  EarlyO<-dplyr::select(earlyO,site1,site2,local_similarity)
  colnames(LateO)[3]<-"lateO"
  colnames(EarlyO)[3]<-"earlyO"
  Occ<-left_join(LateO,EarlyO)
  time<-left_join(Occ,Abn)
  head(time)
  
  
  time$H.A<-time$lateA-time$earlyA
  time$H.O<-time$lateO-time$earlyO
  
  
  cor.test(time$H.A,time$H.O ,method="pearson")
  
  time$quad<-NA
  time$quad[which(time$H.A>0 &time$H.O>0)]<-1
  time$quad[which(time$H.A<=0&time$H.O>0)]<-4  
  time$quad[which(time$H.A<0&time$H.O<0)]<-3 
  time$quad[which(time$H.A>0&time$H.O<=0)]<-2
  time<-filter(time,!is.na(quad))
  
  
  countit<-time %>%group_by(quad) %>% count()    
  countit$perc<-round(countit$n/ sum(countit$n),3)
  sum(countit$n)
  sum(countit$perc)
  
  
  time$dif<-abs(time$H.O-time$H.A)
  round(mean(time$dif,na.rm=TRUE),3)
  round(sd(time$dif,na.rm=TRUE),3)
  time$number<-ifelse(time$dif>=.5,1,0)
  round(table(time$number)[2]/(table(time$number)[1]+table(time$number)[2]),2)
  
  

write.csv(time,"Analyses/AbnOcc/Input/abn_v_occ_temporal.csv")

time<-filter(time,earlyO!=0) 
ggplot()+
  #geom_rect(ymin=0,ymax=1.1,xmin=-1.1,xmax=0, fill="lightgoldenrod1",alpha=.2)+
  #geom_rect(ymin=-1.1,ymax=0,xmin=0,xmax=1.1, fill="lightgoldenrod1",alpha=.2)+
  geom_point(data=time,aes(x=H.O,y=H.A),size=2,alpha=0.3,color="black")+
  geom_hline(yintercept=0,color="gray")+geom_vline(xintercept=0,color="gray")+
  #geom_smooth(data=env,aes(x=H.O,y=H.A),method="lm",color="royalblue3")+
  geom_abline(slope = 1,intercept = 0,color="gray",linetype="dotdash",size=1)+
  #geom_tile(data=heat.test,aes(x=H.O,y=H.A,fill=n),alpha=0.8)+
  #scale_fill_distiller(palette = "RdPu",direction = -1)+  
  ggthemes::theme_few()+
  ylab("abundance")+xlab("occurence")+ylim(-1,1)+xlim(-1,1)  


###
coords<-select(d.end.5,Plot,Lat,Long)
coords<-distinct(coords)
colnames(coords)<-c("site1", "Lat.1", "Lon.1")
coords2<-coords
colnames(coords2)<-c("site2", "Lat.2", "Lon.2")

time<-left_join(time,coords)
time<-left_join(time,coords2)

cord1<-data.frame(Lon=time$Lon.1,Lat=time$Lat.1)
cord2<-data.frame(Lon=time$Lon.2,Lat=time$Lat.2)
library(geodist)
disty1<-geodist(x=cord1,y=cord2,measure="haversine",paired = TRUE)
time$Distance<-disty1*.001
ggplot(time,aes(Distance,dif))+geom_smooth()+geom_point()
time$agree<-ifelse(time$quad %in% c(2,4),1,0)



mod1<-lm(dif~Distance,data=time)
mod<-glm(agree~Distance,data=time,family = binomial(link = 'logit'))
summary(mod1)
time$diffearly<-abs(time$earlyO-time$earlyA)




mod1<-lm(dif~diffearly,data=time)
summary(mod1)
library(brms)
toy.mod.D<-brm(dif~diffearly+(1|mm(site1,site2)),
                     data=time,                                                      
  chains = 4, iter = 2000, warmup = 1000,
  cores = 4, seed = 1234,backend = "cmdstanr")

ggplot(time,aes(diffearly,dif))+geom_point()+geom_smooth(method="lm")


timestarto<-filter(time,diffearly<.1)








invaderdater<-filter(d.end.5, NativeStatus=="I")
max<-invaderdater %>% group_by(Plot) %>% slice_max(PctCov_100)

rich<-invaderdater %>% group_by(Plot) %>% count()



max1<-select(max,Plot,AcceptedTaxonName,PctCov_100)
colnames(max1)[1]<-"site1"
colnames(max1)[2:3]<-paste(colnames(max1)[2:3],"1",sep="_")

rich1<-rich
colnames(rich1)[1]<-"site1"
colnames(rich1)[2]<-paste(colnames(rich1)[2],"1",sep="_")


max2<-select(max,Plot,AcceptedTaxonName,PctCov_100)
colnames(max2)[1]<-"site2"
colnames(max2)[2:3]<-paste(colnames(max2)[2:3],"2",sep="_")

rich2<-rich
colnames(rich2)[1]<-"site2"
colnames(rich2)[2]<-paste(colnames(rich2)[2],"2",sep="_")


timestarto<-left_join(timestarto,max1)
timestarto<-left_join(timestarto,max2)
timestarto<-left_join(timestarto,rich1)
timestarto<-left_join(timestarto,rich2)

timestarto$match<-ifelse(timestarto$AcceptedTaxonName_1==timestarto$AcceptedTaxonName_2,"Same","Different")

timestarto$coverdif<-abs(timestarto$PctCov_100_1-timestarto$PctCov_100_2)
timestarto$richdif<-abs(timestarto$n_1-timestarto$n_2)
#timestartoR<-filter(timestarto, richdif==0)
timestartoR<-filter(timestarto, match=="Same")
#timestartoR<-filter(timestartoR, n_1==1)



ggplot()+
  #geom_rect(ymin=0,ymax=1.1,xmin=-1.1,xmax=0, fill="lightgoldenrod1",alpha=.2)+
  #geom_rect(ymin=-1.1,ymax=0,xmin=0,xmax=1.1, fill="lightgoldenrod1",alpha=.2)+
  geom_point(data=timestarto,aes(x=H.O,y=H.A),size=1,color="black")+
  geom_hline(yintercept=0,color="gray")+geom_vline(xintercept=0,color="gray")+
  #geom_smooth(data=env,aes(x=H.O,y=H.A),method="lm",color="royalblue3")+
  geom_abline(slope = 1,intercept = 0,color="gray",linetype="dotdash",size=1)+
  #geom_tile(data=heat.test,aes(x=H.O,y=H.A,fill=n),alpha=0.8)+
  #scale_fill_distiller(palette = "RdPu",direction = -1)+  
  ggthemes::theme_few()+
  ylab("abundance")+xlab("occurence")+ylim(-1,1)+xlim(-1,1) 


series<-filter(d.res, Plot %in% unique(d.end.5$Plot))
startyear<-dplyr::select(d.ref.early,Plot,Year)
colnames(startyear)[2]<-"start_year"
startyear<-distinct(startyear)
series<-left_join(series,startyear)
series$year<-series$Year-series$start_year






xx0<-filter(series,year==0)
comms.xx0<-matricize.time(xx0)

xx5<-filter(series,year==5)
comms.xx5<-matricize.time(xx5)

xx10<-filter(series,year==10)
comms.xx10<-matricize.time(xx10)

xx15<-filter(series,year==15)
comms.xx15<-matricize.time(xx15)


xx20<-filter(series,year==20)
comms.xx20<-matricize.time(xx20)

xx20.O<-hill_taxa_parti_pairwise(comms.xx20,q = 0)
xx20.O$year<-20

xx15.O<-hill_taxa_parti_pairwise(comms.xx15,q = 0)
xx15.O$year<-15

xx10.O<-hill_taxa_parti_pairwise(comms.xx10,q = 0)
xx10.O$year<-10

xx5.O<-hill_taxa_parti_pairwise(comms.xx5,q = 0)
xx5.O$year<-5

xx0.O<-hill_taxa_parti_pairwise(comms.xx0,q = 0)
xx0.O$year<-0

xxO<-rbind(xx20.O,xx15.O,xx10.O,xx5.O,xx0.O)
ggplot(xxO,aes(year,local_similarity))+geom_point()

xx20.A<-hill_taxa_parti_pairwise(comms.xx20,q = 1)
xx20.A$year<-20

xx15.A<-hill_taxa_parti_pairwise(comms.xx15,q = 1)
xx15.A$year<-15

xx10.A<-hill_taxa_parti_pairwise(comms.xx10,q = 1)
xx10.A$year<-10

xx5.A<-hill_taxa_parti_pairwise(comms.xx5,q = 1)
xx5.A$year<-5

xx0.A<-hill_taxa_parti_pairwise(comms.xx0,q = 1)
xx0.A$year<-0

xxA<-rbind(xx20.A,xx15.A,xx10.O,xx5.A,xx0.A)
ggplot(xxA,aes(year,local_similarity))+geom_point()+geom_line()

xxA$metric<-"Abn"
xxO$metric<-"Occ"
xxT<-rbind(xxA,xxO)

H5.O<-xx5.O$local_similarity-xx0.O$local_similarity
H5.A<-xx5.A$local_similarity-xx0.A$local_similarity

H10.O<-xx10.O$local_similarity-xx0.O$local_similarity
H10.A<-xx10.A$local_similarity-xx0.A$local_similarity

H15.O<-xx15.O$local_similarity-xx0.O$local_similarity
H15.A<-xx15.A$local_similarity-xx0.A$local_similarity

H20.O<-xx20.O$local_similarity-xx0.O$local_similarity
H20.A<-xx20.A$local_similarity-xx0.A$local_similarity

pl<-data.frame(H5A=H5.A,H5O=H5.O,H10A=H10.A,H10O=H10.O,H15A=H15.A,H15O=H15.O,H20A=H20.A,H20O=H20.O)

one<-ggplot()+
  #geom_rect(ymin=0,ymax=1.1,xmin=-1.1,xmax=0, fill="lightgoldenrod1",alpha=.2)+
  #geom_rect(ymin=-1.1,ymax=0,xmin=0,xmax=1.1, fill="lightgoldenrod1",alpha=.2)+
  geom_point(data=pl,aes(x=H5O,y=H5A),size=1,color="black")+
  geom_hline(yintercept=0,color="gray")+geom_vline(xintercept=0,color="gray")+
  #geom_smooth(data=env,aes(x=H.O,y=H.A),method="lm",color="royalblue3")+
  geom_abline(slope = 1,intercept = 0,color="gray",linetype="dotdash",size=1)+
  #geom_tile(data=heat.test,aes(x=H.O,y=H.A,fill=n),alpha=0.8)+
  #scale_fill_distiller(palette = "RdPu",direction = -1)+  
  ggthemes::theme_few()+
  ylab("abundance")+xlab("occurence")+ylim(-1,1)+xlim(-1,1)


two<-ggplot()+
  #geom_rect(ymin=0,ymax=1.1,xmin=-1.1,xmax=0, fill="lightgoldenrod1",alpha=.2)+
  #geom_rect(ymin=-1.1,ymax=0,xmin=0,xmax=1.1, fill="lightgoldenrod1",alpha=.2)+
  geom_point(data=pl,aes(x=H10O,y=H10A),size=1,color="black")+
  geom_hline(yintercept=0,color="gray")+geom_vline(xintercept=0,color="gray")+
  #geom_smooth(data=env,aes(x=H.O,y=H.A),method="lm",color="royalblue3")+
  geom_abline(slope = 1,intercept = 0,color="gray",linetype="dotdash",size=1)+
  #geom_tile(data=heat.test,aes(x=H.O,y=H.A,fill=n),alpha=0.8)+
  #scale_fill_distiller(palette = "RdPu",direction = -1)+  
  ggthemes::theme_few()+
  ylab("abundance")+xlab("occurence")+ylim(-1,1)+xlim(-1,1)

three<-ggplot()+
  #geom_rect(ymin=0,ymax=1.1,xmin=-1.1,xmax=0, fill="lightgoldenrod1",alpha=.2)+
  #geom_rect(ymin=-1.1,ymax=0,xmin=0,xmax=1.1, fill="lightgoldenrod1",alpha=.2)+
  geom_point(data=pl,aes(x=H15O,y=H15A),size=1,color="black")+
  geom_hline(yintercept=0,color="gray")+geom_vline(xintercept=0,color="gray")+
  #geom_smooth(data=env,aes(x=H.O,y=H.A),method="lm",color="royalblue3")+
  geom_abline(slope = 1,intercept = 0,color="gray",linetype="dotdash",size=1)+
  #geom_tile(data=heat.test,aes(x=H.O,y=H.A,fill=n),alpha=0.8)+
  #scale_fill_distiller(palette = "RdPu",direction = -1)+  
  ggthemes::theme_few()+
  ylab("abundance")+xlab("occurence")+ylim(-1,1)+xlim(-1,1)

four<-ggplot()+
  #geom_rect(ymin=0,ymax=1.1,xmin=-1.1,xmax=0, fill="lightgoldenrod1",alpha=.2)+
  #geom_rect(ymin=-1.1,ymax=0,xmin=0,xmax=1.1, fill="lightgoldenrod1",alpha=.2)+
  geom_point(data=pl,aes(x=H20O,y=H20A),size=1,color="black")+
  geom_hline(yintercept=0,color="gray")+geom_vline(xintercept=0,color="gray")+
  #geom_smooth(data=env,aes(x=H.O,y=H.A),method="lm",color="royalblue3")+
  geom_abline(slope = 1,intercept = 0,color="gray",linetype="dotdash",size=1)+
  #geom_tile(data=heat.test,aes(x=H.O,y=H.A,fill=n),alpha=0.8)+
  #scale_fill_distiller(palette = "RdPu",direction = -1)+  
  ggthemes::theme_few()+
  ylab("abundance")+xlab("occurence")+ylim(-1,1)+xlim(-1,1)


ggpubr::ggarrange(one,two,three,four)


cor.test(pl$H5A,pl$H5O ,method="pearson")
cor.test(pl$H10A,pl$H10O ,method="pearson")
cor.test(pl$H15A,pl$H15O ,method="pearson")
cor.test(pl$H20A,pl$H20O ,method="pearson")


quantile(timestarto2$earlyO)
timestarto2A<-filter(timestarto2,earlyO<.2)
timestarto2B<-filter(timestarto2,earlyO>.35)


ggplot()+
  #geom_rect(ymin=0,ymax=1.1,xmin=-1.1,xmax=0, fill="lightgoldenrod1",alpha=.2)+
  #geom_rect(ymin=-1.1,ymax=0,xmin=0,xmax=1.1, fill="lightgoldenrod1",alpha=.2)+
  geom_point(data=timestarto2B,aes(x=H.O,y=H.A),size=1,color="black")+
  geom_hline(yintercept=0,color="gray")+geom_vline(xintercept=0,color="gray")+
  #geom_smooth(data=env,aes(x=H.O,y=H.A),method="lm",color="royalblue3")+
  geom_abline(slope = 1,intercept = 0,color="gray",linetype="dotdash",size=1)+
  #geom_tile(data=heat.test,aes(x=H.O,y=H.A,fill=n),alpha=0.8)+
  #scale_fill_distiller(palette = "RdPu",direction = -1)+  
  ggthemes::theme_few()+
  ylab("abundance")+xlab("occurence")+ylim(-1,1)+xlim(-1,1) 



####lets see if there are any patterns tot he invader
invaderdater<-filter(d.end.5, NativeStatus=="I")
max<-invaderdater %>% group_by(Plot) %>% slice_max(PctCov_100)
table(max$AcceptedTaxonName)


rich<-invaderdater %>% group_by(Plot) %>% count()

max1<-select(max,Plot,AcceptedTaxonName,PctCov_100)
colnames(max1)[1]<-"site1"
colnames(max1)[2:3]<-paste(colnames(max1)[2:3],"1",sep="_")

rich1<-rich
colnames(rich1)[1]<-"site1"
colnames(rich1)[2]<-paste(colnames(rich1)[2],"1",sep="_")


max2<-select(max,Plot,AcceptedTaxonName,PctCov_100)
colnames(max2)[1]<-"site2"
colnames(max2)[2:3]<-paste(colnames(max2)[2:3],"2",sep="_")

rich2<-rich
colnames(rich2)[1]<-"site2"
colnames(rich2)[2]<-paste(colnames(rich2)[2],"2",sep="_")


timestarto<-left_join(time,max1)
timestarto<-left_join(timestarto,max2)

timestarto<-left_join(timestarto,rich1)
timestarto<-left_join(timestarto,rich2)
timestarto<-filter(timestarto,AcceptedTaxonName_1 %in% c("Alliaria petiolata","Lonicera maackii"))
timestarto<-filter(timestarto,AcceptedTaxonName_2 %in% c("Alliaria petiolata","Lonicera maackii"))


timestarto$match<-ifelse(timestarto$AcceptedTaxonName_1==timestarto$AcceptedTaxonName_2,"Same","Different")

timestarto$coverdif<-abs(timestarto$PctCov_100_1-timestarto$PctCov_100_2)
timestarto$richdif<-abs(timestarto$n_1-timestarto$n_2)
#timestartoR<-filter(timestarto, richdif==0)
#timestartoR<-filter(timestarto, match=="Same")
#timestartoR<-filter(timestartoR, n_1==1)
ggplot()+
  #geom_rect(ymin=0,ymax=1.1,xmin=-1.1,xmax=0, fill="lightgoldenrod1",alpha=.2)+
  #geom_rect(ymin=-1.1,ymax=0,xmin=0,xmax=1.1, fill="lightgoldenrod1",alpha=.2)+
  geom_point(data=timestarto,aes(x=H.O,y=H.A,color=match,size=coverdif))+
  geom_hline(yintercept=0,color="gray")+geom_vline(xintercept=0,color="gray")+
  #geom_smooth(data=env,aes(x=H.O,y=H.A),method="lm",color="royalblue3")+
  geom_abline(slope = 1,intercept = 0,color="gray",linetype="dotdash",size=1)+
  #geom_tile(data=heat.test,aes(x=H.O,y=H.A,fill=n),alpha=0.8)+
  #scale_fill_distiller(palette = "RdPu",direction = -1)+  
  ggthemes::theme_few()+facet_wrap(~richdif)
  ylab("abundance")+xlab("occurence")+ylim(-1,1)+xlim(-1,1)  






matricize.YEAR<-function(x){
  temp<-dplyr::select(x,Year,AcceptedTaxonName,PctCov_100) ### select useful columns
  temp<-dplyr::filter(temp,!is.na(temp$AcceptedTaxonName))
  temp<- tidyr::spread(temp,AcceptedTaxonName,PctCov_100)
  temp[is.na(temp)] <- 0 
  temp<-as.data.frame(temp)
  rownames(temp)<-temp$Year# convert plot names to row names
  temp<-dplyr::select(temp,-Year)
  temp<-as.matrix.data.frame(temp)#
  
}
d.res.20<-filter(d.res,Plot %in%unique(d.end.5$Plot))

years<-unique(d.res.20$Plot)
change.dat<-data.frame()
for (i in c(1:length(years))){
  dat.t<-dplyr::filter(d.res.20,Plot==years[i])  
dat.t<-distinct(dat.t)
  comz<-matricize.YEAR(dat.t)
datoO<-hill_taxa_parti_pairwise(comz,q = 0)

datoO$time.dif<-as.numeric(datoO$site2)-as.numeric(datoO$site1)
datoO$Plot<-years[i]

datoA<-hill_taxa_parti_pairwise(comz,q = 1)
datoA$time.dif<-as.numeric(datoA$site2)-as.numeric(datoA$site1)
datoA$Plot<-years[i]
dato<-rbind(datoA,datoO)
change.dat<-rbind(dato,change.dat)
}

length(unique(change.dat$Plot))
ggplot(change.dat,aes(time.dif,local_similarity))+geom_smooth(aes(color=as.factor(q)))
ggplot(change.dat,aes(time.dif))+geom_histogram()
change.dat$metric<-ifelse(change.dat$q==0,"occ","abn")

changeb<-change.dat %>% group_by(Plot,metric) %>% filter(site1 == min(site1))
changeb<-filter(changeb,site1!=site2)

pl<- change.dat %>% group_by(Plot) %>% count()
pl2<-filter(pl,n>19)
goober<-filter(changeb,Plot %in% pl2$Plot)
goober$grouper<-paste(goober$metric,goober$Plot)
goober<-filter(goober,time.dif>4)
ggplot(goober,aes(time.dif,local_similarity))+geom_line(aes(color=metric,group=grouper),size=0.05)+
  geom_smooth(method='lm',aes(color=metric),size=3)

library(brms)
library(lme4)
mod<-lmer(local_similarity~time.dif*metric+1+(time.dif*metric|Plot),data=goober)
summary(mod)

cc<-filter(goober,local_similarity>0)
cc<-filter(cc,local_similarity<1)

modbet<-brm(
  bf(local_similarity ~time.dif*metric+(time.dif*metric|Plot),
     phi ~time.dif*metric+(time.dif*metric|Plot)),
  data = cc,
  family = Beta(),
  chains = 4, iter = 2000, warmup = 1000,
  cores = 4, seed = 1234,backend = "cmdstanr") 

#mod<-brm(local_similarity~time.dif*metric+(time.dif*metric|Plot),data=change.dat.small)

conditional_effects(modbet)

newdat<-data.frame(time.dif=rep(c(5,10,15,20),each=2),metric=rep(c("occ","abn"),4))
library(tidybayes)
bet_pred<- modbet%>% 
  epred_draws(newdata =newdat,ndraws=100,seed = 1234,re_formula = NA)



bet_pred$grouper<-paste(bet_pred$.draw,bet_pred$metric)
ggplot(bet_pred,aes(time.dif,.epred))+stat_pointinterval(aes(shape=metric),.width = .95,size=3)+
  geom_line(aes(color=metric,group=grouper),size=0.01)+scale_shape_manual(values=c(17,1))+
  scale_color_manual(values=c("black","black"))+theme_tidybayes()+
  ylab("similarity to initial survey")+xlab("change in time")

  

closestart<-filter(base)

d.exp<-filter(d.early, Plot %in% example$Plot)
ggplot(d.exp,aes(Plot,PctCov))+geom_bar(aes(fill=AcceptedTaxonName),stat = "identity")+scale_fill_viridis_d()+facet_wrap(~Plot)


    example<-filter(goober,Plot %in% c('069702F_2',	'085502F_1'))
exp.base<-filter(base,site2=='085502F_1')
exp.base<-filter(exp.base,site1=='069702F_2')
exp.base$time<-0

exp.5<-filter(fiver,site2=='085502F_1')
exp.5<-filter(exp.5,site1=='069702F_2')
exp.5$time<-5

exp.10<-filter(tener,site2=='085502F_1')
exp.10<-filter(exp.10,site1=='069702F_2')
exp.10$time<-10

exp.15<-filter(fifteener,site2=='085502F_1')
exp.15<-filter(exp.15,site1=='069702F_2')
exp.15$time<-15

exp.20<-filter(twenter,site2=='085502F_1')
exp.20<-filter(exp.20,site1=='069702F_2')
exp.20$time<-20
exp<-rbind(exp.base,exp.5,exp.10,exp.15,exp.20)

ggpubr::ggarrange(ggplot(exp,aes(time,local_similarity))+geom_point(aes(color=q))+geom_line(aes(color=q,group=q)),
ggplot(example,aes(time.dif,local_similarity))+geom_point(aes(color=metric))+facet_wrap(~Plot)+geom_line(aes(color=metric))+ylim(0,1),ncol=1)

unique(goober$Plot)
goodies<-filter(d.res,Plot %in% unique(goober$Plot))
gooref<-filter(d.early,Plot %in% unique(goober$Plot))
gooref<-select(gooref,Plot, Year)
gooref<-distinct(gooref)###378
colnames(gooref)[2]<-"start"
goodies<-left_join(goodies,gooref)
goodies$year.adj<-goodies$Year-goodies$start
#### do spatial by time
dat.base<-dplyr::filter(goodies,year.adj==0)
dat.5<-dplyr::filter(goodies,year.adj==5)  
dat.10<-dplyr::filter(goodies,year.adj==10)  
dat.15<-dplyr::filter(goodies,year.adj==15)  
dat.20<-dplyr::filter(goodies,year.adj==20)

comms.base<-matricize.time(dat.base)
daty.base.O<-hill_taxa_parti_pairwise(comms.base,q = 0)
daty.base.A<-hill_taxa_parti_pairwise(comms.base,q = 1)
base<-rbind(daty.base.O,daty.base.A)

comms.5<-matricize.time(dat.5)
daty.5.O<-hill_taxa_parti_pairwise(comms.5,q = 0)
daty.5.A<-hill_taxa_parti_pairwise(comms.5,q = 1)
fiver<-rbind(daty.5.O,daty.5.A)

comms.10<-matricize.time(dat.10)
daty.10.O<-hill_taxa_parti_pairwise(comms.10,q = 0)
daty.10.A<-hill_taxa_parti_pairwise(comms.10,q = 1)
tener<-rbind(daty.10.O,daty.10.A)

comms.15<-matricize.time(dat.15)
daty.15.O<-hill_taxa_parti_pairwise(comms.15,q = 0)
daty.15.A<-hill_taxa_parti_pairwise(comms.15,q = 1)
fifteener<-rbind(daty.15.O,daty.15.A)


comms.20<-matricize.time(dat.20)
daty.20.O<-hill_taxa_parti_pairwise(comms.20,q = 0)
daty.20.A<-hill_taxa_parti_pairwise(comms.20,q = 1)
twenter<-rbind(daty.20.O,daty.20.A)

head(twenter)
head(base)

H5<-data.frame(H.O=daty.5.O$local_similarity-daty.base.O$local_similarity,H.A=daty.5.A$local_similarity-daty.base.A$local_similarity)
ggplot(H5,aes(H.O,H.A))+geom_point(size=0.1)+geom_hline(yintercept=0)+
  geom_vline(xintercept=0)


H10<-data.frame(H.O=daty.10.O$local_similarity-daty.base.O$local_similarity,H.A=daty.10.A$local_similarity-daty.base.A$local_similarity)
ggplot(H10,aes(H.O,H.A))+geom_point(size=0.1)+geom_hline(yintercept=0)+
  geom_vline(xintercept=0)

H15<-data.frame(H.O=daty.15.O$local_similarity-daty.base.O$local_similarity,H.A=daty.15.A$local_similarity-daty.base.A$local_similarity)
ggplot(H15,aes(H.O,H.A))+geom_point(size=0.1)+geom_hline(yintercept=0)+
  geom_vline(xintercept=0)



base20.O<-filter(daty.base.O,site1 %in% daty.20.O$site1)
base20.O<-filter(base20.O,site2 %in% daty.20.O$site2)
base20.A<-filter(daty.base.A,site1 %in% daty.20.A$site1)
base20.A<-filter(base20.A,site2 %in% daty.20.A$site2)


H20<-data.frame(H.O=daty.20.O$local_similarity-base20.O$local_similarity,H.A=daty.20.A$local_similarity-base20.A$local_similarity)
ggplot(H20,aes(H.O,H.A))+geom_point(size=0.1)+geom_hline(yintercept=0)+
  geom_vline(xintercept=0)


H5$interval<-5
H10$interval<-10
H15$interval<-15
H20$interval<-20

Hs<-rbind(H5,H10,H15,H20)
summary(lm(Hs$H.O~Hs$interval))
summary(lm(Hs$H.A~Hs$interval))

ggplot(Hs,aes(H.O,H.A))+geom_point(size=0.1,aes(color=as.factor(interval)))+geom_hline(yintercept=0)+
  geom_vline(xintercept=0)+ylim(-1,1)+xlim(-1,1)

Hs %>% group_by(interval) %>% summarize(sdO=sd(H.O),sdA=sd(H.A))

H5$quad<-NA
H5$quad[which(H5$H.A>0 &H5$H.O>0)]<-1
H5$quad[which(H5$H.A<=0&H5$H.O>0)]<-4  
H5$quad[which(H5$H.A<0&H5$H.O<0)]<-3 
H5$quad[which(H5$H.A>0&H5$H.O<=0)]<-2
H5<-filter(H5,!is.na(quad))

H10$quad<-NA
H10$quad[which(H10$H.A>0 &H10$H.O>0)]<-1
H10$quad[which(H10$H.A<=0&H10$H.O>0)]<-4  
H10$quad[which(H10$H.A<0&H10$H.O<0)]<-3 
H10$quad[which(H10$H.A>0&H10$H.O<=0)]<-2
H10<-filter(H10,!is.na(quad))


H15$quad<-NA
H15$quad[which(H15$H.A>0 &H15$H.O>0)]<-1
H15$quad[which(H15$H.A<=0&H15$H.O>0)]<-4  
H15$quad[which(H15$H.A<0&H15$H.O<0)]<-3 
H15$quad[which(H15$H.A>0&H15$H.O<=0)]<-2
H15<-filter(H15,!is.na(quad))

H20$quad<-NA
H20$quad[which(H20$H.A>0 &H20$H.O>0)]<-1
H20$quad[which(H20$H.A<=0&H20$H.O>0)]<-4  
H20$quad[which(H20$H.A<0&H20$H.O<0)]<-3 
H20$quad[which(H20$H.A>0&H20$H.O<=0)]<-2
H20<-filter(H20,!is.na(quad))

H5$dif<-abs(H5$H.O-H5$H.A)
round(mean(H5$dif,na.rm=TRUE),3)

H10$dif<-abs(H10$H.O-H10$H.A)
round(mean(H10$dif,na.rm=TRUE),3)

H15$dif<-abs(H15$H.O-H15$H.A)
round(mean(H15$dif,na.rm=TRUE),3)

H20$dif<-abs(H20$H.O-H20$H.A)
round(mean(H20$dif,na.rm=TRUE),3)


cor.test(H5$H.A,H5$H.O ,method="pearson")
cor.test(H10$H.A,H10$H.O ,method="pearson")
cor.test(H15$H.A,H15$H.O ,method="pearson")
cor.test(H20$H.A,H20$H.O ,method="pearson")

countit<-H5 %>%group_by(quad) %>% count()    
countit$perc<-round(countit$n/ sum(countit$n),3)
sum(countit$n)
sum(countit$perc)


countit10<-H10 %>%group_by(quad) %>% count()    
countit10$perc<-round(countit10$n/ sum(countit10$n),3)
sum(countit$n)
sum(countit$perc)

countit15<-H15 %>%group_by(quad) %>% count()    
countit15$perc<-round(countit15$n/ sum(countit15$n),3)


countit20<-H20 %>%group_by(quad) %>% count()    
countit20$perc<-round(countit20$n/ sum(countit20$n),3)



###Idea do the whole think with these 4 sampled AIM datasets
##look at spatial beta diveristy in year 5, and 10, 15, and 20
#show temporal
### and illustrate with a 2 plot case study