#updarted by datn on July 29 2024 to work on GEB revisions.
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

setwd("~/Documents/git/bioticHogs/")# mac

  d1<-read.csv("Data/FULLDatabase_10272022.csv")
  
  use.list<-read.csv("Analyses/AbnOcc/Input/spatial_matched_plots.csv")

  d<-d1 %>%
    group_by(Plot) %>%
    filter(Year==max(Year)) %>%
    ungroup() 
  
  use.list<-dplyr::select(use.list,-X)
  colnames(use.list)[1]<-"pair"
  
  
  
  d<-left_join(d,use.list)
  
  d<-dplyr::filter(d,!is.na(d$AcceptedTaxonName))
  
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
  
  count.use<-use.list%>%group_by(unit,status) %>% count()

  count.use<-filter(count.use,n>1)

  

  
  list.i<-filter(use.list, status==1)### invaded
  list.p<-filter(use.list, status==0) ### native
  
  
  
  natyA<-data.frame()
  natyO<-data.frame()
  inyA<-data.frame()
  inyO<-data.frame()
  
  ecoregionz<-sort(unique(count.use$unit)) ### 211 ecoregions
  ecoregionz<-ecoregionz[104:211]
  ecoregionz<-ecoregionz[22:108]
     for (z in c(1:length(ecoregionz))){
    
    tlands<-dplyr::filter(d,unit==ecoregionz[z])
    #pairzz<-unique(tlands$pair)
    #use.pairz<-pairzz#sample(pairzz,2, replace=FALSE)
    # d.use<-filter(d,pair %in% tlands$pair)
    
    dat.i<-filter(tlands,Plot %in% list.i$Plot)
    dat.p<-filter(tlands,Plot %in% list.p$Plot) 
    
    
    commsall<-matricize(dat.i)
    commsnative<-matricize(dat.p)
    
    datyall.occ<-hill_taxa_parti_pairwise(commsall,q = 0)
    datyall.occ$unit<-ecoregionz[z]
    datyall.occ$L4_KEY<-unique(tlands$L4_KEY)
    #datyall.occ<-dplyr::select(datyall.occ,site1,site2,local_similarity)
    #colnames(datyall.occ)<-c("site1","site2","beta_I_occ")
    
    datyall.abn<-hill_taxa_parti_pairwise(commsall,q = 1)
    datyall.abn$unit<-ecoregionz[z]
    datyall.abn$L4_KEY<-unique(tlands$L4_KEY)
    #datyall.abn<-dplyr::select(datyall.abn,site1,site2,local_similarity)
    #colnames(datyall.abn)<-c("site1","site2","beta_I_abn")
    
    daty.naty.occ<-hill_taxa_parti_pairwise(commsnative,q = 0)
    daty.naty.occ$unit<-ecoregionz[z]
    daty.naty.occ$L4_KEY<-unique(tlands$L4_KEY)
    #daty.naty.occ<-dplyr::select(daty.naty.occ,site1,site2,local_similarity)
    #colnames(daty.naty.occ)<-c("site1","site2","beta_N_occ")
    
    daty.naty.abn<-hill_taxa_parti_pairwise(commsnative,q = 1)
    daty.naty.abn$unit<-ecoregionz[z]
    daty.naty.abn$L4_KEY<-unique(tlands$L4_KEY)
    #daty.naty.abn<-dplyr::select(daty.naty.abn,site1,site2,local_similarity)
    #colnames(daty.naty.abn)<-c("site1","site2","beta_N_abn")
    
    natyA<-rbind(daty.naty.abn,natyA)
    natyO<-rbind(daty.naty.occ,natyO)
    
    inyA<-rbind(datyall.abn,inyA)
    inyO<-rbind(datyall.occ,inyO)
    
    print(paste("done",which(ecoregionz==ecoregionz[z])))   
  }
  
  InyA<-dplyr::select(inyA,site1,site2,local_similarity,unit)
  NatyA<-dplyr::select(natyA,site1,site2,local_similarity,unit)
  colnames(InyA)[3]<-"invA"
  colnames(NatyA)[3]<-"natA"
  Abn.env<-left_join(InyA,NatyA)
  
  InyO<-dplyr::select(inyO,site1,site2,local_similarity,unit)
  NatyO<-dplyr::select(natyO,site1,site2,local_similarity,unit)
  colnames(InyO)[3]<-"invO"
  colnames(NatyO)[3]<-"natO"
  Occ.env<-left_join(InyO,NatyO)
  SPACE<-left_join(Occ.env,Abn.env)
  
  
  SPACE$H.A<-SPACE$invA-SPACE$natA
  SPACE$H.O<-SPACE$invO-SPACE$natO
  
  SPACE<-dplyr::distinct(SPACE)
  write.csv(SPACE,"Analyses/AbnOcc/Input/abn_v_occ_spatial.csv")  
  space<-read.csv("Analyses/AbnOcc/Input/abn_v_occ_spatial.csv")
  
  space<-SPACE
  space<-distinct(space)
  space$quad<-NA
  space$quad[which(space$H.A>0 &space$H.O>0)]<-1
  space$quad[which(space$H.A<=0&space$H.O>0)]<-4  
  space$quad[which(space$H.A<0&space$H.O<0)]<-3 
  space$quad[which(space$H.A>0&space$H.O<=0)]<-2
  space<-filter(space,!is.na(quad))
  
  space$dif<-abs(space$H.O-space$H.A)
  round(mean(space$dif,na.rm=TRUE),3)
  round(sd(space$dif,na.rm=TRUE),3)
  space$number<-ifelse(space$dif>=.5,1,0)
  round(table(space$number)[2]/(table(space$number)[1]+table(space$number)[2]),2)
  
  
  cor.test(space$H.A,space$H.O ,method="pearson")
  
  countit<-space %>%group_by(quad) %>% count()    
  countit$perc<-round(countit$n/ sum(countit$n),3)
  sum(countit$n)
  sum(countit$perc)

  heat.time<-space
  heat.time$H.A<-round(heat.time$H.A,1)
  heat.time$H.O<-round(heat.time$H.O,1)
  
  

  heat.test<-heat.time %>% group_by(H.O,H.A) %>% count()
  heat.test$frequency<-heat.test$n/sum(heat.test$n)
  heat.test$frequency<-round(heat.test$frequency,2)
  
  heat.test$frequency2<-ifelse(heat.test$frequency==0,NA,heat.test$frequency)
  heat.test$frequency3<-ifelse(is.na(heat.test$frequency2),"",paste(heat.test$frequency2*100,"%",sep=""))
  
  jpeg("Analyses/AbnOcc/spatialheatmap.jpeg",width=9,height=8,units = 'in',res=300)
  ggplot()+
    geom_point(data=space,aes(x=H.O,y=H.A),size=0.01)+
    geom_tile(data=heat.test,aes(x=H.O,y=H.A,fill=frequency),alpha=.8)+
    ggthemes::theme_few(base_size = 14)+
    scale_fill_distiller(palette = "BuGn",direction=1 )+
    geom_hline(yintercept=0,color="darkgray")+geom_vline(xintercept=0,color="darkgray")+
    ylab("abundance")+xlab("occurrence")+
    #geom_text(data=heat.test,aes(x=H.O,y=H.A,label=frequency3),size=2,color="navy")+
    #geom_smooth(method="lm",color="yellow")
    geom_abline(slope = 1,intercept = 0,color="black",linetype="dotdash",size=1)
  dev.off()  
  