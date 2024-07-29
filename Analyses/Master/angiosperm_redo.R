####sensitivity analysis with just angio sperms

rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())
graphics.off()

library(sf)
library(dplyr)
library("mapview")
library("tmap")
library(maps)
library(usmap)
library(phytools)
library(hillR)
library(brms)
library(tidybayes)
library(ggplot2)

setwd("~/Documents/git/bioticHogs/")
load("Analyses/Master/angiosperms.Rda")
d<-read.csv("Data/SPCIS_plant_taxa.csv") ### this is the public data for reproducibility
p<-read.csv("Data/SPCIS_plots.csv") ## plot identities
regionz<-read.csv("Analyses/Master/Input/regionz.csv") ### plots and ecoregion correspondence
regionz<-dplyr::select(regionz,-X)

d<-left_join(d,regionz)
d<-d %>% #select only the most recent survey for resampled plots
  group_by(Plot) %>%
  filter(Year==max(Year)) %>%
  ungroup() 


p<-left_join(p,regionz)
p<-p %>% #select only the most recent survey for resampled plots
  group_by(Plot) %>%
  filter(Year==max(Year)) %>%
  ungroup() 


phy<-read.tree("Analyses/Master/Input/angiosperms.tre") #read in tree


d<-filter(d,SpCode %in% c(phy$tip.label)) ##filter out species from dataset that aren't in phylogeny

##Total # of plots per ecoregion
L4tots<-p%>%group_by(US_L4NAME,NA_L1NAME) %>% count()

## identify number of plots with greater than %5 nonnative cover per ecogregion
L4invs<-d %>% group_by(Plot,NativeStatus,US_L4NAME,NA_L1NAME)%>% summarize(RelCov=sum(PctCov_100,na.rm=TRUE))
L4invs<-filter(L4invs,NativeStatus=="I")
L4invs<-filter(L4invs,RelCov>5)
L4invsmean<- L4invs %>% group_by(US_L4NAME,NA_L1NAME) %>% summarise(meanCov=mean(RelCov,na.rm=TRUE)) ## plot where relcov I is >5
L4invs<- L4invs %>% group_by(US_L4NAME,NA_L1NAME) %>% count() ## plot where relcov I is >5

colnames(L4invs)[3]<-"n_inv"
L4tots<-left_join(L4tots,L4invs)
L4tots$n_inv<-ifelse(is.na(L4tots$n_inv),0,L4tots$n_inv)
L4tots<-filter(L4tots,n>9) ###must have at least 10 total plots

L4tots$regionID<-paste(L4tots$US_L4NAME,L4tots$NA_L1NAME)
L4tots<-filter(L4tots,!is.na(US_L4NAME))
L4tots<-filter(L4tots,!is.na(NA_L1NAME))

d$regionID<-paste(d$US_L4NAME,d$NA_L1NAME)


##function for formatting data
matricize<-function(x){
  temp<-dplyr::select(x,Plot,SpCode,PctCov_100) ### select useful columns
  temp<-dplyr::filter(temp,!is.na(temp$SpCode))
  temp<- tidyr::spread(temp,SpCode,PctCov_100)
  temp[is.na(temp)] <- 0 
  temp<-as.data.frame(temp)
  rownames(temp)<-temp$Plot # convert plot names to row names
  temp<-dplyr::select(temp,-Plot)
  temp<-as.matrix.data.frame(temp)#
}


ecoregionz<-unique(L4tots$regionID) #556 total ecoregiions



###calculate beta diversity ine ach L4 ecoregion
taxdat<-data.frame()
phydat<-data.frame()
for (z in c(1:length(ecoregionz))){
  dat<-dplyr::filter(d,regionID==ecoregionz[z])
  
  
  comms<-matricize(dat)
  
  taxy0<-hill_taxa_parti(comms,q = 0)
  taxy1<-hill_taxa_parti(comms,q = 1)
  taxy2<-hill_taxa_parti(comms,q = 2)
  
  taxy<-rbind(taxy0,taxy1,taxy2)
  taxy$regionID<-ecoregionz[z]
  taxyclass<-"taxonomic"
  
  phy0<-hill_phylo_parti(comms,phy,q=0)
  phy1<-hill_phylo_parti(comms,phy,q=1)
  phy2<-hill_phylo_parti(comms,phy,q=2)
  
  
  phyt<-rbind(phy0,phy1,phy2)
  phyt$regionID<-ecoregionz[z]
  phyt$class<-"phylogenetic"
  
  taxdat<-rbind(taxy,taxdat)
  phydat<-rbind(phyt,phydat)
  
}  

taxdat<-left_join(taxdat,L4tots) 
taxdat$perc_invaded<-taxdat$n_inv/taxdat$n

phydat<-left_join(phydat,L4tots) 
phydat$perc_invaded<-phydat$n_inv/phydat$n


###todo add n as a covariate
taxdat$n_scale<-taxdat$n/max(taxdat$n)
phydat$n_scale<-phydat$n/max(phydat$n)

#pause and take a breath

### next calculate for inv species only and for native species only
d.invy<-filter(d,NativeStatus=="I") ### nonnative
d.native<-filter(d,NativeStatus!="I") ### native


inyonly<-filter(L4tots,n_inv>1)
ecoregionz2<-unique(inyonly$regionID)

taxdat.invy<-data.frame()
phydat.invy<-data.frame()

for (z in c(1:length(ecoregionz2))){
  dat<-dplyr::filter(d.invy,regionID==ecoregionz2[z])
  
  
  comms<-matricize(dat)
  
  taxy0<-hill_taxa_parti(comms,q = 0)
  taxy1<-hill_taxa_parti(comms,q = 1)
  taxy2<-hill_taxa_parti(comms,q = 2)
  
  taxy<-rbind(taxy0,taxy1,taxy2)
  taxy$regionID<-ecoregionz2[z]
  taxyclass<-"taxonomic"
  
  phy0<-hill_phylo_parti(comms,phy,q=0)
  phy1<-hill_phylo_parti(comms,phy,q=1)
  phy2<-hill_phylo_parti(comms,phy,q=2)
  
  
  phyt<-rbind(phy0,phy1,phy2)
  phyt$regionID<-ecoregionz2[z]
  phyt$class<-"phylogenetic"
  
  taxdat.invy<-rbind(taxy,taxdat.invy)
  phydat.invy<-rbind(phyt,phydat.invy)
  
}  

taxdat.invy<-left_join(taxdat.invy,L4tots)
taxdat.invy$perc_invaded<-taxdat.invy$n_inv/taxdat.invy$n

phydat.invy<-left_join(phydat.invy,L4tots)
phydat.invy$perc_invaded<- phydat.invy$n_inv/ phydat.invy$n

taxdat.invy$n_scale<-taxdat.invy$n/max(taxdat.invy$n)
phydat.invy$n_scale<-phydat.invy$n/max(phydat.invy$n)

###native

taxdat.nat<-data.frame()
phydat.nat<-data.frame()
for (z in c(1:length(ecoregionz))){
  dat<-dplyr::filter(d.native,regionID==ecoregionz[z])
  
  
  comms<-matricize(dat)
  
  taxy0<-hill_taxa_parti(comms,q = 0)
  taxy1<-hill_taxa_parti(comms,q = 1)
  taxy2<-hill_taxa_parti(comms,q = 2)
  
  taxy<-rbind(taxy0,taxy1,taxy2)
  taxy$regionID<-ecoregionz[z]
  taxyclass<-"taxonomic"
  
  phy0<-hill_phylo_parti(comms,phy,q=0)
  phy1<-hill_phylo_parti(comms,phy,q=1)
  phy2<-hill_phylo_parti(comms,phy,q=2)
  
  
  phyt<-rbind(phy0,phy1,phy2)
  phyt$regionID<-ecoregionz[z]
  phyt$class<-"phylogenetic"
  
  taxdat.nat<-rbind(taxy,taxdat.nat)
  phydat.nat<-rbind(phyt,phydat.nat)
  
}

taxdat.nat<-left_join(taxdat.nat,L4tots)
taxdat.nat$perc_invaded<-taxdat.nat$n_inv/taxdat.nat$n

phydat.nat<-left_join(phydat.nat,L4tots)
phydat.nat$perc_invaded<- phydat.nat$n_inv/ phydat.nat$n

taxdat.nat$n_scale<-taxdat.nat$n/max(taxdat.nat$n)
phydat.nat$n_scale<-phydat.nat$n/max(phydat.nat$n)

taxdatq2<-filter(taxdat,q==2)
phydatq2<-filter(phydat,q==2)

q2.tax<- brm(
  bf(local_similarity ~n_scale+ perc_invaded+(perc_invaded|NA_L1NAME),
     phi ~n_scale+ perc_invaded+(perc_invaded|NA_L1NAME)),
  data = taxdatq2,
  family = Beta(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr")

q2.phy<- brm(
  bf(local_similarity ~n_scale+ perc_invaded+(perc_invaded|NA_L1NAME),
     phi ~ n_scale+perc_invaded+(perc_invaded|NA_L1NAME)),
  data = phydatq2,
  family = Beta(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr")

round(fixef(q2.tax,probs = c(.25,.75,.05,.95)),2)#*
round(fixef(q2.phy,probs = c(.25,.75,.05,.95)),2)#**

coef(q2.tax,probs = c(.05,.95))#*
coef(q2.phy,probs = c(.05,.95))#*

new.data<-data.frame(NA_L1NAME=rep(unique(taxdat$NA_L1NAME),each=11),perc_invaded=rep(c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1),10))
new.data$n_scale<-median(taxdat$n_scale)

library(tidybayes)


toy_pred.l4.phy2<- q2.phy %>% 
  epred_draws(newdata =new.data,ndraws = 4000)
toy_pred.l4.tax2<- q2.tax %>% 
  epred_draws(newdata =new.data,ndraws = 4000)

toy_pred.l4.tax2$class<-"taxonomic"
toy_pred.l4.phy2$class<-"phylogenetic"

toy_pred.l4.phy2$grouper<-paste(toy_pred.l4.phy2$NA_L1NAME,toy_pred.l4.phy2$class,toy_pred.l4.phy2$.draw)
toy_pred.l4.tax2$grouper<-paste(toy_pred.l4.tax2$NA_L1NAME,toy_pred.l4.tax2$class,toy_pred.l4.tax2$.draw)

mainq2tax<-filter(taxdatq2,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))
mainq2phy<-filter(phydatq2,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))

toy_mainq2tax<-filter(toy_pred.l4.tax2,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))
toy_mainq2phy<-filter(toy_pred.l4.phy2,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))

jpeg("Analyses/Master/Plots/angiocomps.jpeg",width = 9,height = 7, unit="in",res=300)
ggpubr::ggarrange(ggplot()+
  geom_point(data=mainq2tax,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
  
  geom_line(data=toy_mainq2tax,aes(x=perc_invaded*100,y=.epred, group=grouper),size=.01,alpha=.5)+facet_wrap(~NA_L1NAME,nrow=1)+
  geom_smooth(data=toy_mainq2tax,aes(x=perc_invaded*100,y=.epred),size=1,color="gray")+
  coord_cartesian(ylim=c(.05,.9))+ggthemes::theme_few(base_size = 8)+ylab("taxonomic similarity")+xlab("% invaded"),



ggplot()+
  geom_point(data=mainq2phy,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
  
  geom_line(data=toy_mainq2phy,aes(x=perc_invaded*100,y=.epred, group=grouper),size=.01,alpha=.5)+facet_wrap(~NA_L1NAME,ncol=5)+
  geom_smooth(data=toy_mainq2phy,aes(x=perc_invaded*100,y=.epred),size=1,color="gray")+
  coord_cartesian(ylim=c(.04,.9))+ggthemes::theme_few(base_size = 8)+ylab("phylogenetic similarity")+xlab("% invaded"),nrow=2)
dev.off()

save.image("Analyses/Master/angiosperms.Rda")
