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
dplot<-read.csv("Data/FULLDatabase_10272022.csv")
dplot<-dplyr::select(dplot,Dataset,Plot)
dplot<-distinct(dplot)
d<-read.csv("Data/SPCIS_plant_taxa.csv") ### this is the public data for reproducibility
d<-left_join(d,dplot)
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

if(FALSE){# plot gymnosperms ~ invasoive
d$gymnosperm<-ifelse(d$SpCode %in% c(phy$tip.label),0,1)  

L4gym<-d %>% group_by(Plot,gymnosperm,US_L4NAME,NA_L1NAME)%>% summarize(RelCov=sum(PctCov_100,na.rm=TRUE))
L4gym<-filter(L4gym,gymnosperm==1)
#L4gym<-filter(L4gym,NA_L1NAME %in% c("EASTERN TEMPERATE FORESTS","NORTHERN FORESTS","NORTHWESTERN FORESTED MOUNTAINS" ))
#ggplot(L4gym, aes(RelCov))+geom_density()+facet_wrap(~NA_L1NAME)+xlim(0,100)


#L4gym<-filter(L4gym,RelCov>5)
L4gymsmean<- L4gym %>% group_by(US_L4NAME,NA_L1NAME) %>% summarise(meanCov=mean(RelCov,na.rm=TRUE)) ## plot where relcov I is >5
L4gym<- L4gym %>% group_by(US_L4NAME,NA_L1NAME) %>% count() 


  


##Total # of plots per ecoregion
L4tots<-p%>%group_by(US_L4NAME,NA_L1NAME) %>% count()

## identify number of plots with greater than %5 nonnative cover per ecogregion
L4invs<-d %>% group_by(Plot,NativeStatus,US_L4NAME,NA_L1NAME)%>% summarize(RelCov=sum(PctCov_100,na.rm=TRUE))
L4invs<-filter(L4invs,NativeStatus=="I")
L4invs<-filter(L4invs,RelCov>5)
L4invsmean<- L4invs %>% group_by(US_L4NAME,NA_L1NAME) %>% summarise(meanCov=mean(RelCov,na.rm=TRUE)) ## plot where relcov I is >5
L4invs<- L4invs %>% group_by(US_L4NAME,NA_L1NAME) %>% count() ## plot where relcov I is >5

colnames(L4invs)[3]<-"n_inv"
colnames(L4gym)[3]<-"n_gym"
L4tots<-left_join(L4tots,L4invs)
L4tots<-left_join(L4tots,L4gym)

L4tots$n_inv<-ifelse(is.na(L4tots$n_inv),0,L4tots$n_inv)
L4tots$n_gym<-ifelse(is.na(L4tots$n_gym),0,L4tots$n_gym)
L4tots$perc_inv<-L4tots$n_inv/L4tots$n
L4tots$perc_gym<-L4tots$n_gym/L4tots$n

unique(L4tots$NA_L1NAME)
L4tots<-left_join(L4tots,L4gymsmean)
L4tots$meanCov<-ifelse(is.na(L4tots$meanCov),0,L4tots$meanCov)
L4tots<-filter(L4tots,n>9) ###must have at least 10 total plots
L4tots$gymno<-L4tots$perc_gym*L4tots$meanCov
gymplot<-filter(L4tots,NA_L1NAME %in% c("EASTERN TEMPERATE FORESTS","NORTHWESTERN FORESTED MOUNTAINS", "NORTH AMERICAN DESERTS" ,"GREAT PLAINS", "NORTHERN FORESTS"))
  gymplot2<-filter(L4tots,NA_L1NAME %in% c("NORTHWESTERN FORESTED MOUNTAINS" ))

  cor(gymplot2$gymno,gymplot2$perc_inv)
  cor(gymplot$gymno,gymplot$perc_inv)
  
gympy2<-ggplot(gymplot,aes(gymno,perc_inv))+geom_point(size=0.1,aes(color=NA_L1NAME))+
  geom_smooth(method="glm",method.args = list(family = "quasibinomial"),aes(color=NA_L1NAME,fill=NA_L1NAME),level=.5)+
  ggthemes::theme_few(base_size = 10)+
theme(legend.position = "bottom")+scale_color_viridis_d(name="")+scale_fill_viridis_d(name="")+
  guides(fill=guide_legend(ncol=2),color=guide_legend(ncol=2))+ylab("% invaded")+xlab("% non-angiosperm x relative abundance")

nonnat<-filter(d,NativeStatus=="I")
nonnat<-select(nonnat,SpCode,gymnosperm)
nonnat<-distinct(nonnat)
gympy1<-ggplot(nonnat,aes(as.factor(gymnosperm)))+geom_bar()+
  scale_x_discrete(name="taxonomic position",labels=c("angiosperm","non-angiosperm"))+
  ggthemes::theme_few(base_size = 10)+ylab("number of non-native species")

jpeg("Analyses/Master/Plots/gym4sup.jpeg",width=9,height=5,unit='in',res=200)
ggpubr::ggarrange(gympy1,gympy2,nrow=1,widths=c(1,2),labels = c("a)","b)"))
dev.off()


L4tots$regionID<-paste(L4tots$US_L4NAME,L4tots$NA_L1NAME)
L4tots<-filter(L4tots,!is.na(US_L4NAME))
L4tots<-filter(L4tots,!is.na(NA_L1NAME))

}

d<-filter(d,SpCode %in% c(phy$tip.label)) ##filter out species from dataset that aren't in phylogeny

L4tots<-p%>%group_by(US_L4NAME,NA_L1NAME,Dataset) %>% count()

## identify number of plots with greater than %5 nonnative cover per ecogregion
L4invs<-d %>% group_by(Plot,NativeStatus,US_L4NAME,NA_L1NAME,Dataset)%>% summarize(RelCov=sum(PctCov_100,na.rm=TRUE))
L4invs<-filter(L4invs,NativeStatus=="I")
L4invs<-filter(L4invs,RelCov>5)
L4invsmean<- L4invs %>% group_by(US_L4NAME,NA_L1NAME,Dataset) %>% summarise(meanCov=mean(RelCov,na.rm=TRUE)) ## plot where relcov I is >5
L4invs<- L4invs %>% group_by(US_L4NAME,NA_L1NAME,Dataset) %>% count() ## plot where relcov I is >5

colnames(L4invs)[4]<-"n_inv"
L4tots<-left_join(L4tots,L4invs)
L4tots$n_inv<-ifelse(is.na(L4tots$n_inv),0,L4tots$n_inv)
L4tots<-filter(L4tots,n>2) ###must have at least 2 total plots

L4tots$regionID<-paste(L4tots$US_L4NAME,L4tots$NA_L1NAME,L4tots$Dataset)
L4tots<-filter(L4tots,!is.na(US_L4NAME))
L4tots<-filter(L4tots,!is.na(NA_L1NAME))

d$regionID<-paste(d$US_L4NAME,d$NA_L1NAME,d$Dataset)


L4tots$regionID<-paste(L4tots$US_L4NAME,L4tots$NA_L1NAME,L4tots$Dataset)


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


ecoregionz<-unique(L4tots$regionID)[-945]
#"Shadscale-Dominated Saline Basins NORTH AMERICAN DESERTS NWCA"
length(unique(L4tots$US_L4NAME))#683 total ecoregiions



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

if(FALSE){
ecoregionz2<-unique(inyonly$regionID)[-4]


taxdat.invy<-data.frame()
phydat.invy<-data.frame()

ecoregionz2<-ecoregionz2[234:806]
ecoregionz2<-ecoregionz2[203:573]
ecoregionz2<-ecoregionz2[114:371]


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
}
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


save.image("Analyses/Master/angiosperms.Rda")

taxdatq2<-filter(taxdat,q==2)
phydatq2<-filter(phydat,q==2)

nattaxdatq2<-filter(taxdat.nat,q==2)
natphydatq2<-filter(phydat.nat,q==2)

q2.tax<- brm(
  bf(local_similarity ~n_scale+ perc_invaded+(perc_invaded|NA_L1NAME)+(1|Dataset),
     phi ~n_scale+ perc_invaded+(perc_invaded|NA_L1NAME)+(1|Dataset)),
  data = taxdatq2,
  family = Beta(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr")

q2.phy<- brm(
  bf(local_similarity ~n_scale+ perc_invaded+(perc_invaded|NA_L1NAME)+(1|Dataset),
     phi ~ n_scale+perc_invaded+(perc_invaded|NA_L1NAME)+(1|Dataset)),
  data = phydatq2,
  family = Beta(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr")


nattaxdatq2<-dplyr::filter(nattaxdatq2,local_similarity>0)

q2.tax.nat<- brm(
  bf(local_similarity ~n_scale+ perc_invaded+(perc_invaded|NA_L1NAME)+(1|Dataset),
     phi ~n_scale+ perc_invaded+(perc_invaded|NA_L1NAME)+(1|Dataset)),
  data = nattaxdatq2,
  family = Beta(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr")


q2.phy.nat<- brm(
  bf(local_similarity ~n_scale+ perc_invaded+(perc_invaded|NA_L1NAME)+(1|Dataset),
     phi ~ n_scale+perc_invaded+(perc_invaded|NA_L1NAME)+(1|Dataset)),
  data = natphydatq2,
  family = Beta(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr")




###########
new.data<-data.frame(NA_L1NAME=rep(unique(taxdat$NA_L1NAME),each=11),perc_invaded=rep(c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1),10))
new.data$n_scale<-median(taxdat$n_scale)

new.data2<-data.frame(perc_invaded=c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1))
new.data2$n_scale<-median(taxdat$n_scale)

toy_pred.l4.phy2<- q2.phy %>% 
  epred_draws(newdata =new.data,ndraws = 4000,re_formula = ~(perc_invaded|NA_L1NAME))
toy_pred.l4.tax2<- q2.tax %>% 
  epred_draws(newdata =new.data,ndraws = 4000,re_formula = ~(perc_invaded|NA_L1NAME))

toy_pred.l4.phy.noreg2<- q2.phy %>% 
  epred_draws(newdata =new.data2,ndraws = 4000,re_formula = NA)
toy_pred.l4.tax.noreg2<- q2.tax %>% 
  epred_draws(newdata =new.data2,ndraws = 4000,re_formula = NA)

toy_pred.l4.tax2$class<-"taxonomic"
toy_pred.l4.phy2$class<-"phylogenetic"

toy_pred.l4.tax.noreg2$class<-"taxonomic"
toy_pred.l4.phy.noreg2$class<-"phylogenetic"

toy_pred.l4.phy2$grouper<-paste(toy_pred.l4.phy2$NA_L1NAME,toy_pred.l4.phy2$class,toy_pred.l4.phy2$.draw)
toy_pred.l4.tax2$grouper<-paste(toy_pred.l4.tax2$NA_L1NAME,toy_pred.l4.tax2$class,toy_pred.l4.tax2$.draw)


toy_pred.l4.phy.noreg2$grouper<-paste(toy_pred.l4.phy.noreg2$class,toy_pred.l4.phy.noreg2$.draw)
toy_pred.l4.tax.noreg2$grouper<-paste(toy_pred.l4.tax.noreg2$class,toy_pred.l4.tax.noreg2$.draw)


nattoy_pred.l4.phy2<- q2.phy.nat %>% 
  epred_draws(newdata =new.data,ndraws = 4000,re_formula = ~(perc_invaded|NA_L1NAME))
nattoy_pred.l4.tax2<- q2.tax.nat %>% 
  epred_draws(newdata =new.data,ndraws = 4000,re_formula = ~(perc_invaded|NA_L1NAME))


nattoy_pred.l4.phy.noreg2<- q2.phy.nat %>% 
  epred_draws(newdata =new.data2,ndraws = 4000,re_formula = NA)
nattoy_pred.l4.tax.noreg2<- q2.tax.nat %>% 
  epred_draws(newdata =new.data2,ndraws = 4000,re_formula = NA)


nattoy_pred.l4.tax2$class<-"taxonomic"
nattoy_pred.l4.phy2$class<-"phylogenetic"
nattoy_pred.l4.tax.noreg2$class<-"taxonomic"
nattoy_pred.l4.phy.noreg2$class<-"phylogenetic"

nattoy_pred.l4.phy2$grouper<-paste(nattoy_pred.l4.phy2$NA_L1NAME,nattoy_pred.l4.phy2$class,nattoy_pred.l4.phy2$.draw)
nattoy_pred.l4.tax2$grouper<-paste(nattoy_pred.l4.tax2$NA_L1NAME,nattoy_pred.l4.tax2$class,nattoy_pred.l4.tax2$.draw)

nattoy_pred.l4.phy.noreg2$grouper<-paste(nattoy_pred.l4.phy.noreg2$class,nattoy_pred.l4.phy.noreg2$.draw)
nattoy_pred.l4.tax.noreg2$grouper<-paste(nattoy_pred.l4.tax.noreg2$class,nattoy_pred.l4.tax.noreg2$.draw)


###regions

mainq2tax<-filter(taxdatq2,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))
mainq2phy<-filter(phydatq2,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))

toy_mainq2tax<-filter(toy_pred.l4.tax2,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))
toy_mainq2phy<-filter(toy_pred.l4.phy2,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))





ccc<-ggpubr::ggarrange(ggplot()+
                         geom_point(data=taxdatq2,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
                         
                         geom_line(data=toy_pred.l4.tax.noreg2,aes(x=perc_invaded*100,y=.epred, group=grouper),size=.01,alpha=.5)+
                         geom_smooth(data=toy_pred.l4.tax.noreg2,aes(x=perc_invaded*100,y=.epred),size=1,color="darkgray")+
                         coord_cartesian(ylim=c(0,.5))+ggthemes::theme_few(base_size = 10)+ylab("taxonomic similarity")+xlab("% invaded"),
                       
                       
                       ggplot()+
                         geom_point(data=phydatq2,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
                         
                         geom_line(data=toy_pred.l4.phy.noreg2,aes(x=perc_invaded*100,y=.epred, group=grouper),size=.01,alpha=.5)+
                         geom_smooth(data=toy_pred.l4.phy.noreg2,aes(x=perc_invaded*100,y=.epred),size=1,color="darkgray")+
                         coord_cartesian(ylim=c(.5,1))+ggthemes::theme_few(base_size = 10)+ylab("phylogenetic similarity")+xlab("% invaded"),ncol=2,labels=c("a)","b)"))


aa2<-ggpubr::ggarrange(ggplot()+
                         geom_point(data=mainq2tax,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
                         
                         geom_line(data=toy_mainq2tax,aes(x=perc_invaded*100,y=.epred, group=grouper),size=.01,alpha=.5)+facet_wrap(~NA_L1NAME,ncol=5)+
                         geom_smooth(data=toy_mainq2tax,aes(x=perc_invaded*100,y=.epred),size=1,color="gray")+
                         coord_cartesian(ylim=c(0,.6))+ggthemes::theme_few(base_size = 9)+ylab("taxonomic similarity")+xlab("% invaded"),
                       
                       
                       ggplot()+
                         geom_point(data=mainq2phy,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
                         
                         geom_line(data=toy_mainq2phy,aes(x=perc_invaded*100,y=.epred, group=grouper),size=.01,alpha=.5)+facet_wrap(~NA_L1NAME,ncol=5)+
                         geom_smooth(data=toy_mainq2phy,aes(x=perc_invaded*100,y=.epred),size=1,color="gray")+
                         coord_cartesian(ylim=c(.5,1))+ggthemes::theme_few(base_size = 9)+ylab("phylogenetic similarity")+xlab("% invaded"),nrow=2,labels=c("c)","d)"))

ggpubr::ggarrange(ccc,aa2,ncol=1,heights=c(1,1.5))


jpeg("Analyses/Master/Plots/phyandtax_angioonly.jpeg",width = 11,height = 8, unit="in",res=300)
ggpubr::ggarrange(ccc,aa2,ncol=1,heights=c(1,1.5))
dev.off()


natmainq2tax<-filter(nattaxdatq2,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))
natmainq2phy<-filter(natphydatq2,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))


nattoy_mainq2tax<-filter(nattoy_pred.l4.tax2,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))
nattoy_mainq2phy<-filter(nattoy_pred.l4.phy2,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))


nattoy_mainq2tax$status<-"native"
toy_mainq2tax$status<-"all"


nattoy_mainq2phy$status<-"native"
toy_mainq2phy$status<-"all"

toyphy<-rbind(nattoy_mainq2phy,toy_mainq2phy)

natmainq2phy$status<-"native"
mainq2phy$status<-"all"
natphy<-rbind(natmainq2phy,mainq2phy)
toyphy$grouper<-paste(toyphy$grouper,toyphy$status)



toytax<-rbind(nattoy_mainq2tax,toy_mainq2tax)

natmainq2tax$status<-"native"
mainq2tax$status<-"all"
nattax<-rbind(natmainq2tax,mainq2tax)
toytax$grouper<-paste(toytax$grouper,toytax$status)


jpeg("Analyses/Master/Plots/q2comps_angios-only.jpeg",width = 11,height = 6, unit="in",res=300)

ggpubr::ggarrange(ggplot()+
                    #geom_point(data=nattax,aes(x=perc_invaded*100,y=local_similarity,color=status),size=0.1)+
                    geom_line(data=toytax,aes(x=perc_invaded*100,y=.epred,color=status,group=grouper),alpha=0.8,size=0.002)+
                    
                    geom_smooth(data=toytax,aes(x=perc_invaded*100,y=.epred,color=status),size=2)+
                    facet_wrap(~NA_L1NAME,nrow=1,scales="free_y")+
                    #coord_cartesian(ylim=c(.1,.6))+
                    ggthemes::theme_few(base_size = 9)+ylab("taxonomic similarity")+xlab("% invaded")+scale_color_manual(values=c("black","forestgreen")),
                  
                  ggplot()+
                    #geom_point(data=natphy,aes(x=perc_invaded*100,y=local_similarity,color=status),size=0.1)+
                    geom_line(data=toyphy,aes(x=perc_invaded*100,y=.epred,color=status,group=grouper),alpha=0.8,size=0.002)+
                    
                    geom_smooth(data=toyphy,aes(x=perc_invaded*100,y=.epred,color=status),size=2)+
                    facet_wrap(~NA_L1NAME,nrow=1,scales="free_y")+
                    #coord_cartesian(ylim=c(.7,1))+
                    ggthemes::theme_few(base_size = 9)+
                    ylab("phylogenetic similarity")+xlab("% invaded")+
                    scale_color_manual(values=c("black","forestgreen")),ncol=1,common.legend = TRUE)

dev.off()


save.image("Analyses/Master/angiosperms.Rda")
