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
#load("Analyses/Master/hogzmain.Rda")
load("Analyses/Master/hogzmainwDataset.Rda")
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


phy<-read.tree("Analyses/Master/Input/FINAL_PHYLOGENY_SpCodes_Apr2023.tr") #read in tree


d<-filter(d,SpCode %in% c(phy$tip.label)) ##filter out species from dataset that aren't in phylogeny

###July adding dataset to account for different sizes and sampling methods
##Total # of plots per ecoregion
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

L4tots$regionID<-paste(L4tots$US_L4NAME,L4tots$NA_L1NAME,Dataset)
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
ecoregionz2<-unique(inyonly$regionID)[-4]
ecoregionz2<-ecoregionz2[-234]
ecoregionz2<-ecoregionz2[-437]
taxdat.invy<-data.frame()
phydat.invy<-data.frame()

ecoregionz2<-ecoregionz2[550:808]
ecoregionz2<-ecoregionz2[175:259]
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



###models for q1


taxdatq2<-filter(taxdat,q==2)
phydatq2<-filter(phydat,q==2)

nattaxdatq2<-filter(taxdat.nat,q==2)
natphydatq2<-filter(phydat.nat,q==2)

if(FALSE){
nattaxdatq1<-filter(taxdat.nat,q==1)
natphydatq1<-filter(phydat.nat,q==1)

nattaxdatq0<-filter(taxdat.nat,q==0)
natphydatq0<-filter(phydat.nat,q==0)


taxdatq1<-filter(taxdat,q==1)
phydatq1<-filter(phydat,q==1)

taxdatq0<-filter(taxdat,q==0)
phydatq0<-filter(phydat,q==0)
}

if(FALSE){
q1.tax<- brm(
  bf(local_similarity ~n_scale+ perc_invaded+(perc_invaded|NA_L1NAME),
     phi ~n_scale+ perc_invaded+(perc_invaded|NA_L1NAME)),
  data = taxdatq1,
  family = Beta(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr")


q0.tax<- brm(
  bf(local_similarity ~n_scale+ perc_invaded+(perc_invaded|NA_L1NAME),
     phi ~n_scale+ perc_invaded+(perc_invaded|NA_L1NAME)),
  data = taxdatq0,
  family = Beta(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr")
}
q2.tax<- brm(
  bf(local_similarity ~n_scale+ perc_invaded+(perc_invaded|NA_L1NAME)+(1|Dataset),
     phi ~n_scale+ perc_invaded+(perc_invaded|NA_L1NAME)+(1|Dataset)),
  data = taxdatq2,
  family = Beta(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr")

######phlo
if(FALSE){
q1.phy<- brm(
  bf(local_similarity ~n_scale+ perc_invaded+(perc_invaded|NA_L1NAME),
     phi ~ n_scale+perc_invaded+(perc_invaded|NA_L1NAME)),
  data = phydatq1,
  family = Beta(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr")

q0.phy<- brm(
  bf(local_similarity ~n_scale+ perc_invaded+(perc_invaded|NA_L1NAME),
     phi ~ n_scale+perc_invaded+(perc_invaded|NA_L1NAME)),
  data = phydatq0,
  family = Beta(),
  control=list(adapt_delta=.995),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr")
}
q2.phy<- brm(
  bf(local_similarity ~n_scale+ perc_invaded+(perc_invaded|NA_L1NAME)+(1|Dataset),
     phi ~ n_scale+perc_invaded+(perc_invaded|NA_L1NAME)+(1|Dataset)),
  data = phydatq2,
  family = Beta(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr")


###native species only
if(FALSE){
q1.tax.nat<- brm(
  bf(local_similarity ~n_scale+ perc_invaded+(perc_invaded|NA_L1NAME),
     phi ~n_scale+ perc_invaded+(perc_invaded|NA_L1NAME)),
  data = nattaxdatq1,
  family = Beta(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr")

q0.tax.nat<- brm(
  bf(local_similarity ~n_scale+ perc_invaded+(perc_invaded|NA_L1NAME),
     phi ~n_scale+ perc_invaded+(perc_invaded|NA_L1NAME)),
  data = nattaxdatq0,
  family = Beta(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr")
}

nattaxdatq2<-dplyr::filter(nattaxdatq2,local_similarity>0)

q2.tax.nat<- brm(
  bf(local_similarity ~n_scale+ perc_invaded+(perc_invaded|NA_L1NAME)+(1|Dataset),
     phi ~n_scale+ perc_invaded+(perc_invaded|NA_L1NAME)+(1|Dataset)),
  data = nattaxdatq2,
  family = Beta(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr")

if(FALSE){
q1.phy.nat<- brm(
  bf(local_similarity ~n_scale+ perc_invaded+(perc_invaded|NA_L1NAME),
     phi ~ n_scale+perc_invaded+(perc_invaded|NA_L1NAME)),
  data = natphydatq1,
  family = Beta(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr")

q0.phy.nat<- brm(
  bf(local_similarity ~n_scale+ perc_invaded+(perc_invaded|NA_L1NAME),
     phi ~ n_scale+perc_invaded+(perc_invaded|NA_L1NAME)),
  data = natphydatq0,
  family = Beta(),
  control=list(adapt_delta=.995),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr")
}
q2.phy.nat<- brm(
  bf(local_similarity ~n_scale+ perc_invaded+(perc_invaded|NA_L1NAME)+(1|Dataset),
     phi ~ n_scale+perc_invaded+(perc_invaded|NA_L1NAME)+(1|Dataset)),
  data = natphydatq2,
  family = Beta(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr")



###########################################################

new.data<-data.frame(NA_L1NAME=rep(unique(taxdat$NA_L1NAME),each=11),perc_invaded=rep(c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1),10))
new.data$n_scale<-median(taxdat$n_scale)

new.data2<-data.frame(perc_invaded=c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1))
new.data2$n_scale<-median(taxdat$n_scale)

#### predictions
library(tidybayes)
if(FALSE){
toy_pred.l4.phy<- q1.phy %>% 
  epred_draws(newdata =new.data,ndraws = 4000)
toy_pred.l4.tax<- q1.tax %>% 
  epred_draws(newdata =new.data,ndraws = 4000)
}

toy_pred.l4.phy2<- q2.phy %>% 
  epred_draws(newdata =new.data,ndraws = 4000,re_formula = ~(perc_invaded|NA_L1NAME))
toy_pred.l4.tax2<- q2.tax %>% 
  epred_draws(newdata =new.data,ndraws = 4000,re_formula = ~(perc_invaded|NA_L1NAME))
if(FALSE){
toy_pred.l4.phy.noreg<- q1.phy %>% 
  epred_draws(newdata =new.data2,ndraws = 4000,re_formula = NA)
toy_pred.l4.tax.noreg<- q1.tax %>% 
  epred_draws(newdata =new.data2,ndraws = 4000,re_formula = NA)

toy_pred.l4.phy.noreg0<- q0.phy %>% 
  epred_draws(newdata =new.data2,ndraws = 4000,re_formula = NA)
toy_pred.l4.tax.noreg0<- q0.tax %>% 
  epred_draws(newdata =new.data2,ndraws = 4000,re_formula = NA)
}
toy_pred.l4.phy.noreg2<- q2.phy %>% 
  epred_draws(newdata =new.data2,ndraws = 4000,re_formula = NA)
toy_pred.l4.tax.noreg2<- q2.tax %>% 
  epred_draws(newdata =new.data2,ndraws = 4000,re_formula = NA)




if(FALSE){
toy_pred.l4.tax$class<-"taxonomic"
toy_pred.l4.phy$class<-"phylogenetic"
}
toy_pred.l4.tax2$class<-"taxonomic"
toy_pred.l4.phy2$class<-"phylogenetic"

if(FALSE){
toy_pred.l4.tax.noreg$class<-"taxonomic"
toy_pred.l4.phy.noreg$class<-"phylogenetic"

toy_pred.l4.tax.noreg0$class<-"taxonomic"
toy_pred.l4.phy.noreg0$class<-"phylogenetic"
}
toy_pred.l4.tax.noreg2$class<-"taxonomic"
toy_pred.l4.phy.noreg2$class<-"phylogenetic"

if(FALSE){
toy_pred.l4.phy$grouper<-paste(toy_pred.l4.phy$NA_L1NAME,toy_pred.l4.phy$class,toy_pred.l4.phy$.draw)
toy_pred.l4.tax$grouper<-paste(toy_pred.l4.tax$NA_L1NAME,toy_pred.l4.tax$class,toy_pred.l4.tax$.draw)
}
toy_pred.l4.phy2$grouper<-paste(toy_pred.l4.phy2$NA_L1NAME,toy_pred.l4.phy2$class,toy_pred.l4.phy2$.draw)
toy_pred.l4.tax2$grouper<-paste(toy_pred.l4.tax2$NA_L1NAME,toy_pred.l4.tax2$class,toy_pred.l4.tax2$.draw)

if(FALSE){
toy_pred.l4.phy.noreg$grouper<-paste(toy_pred.l4.phy.noreg$class,toy_pred.l4.phy.noreg$.draw)
toy_pred.l4.tax.noreg$grouper<-paste(toy_pred.l4.tax.noreg$class,toy_pred.l4.tax.noreg$.draw)
}





if(FALSE){
aaa<-ggpubr::ggarrange(ggplot()+
                         geom_point(data=taxdatq1,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
                         
                         geom_line(data=toy_pred.l4.tax.noreg,aes(x=perc_invaded*100,y=.epred, group=grouper),size=.01,alpha=.5)+
                         coord_cartesian(ylim=c(0,1))+ggthemes::theme_few(base_size = 10)+geom_smooth(method="glm",color="black")+ylab("taxonomic similarity")+xlab("% invaded"),
                       
                       
                       ggplot()+
                         geom_point(data=phydatq1,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
                         
                         geom_line(data=toy_pred.l4.phy.noreg,aes(x=perc_invaded*100,y=.epred, group=grouper),size=.01,alpha=.5)+
                         coord_cartesian(ylim=c(0.7,1))+ggthemes::theme_few(base_size = 10)+geom_smooth(method="glm",color="black")+ylab("phylogenetic similarity")+xlab("% invaded"),ncol=2)

bbb<-ggpubr::ggarrange(ggplot()+
                         geom_point(data=taxdatq0,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
                         
                         geom_line(data=toy_pred.l4.tax.noreg0,aes(x=perc_invaded*100,y=.epred, group=grouper),size=.01,alpha=.5)+
                         coord_cartesian(ylim=c(0,1))+ggthemes::theme_few(base_size = 10)+geom_smooth(method="glm",color="black")+ylab("taxonomic similarity")+xlab("% invaded"),
                       
                       
                       ggplot()+
                         geom_point(data=phydatq0,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
                         
                         geom_line(data=toy_pred.l4.phy.noreg0,aes(x=perc_invaded*100,y=.epred, group=grouper),size=.01,alpha=.5)+
                         coord_cartesian(ylim=c(0.7,1))+ggthemes::theme_few(base_size = 10)+geom_smooth(method="glm",color="black")+ylab("phylogenetic similarity")+xlab("% invaded"),ncol=2)


}
ccc<-ggpubr::ggarrange(ggplot()+
                    geom_point(data=taxdatq2,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
                    
                    geom_line(data=toy_pred.l4.tax.noreg2,aes(x=perc_invaded*100,y=.epred, group=grouper),size=.01,alpha=.5)+
                      geom_smooth(data=toy_pred.l4.tax.noreg2,aes(x=perc_invaded*100,y=.epred),size=1,color="darkgray")+
                      coord_cartesian(ylim=c(0,.5))+ggthemes::theme_few(base_size = 10)+ylab("taxonomic similarity")+xlab("% invaded"),
                  
                  
                  ggplot()+
                    geom_point(data=phydatq2,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
                    
                    geom_line(data=toy_pred.l4.phy.noreg2,aes(x=perc_invaded*100,y=.epred, group=grouper),size=.01,alpha=.5)+
                    geom_smooth(data=toy_pred.l4.phy.noreg2,aes(x=perc_invaded*100,y=.epred),size=1,color="darkgray")+
                     coord_cartesian(ylim=c(.6,1))+ggthemes::theme_few(base_size = 10)+ylab("phylogenetic similarity")+xlab("% invaded"),ncol=2,labels=c("a)","b)"))

phydatq2$metric<-"inverse simpson"
phydatq0$metric<-"occurence"
phydatq1$metric<-"shannon"
toy_pred.l4.phy.noreg2$metric<-"inverse simpson"
toy_pred.l4.phy.noreg$metric<-"shannon"
toy_pred.l4.phy.noreg0$metric<-"occurence"


taxdatq2$metric<-"inverse simpson"
taxdatq0$metric<-"occurence"
taxdatq1$metric<-"shannon"
colnames(taxdatq0)
taxdatq11<-select(taxdatq1,-status)
toy_pred.l4.tax.noreg2$metric<-"inverse simpson"
toy_pred.l4.tax.noreg$metric<-"shannon"
toy_pred.l4.tax.noreg0$metric<-"occurence"


toy_pred.l4.phy.noreg0$grouper<-paste(toy_pred.l4.phy.noreg0$class,toy_pred.l4.phy.noreg0$.draw,toy_pred.l4.phy.noreg0$metric)
toy_pred.l4.phy.noreg$grouper<-paste(toy_pred.l4.phy.noreg$class,toy_pred.l4.phy.noreg$.draw,toy_pred.l4.phy.noreg$metric)

toy_pred.l4.tax.noreg$grouper<-paste(toy_pred.l4.tax.noreg$class,toy_pred.l4.tax.noreg$.draw,toy_pred.l4.tax.noreg$metric)

toy_pred.l4.tax.noreg0$grouper<-paste(toy_pred.l4.tax.noreg0$class,toy_pred.l4.tax.noreg0$.draw,toy_pred.l4.tax.noreg0$metric)

toy_pred.l4.phy.noreg2$grouper<-paste(toy_pred.l4.phy.noreg2$class,toy_pred.l4.phy.noreg2$.draw,toy_pred.l4.phy.noreg2$metric)
toy_pred.l4.tax.noreg2$grouper<-paste(toy_pred.l4.tax.noreg2$class,toy_pred.l4.tax.noreg2$.draw,toy_pred.l4.tax.noreg2$metric)

toy_pred.l4.phy.noregcomp<-rbind(toy_pred.l4.phy.noreg2,toy_pred.l4.phy.noreg0,toy_pred.l4.phy.noreg)
toy_pred.l4.tax.noregcomp<-rbind(toy_pred.l4.tax.noreg2,toy_pred.l4.tax.noreg0,toy_pred.l4.tax.noreg)

recap<-toy_pred.l4.phy.noreg0%>%group_by(perc_invaded) %>%summarise(delta=mean(.epred))


phydatcomp<-rbind(phydatq0,phydatq1,phydatq2)
taxdatcomp<-rbind(taxdatq0,taxdatq11,taxdatq2)

jpeg("Analyses/Master/Plots/metriccomp.jpeg",width = 6,height = 3, unit="in",res=300)
ggpubr::ggarrange(ggplot()+
  geom_point(data=taxdatcomp,aes(x=perc_invaded*100,y=local_similarity,color=metric),size=.1)+
  geom_line(data=toy_pred.l4.tax.noregcomp,aes(x=perc_invaded*100,y=.epred, group=grouper,color=metric),size=.01,alpha=.2)+
  geom_smooth(data=toy_pred.l4.tax.noregcomp,aes(x=perc_invaded*100,y=.epred,color=metric),size=1)+
  coord_cartesian(ylim=c(0.1,.9))+ggthemes::theme_few(base_size = 10)+ylab("taxonomic similarity")+xlab("% invaded")+scale_color_viridis_d(),


ggplot()+
  geom_point(data=phydatcomp,aes(x=perc_invaded*100,y=local_similarity,color=metric),size=.1)+
  geom_line(data=toy_pred.l4.phy.noregcomp,aes(x=perc_invaded*100,y=.epred, group=grouper,color=metric),size=.01,alpha=.2)+
  geom_smooth(data=toy_pred.l4.phy.noregcomp,aes(x=perc_invaded*100,y=.epred,color=metric),size=1)+
  coord_cartesian(ylim=c(0.7,.95))+ggthemes::theme_few(base_size = 10)+ylab("phylogenetic similarity")+xlab("% invaded")+scale_color_viridis_d(),ncol=2,common.legend = TRUE)

dev.off()


jpeg("Analyses/Master/Plots/q2allsps.jpeg",width = 6,height = 3, unit="in",res=300)
ccc
dev.off()
ggpubr::ggarrange(bbb,aaa,ccc,ncol=1,labels = c("a)","b)","c)"))

round(fixef(q1.tax,probs = c(.25,.75,.05,.95)),2)#-
round(fixef(q1.phy,probs = c(.25,.75,.05,.95)),2)#**

round(fixef(q0.tax,probs = c(.25,.75,.05,.95)),2)#-
round(fixef(q0.phy,probs = c(.25,.75,.05,.95)),2)#* negative

round(fixef(q2.tax,probs = c(.25,.75,.05,.95)),2)#*
round(fixef(q2.phy,probs = c(.25,.75,.05,.95)),2)#**

coef(q2.tax,probs = c(.05,.95))#*
coef(q2.phy,probs = c(.05,.95))#*

#q2tax mu
round(c(-0.5632522, 0.2313231, -0.9425615, -0.1734276),2)
round(c( 0.4403011, 0.1956091,  0.1211382 , 0.7555733),2)
round(c( 0.9398839, 0.1756180,  0.6341130, 1.2165341),2)
round(c(-0.3271651, 0.4698648, -1.1117038,  0.4392207),2)
round(c( 0.6184956, 0.2023728,  0.2933694,  0.9541079),2)
round(c( 0.8030746, 0.5468270, -0.1086974,  1.6848628),2)
round(c( 0.1800381, 0.3447946, -0.4087163,  0.7002513),2)
round(c( -0.1020124, 0.8605675, -1.5901859,  1.1509557),2)
round(c( -0.1237078, 0.6918279, -1.3218815 , 0.9123928),2)
round(c( 0.1285131 ,0.9976464,-1.5455668,  1.6501564),2)

#q2tax phi
round(c(0.8640548, 0.4720603,  0.14157172, 1.6723029),2)
round(c( -0.1005729, 0.4822938, -0.93843762, 0.6114199),2)
round(c( 0.1694626, 0.3626762, -0.47468990, 0.7449577),2)
round(c(0.6250187, 1.0183465, -0.69571920, 2.6630189),2)
round(c(0.6488497, 0.4462245, -0.02803795, 1.4373219),2)
round(c(0.5395969, 0.8179024, -0.72677690, 1.9750389),2)
round(c( 0.2504195, 0.6498087, -0.88796785, 1.2255285),2)
round(c( 0.4256878, 0.9771595, -1.07823200, 2.0405426),2)
round(c(0.2242332, 1.0301006, -1.47974140, 1.7767599),2)
round(c(0.6219327, 0.9658447, -0.63637585, 2.2691923),2)

coef(q2.phy,probs = c(.05,.95))#*
round(c( 0.6405518, 0.1663682,  0.35515650, 0.9045313),2)
round(c( 1.0606976, 0.2295371,  0.71936996, 1.4611367),2)
round(c( 0.8127987, 0.1728507,  0.52877195, 1.0883322),2)
round(c(0.7200951, 0.2684366,  0.23243185, 1.1244768),2)
round(c(0.9500321, 0.1703042,  0.68508032, 1.2443046),2)
round(c(0.7470003, 0.3203231,  0.17650175, 1.2156824),2)
round(c(0.6485267, 0.4047323, -0.08613385, 1.2208505),2)
round(c(0.8078752, 0.4621890,  0.10055265, 1.5040407),2)
round(c( 0.9814089, 0.4034235,  0.38735835, 1.6739195),2)
round(c(0.9427012, 0.5046454,  0.25691905, 1.7677540),2)
            
            
round(fixef(q2.tax.nat,probs = c(.05,.95)),2)#*
coef(q2.tax.nat,probs = c(.05,.95))
coef(q1.tax.nat,probs = c(.05,.95))
coef(q0.tax.nat,probs = c(.05,.95))

round(fixef(q1.tax.nat,probs = c(.05,.95)),2)#*
round(fixef(q0.tax.nat,probs = c(.05,.95)),2)#*


round(fixef(q2.phy.nat,probs = c(.05,.95)),2)#**
coef(q2.phy.nat,probs = c(.05,.95))
coef(q1.phy.nat,probs = c(.05,.95))
coef(q0.phy.nat,probs = c(.05,.95))

round(fixef(q1.phy.nat,probs = c(.05,.95)),2)#**
round(fixef(q0.phy.nat,probs = c(.05,.95)),2)#**

nattoy_pred.l4.phy<- q1.phy.nat %>% 
  epred_draws(newdata =new.data,ndraws = 4000)
nattoy_pred.l4.tax<- q1.tax.nat %>% 
  epred_draws(newdata =new.data,ndraws = 4000)


nattoy_pred.l4.phy2<- q2.phy.nat %>% 
  epred_draws(newdata =new.data,ndraws = 4000,re_formula = ~(perc_invaded|NA_L1NAME))
nattoy_pred.l4.tax2<- q2.tax.nat %>% 
  epred_draws(newdata =new.data,ndraws = 4000,re_formula = ~(perc_invaded|NA_L1NAME))

if(FALSE){
nattoy_pred.l4.phy.noreg<- q1.phy.nat %>% 
  epred_draws(newdata =new.data2,ndraws = 4000,re_formula = NA)
nattoy_pred.l4.tax.noreg<- q1.tax.nat %>% 
  epred_draws(newdata =new.data2,ndraws = 4000,re_formula = NA)


nattoy_pred.l4.phy.noreg0<- q0.phy.nat %>% 
  epred_draws(newdata =new.data2,ndraws = 4000,re_formula = NA)
nattoy_pred.l4.tax.noreg0<- q0.tax.nat %>% 
  epred_draws(newdata =new.data2,ndraws = 4000,re_formula = NA)

}
nattoy_pred.l4.phy.noreg2<- q2.phy.nat %>% 
  epred_draws(newdata =new.data2,ndraws = 4000,re_formula = NA)
nattoy_pred.l4.tax.noreg2<- q2.tax.nat %>% 
  epred_draws(newdata =new.data2,ndraws = 4000,re_formula = NA)


if (FALSE){
  nattoy_pred.l4.tax$class<-"taxonomic"
nattoy_pred.l4.phy$class<-"phylogenetic"
nattoy_pred.l4.tax.noreg$class<-"taxonomic"
nattoy_pred.l4.phy.noreg$class<-"phylogenetic"

}
nattoy_pred.l4.tax2$class<-"taxonomic"
nattoy_pred.l4.phy2$class<-"phylogenetic"
nattoy_pred.l4.tax.noreg2$class<-"taxonomic"
nattoy_pred.l4.phy.noreg2$class<-"phylogenetic"

#nattoy_pred.l4.phy$grouper<-paste(nattoy_pred.l4.phy$NA_L1NAME,nattoy_pred.l4.phy$class,nattoy_pred.l4.phy$.draw)
#nattoy_pred.l4.tax$grouper<-paste(nattoy_pred.l4.tax$NA_L1NAME,nattoy_pred.l4.tax$class,nattoy_pred.l4.tax$.draw)

nattoy_pred.l4.phy2$grouper<-paste(nattoy_pred.l4.phy2$NA_L1NAME,nattoy_pred.l4.phy2$class,nattoy_pred.l4.phy2$.draw)
nattoy_pred.l4.tax2$grouper<-paste(nattoy_pred.l4.tax2$NA_L1NAME,nattoy_pred.l4.tax2$class,nattoy_pred.l4.tax2$.draw)


#nattoy_pred.l4.phy.noreg$grouper<-paste(nattoy_pred.l4.phy.noreg$class,nattoy_pred.l4.phy.noreg$.draw)
#nattoy_pred.l4.tax.noreg$grouper<-paste(nattoy_pred.l4.tax.noreg$class,nattoy_pred.l4.tax.noreg$.draw)

nattoy_pred.l4.phy.noreg2$grouper<-paste(nattoy_pred.l4.phy.noreg2$class,nattoy_pred.l4.phy.noreg2$.draw)
nattoy_pred.l4.tax.noreg2$grouper<-paste(nattoy_pred.l4.tax.noreg2$class,nattoy_pred.l4.tax.noreg2$.draw)





aab<-ggpubr::ggarrange(ggplot()+
                         geom_point(data=nattaxdatq1,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
                         
                         geom_line(data=nattoy_pred.l4.tax.noreg,aes(x=perc_invaded*100,y=.epred, group=grouper),size=.01,alpha=.5,color="darkgreen")+
                         coord_cartesian(ylim=c(.3,.8))+ggthemes::theme_few(base_size = 10)+geom_smooth(method="glm",color="black")+ylab("taxonomic similarity")+xlab("% invaded"),
                       
                       
                       ggplot()+
                         geom_point(data=natphydatq1,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
                         
                         geom_line(data=nattoy_pred.l4.phy.noreg,aes(x=perc_invaded*100,y=.epred, group=grouper),size=.01,alpha=.5,color="darkgreen")+
                         coord_cartesian(ylim=c(.8,.95))+ggthemes::theme_few(base_size = 10)+geom_smooth(method="glm",color="darkgreen")+ylab("phylogenetic similarity")+xlab("% invaded"),ncol=2)


aab2<-ggpubr::ggarrange(ggplot()+
                         geom_point(data=nattaxdatq2,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
                         
                         geom_line(data=nattoy_pred.l4.tax.noreg2,aes(x=perc_invaded*100,y=.epred, group=grouper),size=.01,alpha=.5,color="darkgreen")+
                          geom_smooth(data=nattoy_pred.l4.tax.noreg2,aes(x=perc_invaded*100,y=.epred),size=1,color="darkseagreen3")+
                          coord_cartesian(ylim=c(.1,.4))+ggthemes::theme_few(base_size = 10)+geom_smooth(method="glm",color="black")+ylab("taxonomic similarity")+xlab("% invaded"),
                       
                       
                       ggplot()+
                         geom_point(data=natphydatq2,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
                         
                         geom_line(data=nattoy_pred.l4.phy.noreg2,aes(x=perc_invaded*100,y=.epred, group=grouper),size=.01,alpha=.5,color="darkgreen")+
                         geom_smooth(data=nattoy_pred.l4.phy.noreg2,aes(x=perc_invaded*100,y=.epred),size=1,color="darkseagreen3")+
                         coord_cartesian(ylim=c(.7,.95))+ggthemes::theme_few(base_size = 10)+ylab("phylogenetic similarity")+xlab("% invaded"),ncol=2)


jpeg("Analyses/Master/Plots/q2allnatives.jpeg",width = 6,height = 3, unit="in",res=300)
aab2
dev.off()



###regions
mainq1tax<-filter(taxdatq1,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))
mainq1phy<-filter(phydatq1,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))

mainq2tax<-filter(taxdatq2,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))
mainq2phy<-filter(phydatq2,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))


toy_mainq1tax<-filter(toy_pred.l4.tax,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))
toy_mainq1phy<-filter(toy_pred.l4.phy,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))

toy_mainq2tax<-filter(toy_pred.l4.tax2,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))
toy_mainq2phy<-filter(toy_pred.l4.phy2,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))


aa<-ggpubr::ggarrange(ggplot()+
                        geom_point(data=mainq1tax,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
                        
                        geom_line(data=toy_mainq1tax,aes(x=perc_invaded*100,y=.epred, group=grouper),size=.01,alpha=.5)+facet_wrap(~NA_L1NAME,nrow=5)+
                        coord_cartesian(ylim=c(.1,.8))+ggthemes::theme_few(base_size = 9)+geom_smooth(method="glm",color="black")+ylab("taxonomic similarity")+xlab("% invaded"),
                      
                      
                      ggplot()+
                        geom_point(data=mainq1phy,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
                        
                        geom_line(data=toy_mainq1phy,aes(x=perc_invaded*100,y=.epred, group=grouper),size=.01,alpha=.5)+facet_wrap(~NA_L1NAME,nrow=5)+
                        coord_cartesian(ylim=c(.7,1))+ggthemes::theme_few(base_size = 9)+geom_smooth(method="glm",color="black")+ylab("phylogenetic similarity")+xlab("% invaded"),ncol=2)

jpeg("Analyses/Master/Plots/taxoplots.jpeg",width = 9,height = 5.5, unit="in",res=300)
ggpubr::ggarrange(ggplot()+
  geom_point(data=taxdatq2,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
  
  geom_line(data=toy_pred.l4.tax.noreg2,aes(x=perc_invaded*100,y=.epred, group=grouper),size=.01,alpha=.5)+
  geom_smooth(data=toy_pred.l4.tax.noreg2,aes(x=perc_invaded*100,y=.epred),size=1,color="darkgray")+
  coord_cartesian(ylim=c(0,.6))+ggthemes::theme_few(base_size = 10)+ylab("taxonomic similarity")+xlab("% invaded"),


ggplot()+
  geom_point(data=mainq2tax,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
  
  geom_line(data=toy_mainq2tax,aes(x=perc_invaded*100,y=.epred, group=grouper),size=.01,alpha=.5)+facet_wrap(~NA_L1NAME,nrow=1)+
  geom_smooth(data=toy_mainq2tax,aes(x=perc_invaded*100,y=.epred),size=1,color="gray")+
  coord_cartesian(ylim=c(.05,.6))+ggthemes::theme_few(base_size = 8)+ylab("taxonomic similarity")+xlab("% invaded"),nrow=2,heights=c(1.5,2),labels=c("a)","b)"))
dev.off()

jpeg("Analyses/Master/Plots/phyloplots.jpeg",width = 9,height = 5.5, unit="in",res=300)
ggpubr::ggarrange(
ggplot()+
  geom_point(data=phydatq2,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
  
  geom_line(data=toy_pred.l4.phy.noreg2,aes(x=perc_invaded*100,y=.epred, group=grouper),size=.01,alpha=.5)+
  geom_smooth(data=toy_pred.l4.phy.noreg2,aes(x=perc_invaded*100,y=.epred),size=1,color="darkgray")+
  coord_cartesian(ylim=c(0.7,1))+ggthemes::theme_few(base_size = 10)+ylab("phylogenetic similarity")+xlab("% invaded"),

ggplot()+
  geom_point(data=mainq2phy,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
  
  geom_line(data=toy_mainq2phy,aes(x=perc_invaded*100,y=.epred, group=grouper),size=.01,alpha=.5)+facet_wrap(~NA_L1NAME,ncol=5)+
  geom_smooth(data=toy_mainq2phy,aes(x=perc_invaded*100,y=.epred),size=1,color="gray")+
  coord_cartesian(ylim=c(.7,1))+ggthemes::theme_few(base_size = 8)+ylab("phylogenetic similarity")+xlab("% invaded"),nrow=2,heights=c(1.5,2),labels=c("a)","b)"))
dev.off()


aa2<-ggpubr::ggarrange(ggplot()+
                        geom_point(data=mainq2tax,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
                        
                        geom_line(data=toy_mainq2tax,aes(x=perc_invaded*100,y=.epred, group=grouper),size=.01,alpha=.5)+facet_wrap(~NA_L1NAME,ncol=5)+
                         geom_smooth(data=toy_mainq2tax,aes(x=perc_invaded*100,y=.epred),size=1,color="gray")+
                         coord_cartesian(ylim=c(0,.6))+ggthemes::theme_few(base_size = 9)+ylab("taxonomic similarity")+xlab("% invaded"),
                      
                      
                      ggplot()+
                        geom_point(data=mainq2phy,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
                        
                        geom_line(data=toy_mainq2phy,aes(x=perc_invaded*100,y=.epred, group=grouper),size=.01,alpha=.5)+facet_wrap(~NA_L1NAME,ncol=5)+
                        geom_smooth(data=toy_mainq2phy,aes(x=perc_invaded*100,y=.epred),size=1,color="gray")+
                        coord_cartesian(ylim=c(.7,1))+ggthemes::theme_few(base_size = 9)+ylab("phylogenetic similarity")+xlab("% invaded"),nrow=2,labels=c("c)","d)"))

jpeg("Analyses/Master/Plots/q2big5all.jpeg",width = 5,height = 6, unit="in",res=300)
aa2
dev.off()


jpeg("Analyses/Master/Plots/phyandtax.jpeg",width = 11,height = 8, unit="in",res=300)
ggpubr::ggarrange(ccc,aa2,ncol=1,heights=c(1,1.5))
dev.off()










natmainq1tax<-filter(nattaxdatq1,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))
natmainq1phy<-filter(natphydatq1,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))

nattoy_mainq1tax<-filter(nattoy_pred.l4.tax,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))
nattoy_mainq1phy<-filter(nattoy_pred.l4.phy,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))



natmainq2tax<-filter(nattaxdatq2,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))
natmainq2phy<-filter(natphydatq2,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))


nattoy_mainq2tax<-filter(nattoy_pred.l4.tax2,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))
nattoy_mainq2phy<-filter(nattoy_pred.l4.phy2,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))



ab<-ggpubr::ggarrange(ggplot()+
                        geom_point(data=natmainq1tax,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
                        
                        geom_line(data=nattoy_mainq1tax,color="darkgreen",aes(x=perc_invaded*100,y=.epred, group=grouper),size=.01,alpha=.5)+facet_wrap(~NA_L1NAME,nrow=5)+
                        coord_cartesian(ylim=c(.25,.8))+ggthemes::theme_few(base_size = 9)+geom_smooth(method="glm",color="darkgreen")+ylab("taxonomic similarity")+xlab("% invaded"),
                      
                      
                      ggplot()+
                        geom_point(data=natmainq1phy,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
                        
                        geom_line(data=nattoy_mainq1phy,color="darkgreen",aes(x=perc_invaded*100,y=.epred, group=grouper),size=.01,alpha=.5)+facet_wrap(~NA_L1NAME,nrow=5)+
                        coord_cartesian(ylim=c(.8,1))+ggthemes::theme_few(base_size = 9)+geom_smooth(method="glm",color="darkgreen")+ylab("phylogenetic similarity")+xlab("% invaded"),ncol=2)

save.image("Analyses/Master/hogzmain.Rda")
save.image("Analyses/Master/hogzmainwDataset.Rda")
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


jpeg("Analyses/Master/Plots/q2comps.jpeg",width = 11,height = 6, unit="in",res=300)

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
coef(q2.tax.nat,probs = c(.25,.75,.05,.95))
coef(q2.phy.nat,probs = c(.25,.75,.05,.95))


taxinvyq2<-filter(taxdat.invy,q==2)

#taxinvyq2<-filter(taxinvyq2,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))



phyinvyq2<-filter(phydat.invy,q==2)
#phyinvyq2<-filter(phyinvyq2,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))
q2.phy.invy<- brm(
  bf(local_similarity ~n_scale+ perc_invaded+(perc_invaded|NA_L1NAME)+(1|Dataset),
     phi ~ n_scale+perc_invaded+(perc_invaded|NA_L1NAME)+(1|Dataset)),
  data = phyinvyq2,
  family = Beta(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr")

taxinvyq2<-filter(taxinvyq2,local_similarity>0)

q2.tax.invy<- brm(
  bf(local_similarity ~n_scale+ perc_invaded+(perc_invaded|NA_L1NAME)+(1|Dataset),
     phi ~ n_scale+perc_invaded+(perc_invaded|NA_L1NAME)+(1|Dataset)),
  data = taxinvyq2,
  family = Beta(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr")


pred.inv.phy2noran<- q2.phy.invy %>% 
  epred_draws(newdata =new.data,ndraws = 4000,re_formula = NA)

pred.inv.tax2.noran<- q2.tax.invy %>% 
  epred_draws(newdata =new.data,ndraws = 4000,re_formula = NA)

pred.inv.tax2.noran$class<-"taxonomic"
pred.inv.phy2noran$class<-"phylogenetic"




pred.inv.phy2noran$grouper<-paste(pred.inv.phy2noran$class,pred.inv.phy2noran$.draw)
pred.inv.tax2.noran$grouper<-paste(pred.inv.tax2.noran$class,pred.inv.tax2.noran$.draw)

pred.inv.phy2<- q2.phy.invy %>% 
  epred_draws(newdata =new.data,ndraws = 4000,re_formula = ~(perc_invaded|NA_L1NAME))

pred.inv.tax2<- q2.tax.invy %>% 
  epred_draws(newdata =new.data,ndraws = 4000,re_formula = ~(perc_invaded|NA_L1NAME))

pred.inv.tax2$class<-"taxonomic"
pred.inv.phy2$class<-"phylogenetic"

pred.inv.phy2$grouper<-paste(pred.inv.phy2$NA_L1NAME,pred.inv.phy2$class,pred.inv.phy2$.draw)
pred.inv.tax2$grouper<-paste(pred.inv.tax2$NA_L1NAME,pred.inv.tax2$class,pred.inv.tax2$.draw)

invasivesslopes<-rbind(pred.inv.phy2noran,pred.inv.tax2.noran)

ggplot()+
  #geom_point(data=taxinvyq2,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
 
  geom_line(data=invasivesslopes,aes(x=perc_invaded*100,y=.epred,group=grouper,color=class),size=.01,alpha=.5)+#facet_wrap(~NA_L1NAME,ncol=5)+
  geom_smooth(data=invasivesslopes,aes(x=perc_invaded*100,y=.epred, color=class),size=1)+scale_color_viridis_d()+
  coord_cartesian(ylim=c(0,.8))+ggthemes::theme_few(base_size = 9)+ylab("biotic similarity")+xlab("% invaded")


ggpubr::ggarrange(ggplot()+
                    geom_point(data=taxinvyq2,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
                    
                    geom_line(data=pred.inv.tax2.noran,aes(x=perc_invaded*100,y=.epred,group=grouper),size=.01,alpha=.5)+#facet_wrap(~NA_L1NAME,ncol=5)+
                    geom_smooth(data=pred.inv.tax2.noran,aes(x=perc_invaded*100,y=.epred),size=1,color="firebrick4")+
                    coord_cartesian(ylim=c(0,.8))+ggthemes::theme_few(base_size = 9)+ylab("taxonomic similarity")+xlab("% invaded"),
                  
                  
                  ggplot()+
                    geom_point(data=phyinvyq2,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
                    
                    geom_line(data=pred.inv.phy2noran,aes(x=perc_invaded*100,y=.epred,group=grouper),size=.01,alpha=.5)+#facet_wrap(~NA_L1NAME,ncol=5)+
                    geom_smooth(data=pred.inv.phy2noran,aes(x=perc_invaded*100,y=.epred),size=1,color="firebrick4")+
                    coord_cartesian(ylim=c(0,.9))+ggthemes::theme_few(base_size = 9)+ylab("phylogenetic similarity")+xlab("% invaded"),nrow=2,labels=c("c)","d)"))






pred.inv.phy2$grouper<-paste(pred.inv.phy2$NA_L1NAME,pred.inv.phy2$class,pred.inv.phy2$.draw)
pred.inv.tax2$grouper<-paste(pred.inv.tax2$NA_L1NAME,pred.inv.tax2$class,pred.inv.tax2$.draw)


pred.inv.tax2<-filter(pred.inv.tax2, NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))
pred.inv.phy2<-filter(pred.inv.phy2, NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))

short.tax2<-filter(taxinvyq2, NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))
short.inv.phy2<-filter(phyinvyq2, NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))



invytrends<-ggpubr::ggarrange(ggplot()+
                    geom_point(data=short.tax2,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
                    
                    geom_line(data=pred.inv.tax2,aes(x=perc_invaded*100,y=.epred, group=grouper),size=.002,alpha=.9,color="firebrick4")+facet_wrap(~NA_L1NAME,ncol=5)+
                    geom_smooth(data=pred.inv.tax2,aes(x=perc_invaded*100,y=.epred),size=1,color="firebrick4")+
                    coord_cartesian(ylim=c(0,.8))+ggthemes::theme_few(base_size = 9)+ylab("taxonomic similarity")+xlab("% invaded"),
                  
                  
                  ggplot()+
                    geom_point(data=short.inv.phy2,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
                    
                    geom_line(data=pred.inv.phy2,aes(x=perc_invaded*100,y=.epred, group=grouper),size=.002,alpha=.8,color="firebrick4")+facet_wrap(~NA_L1NAME,ncol=5)+
                    geom_smooth(data=pred.inv.phy2,aes(x=perc_invaded*100,y=.epred),size=1,color="firebrick4")+
                    coord_cartesian(ylim=c(0,.9))+ggthemes::theme_few(base_size = 9)+ylab("phylogenetic similarity")+xlab("% invaded"),nrow=2,labels=c("c)","d)"))






gammamod<-brm(
  bf(TD_gamma ~n_scale+(1|NA_L1NAME)+(1|Dataset)),
  data = taxinvyq2,
  family =  gaussian(),
  control=list(adapt_delta=.95),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr")

alphamod<-brm(
  bf(TD_alpha ~n_scale+(1|NA_L1NAME)+(1|Dataset)),
  data = taxinvyq2,
  family =  gaussian(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr")

#gammamodphy<-brm(
  bf(PD_gamma ~n_scale+(1|NA_L1NAME)+(1|Dataset)),
  data = phyinvyq2,
  family =  gaussian(),
  control=list(adapt_delta=.95),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr")

alphamodphy<-brm(
  bf(PD_alpha ~n_scale+(1|NA_L1NAME)+(1|Dataset)),
  data = phyinvyq2,
  family =  gaussian(),
  control=list(adapt_delta=.95),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr")


expdat<-data.frame(NA_L1NAME=c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))
expdat$n_scale<-median(taxinvyq2$n_scale)

predgamma<-gammamod %>% 
  epred_draws(newdata =expdat,ndraws = 1000,re_formula =~(1|NA_L1NAME) )

predalpha<-alphamod %>% 
  epred_draws(newdata =expdat,ndraws = 1000,re_formula =~(1|NA_L1NAME))

predgammaphy<-gammamodphy %>% 
  epred_draws(newdata =expdat,ndraws = 1000,re_formula =~(1|NA_L1NAME))

ggplot(predgamma,aes(.epred,reorder(NA_L1NAME,.epred)))+stat_pointinterval(.width = c(.5,.9))+ggthemes::theme_few(base_size = 9)+xlab("non-native gamma diversity")+ylab("")


taxyinvq22<-filter(taxinvyq2, local_similarity>0)
simmod<-brm(
  bf(local_similarity ~n_scale+(1|NA_L1NAME)+(1|Dataset),
     phi ~ n_scale+(1|NA_L1NAME)+(1|Dataset)),
  data = taxyinvq22,
  control=list(adapt_delta=.99),
  family =  Beta(),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr")

#simmodphy<-brm(
#  bf(local_similarity ~n_scale+(1|NA_L1NAME)+(1|Dataset),
 #    phi ~ n_scale+(1|NA_L1NAME)+(1|Dataset)),
#  data = phyinvyq2,
#  control=list(adapt_delta=.99),
#  family =  Beta(),
#  chains = 4, iter = 4000, warmup = 3000,
 # cores = 4, seed = 1234,backend = "cmdstanr")

predsim<-simmod %>% 
  epred_draws(newdata =expdat,ndraws = 1000,re_formula =~(1|NA_L1NAME))

#predsimphy<-simmodphy %>% 
 # epred_draws(newdata =expdat,ndraws = 1000,re_formula =~(1|NA_L1NAME))

predalpha$guide<-predgamma$`.epred`
pbeta<-ggplot(predsim,aes(.epred,reorder(NA_L1NAME,-.epred)))+stat_pointinterval(.width = c(.5,.9),color="firebrick4")+ggthemes::theme_few(base_size = 9)+xlab("non-native taxonomic similarity")+ylab("")
#pbetaphy<-ggplot(predsimphy,aes(.epred,reorder(NA_L1NAME,-.epred)))+stat_pointinterval(.width = c(.5,.9))+ggthemes::theme_few(base_size = 9)+xlab("non-native taxonomic similarity")+ylab("")


pgamma<-ggplot(predgamma,aes(.epred,reorder(NA_L1NAME,.epred)))+stat_pointinterval(.width = c(.5,.9),color="firebrick4")+ggthemes::theme_few(base_size = 9)+xlab("non-native gamma diversity")+ylab("")+
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank())
palpha<-ggplot(predalpha,aes(.epred,reorder(NA_L1NAME,guide)))+stat_pointinterval(.width = c(.5,.9),color="firebrick4")+ggthemes::theme_few(base_size = 9)+
  xlab("non-native alpha diversity")+ylab("")+theme(axis.text.y = element_blank(),axis.ticks.y = element_blank())

#jpeg("Analyses/Master/Plots/q2explaintheinvies.jpeg",width = 8,height = 3, unit="in",res=300)

placeholder<-ggplot()+ggthemes::theme_few()


#jpeg("Analyses/Master/Plots/q2explaintheinvies.jpeg",width = 10,height = 11, unit="in",res=300)
pdf("Analyses/Master/Plots/q2explaintheinvies.pdf",width = 10,height = 11)
ggpubr::ggarrange(
ggpubr::ggarrange(pbeta,pgamma,palpha,ncol=3,widths = c(.35,.2,.2)),placeholder,invytrends,ncol=1,labels=c("a)","",""))
dev.off()

###map maker

library(sf)
#demo(nc, ask = FALSE, echo = FALSE)
#plot(nc["AREA"])
#View(nc)
regie<-st_read("Data/us_eco_l4 (1)/us_eco_l4_no_st.shx")
regie2<-left_join(regie,taxdatq2)

regiea <- regie2 %>% select(perc_invaded,local_similarity, geometry)
colnames(regiea)[2]<-"tax_similarity"

regie3<-left_join(regie,phydatq2)
regieb <- regie3 %>% select(perc_invaded,local_similarity, geometry)
colnames(regieb)[2]<-"phy_similarity"

regiec<-st_join(regiea, regieb,largest=TRUE)
regied<-regiec %>% gather(var,scale , -geometry)
regied<-filter(regied,var!="perc_invaded.y")
jpeg("Analyses/Master/Plots/maps1.jpeg",width = 10,height=9,unit='in',res=200)
ggplot() + 
  geom_sf(data = regied, aes(fill = scale),lwd = 0.01)+scale_fill_viridis_c(option="plasma")+
  ggthemes::theme_few()+theme(legend.position = "top")+
  facet_wrap(~var,ncol=3,labeller = labeller(var = c(`perc_invaded.x` = "% invaded",`phy_similarity`=" phylogenetic\nsimilarity",`tax_similarity`="taxonomic\nsimilarity")))
dev.off()

pp1<-ggplot() + 
  geom_sf(data = regiea, aes(fill = perc_invaded),lwd = 0.01)+scale_fill_viridis_c()+ggthemes::theme_few()+theme(legend.position = "top")+labs(fill = "% invaded") 

pp2<-    ggplot() + 
  geom_sf(data = regiea, aes(fill = tax_similarity),lwd = 0.01)+scale_fill_viridis_c(option = "A")+ggthemes::theme_few()+theme(legend.position = "top")+labs(fill = "taxonomic similarity") 

pp3<-    ggplot() + 
  geom_sf(data = regieb, aes(fill = phy_similarity),lwd = 0.01)+scale_fill_viridis_c(option = "C")+ggthemes::theme_few()+theme(legend.position = "top")+labs(fill = "phylogenetic similarity") 

jpeg("Analyses/Master/Plots/maps2.jpeg",width = 10,height=9,unit='in',res=200)
ggpubr::ggarrange(pp1,pp2,pp3,nrow=1,ncol=3)
dev.off()

invasivelist<-d.invy %>% group_by(NA_L1NAME,SpCode, AcceptedTaxonName) %>%summarise(ave_cov=mean(PctCov_100))

family<-read.csv("~/OneDrive - University of Massachusetts/Desktop/SpatialEcology/data/L1/data/NorthernAmerica2_G_F.csv")
family<-select(family,Taxon, Accepted_family)
family2<-read.csv("~/OneDrive - University of Massachusetts/Desktop/SpatialEcology/data/L1/data/AsiaTemperate2_G_F.csv")
family2<-select(family2,Taxon, Accepted_family)
family3<-read.csv("~/OneDrive - University of Massachusetts/Desktop/SpatialEcology/data/L1/data/Europe2_G_F.csv")
family3<-select(family3,Taxon, Accepted_family)

family4<-read.csv("~/OneDrive - University of Massachusetts/Desktop/SpatialEcology/data/L1/data/AsiaTropical2_G_F.csv")
family4<-select(family4,Taxon, Accepted_family)

family5<-read.csv("~/OneDrive - University of Massachusetts/Desktop/SpatialEcology/data/L1/data/Pacific2_G_F.csv")
family5<-select(family5,Taxon, Accepted_family)

family6<-read.csv("~/OneDrive - University of Massachusetts/Desktop/SpatialEcology/data/L1/data/AsiaTropical2_G_F.csv")
family6<-select(family6,Taxon, Accepted_family)

family7<-read.csv("~/OneDrive - University of Massachusetts/Desktop/SpatialEcology/data/L1/data/Australasia2_G_F.csv")
family7<-select(family7,Taxon, Accepted_family)

family8<-read.csv("~/OneDrive - University of Massachusetts/Desktop/SpatialEcology/data/L1/data/SouthernAmerica2_G_F.csv")
family8<-select(family8,Taxon, Accepted_family)


colnames(family)[1]<-"AcceptedTaxonName"
colnames(family2)[1]<-"AcceptedTaxonName"
colnames(family3)[1]<-"AcceptedTaxonName"
colnames(family4)[1]<-"AcceptedTaxonName"
colnames(family5)[1]<-"AcceptedTaxonName"
colnames(family6)[1]<-"AcceptedTaxonName"
colnames(family7)[1]<-"AcceptedTaxonName"
colnames(family8)[1]<-"AcceptedTaxonName"

family<-rbind(family,family2,family3,family4,family5,family6,family7,family8)


family<-distinct(family)


invasivelist<-left_join(invasivelist,family)
table(is.na(invasivelist$Accepted_family))
invasivelist<-filter(invasivelist,!is.na(Accepted_family))
rank<-invasivelist %>% group_by(NA_L1NAME,Accepted_family) %>% count()
rankETF<-filter(rank,NA_L1NAME=="EASTERN TEMPERATE FORESTS")
rankETF$perc<-rankETF$n/sum(rankETF$n)

rankET<-rankETF %>%arrange(desc(perc)) 
rankET<-rankET[1:4,]

rankNAD<-filter(rank,NA_L1NAME=="NORTH AMERICAN DESERTS")
rankNAD$perc<-rankNAD$n/sum(rankNAD$n)
rankNAD$perc<-rankNAD$n/sum(rankNAD$n)

rankNA<-rankNAD %>%arrange(desc(perc)) 
rankNA<-rankNA[1:4,]
sum(rankNA$perc)
ggplot(invasivelist,aes(Accepted_family))+geom_histogram(stat="count")+facet_wrap(~NA_L1NAME)


type<-d %>% select(SpCode,Plot,NA_L1NAME,GrowthHabit,NativeStatus, PctCov_100)
type<-distinct(type)
type<-type%>% group_by(Plot,NA_L1NAME,GrowthHabit,NativeStatus) %>%summarise(cover=mean(PctCov_100))
type<-type%>% group_by(NA_L1NAME,GrowthHabit,NativeStatus) %>%summarise(cover=mean(cover))
type<-filter(type,GrowthHabit %in% c("Tree","Graminoid","Shrub","Forb/herb"))
type$class<-ifelse(type$GrowthHabit %in% c("Tree","Shrub"),"canopy","groundlayer")
type<-filter(type,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))
type<-filter(type, NativeStatus %in% c("I","N"))
unique(type$GrowthHabit)
ggplot(type,aes(class,cover))+geom_bar(aes(fill=NativeStatus),stat="identity",position="dodge")+facet_wrap(~NA_L1NAME,scales="free_y")

                                             