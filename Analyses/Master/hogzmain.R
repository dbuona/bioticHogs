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

setwd("~/Documents/git/bioticHogs/")

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


phy<-read.tree("Analyses/Master/Input/FINAL_PHYLOGENY_SpCodes_Apr2023.tr") #read in tree


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



###models for q1

taxdatq1<-filter(taxdat,q==1)
phydatq1<-filter(phydat,q==1)

taxdatq0<-filter(taxdat,q==0)
phydatq0<-filter(phydat,q==0)

taxdatq2<-filter(taxdat,q==2)
phydatq2<-filter(phydat,q==2)

nattaxdatq1<-filter(taxdat.nat,q==1)
natphydatq1<-filter(phydat.nat,q==1)

nattaxdatq0<-filter(taxdat.nat,q==0)
natphydatq0<-filter(phydat.nat,q==0)

nattaxdatq2<-filter(taxdat.nat,q==2)
natphydatq2<-filter(phydat.nat,q==2)


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

q2.tax<- brm(
  bf(local_similarity ~n_scale+ perc_invaded+(perc_invaded|NA_L1NAME),
     phi ~n_scale+ perc_invaded+(perc_invaded|NA_L1NAME)),
  data = taxdatq2,
  family = Beta(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr")

######phlo

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

q2.phy<- brm(
  bf(local_similarity ~n_scale+ perc_invaded+(perc_invaded|NA_L1NAME),
     phi ~ n_scale+perc_invaded+(perc_invaded|NA_L1NAME)),
  data = phydatq2,
  family = Beta(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr")


###native species only
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

q2.tax.nat<- brm(
  bf(local_similarity ~n_scale+ perc_invaded+(perc_invaded|NA_L1NAME),
     phi ~n_scale+ perc_invaded+(perc_invaded|NA_L1NAME)),
  data = nattaxdatq2,
  family = Beta(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr")


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

q2.phy.nat<- brm(
  bf(local_similarity ~n_scale+ perc_invaded+(perc_invaded|NA_L1NAME),
     phi ~ n_scale+perc_invaded+(perc_invaded|NA_L1NAME)),
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
toy_pred.l4.phy<- q1.phy %>% 
  epred_draws(newdata =new.data,ndraws = 4000)
toy_pred.l4.tax<- q1.tax %>% 
  epred_draws(newdata =new.data,ndraws = 4000)

toy_pred.l4.phy.noreg<- q1.phy %>% 
  epred_draws(newdata =new.data2,ndraws = 4000,re_formula = NA)
toy_pred.l4.tax.noreg<- q1.tax %>% 
  epred_draws(newdata =new.data2,ndraws = 4000,re_formula = NA)

toy_pred.l4.phy.noreg0<- q0.phy %>% 
  epred_draws(newdata =new.data2,ndraws = 4000,re_formula = NA)
toy_pred.l4.tax.noreg0<- q0.tax %>% 
  epred_draws(newdata =new.data2,ndraws = 4000,re_formula = NA)

toy_pred.l4.phy.noreg2<- q2.phy %>% 
  epred_draws(newdata =new.data2,ndraws = 4000,re_formula = NA)
toy_pred.l4.tax.noreg2<- q2.tax %>% 
  epred_draws(newdata =new.data2,ndraws = 4000,re_formula = NA)





toy_pred.l4.tax$class<-"taxonomic"
toy_pred.l4.phy$class<-"phylogenetic"

toy_pred.l4.tax.noreg$class<-"taxonomic"
toy_pred.l4.phy.noreg$class<-"phylogenetic"

toy_pred.l4.tax.noreg0$class<-"taxonomic"
toy_pred.l4.phy.noreg0$class<-"phylogenetic"

toy_pred.l4.tax.noreg2$class<-"taxonomic"
toy_pred.l4.phy.noreg2$class<-"phylogenetic"

toy_pred.l4.phy$grouper<-paste(toy_pred.l4.phy$NA_L1NAME,toy_pred.l4.phy$class,toy_pred.l4.phy$.draw)
toy_pred.l4.tax$grouper<-paste(toy_pred.l4.tax$NA_L1NAME,toy_pred.l4.tax$class,toy_pred.l4.tax$.draw)


toy_pred.l4.phy.noreg$grouper<-paste(toy_pred.l4.phy.noreg$class,toy_pred.l4.phy.noreg$.draw)
toy_pred.l4.tax.noreg$grouper<-paste(toy_pred.l4.tax.noreg$class,toy_pred.l4.tax.noreg$.draw)

toy_pred.l4.phy.noreg0$grouper<-paste(toy_pred.l4.phy.noreg0$class,toy_pred.l4.phy.noreg0$.draw)
toy_pred.l4.tax.noreg0$grouper<-paste(toy_pred.l4.tax.noreg0$class,toy_pred.l4.tax.noreg0$.draw)

toy_pred.l4.phy.noreg2$grouper<-paste(toy_pred.l4.phy.noreg2$class,toy_pred.l4.phy.noreg2$.draw)
toy_pred.l4.tax.noreg2$grouper<-paste(toy_pred.l4.tax.noreg2$class,toy_pred.l4.tax.noreg2$.draw)





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



ccc<-ggpubr::ggarrange(ggplot()+
                    geom_point(data=taxdatq2,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
                    
                    geom_line(data=toy_pred.l4.tax.noreg2,aes(x=perc_invaded*100,y=.epred, group=grouper),size=.01,alpha=.5)+
                    coord_cartesian(ylim=c(0,1))+ggthemes::theme_few(base_size = 10)+geom_smooth(method="glm",color="black")+ylab("taxonomic similarity")+xlab("% invaded"),
                  
                  
                  ggplot()+
                    geom_point(data=phydatq2,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
                    
                    geom_line(data=toy_pred.l4.phy.noreg2,aes(x=perc_invaded*100,y=.epred, group=grouper),size=.01,alpha=.5)+
                    coord_cartesian(ylim=c(0.7,1))+ggthemes::theme_few(base_size = 10)+geom_smooth(method="glm",color="black")+ylab("phylogenetic similarity")+xlab("% invaded"),ncol=2)



ggpubr::ggarrange(bbb,aaa,ccc,ncol=1,labels = c("a)","b)","c)"))

round(fixef(q1.tax,probs = c(.25,.75,.05,.95)),2)#-
round(fixef(q1.phy,probs = c(.25,.75,.05,.95)),2)#**

round(fixef(q0.tax,probs = c(.25,.75,.05,.95)),2)#-
round(fixef(q0.phy,probs = c(.25,.75,.05,.95)),2)#* negative

round(fixef(q2.tax,probs = c(.25,.75,.05,.95)),2)#*
round(fixef(q2.phy,probs = c(.25,.75,.05,.95)),2)#**

nattoy_pred.l4.phy<- q1.phy.nat %>% 
  epred_draws(newdata =new.data,ndraws = 4000)
nattoy_pred.l4.tax<- q1.tax.nat %>% 
  epred_draws(newdata =new.data,ndraws = 4000)

nattoy_pred.l4.phy.noreg<- q1.phy.nat %>% 
  epred_draws(newdata =new.data2,ndraws = 4000,re_formula = NA)
nattoy_pred.l4.tax.noreg<- q1.tax.nat %>% 
  epred_draws(newdata =new.data2,ndraws = 4000,re_formula = NA)


nattoy_pred.l4.phy.noreg0<- q0.phy.nat %>% 
  epred_draws(newdata =new.data2,ndraws = 4000,re_formula = NA)
nattoy_pred.l4.tax.noreg0<- q0.tax.nat %>% 
  epred_draws(newdata =new.data2,ndraws = 4000,re_formula = NA)


nattoy_pred.l4.phy.noreg2<- q2.phy.nat %>% 
  epred_draws(newdata =new.data2,ndraws = 4000,re_formula = NA)
nattoy_pred.l4.tax.noreg2<- q2.tax.nat %>% 
  epred_draws(newdata =new.data2,ndraws = 4000,re_formula = NA)



nattoy_pred.l4.tax$class<-"taxonomic"
nattoy_pred.l4.phy$class<-"phylogenetic"

nattoy_pred.l4.tax.noreg$class<-"taxonomic"
nattoy_pred.l4.phy.noreg$class<-"phylogenetic"

nattoy_pred.l4.phy$grouper<-paste(nattoy_pred.l4.phy$NA_L1NAME,nattoy_pred.l4.phy$class,nattoy_pred.l4.phy$.draw)
nattoy_pred.l4.tax$grouper<-paste(nattoy_pred.l4.tax$NA_L1NAME,nattoy_pred.l4.tax$class,nattoy_pred.l4.tax$.draw)

nattoy_pred.l4.phy.noreg$grouper<-paste(nattoy_pred.l4.phy.noreg$class,nattoy_pred.l4.phy.noreg$.draw)
nattoy_pred.l4.tax.noreg$grouper<-paste(nattoy_pred.l4.tax.noreg$class,nattoy_pred.l4.tax.noreg$.draw)



aab<-ggpubr::ggarrange(ggplot()+
                         geom_point(data=nattaxdatq1,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
                         
                         geom_line(data=nattoy_pred.l4.tax.noreg,aes(x=perc_invaded*100,y=.epred, group=grouper),size=.01,alpha=.5,color="darkgreen")+
                         coord_cartesian(ylim=c(.3,.8))+ggthemes::theme_few(base_size = 10)+geom_smooth(method="glm",color="black")+ylab("taxonomic similarity")+xlab("% invaded"),
                       
                       
                       ggplot()+
                         geom_point(data=natphydatq1,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
                         
                         geom_line(data=nattoy_pred.l4.phy.noreg,aes(x=perc_invaded*100,y=.epred, group=grouper),size=.01,alpha=.5,color="darkgreen")+
                         coord_cartesian(ylim=c(.8,.95))+ggthemes::theme_few(base_size = 10)+geom_smooth(method="glm",color="darkgreen")+ylab("phylogenetic similarity")+xlab("% invaded"),ncol=2)


###regions
mainq1tax<-filter(taxdatq1,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))
mainq1phy<-filter(phydatq1,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))

toy_mainq1tax<-filter(toy_pred.l4.tax,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))
toy_mainq1phy<-filter(toy_pred.l4.phy,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))

aa<-ggpubr::ggarrange(ggplot()+
                        geom_point(data=mainq1tax,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
                        
                        geom_line(data=toy_mainq1tax,aes(x=perc_invaded*100,y=.epred, group=grouper),size=.01,alpha=.5)+facet_wrap(~NA_L1NAME,nrow=5)+
                        coord_cartesian(ylim=c(.1,.8))+ggthemes::theme_few(base_size = 9)+geom_smooth(method="glm",color="black")+ylab("taxonomic similarity")+xlab("% invaded"),
                      
                      
                      ggplot()+
                        geom_point(data=mainq1phy,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
                        
                        geom_line(data=toy_mainq1phy,aes(x=perc_invaded*100,y=.epred, group=grouper),size=.01,alpha=.5)+facet_wrap(~NA_L1NAME,nrow=5)+
                        coord_cartesian(ylim=c(.7,1))+ggthemes::theme_few(base_size = 9)+geom_smooth(method="glm",color="black")+ylab("phylogenetic similarity")+xlab("% invaded"),ncol=2)




natmainq1tax<-filter(nattaxdatq1,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))
natmainq1phy<-filter(natphydatq1,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))

nattoy_mainq1tax<-filter(nattoy_pred.l4.tax,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))
nattoy_mainq1phy<-filter(nattoy_pred.l4.phy,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))


ab<-ggpubr::ggarrange(ggplot()+
                        geom_point(data=natmainq1tax,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
                        
                        geom_line(data=nattoy_mainq1tax,color="darkgreen",aes(x=perc_invaded*100,y=.epred, group=grouper),size=.01,alpha=.5)+facet_wrap(~NA_L1NAME,nrow=5)+
                        coord_cartesian(ylim=c(.25,.8))+ggthemes::theme_few(base_size = 9)+geom_smooth(method="glm",color="darkgreen")+ylab("taxonomic similarity")+xlab("% invaded"),
                      
                      
                      ggplot()+
                        geom_point(data=natmainq1phy,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
                        
                        geom_line(data=nattoy_mainq1phy,color="darkgreen",aes(x=perc_invaded*100,y=.epred, group=grouper),size=.01,alpha=.5)+facet_wrap(~NA_L1NAME,nrow=5)+
                        coord_cartesian(ylim=c(.8,1))+ggthemes::theme_few(base_size = 9)+geom_smooth(method="glm",color="darkgreen")+ylab("phylogenetic similarity")+xlab("% invaded"),ncol=2)

save.image("Analyses/Master/hogzmain.Rda")

nattoy_mainq1tax$status<-"native"
toy_mainq1tax$status<-"all"


nattoy_mainq1phy$status<-"native"
toy_mainq1phy$status<-"all"

toyphy<-rbind(nattoy_mainq1phy,toy_mainq1phy)

natmainq1phy$status<-"native"
mainq1phy$status<-"all"
natphy<-rbind(natmainq1phy,mainq1phy)
toyphy$grouper<-paste(toyphy$grouper,toyphy$status)



toytax<-rbind(nattoy_mainq1tax,toy_mainq1tax)

natmainq1tax$status<-"native"
mainq1tax$status<-"all"
nattax<-rbind(natmainq1tax,mainq1tax)
toytax$grouper<-paste(toytax$grouper,toytax$status)

ggpubr::ggarrange(ggplot()+
  geom_point(data=nattax,aes(x=perc_invaded*100,y=local_similarity,color=status),size=0.1)+
  geom_line(data=toytax,aes(x=perc_invaded*100,y=.epred,color=status,group=grouper),alpha=.4,size=0.01)+

  #geom_smooth(data=toytax,aes(x=perc_invaded*100,y=.epred,color=status))+
  facet_wrap(~NA_L1NAME,nrow=5)+
  coord_cartesian(ylim=c(.2,.7))+
  ggthemes::theme_few(base_size = 9)+ylab("taxonomic similarity")+xlab("% invaded")+scale_color_manual(values=c("black","forestgreen")),

ggplot()+
  geom_point(data=natphy,aes(x=perc_invaded*100,y=local_similarity,color=status),size=0.1)+
  geom_line(data=toyphy,aes(x=perc_invaded*100,y=.epred,color=status,group=grouper),alpha=.4,size=0.01)+
  
  #geom_smooth(data=toyphy,aes(x=perc_invaded*100,y=.epred,color=status))+
  facet_wrap(~NA_L1NAME,nrow=5)+
  coord_cartesian(ylim=c(.8,.95))+
  ggthemes::theme_few(base_size = 9)+
  ylab("taxonomic similarity")+xlab("% invaded")+
  scale_color_manual(values=c("black","forestgreen")),ncol=2,common.legend = TRUE)
