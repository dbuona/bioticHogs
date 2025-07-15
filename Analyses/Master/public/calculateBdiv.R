rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())
graphics.off()

setwd("~/Documents/git/bioticHogs/")

library(dplyr)
library(phytools)
library(hillR)
###read in data and phylogeny#########
d<-read.csv("Data/SPCIS_plant_taxa.csv") ### this is the public data for reproducibility
p<-read.csv("Data/SPCIS_plots.csv") ## plot identities
phy<-read.tree("Analyses/Master/Input/FINAL_PHYLOGENY_SpCodes_Apr2023.tr") #read in tree
regionz<-read.csv("Analyses/Master/Input/regionz.csv") ### plots and ecoregion correspondence
##################################

regionz<-dplyr::select(regionz,-X) #get rid of extra column

d<-left_join(d,regionz) ###merge all ecoregions to data
  d<-d %>% #select only the most recent survey for resampled plots
  group_by(Plot) %>%
  filter(Year==max(Year)) %>%
  ungroup() 


p<-left_join(p,regionz)  ###merge all ecoregions to data
  p<-p %>% #select only the most recent survey for resampled plots
  group_by(Plot) %>%
  filter(Year==max(Year)) %>%
  ungroup() 


d<-filter(d,SpCode %in% c(phy$tip.label)) ##filter out species from dataset that aren't in phylogeny
d<-left_join(d,p)

#######Below will before calculating invasion prevelence#####
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

L4tots$regionID<-paste(L4tots$US_L4NAME,L4tots$NA_L1NAME,L4tots$Dataset)
L4tots<-filter(L4tots,!is.na(US_L4NAME))
L4tots<-filter(L4tots,!is.na(NA_L1NAME))

d$regionID<-paste(d$US_L4NAME,d$NA_L1NAME,d$Dataset)

matricize<-function(x){
  temp<-dplyr::select(x,Plot,SpCode,PctCov_100) ### select useful columns
  temp<-dplyr::filter(temp,!is.na(temp$SpCode))
  temp<- tidyr::spread(temp,SpCode,PctCov_100)
  temp[is.na(temp)] <- 0 
  temp<-as.data.frame(temp)
  rownames(temp)<-temp$Plot # convert plot names to row names
  temp<-dplyr::select(temp,-Plot)
  temp<-as.matrix.data.frame(temp)#
} #function for formating data

if(FALSE){
taxdat<-data.frame()
phydat<-data.frame()


ecoregionz<-unique(L4tots$regionID)
ecoregionz<-ecoregionz[946:1213]

for (z in c(1:length(ecoregionz))){
  dat<-dplyr::filter(d,regionID==ecoregionz[z])
  
  
  comms<-matricize(dat)
  

  taxy<-hill_taxa_parti_pairwise(comms,q = 2)
  taxy$regionID<-ecoregionz[z]
  taxyclass<-"taxonomic"

  phyt<-hill_phylo_parti_pairwise(comms,phy,q=2)
  phyt$regionID<-ecoregionz[z]
  phyt$class<-"phylogenetic"
  
  taxdat<-rbind(taxy,taxdat)
  phydat<-rbind(phyt,phydat)
  
}  

raw<-list(taxdat,phydat)
saveRDS(raw,"Analyses/Master/pairwiseBdiv_allsps.Rda") #19,558,091

tax.data<-taxdat %>% group_by(regionID) %>%summarise(
  similarity= mean(local_similarity, na.rm=TRUE),
  sd = sd(local_similarity,na.rm=TRUE),
  n = n(),
  se = sd / sqrt(n))

phy.data<-phydat %>% group_by(regionID) %>%summarise(
  similarity= mean(local_similarity, na.rm=TRUE),
  sd = sd(local_similarity,na.rm=TRUE),
  n = n(),
  se = sd / sqrt(n))

tax.data$metric<-"taxonomic"
phy.data$metric<-"phylogenetic"

full.dat<-rbind(tax.data,phy.data)
full.dat$set<-"all species"
write.csv(full.dat,"Analyses/Master/public/allspeciesBdata.csv",row.names = FALSE)
}

d.invy<-filter(d,NativeStatus=="I") ### nonnative
d.native<-filter(d,NativeStatus!="I") ### native

###native species

taxdat.nat<-data.frame()
phydat.nat<-data.frame()

ecoregionz<-unique(L4tots$regionID)[-945]

for (z in c(1:length(ecoregionz))){
  dat<-dplyr::filter(d.native,regionID==ecoregionz[z])
  
  
  comms<-matricize(dat)
  
  
  taxy<-hill_taxa_parti_pairwise(comms,q = 2)
  taxy$regionID<-ecoregionz[z]
  taxy$class<-"taxonomic"
  
  phyt<-hill_phylo_parti_pairwise(comms,phy,q=2)
  phyt$regionID<-ecoregionz[z]
  phyt$class<-"phylogenetic"
  
  taxdat.nat<-rbind(taxy,taxdat.nat)
  phydat.nat<-rbind(phyt,phydat.nat)
  
}  

raw2<-list(taxdat.nat,phydat.nat)
saveRDS(raw2,"Analyses/Master/pairwiseBdiv_nativeonly.Rda") #19,558,091

tax.data.nat<-taxdat.nat %>% group_by(regionID) %>%summarise(
  similarity= mean(local_similarity, na.rm=TRUE),
  sd = sd(local_similarity,na.rm=TRUE),
  n = n(),
  se = sd / sqrt(n))

phy.data.nat<-phydat.nat %>% group_by(regionID) %>%summarise(
  similarity= mean(local_similarity, na.rm=TRUE),
  sd = sd(local_similarity,na.rm=TRUE),
  n = n(),
  se = sd / sqrt(n))

tax.data.nat$metric<-"taxonomic"
phy.data.nat$metric<-"phylogenetic"

full.dat.nat<-rbind(tax.data.nat,phy.data.nat)
full.dat.nat$set<-"natives"

mydat<-rbind(full.dat,full.dat.nat)

predo<-L4tots
predo$inv_prev<-predo$n_inv/predo$n

predo<-select(predo,-n,-n_inv)

mydat<-left_join(mydat,predo)

inv_rich<-d.invy %>% dplyr::select(US_L4NAME,AcceptedTaxonName)
inv_rich<-distinct(inv_rich)

gamma<-inv_rich %>% group_by(US_L4NAME) %>% count() 
gamma<-filter(gamma,!is.na(US_L4NAME))

colnames(gamma)[2]<-"inv_gamma"

inv_alpha<-d.invy %>% dplyr::select(Plot,regionID,AcceptedTaxonName)
inv_alpha<-distinct(inv_alpha)

inv_alpha<- inv_alpha %>% group_by(Plot,regionID) %>% count()

alpha<-inv_alpha %>% group_by(regionID) %>%summarise(
  meanI_alpha= mean(n, na.rm=TRUE),
  sdI_alpha = sd(n,na.rm=TRUE),
  nI_alpha = n(),
  seI_alpha = sdI_alpha / sqrt(nI_alpha))

alpha[is.na(alpha)] <- 0

mydat<-left_join(mydat,gamma)
mydat<-left_join(mydat,alpha)

write.csv(mydat,"Analyses/Master/public/biohogsData.csv",row.names = FALSE)



