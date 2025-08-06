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
phy<-read.tree("Analyses/Master/Input/angiosperms.tre") #read in tree
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
L4invs<-d %>% group_by(Plot,NativeStatus,US_L4NAME,NA_L1NAME,Dataset)%>% dplyr::summarize(RelCov=sum(PctCov_100,na.rm=TRUE))
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
  saveRDS(raw,"Analyses/Master/pairwiseBdiv_allsps_angio.Rda") #19,558,091
  
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
  write.csv(full.dat,"Analyses/Master/public/allspeciesBdata_angio.csv",row.names = FALSE)
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
saveRDS(raw2,"Analyses/Master/pairwiseBdiv_nativeonly_angio.Rda") #19,558,091

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


L4tots$inv_prev<-L4tots$n_inv/L4tots$n
regio<-dplyr::select(L4tots,US_L4NAME,NA_L1NAME,Dataset,regionID,inv_prev)
regio<-distinct(regio)

mydat<-left_join(mydat,regio)
write.csv(mydat,"Analyses/Master/public/bothspeciesBdata_angio.csv",row.names = FALSE)


p$regionID<-paste(p$US_L4NAME,p$NA_L1NAME,p$Dataset)


invcount<-p %>% group_by(regionID) %>% count()
invcount<-filter(invcount,n>1)



invcount2<-d.invy%>% group_by(regionID,AcceptedTaxonName) %>% count()
invcount2<-invcount2%>% group_by(regionID) %>% count()
invcount2<-filter(invcount2,n>1)

taxdat.inv<-data.frame()
phydat.inv<-data.frame()

ecoregionz<-intersect(invcount$regionID,invcount2$regionID)
ecoregionz<-ecoregionz[17:1151]
ecoregionz<-ecoregionz[41:1135]
ecoregionz<-ecoregionz[107:1095]
ecoregionz<-ecoregionz[83:989]
ecoregionz<-ecoregionz[35:907]
ecoregionz<-ecoregionz[2:873]
ecoregionz<-ecoregionz[56:872]
ecoregionz<-ecoregionz[166:817]
ecoregionz<-ecoregionz[18:652]
ecoregionz<-ecoregionz[37:635]
ecoregionz<-ecoregionz[68:599]
ecoregionz<-ecoregionz[27:532]
ecoregionz<-ecoregionz[61:506]
ecoregionz<-ecoregionz[12:446]
ecoregionz<-ecoregionz[4:435]
ecoregionz<-ecoregionz[38:432]
ecoregionz<-ecoregionz[2:395]
ecoregionz<-ecoregionz[28:394]
ecoregionz<-ecoregionz[49:367]
ecoregionz<-ecoregionz[172:319]
ecoregionz<-ecoregionz[2:148]
ecoregionz<-ecoregionz[84:147]
ecoregionz<-ecoregionz[49:64]

for (z in c(1:length(ecoregionz))){
  dat<-dplyr::filter(d.invy,regionID==ecoregionz[z])
  
  
  comms<-matricize(dat)
  
  
  taxy<-hill_taxa_parti_pairwise(comms,q = 2)
  taxy$regionID<-ecoregionz[z]
  taxy$class<-"taxonomic"
  
  phyt<-hill_phylo_parti_pairwise(comms,phy,q=2)
  phyt$regionID<-ecoregionz[z]
  phyt$class<-"phylogenetic"
  
  taxdat.inv<-rbind(taxy,taxdat.inv)
  phydat.inv<-rbind(phyt,phydat.inv)
  
}  

raw3<-list(taxdat.inv,phydat.inv)
saveRDS(raw2,"Analyses/Master/pairwiseBdiv_invonly_angio.Rda") #19,558,091

colnames(taxdat.inv)

tax.data.invy<-taxdat.inv %>% group_by(regionID) %>%summarise(
  similarity= mean(local_similarity, na.rm=TRUE),
  gamma= mean(TD_gamma, na.rm=TRUE),
  alpha= mean(TD_alpha, na.rm=TRUE),
  sd_sim = sd(local_similarity,na.rm=TRUE),
  sd_gam = sd(TD_gamma,na.rm=TRUE),
  sd_alph = sd(TD_alpha,na.rm=TRUE),
  n = n(),
  se_sim = sd_sim / sqrt(n),
  se_gam = sd_gam / sqrt(n),
  se_alph = sd_alph / sqrt(n),
  )

head(tax.data.invy)


phy.data.invy<-phydat.inv %>% group_by(regionID) %>%summarise(
  similarity= mean(local_similarity, na.rm=TRUE),
  gamma= mean(PD_gamma, na.rm=TRUE),
  alpha= mean(PD_alpha, na.rm=TRUE),
  sd_sim = sd(local_similarity,na.rm=TRUE),
  sd_gam = sd(PD_gamma,na.rm=TRUE),
  sd_alph = sd(PD_alpha,na.rm=TRUE),
  n = n(),
  se_sim = sd_sim / sqrt(n),
  se_gam = sd_gam / sqrt(n),
  se_alph = sd_alph / sqrt(n),
)


tax.data.invy$metric<-"taxonomic"
phy.data.invy$metric<-"phylogenetic"


full.dat.invy<-rbind(tax.data.invy,phy.data.invy)

mydat.invy<-left_join(full.dat.invy,regio)
write.csv(mydat.invy,"Analyses/Master/public/invaderBdata_angio.csv",row.names = FALSE)


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



