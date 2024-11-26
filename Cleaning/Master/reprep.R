###start cleaning concept papers, match plots by characteristics

rm(list=ls())
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())
#rstan_options(auto_write = TRUE)

graphics.off()
library(dplyr)
library(tidyr)

library(ggplot2)

library(ape)
library(geodist)
library("MatchIt")

setwd("/Users/danielbuonaiuto/Documents/git/bioticHogs/")
p<-read.csv("Data/FULLDatabase_10272022_diversity_PctCov_100_10272022_INHABIT.csv")
d<-read.csv("Data/FULLDatabase_10272022.csv")

p<-p %>% #select only the most recent survey for resampled plots
  group_by(Plot) %>%
  filter(Year==max(Year)) %>%
  ungroup() 

d<-d %>% #select only the most recent survey for resampled plots
  group_by(Plot) %>%
  filter(Year==max(Year)) %>%
  ungroup() 

phy<-read.tree("Analyses/Master/Input/FINAL_PHYLOGENY_SpCodes_Apr2023.tr") #read in tree
d<-filter(d,SpCode %in% c(phy$tip.label)) ##filter out species from dataset that aren't in phylogeny


p.i<-filter(p,RelCov_I>=5)
p.p<-filter(p,RelCov_I<5)

d.i<-filter(d,Plot %in% c(p.i$Plot))
d.p<-filter(d,Plot %in% c(p.p$Plot))
d.r<-filter(d.i,NativeStatus=="N")
d.ro<-filter(d.i,NativeStatus=="I")

counters<- p%>% group_by(L4_KEY,Dataset) %>% count() %>% arrange(n)
colnames(counters)[3]<-"n_all"

counters.i<- p.i%>% group_by(L4_KEY,Dataset) %>% count() %>% arrange(n)
colnames(counters.i)[3]<-"n_i"
counters.p<- p.p%>% group_by(L4_KEY,Dataset) %>% count() %>% arrange(n)
colnames(counters.p)[3]<-"n_p"

counters<-left_join(counters,counters.i)
counters<-left_join(counters,counters.p)
counters<-filter(counters,n_i>1)
counters<-filter(counters,n_p>1)            
counters<-filter(counters,!is.na(L4_KEY))


p.i<-filter(p.i, L4_KEY %in% counters$L4_KEY)
p.i<-filter(p.i, Dataset %in% counters$Dataset)

p.p<-filter(p.p, L4_KEY %in% counters$L4_KEY)
p.p<-filter(p.p, Dataset %in% counters$Dataset)

d.i<-filter(d.i,Plot %in% c(p.i$Plot))
d.p<-filter(d.p,Plot %in% c(p.p$Plot))
d.r<-filter(d.r,Plot %in% c(p.i$Plot))
d.ro<-filter(d.ro,Plot %in% c(p.i$Plot))

###Next,make sure there are multiple species per site\

p.enough<-d.p%>%group_by(Plot) %>%count()
p.enough<-filter(p.enough,n>2)
r.enough<-d.r%>%group_by(Plot) %>%count()
r.enough<-filter(r.enough,n>2)


d.p<-filter(d.p,Plot %in% p.enough$Plot)
d.r<-filter(d.r,Plot %in% r.enough$Plot)


i.enough<-d.i%>%group_by(Plot) %>%count()
i.enough<-filter(i.enough,n>2)

d.i<-filter(d.i,Plot %in% i.enough$Plot)

p.r<-p.i
p.r<-filter(p.r,Plot %in% c(d.r$Plot))
p.p<-filter(p.p,Plot %in% c(d.p$Plot))
p.i<-filter(p.i,Plot %in% c(d.i$Plot))
p.ro<-p.i
p.ro<-filter(p.ro,Plot %in% c(d.ro$Plot))

#### Prepare formating for diversity
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


counters$group<-paste(counters$Dataset,counters$L4_KEY,sep="_")
######

p.i$group<-paste(p.i$Dataset,p.i$L4_KEY,sep="_")
plotref.i<-dplyr::select(p.i,Plot,Dataset,L4_KEY,group)
d.i<-left_join(d.i,plotref.i)



taxdat_i<-data.frame()
phydat_i<-data.frame()
ecoregionz<-sort(counters$group)
#NWCA_15o  Coeur d Alene Metasedimentary Zone
ecoregionz<-ecoregionz[458:692]
for (z in c(1:length(ecoregionz))){
  dat<-dplyr::filter(d.i,group==ecoregionz[z])
  
  
  comms<-matricize(dat)
  
  taxy2<-hill_taxa_parti_pairwise(comms,q = 2)
  
  taxy2$group<-ecoregionz[z]
  taxy2$class<-"taxonomic"
  

  phyt<-hill_phylo_parti_pairwise(comms,phy,q=2)
  
   phyt$group<-ecoregionz[z]
  phyt$class<-"phylogenetic"
  
  taxdat_i<-rbind(taxy2,taxdat_i)
  phydat_i<-rbind(phyt,phydat_i)
  
}  

write.csv(taxdat_i,"Analyses/Master/Input/pairwise_invaded_taxo.csv")
write.csv(phydat_i,"Analyses/Master/Input/pairwise_invaded_phylo.csv")

p.p$group<-paste(p.p$Dataset,p.p$L4_KEY,sep="_")
plotref.p<-dplyr::select(p.p,Plot,Dataset,L4_KEY,group)
d.p<-left_join(d.p,plotref.p)



taxdat_p<-data.frame()
phydat_p<-data.frame()
ecoregionz<-sort(counters$group)
#NWCA_15o  Coeur d Alene Metasedimentary Zone
ecoregionz<-ecoregionz[-457]
ecoregionz<-ecoregionz[40:691]
#AIM_9j  Klamath Juniper Woodland/Devils Garden
#CVS_66i  High Mountains

ecoregionz<-ecoregionz[109:652]
ecoregionz<-ecoregionz[25:544]
#"NPS_3a  Portland/Vancouver Basin"
ecoregionz<-ecoregionz[189:520]
#NWCA_21j  Grassland Parks
ecoregionz<-ecoregionz[126:332]

for (z in c(1:length(ecoregionz))){
  dat<-dplyr::filter(d.p,group==ecoregionz[z])
  
  
  comms<-matricize(dat)
  
  taxy2<-hill_taxa_parti_pairwise(comms,q = 2)
  
  taxy2$group<-ecoregionz[z]
  taxy2$class<-"taxonomic"
  
  
  phyt<-hill_phylo_parti_pairwise(comms,phy,q=2)
  
  phyt$group<-ecoregionz[z]
  phyt$class<-"phylogenetic"
  
  taxdat_p<-rbind(taxy2,taxdat_p)
  phydat_p<-rbind(phyt,phydat_p)
  
} 

write.csv(taxdat_p,"Analyses/Master/Input/pairwise_pristine_taxo.csv")
write.csv(phydat_p,"Analyses/Master/Input/pairwise_pristine_phylo.csv")



p.r$group<-paste(p.r$Dataset,p.r$L4_KEY,sep="_")
plotref.r<-dplyr::select(p.r,Plot,Dataset,L4_KEY,group)
d.r<-left_join(d.r,plotref.r)



taxdat_r<-data.frame()
phydat_r<-data.frame()
ecoregionz<-sort(counters$group)
#NWCA_15o  Coeur d Alene Metasedimentary Zone
# AIM_14f  Mojave Playas
ecoregionz<-ecoregionz[-457]
ecoregionz<-ecoregionz[40:691]
#NPS_26d  Semiarid Canadian Breaks
ecoregionz<-ecoregionz[311:652]
#NPS_78d  Serpentine Siskiyous
ecoregionz<-ecoregionz[87:342]
#NWCA_15l  Salish Mountains
ecoregionz<-ecoregionz[22:256]
#NWCA_63g  Carolinian Barrier Islands and Coastal Mars
ecoregionz<-ecoregionz[126:235]
#"NWCA_80g  High Lava Plains"
ecoregionz<-ecoregionz[63:110]


for (z in c(1:length(ecoregionz))){
  dat<-dplyr::filter(d.r,group==ecoregionz[z])
  
  
  comms<-matricize(dat)
  
  taxy2<-hill_taxa_parti_pairwise(comms,q = 2)
  
  taxy2$group<-ecoregionz[z]
  taxy2$class<-"taxonomic"
  
  
  phyt<-hill_phylo_parti_pairwise(comms,phy,q=2)
  
  phyt$group<-ecoregionz[z]
  phyt$class<-"phylogenetic"
  
  taxdat_r<-rbind(taxy2,taxdat_r)
  phydat_r<-rbind(phyt,phydat_r)
  
}  

write.csv(taxdat_r,"Analyses/Master/Input/pairwise_remnant_taxo.csv")
write.csv(phydat_r,"Analyses/Master/Input/pairwise_remnant_phylo.csv")



p.ro$group<-paste(p.ro$Dataset,p.ro$L4_KEY,sep="_")
plotref.ro<-dplyr::select(p.ro,Plot,Dataset,L4_KEY,group)
d.ro<-left_join(d.ro,plotref.ro)


taxdat_io<-data.frame()
phydat_io<-data.frame()

ecoregionz<-sort(counters$group)
#AIM_13w  Tonopah Uplands
ecoregionz<-ecoregionz[33:692]
#AIM_16d  Dry Partly Wooded Mountains
ecoregionz<-ecoregionz[13:660]
#AIM_21a  Alpine Zone
ecoregionz<-ecoregionz[34:648]
#AIM_22m  Albuquerque Basin
ecoregionz<-ecoregionz[20:615]
#AIM_9e  Pumice Plateau
ecoregionz<-ecoregionz[51:596]
#FIA_50n  Boundary Lakes and Hills
ecoregionz<-ecoregionz[72:546]
ecoregionz<-ecoregionz[2:475]
#NWCA_14f  Mojave Playas
ecoregionz<-ecoregionz[237:474]
#NWCA_15o  Coeur d Alene Metasedimentary Zone
ecoregionz<-ecoregionz[4:238]
#NWCA_22f  Taos Plateau
ecoregionz<-ecoregionz[30:235]
#NWCA_50x  Grand Marais Lakeshore
ecoregionz<-ecoregionz[70:206]
#NWCA_76b  Big Cypress
ecoregionz<-ecoregionz[83:137]
#VNHP_63d  Virginian Barrier Islands and Coastal Mars
ecoregionz<-ecoregionz[27:55]

for (z in c(1:length(ecoregionz))){
  dat<-dplyr::filter(d.ro,group==ecoregionz[z])
  
  
  comms<-matricize(dat)
  
  taxy2<-hill_taxa_parti_pairwise(comms,q = 2)
  
  taxy2$group<-ecoregionz[z]
  taxy2$class<-"taxonomic"
  
  
  phyt<-hill_phylo_parti_pairwise(comms,phy,q=2)
  
  phyt$group<-ecoregionz[z]
  phyt$class<-"phylogenetic"
  
  taxdat_io<-rbind(taxy2,taxdat_io)
  phydat_io<-rbind(phyt,phydat_io)
  
}  

write.csv(taxdat_io,"Analyses/Master/Input/pairwise_invaderonly_taxo.csv")
write.csv(phydat_io,"Analyses/Master/Input/pairwise_invaderonly_phylo.csv")


coords<-dplyr::select(d,Plot,Lat,Long)
coords<-distinct(coords)
colnames(coords)<-c("site1", "Lat.1", "Lon.1")
coords2<-coords
colnames(coords2)<-c("site2", "Lat.2", "Lon.2")



taxdat_i<-left_join(taxdat_i,coords)
taxdat_i<-left_join(taxdat_i,coords2)

phydat_i<-left_join(phydat_i,coords)
phydat_i<-left_join(phydat_i,coords2)

cord1<-data.frame(Lon=taxdat_i$Lon.1,Lat=taxdat_i$Lat.1)
cord2<-data.frame(Lon=taxdat_i$Lon.2,Lat=taxdat_i$Lat.2)

library(geodist)
disty1<-geodist(x=cord1,y=cord2,measure="haversine",paired = TRUE)
head(disty1)
taxdat_i$Distance<-disty1/100
phydat_i$Distance<-disty1/100

write.csv(taxdat_i,"Analyses/Master/Input/pairwise_invaded_taxo.csv")
write.csv(phydat_i,"Analyses/Master/Input/pairwise_invaded_phylo.csv")

taxdat_p<-left_join(taxdat_p,coords)
taxdat_p<-left_join(taxdat_p,coords2)

phydat_p<-left_join(phydat_p,coords)
phydat_p<-left_join(phydat_p,coords2)

cord1<-data.frame(Lon=taxdat_p$Lon.1,Lat=taxdat_p$Lat.1)
cord2<-data.frame(Lon=taxdat_p$Lon.2,Lat=taxdat_p$Lat.2)

disty1<-geodist(x=cord1,y=cord2,measure="haversine",paired = TRUE)

taxdat_p$Distance<-disty1/100
phydat_p$Distance<-disty1/100

write.csv(taxdat_p,"Analyses/Master/Input/pairwise_pristine_taxo.csv")
write.csv(phydat_p,"Analyses/Master/Input/pairwise_pristine_phylo.csv")

taxdat_r<-left_join(taxdat_r,coords)
taxdat_r<-left_join(taxdat_r,coords2)

phydat_r<-left_join(phydat_r,coords)
phydat_r<-left_join(phydat_r,coords2)

cord1<-data.frame(Lon=taxdat_r$Lon.1,Lat=taxdat_r$Lat.1)
cord2<-data.frame(Lon=taxdat_r$Lon.2,Lat=taxdat_r$Lat.2)

disty1<-geodist(x=cord1,y=cord2,measure="haversine",paired = TRUE)

taxdat_r$Distance<-disty1/100
phydat_r$Distance<-disty1/100

write.csv(taxdat_r,"Analyses/Master/Input/pairwise_remnant_taxo.csv")
write.csv(phydat_r,"Analyses/Master/Input/pairwise_remnant_phylo.csv")

taxdat_io<-left_join(taxdat_io,coords)
taxdat_io<-left_join(taxdat_io,coords2)

phydat_io<-left_join(phydat_io,coords)
phydat_io<-left_join(phydat_io,coords2)

cord1<-data.frame(Lon=taxdat_io$Lon.1,Lat=taxdat_io$Lat.1)
cord2<-data.frame(Lon=taxdat_io$Lon.2,Lat=taxdat_io$Lat.2)

disty1<-geodist(x=cord1,y=cord2,measure="haversine",paired = TRUE)

taxdat_io$Distance<-disty1/100
phydat_io$Distance<-disty1/100

write.csv(taxdat_io,"Analyses/Master/Input/pairwise_invaderonly_taxo.csv")
write.csv(phydat_io,"Analyses/Master/Input/pairwise_invaderonly_phylo.csv")

