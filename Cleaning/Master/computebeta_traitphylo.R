##calculate beta diversity
###get ready
rm(list=ls())
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())
setwd("/Users/danielbuonaiuto/Documents/git/bioticHogs")

d<-read.csv("Analyses/Master/Input/specis_wtraitphylo.csv")
d.refL3<-read.csv("Analyses/Master/Input/balancedL3_guide.csv")
phy<-read.tree("Analyses/Master/Input/FINAL_PHYLOGENY_SpCodes_Apr2023.tr")


##get trait data alone

d.trait<-select(d,SpCode,heightveg_m,SLA_mm2.mg)
d.trait<-distinct(d.trait)
d.trait<-column_to_rownames(d.trait, var = "SpCode")
trait.height<-dplyr::select(d.trait,heightveg_m)
trait.SLA<-dplyr::select(d.trait,-heightveg_m)


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

d.cov.p<-filter(d.refL3,Status==0)
d.cov.i<-filter(d.refL3,Status!=0)

###add distance between plots
coords<-dplyr::select(d,Plot,Lat,Long)
coords<-distinct(coords)
colnames(coords)<-c("site1", "Lat.1", "Lon.1")
coords2<-coords
colnames(coords2)<-c("site2", "Lat.2", "Lon.2")
library(geodist)


ecoregionz<-unique(d.refL3$NA_L3NAME) ##81
ecoregionz<-ecoregionz[21:81] #(skip 20 Arizona new nexico mnts for mow)
ecoregionz<-ecoregionz[8:61] ## you will need to go back and fdi arizon new mexico
ecoregionz<-ecoregionz[2:54]
ecoregionz<-ecoregionz[11:53]
ecoregionz<-ecoregionz[2:43]
ecoregionz<-ecoregionz[18:42]
ecoregionz<-ecoregionz[2:25]
ecoregionz<-ecoregionz[7:24]
ecoregionz<-ecoregionz[6:18]
ecoregionz<-ecoregionz[2:13]
ecoregionz<-ecoregionz[3:13]
ecoregionz<-ecoregionz[10]
#73 ecoregions works
d.cov.i$NA_L3NAME<-gsub("/","_", d.cov.i$NA_L3NAME) #### try this again 
ecoregionz<-gsub("/","_", ecoregionz) 

setwd("Analyses/Master/Input/L3/invaded/")
library(hillR)
for (z in c(1:length(ecoregionz))){
  mydat<-data.frame()
  tlands<-dplyr::filter(d.cov.i,NA_L3NAME==ecoregionz[z])
  dat<-filter(d, Plot %in% c(unique(tlands$Plot)))
  comms<-matricize(dat)
  
  taxy0<-hill_taxa_parti_pairwise(comms,q = 0)
  taxy1<-hill_taxa_parti_pairwise(comms,q = 1)
  taxy2<-hill_taxa_parti_pairwise(comms,q = 2)
  taxy<-rbind(taxy0,taxy1,taxy2)
  taxy$class<-"taxonomic"
  
  phylo0<-hill_phylo_parti_pairwise(comms,phy,q = 0)
  phylo1<-hill_phylo_parti_pairwise(comms,phy,q = 1)
  phylo2<-hill_phylo_parti_pairwise(comms,phy,q = 2)
  phylo<-rbind(phylo0,phylo1,phylo2)
  phylo$class<-"phylogenetic"
  
  #t.h<-filter(trait.height,rownames(trait.height)%in% unique(dat$SpCode))
  
  height0<-hill_func_parti_pairwise(comms,traits =trait.height,q = 0)
  height1<-hill_func_parti_pairwise(comms,traits =trait.height,q = 1)
  height2<-hill_func_parti_pairwise(comms,traits =trait.height,q = 2)
  height<-rbind(height0,height1,height2)
  height$class<-"height"
  
  SLA0<-hill_func_parti_pairwise(comms,traits =trait.SLA ,q = 0)
  SLA1<-hill_func_parti_pairwise(comms,traits =trait.SLA ,q = 1)
  SLA2<-hill_func_parti_pairwise(comms,traits =trait.SLA ,q = 2)
  SLA<-rbind(SLA0,SLA1,SLA2)
  SLA$class<-"SLA"
  
  
  
  colnames(phylo)[4:6]<-c("gamma","alpha","beta")
  colnames(taxy)[4:6]<-c("gamma","alpha","beta")
  colnames(height)[4:6]<-c("gamma","alpha","beta")
  colnames(SLA)[4:6]<-c("gamma","alpha","beta")
  
  dathere<-rbind(phylo,taxy,height,SLA)
  dathere$NA_L3NAME<-ecoregionz[z]
  mydat<-dathere
  mydat<-left_join(mydat,coords)
  mydat<-left_join(mydat,coords2)
  
  cord1<-data.frame(Lon=mydat$Lon.1,Lat=mydat$Lat.1)
  cord2<-data.frame(Lon=mydat$Lon.2,Lat=mydat$Lat.2)
  
  
  disty1<-geodist(x=cord1,y=cord2,measure="haversine",paired = TRUE)
  head(disty1)
  mydat$Distance<-disty1/100   
  mydat$Distance<-ifelse(mydat$Distance==0,runif(length(mydat$Distance==0),0,1),mydat$Distance)
  mydat$logDist<-log(mydat$Distance)
  
  write.csv(mydat, paste0(ecoregionz[z], "invaded",".csv"))
  
  rm(mydat,dathere,dat)
  print(ecoregionz[z])
}

####now do uninvaded
### read back in datalist to make ecoregionz the same
temp<- list.files(pattern="*.csv")

library(data.table)
mydat<- rbindlist(sapply(temp, fread,simplify = FALSE))##load data 60372636 row

ecoregionz<-unique(mydat$NA_L3NAME)
rm(mydat)

setwd("/Users/danielbuonaiuto/Documents/git/bioticHogs/Analyses/Master/Input/L3/uninvaded")
d.cov.p$NA_L3NAME<-gsub("/","_", d.cov.p$NA_L3NAME)


for (z in c(1:length(ecoregionz))){
  mydat<-data.frame()
  tlands<-dplyr::filter(d.cov.p,NA_L3NAME==ecoregionz[z])
  dat<-filter(d, Plot %in% c(unique(tlands$Plot)))
  comms<-matricize(dat)
  
  taxy0<-hill_taxa_parti_pairwise(comms,q = 0)
  taxy1<-hill_taxa_parti_pairwise(comms,q = 1)
  taxy2<-hill_taxa_parti_pairwise(comms,q = 2)
  taxy<-rbind(taxy0,taxy1,taxy2)
  taxy$class<-"taxonomic"
  
  phylo0<-hill_phylo_parti_pairwise(comms,phy,q = 0)
  phylo1<-hill_phylo_parti_pairwise(comms,phy,q = 1)
  phylo2<-hill_phylo_parti_pairwise(comms,phy,q = 2)
  phylo<-rbind(phylo0,phylo1,phylo2)
  phylo$class<-"phylogenetic"
  
  #t.h<-filter(trait.height,rownames(trait.height)%in% unique(dat$SpCode))
  
  height0<-hill_func_parti_pairwise(comms,traits =trait.height,q = 0)
  height1<-hill_func_parti_pairwise(comms,traits =trait.height,q = 1)
  height2<-hill_func_parti_pairwise(comms,traits =trait.height,q = 2)
  height<-rbind(height0,height1,height2)
  height$class<-"height"
  
  SLA0<-hill_func_parti_pairwise(comms,traits =trait.SLA ,q = 0)
  SLA1<-hill_func_parti_pairwise(comms,traits =trait.SLA ,q = 1)
  SLA2<-hill_func_parti_pairwise(comms,traits =trait.SLA ,q = 2)
  SLA<-rbind(SLA0,SLA1,SLA2)
  SLA$class<-"SLA"
  
  
  
  colnames(phylo)[4:6]<-c("gamma","alpha","beta")
  colnames(taxy)[4:6]<-c("gamma","alpha","beta")
  colnames(height)[4:6]<-c("gamma","alpha","beta")
  colnames(SLA)[4:6]<-c("gamma","alpha","beta")
  
  dathere<-rbind(phylo,taxy,height,SLA)
  dathere$NA_L3NAME<-ecoregionz[z]
  mydat<-dathere
  mydat<-left_join(mydat,coords)
  mydat<-left_join(mydat,coords2)
  
  cord1<-data.frame(Lon=mydat$Lon.1,Lat=mydat$Lat.1)
  cord2<-data.frame(Lon=mydat$Lon.2,Lat=mydat$Lat.2)
  
  
  disty1<-geodist(x=cord1,y=cord2,measure="haversine",paired = TRUE)
  head(disty1)
  mydat$Distance<-disty1/100   
  mydat$Distance<-ifelse(mydat$Distance==0,runif(length(mydat$Distance==0),0,1),mydat$Distance)
  mydat$logDist<-log(mydat$Distance)
  
  write.csv(mydat, paste0(ecoregionz[z], ".uninvaded",".csv"))
  
  rm(mydat,dathere,dat)
  print(ecoregionz[z])
}



