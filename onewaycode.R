rm(list=ls())
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

graphics.off()

library(dplyr)
library(vegan)
library("geodist")
library(ggplot2)
library(hillR)
library(betapart)
library(reshape2)


setwd("~/Desktop/Powell/")# for mac
set.seed(3)

###read in data
d<-read.csv("Data/FULLDatabase_05272022.csv")
d.env<-read.csv("Data/FULLDatabase_05272022_plotsenv4Jul2022.csv")
d.ref<-read.csv("Data/FullDatabase_diversity_July4.csv")

###how many plots are in each ecoregion
counters<- d.env %>% group_by(US_L3NAME) %>% count() %>% arrange(n)
### choose 3 forest and 3 grassland desert be for now, one

d.env<-filter(d.env, US_L3NAME %in% c("Northern Basin and Range"))
d<-filter(d, Plot %in% unique(d.env$Plot))

d<-filter(d, Plot %in% unique(d.env$Plot))

table(d$Resampled)

d<- d%>% group_by(Plot) %>% filter(Year==max(Year)) %>% ungroup()

d.ref2<-filter(d.ref,Plot %in% c(unique(d$Plot)))

d.ref2<-d.ref2 %>% group_by(Plot) %>% filter(Year==max(Year)) %>% ungroup()

table(d.ref2$Dataset)
#### filter to one dataset
d<-filter(d,Dataset=="AIM")
d.ref2<-filter(d.ref2,Dataset=="AIM")
d.env<-filter(d.env,Dataset=="AIM")
######now you have one dataset per ecoregion
d.ref2<-select(d.ref2,Plot,RelCov_I,Long,Lat)

distm<-select(d.ref2,Long,Lat)
colnames(distm)<-c("lon","lat")


distn<-c(-117.7458, 42.02338) #for nothern basin and range

#distn<-c(-76.64072,34.68239) ## for soutehr atlantic coastal plain
#distn<-c(-106.3050, 35.72900) #Souterh Rockies
#distn<-c(-107.9667, 43.77188) # Wyoming
#distn<-c(-105.2695,46.39615)# NW great plains
#distn<-c(-114.2226, 36.59163) #mojave
disty<-geodist(distn , distm,measure="haversine")
disty<-t(disty)
d.ref2<-cbind(d.ref2,disty)

d.ref3<-filter(d.ref2,disty<30000) #30 or 50  km dpending on how many point

d<-filter(d,Plot %in% c(d.ref3$Plot))
d.native<-filter(d,NativeStatus!="I")



###calcularte disim
dater.high<-dplyr::select(d,Plot,AcceptedTaxonName,PctCov_100) ## rid of complilcating columns
dater.high<-dplyr::filter(dater.high,!is.na(dater.high$AcceptedTaxonName)) ## spreading the data to long for (matrix) won't work if there are NA's in taxon name so get rid of them here
range(dater.high$PctCov_100) 
high<- tidyr::spread(dater.high,AcceptedTaxonName,PctCov_100)  ## spread to long       
high[is.na(high)] <- 0 ## make na's zeros
high<-as.data.frame(high)
rownames(high)<-high$Plot # convert plot names to row names
high<-dplyr::select(high,-Plot) # remove plot column
high<-as.matrix.data.frame(high) # make it a matrix

high.hill<-hill_taxa_parti_pairwise(high,q=2,output = "matrix")


betapart::beta.multi.abund()

######### now native 
dater.native<-dplyr::select(d.native,Plot,AcceptedTaxonName,PctCov_100) ## rid of complilcating columns
dater.native<-dplyr::filter(dater.native,!is.na(dater.native$AcceptedTaxonName)) ## spreading the data to long for (matrix) won't work if there are NA's in taxon name so get rid of them here
range(dater.native$PctCov_100) 
native<- tidyr::spread(dater.native,AcceptedTaxonName,PctCov_100)  ## spread to long       
native[is.na(native)] <- 0 ## make na's zeros
native<-as.data.frame(native)
rownames(native)<-native$Plot # convert plot names to row names
native<-dplyr::select(native,-Plot) # remove plot column
native<-as.matrix.data.frame(native) # make it a matrix

native.hill<-hill_taxa_parti_pairwise(native,q=2,output= "matrix")


###  
dater2<-dplyr::select(d,Plot,Original.Long,Original.Lat) ## subset to lat longs
dater2<-distinct(dater2) ## remove dplicates casue each species was a row
remove<-setdiff(dater2$Plot,rownames(high))
dater2<-filter(dater2,!Plot %in% c(remove))
colnames(dater2)<-c("Plot","Long","Lat")
b<-geodist(dater2,paired = TRUE,measure = "haversine")
max(b)
dater2.native<-dplyr::select(d.native,Plot,Original.Long,Original.Lat) ## subset to lat longs
dater2.native<-distinct(dater2.native) ## remove dplicates casue each species was a row
remove.native<-setdiff(dater2.native$Plot,rownames(native))
dater2.native<-filter(dater2.native,!Plot %in% c(remove.native))
colnames(dater2.native)<-c("Plot","Long","Lat")
b.native<-geodist(dater2.native,paired = TRUE,measure = "haversine")


mydat<-data.frame(local_sim=numeric(),Beta=numeric(),regional_sim=numeric(), Distance=numeric(),RowID=character(),ColumnID=character())

for(i in 1:nrow(b)) {
  for(j in 1:ncol(b)) {
    if(i<j) { #only include the lower diagonal of the distance matrix
      mydat<-rbind(mydat,data.frame(
        Distance=b[i,j],# distance
        local_sim=high.hill$local_similarity[i,j],
        regional_sim=high.hill$region_similarity[i,j],
        Beta=high.hill$TD_beta[i,j],
        RowID=as.character(i), # row ID
        ColumnID=as.character(j))) # column ID
    }
  }
}

mydatnative<-data.frame(local_sim=numeric(),Beta=numeric(),regional_sim=numeric(), Distance=numeric(),RowID=character(),ColumnID=character())

for(i in 1:nrow(b)) {
  for(j in 1:ncol(b)) {
    if(i<j) { #only include the lower diagonal of the distance matrix
      mydatnative<-rbind(mydatnative,data.frame(
        Distance=b.native[i,j],# distance
        local_sim=native.hill$local_similarity[i,j],
        regional_sim=native.hill$region_similarity[i,j],
        Beta=native.hill$TD_beta[i,j],
        RowID=as.character(i), # row ID
        ColumnID=as.character(j))) # column ID
    }
  }
}


mydat$invasion<-"invaded"
mydatnative$invasion<-"uninvaded"




mydat2<-rbind(mydat,mydatnative)

ggplot(mydat2,aes(Distance,local_sim))+geom_smooth(aes(color=invasion))
  
ggplot(mydat2,aes(Distance,local_sim))+  geom_smooth(method="glm",method.args = list(family = "quasibinomial"),aes(color=invasion),size=1.2)


summariz<-high.hill %>%group_by(site1) %>%summarise(mean_sim=mean(local_similarity))
colnames(summariz)[1]<-c("Plot")
goober<-left_join(d.ref3,summariz)
goober2<-filter(goober,RelCov_I<50)

stop()
jpeg(filename = "~/Documents/git/bioticHogs/plots/MOHJAVEoneway.jpeg")
ggplot(goober2,aes(RelCov_I,mean_sim))+geom_point()+geom_smooth(method="lm",se = FALSE)+ggthemes::theme_few()+
  ggtitle("Mojave")
dev.off()

cor(goober$RelCov_I,goober$mean_sim,use = "complete.obs")
