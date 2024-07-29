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
### atempted plot matching
##data prep

##3read in plot summeries
d<-read.csv("Data/FULLDatabase_10272022_diversity_PctCov_100_10272022_INHABIT.csv")
d2<-dplyr::select(d,gHM.tif,Dataset,Plot,Year,RelCov_I,NDMI_median_Apr1_Sept30_1985_2020_90m.tif,SG_n_tot_M_sl2_100m.tif,TMIN_Dec_Feb_1981_2018_mean.tif,VCF_percent_treeCover_2000_2016_mean_integer.tif,gHM.tif,L4_KEY)
write.csv(d2,"Analyses/AbnOcc/GEB/envplotdata.csv",row.names = FALSE)
d<-d %>%
  group_by(Plot) %>%
  filter(Year==max(Year)) %>%
  ungroup() 

d$status<-ifelse(d$RelCov_I>0,1,0)

d.n<-filter(d,status==0)
d.i<-filter(d,status==1)
d.i<-filter(d.i,RelCov_I>5)
d<-rbind(d.i,d.n)


d$Richness_tot<-d$Richness_I+d$Richness_N+d$Richness_Unk


#qants<-d %>%group_by (L4_KEY) %>% summarise(Qs=round(quantile(RelCov_I,na.rm=TRUE,prob=c(0,.25)),1))
#qants$Q<-rep(c("Q25","Q75"),length(unique(d$L4_KEY)))
#qants<-tidyr::spread(qants,"Q","Qs")
#qants<-filter(qants,Q25!=Q75)

#d<-left_join(d,qants)
#d<-filter(d,!is.na(Q75))
#d$status<-NA
#d$status[which(d$RelCov_I<=d$Q25)]<- 0
#d$status[which(d$RelCov_I>=d$Q75)]<- 1

d2<-dplyr::select(d,Site,Lat,Long,gHM.tif,Dataset,Plot,Year,status,Richness_tot,NDMI_median_Apr1_Sept30_1985_2020_90m.tif,SG_n_tot_M_sl2_100m.tif,TMIN_Dec_Feb_1981_2018_mean.tif,VCF_percent_treeCover_2000_2016_mean_integer.tif,gHM.tif,L4_KEY)
#d2<-dplyr::select(d,Site,Lat,Long,Dataset,Plot,Year,status,Richness_tot,NDMI_median_Apr1_Sept30_1985_2020_90m.tif,TMIN_Dec_Feb_1981_2018_mean.tif,gHM.tif,L4_KEY)
d2<- d2[complete.cases(d2), ] ##55809

counters<-d2 %>% group_by(L4_KEY)%>%count()
counterss<-filter(counters,n>1)
d2<-filter(d2,L4_KEY %in% counterss$L4_KEY)

set.seed(111)
m2<-matchit(status~L4_KEY+Dataset+
                  NDMI_median_Apr1_Sept30_1985_2020_90m.tif+
                  TMIN_Dec_Feb_1981_2018_mean.tif+
                  SG_n_tot_M_sl2_100m.tif+
                  VCF_percent_treeCover_2000_2016_mean_integer.tif+
                  gHM.tif,data=d2,method = "nearest",distance="glm",exact=~L4_KEY+Dataset)





m3<-matchit(status~L4_KEY+Dataset+
            Lat+  
            Long,data=d2,method = "nearest",distance="glm",exact=~L4_KEY+Dataset)

mdat2<-match.data(m2)###20900 plots
mdat3<-match.data(m3)### plots

nater<-filter(mdat2,status==0)
iner<-filter(mdat2,status==1)
iner<-arrange(iner,subclass)
nater<-arrange(nater,subclass)

head(nater$subclass)
head(iner$subclass)
head(iner)

cordziner<-data.frame(Lon=iner$Long,Lat=iner$Lat)
cordznater<-data.frame(Lon=nater$Long,Lat=nater$Lat)
disty1<-geodist(x=cordziner,y=cordznater,measure="haversine",paired = TRUE)
iner$geodist<-disty1
iner$geodist<-iner$geodist/1000


inerclose<-filter(iner,geodist<10000)

quantile(iner$geodist)
?matchit()
pdf("Analyses/AbnOcc/propplot.pdf")
plot(m2, type = "density", interactive = FALSE,
     which.xs = ~NDMI_median_Apr1_Sept30_1985_2020_90m.tif+ TMIN_Dec_Feb_1981_2018_mean.tif + SG_n_tot_M_sl2_100m.tif+VCF_percent_treeCover_2000_2016_mean_integer.tif+ gHM.tif)
dev.off()

envpair<-dplyr::select(mdat2,Plot,L4_KEY,status,subclass)
write.csv(envpair,"Analyses/AbnOcc/Input/env_matched_plots.csv")

library(reshape2)

pairz.data<-data.frame()
ecoregionz<-sort(unique(counterss$L4_KEY))
for (z in c(1:length(ecoregionz))){
  goo<-filter(d2,L4_KEY==ecoregionz[z]) 
  cordz<-data.frame(Lon=goo$Long,Lat=goo$Lat)
  disty1<-geodist(cordz,measure="haversine")
  rownames(disty1)<-goo$Plot
  colnames(disty1)<-goo$Plot
  diss<-melt(disty1)
  diss<-filter(diss,Var1!=Var2)
  
  
  
  colnames(diss)<-c("plot1","plot2","distance")
  
  #diss<-filter(diss,plot.native %in% d.n$Plot) ### fitler this column to only native
  #diss<-filter(diss,plot.invaded %in% d.i$Plot) ### filter this column to only invaded
  
  
  diss$Distance<-diss$distance/1000
  
  
  #check<-diss %>% dplyr::group_by(plot.native)%>%slice(which.min(Distance))
  
  diss$L4_KEY<-ecoregionz[z]
  pairz.data<-rbind(diss,pairz.data)
  print(paste("done",which(ecoregionz==ecoregionz[z]))) 
}

pairz<-filter(pairz.data,plot1 %in% d.n$Plot)
pairz<-filter(pairz,plot2 %in% d.i$Plot)

pairz2<-pairz

pairz<-pairz %>% dplyr::group_by(plot1) %>%slice(which.min(Distance)) %>%ungroup()
pairz<-pairz %>% dplyr::group_by(plot2) %>%slice(which.min(Distance)) %>%ungroup()

pairz<-tibble::rownames_to_column(pairz,"match2")
pairz$match<-paste("match",pairz$match2,sep="_")

pairz.t<-dplyr::select(pairz,match,L4_KEY,plot1)
pairz.t$status<-0
pairz.t2<-dplyr::select(pairz,match,L4_KEY,plot2)
pairz.t2$status<-1

colnames(pairz.t)[3]<-"Plot"
colnames(pairz.t2)[3]<-"Plot"

distpair<-rbind(pairz.t,pairz.t2)
nrow(distpair)/2
nrow(envpair)/2

pairzclose<-filter(pairz2,Distance<10)
plots<-c(unique(pairzclose$plot1),unique(pairzclose$plot2))

d4<-filter(d2,Plot %in% plots)


set.seed(111)
m.short<-matchit(status~L4_KEY+Dataset+
              NDMI_median_Apr1_Sept30_1985_2020_90m.tif+
              TMIN_Dec_Feb_1981_2018_mean.tif+
              SG_n_tot_M_sl2_100m.tif+
              VCF_percent_treeCover_2000_2016_mean_integer.tif+
              gHM.tif,data=d4,method = "nearest",distance="glm",exact=~L4_KEY+Dataset)



mdatshort<-match.data(m.short)###9152 plots


nater<-filter(mdatshort,status==0)
iner<-filter(mdatshort,status==1)

iner<-arrange(iner,subclass)
nater<-arrange(nater,subclass)

head(nater$subclass)
head(iner$subclass)
head(iner)

cordziner<-data.frame(Lon=iner$Long,Lat=iner$Lat)
cordznater<-data.frame(Lon=nater$Long,Lat=nater$Lat)
disty1<-geodist(x=cordziner,y=cordznater,measure="haversine",paired = TRUE)
iner$geodist<-disty1
iner$geodist<-iner$geodist/1000










