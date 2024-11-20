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
load("Analyses/AbnOcc/spacematches.Rda")
### atempted plot matching
##data prep

##3read in plot summeries
d<-read.csv("Data/FULLDatabase_10272022_diversity_PctCov_100_10272022_INHABIT.csv")
#d2<-dplyr::select(d,gHM.tif,Dataset,Plot,Year,RelCov_I,NDMI_median_Apr1_Sept30_1985_2020_90m.tif,SG_n_tot_M_sl2_100m.tif,TMIN_Dec_Feb_1981_2018_mean.tif,VCF_percent_treeCover_2000_2016_mean_integer.tif,gHM.tif,L4_KEY)
#write.csv(d2,"Analyses/AbnOcc/GEB/envplotdata.csv",row.names = FALSE)
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
mdat3<-match.data(m3)###20900 plots

nater<-dplyr::filter(mdat2,status==0)
iner<-dplyr::filter(mdat2,status==1)
iner<-dplyr::arrange(iner,subclass)
nater<-dplyr::arrange(nater,subclass)

head(nater$subclass)
head(iner$subclass)
head(iner)

cordziner<-data.frame(Lon=iner$Long,Lat=iner$Lat)
cordznater<-data.frame(Lon=nater$Long,Lat=nater$Lat)
head(cordziner)
head(cordznater)
disty1<-geodist(x=cordznater,y=cordziner,measure="haversine",paired = TRUE)
iner$geodist<-disty1
iner$geodist<-iner$geodist/1000

nater$geodist<-disty1
nater$geodist<-nater$geodist/1000

median(iner$geodist)
round(mean(iner$geodist),1) #9545.9
round(sd(iner$geodist),1) # 5118.1
ggplot(iner,aes(geodist))+geom_histogram(bins=100)


inerclose<-filter(iner,geodist<9769.892)

naterclose<-filter(nater,subclass %in% c(inerclose$subclass))


inerclose<-dplyr::select(inerclose,Dataset,Plot,L4_KEY,status,subclass)
naterclose<-dplyr::select(naterclose,Dataset,Plot,L4_KEY,status,subclass)
#inerclose<-dplyr::select(inerclose,-geodist)
#naterclose<-dplyr::select(naterclose,-geodist)
#naterclose$subclass<-as.character(naterclose$subclass)
#inerclose$subclass<-as.character(inerclose$subclass)

library(data.table)

close<- rbindlist(list(naterclose,inerclose))
head(close)


quantile(iner$geodist)
?matchit()
pdf("Analyses/AbnOcc/propplot.pdf")
plot(m2, type = "density", interactive = FALSE,
     which.xs = ~NDMI_median_Apr1_Sept30_1985_2020_90m.tif+ TMIN_Dec_Feb_1981_2018_mean.tif + SG_n_tot_M_sl2_100m.tif+VCF_percent_treeCover_2000_2016_mean_integer.tif+ gHM.tif)
dev.off()

envpair<-dplyr::select(mdat2,Plot,Dataset,L4_KEY,status,subclass)
write.csv(envpair,"Analyses/AbnOcc/Input/env_matched_plots.csv")

write.csv(close,"Analyses/AbnOcc/Input/env_matched_plots_close.csv")



###minimum distance
library(reshape2)

d$unit<-paste(d$L4_KEY,d$Dataset,sep="_")
counters2<-d %>% group_by(L4_KEY,Dataset)%>%count()
counterss2<-filter(counters2,n>1)
counterss2$unit<-paste(counterss2$L4_KEY,counterss2$Dataset,sep="_")
d3<-filter(d,unit %in% counterss2$unit)


pairz.data<-data.frame()

ecoregionz<-sort(unique(counterss2$unit),decreasing = TRUE)
head(ecoregionz,10)
ecoregionz<-ecoregionz[7:1199]


for (z in c(1:length(ecoregionz))){
  goo<-filter(d3,unit==ecoregionz[z]) 
  cordz<-data.frame(Lon=goo$Long,Lat=goo$Lat)
  disty1<-geodist(cordz,measure="haversine")
  rownames(disty1)<-goo$Plot
  colnames(disty1)<-goo$Plot
  diss<-melt(disty1)
  diss<-filter(diss,Var1!=Var2)
  diss<-filter(diss,Var2!=Var1)
  colnames(diss)<-c("plot1","plot2","distance")
  diss<-filter(diss,plot1 %in% unique(d.n$Plot))###make plot1 native
  diss<-filter(diss,plot2 %in% unique(d.i$Plot)) 
  diss<-diss %>% dplyr::group_by(plot2) %>% dplyr::slice(which.min(distance)) %>%ungroup() ###
  #diss<-filter(diss,plot1 %in% d.n$Plot) ### fitler this column to only native
  #diss<-filter(diss,plot2%in% d.i$Plot) ### filter this column to only invaded
  
  
  diss$Distance<-diss$distance/1000
  
  
  #check<-diss %>% dplyr::group_by(plot.native)%>%slice(which.min(Distance))
  
  diss$unit<-ecoregionz[z]
  #diss<-diss %>% dplyr::group_by(plot2) %>% dplyr::slice(which.min(Distance)) %>%ungroup()
  #isna<-((grepl("NA_",diss$unit)))
  #diss<-cbind(diss,isna)
  #pairz2<-filter(pairz2,isna==FALSE)

  
  pairz.data<-rbind(diss,pairz.data)
  print(paste("done",which(ecoregionz==ecoregionz[z]))) 
}

mean(pairz.data$Distance) #2006.km

write.csv(pairz.data,"spatialclose_all_plots.csv")


#so pairz in now a list of all native plots and their distance to invaded ones




pairz2<-pairz.data %>% dplyr::filter(Distance<200) #se the min distance to every invaded
mean(pairz2$Distance)
median(pairz2$Distance)
sd(pairz2$Distance)


ggplot(pairz2,aes(Distance))+geom_histogram(bins=100)+xlim(0,200)





#pairz<-pairz %>% dplyr::group_by(plot2) %>%slice(which.min(Distance)) %>%ungroup()

pairz2<-tibble::rownames_to_column(pairz2,"match2")
pairz2$match<-paste("match",pairz2$match2,sep="_")


pairzclose<-pairz2
median(pairzclose$Distance)

pairz.t<-dplyr::select(pairzclose,match,unit,plot1)
pairz.t$status<-0
pairz.t2<-dplyr::select(pairzclose,match,unit,plot2)
pairz.t2$status<-1

colnames(pairz.t)[3]<-"Plot"
colnames(pairz.t2)[3]<-"Plot"

distpair<-rbind(pairz.t,pairz.t2)
nrow(distpair)/2



save.image("Analyses/AbnOcc/spacematches.Rda")
write.csv(distpair,"Analyses/AbnOcc/Input/spatial_matched_plots.csv")








stop()##################






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










