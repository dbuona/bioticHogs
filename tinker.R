####see highplainsdrifter.R for most uptodat version of these loops.

rm(list=ls())
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#graphics.off()

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
d1<-read.csv("Data/FULLDatabase_05272022.csv")
d.env1<-read.csv("Data/FULLDatabase_05272022_plotsenv4Jul2022.csv")
d.ref<-read.csv("Data/FullDatabase_diversity_July4.csv")

###how many plots are in each ecoregion
counters<- d.env1 %>% group_by(US_L3NAME) %>% count() %>% arrange(n)
twohundo<-filter(counters,n<600)
twohundo<-filter(twohundo,n>=300)

###choose you sandbox. High Plains is good. 334 Plots so 
d.env<-filter(d.env1, US_L3NAME %in% c(twohundo$US_L3NAME))
d<-filter(d1, Plot %in% unique(d.env$Plot)) ### filter the main dataset

#### deal with resampled plots
table(d$Resampled)

d<- d%>% group_by(Plot) %>% 
  filter(Year==max(Year)) %>%
  ungroup() ### this take resapled plots and chooses the most recent

d.ref2<-filter(d.ref,Plot %in% c(unique(d$Plot))) ### this filter the invasion levels data

d.ref2<-d.ref2 %>%
  group_by(Plot) %>%
  filter(Year==max(Year)) %>%
  ungroup() #  Same as above, only one on observation per plot
#################################

table(d.ref2$Dataset) ### here in theory you could select for a dataset though i dont

#######make abundace maxtrix

matricize<-function(x){
  temp<-dplyr::select(x,Plot,AcceptedTaxonName,PctCov_100) ### select useful columns
  temp<-dplyr::filter(temp,!is.na(temp$AcceptedTaxonName))
  temp<- tidyr::spread(temp,AcceptedTaxonName,PctCov_100)
  temp[is.na(temp)] <- 0 
  temp<-as.data.frame(temp)
  rownames(temp)<-temp$Plot # convert plot names to row names
  temp<-dplyr::select(temp,-Plot)
  temp<-as.matrix.data.frame(temp)# remove plot column
  }

distancize<-function(x,y){
  temp2<-dplyr::select(x,Plot,Original.Long,Original.Lat) ## subset to lat longs
  temp2<-dplyr::distinct(temp2) ## remove dplicates casue each species was a row
  remove<-setdiff(temp2$Plot,rownames(y))
  temp2<-filter(temp2,!Plot %in% c(remove))
  setNames(data.frame(temp2),c("Plot","Long","Lat"))
  
}
datacize<-function(x,y,z){
  test<-setNames(melt(as.matrix(x$beta.bray)),c("site1","site2", "beta_div"))
  test2<-setNames(melt(as.matrix(x$beta.bray.bal)),c("site1","site2", "turnover"))
  test3<-setNames(melt(as.matrix(x$beta.bray.gra)),c("site1","site2", "nestedness"))
  test<-left_join(test,test2)
  test<-left_join(test,test3)
  testD<-setNames(melt(as.matrix(y)),c("site1","site2", "Distance"))
  test$Distance<-testD$Distance
  test$metric<-"abundance"
  
  testo<-setNames(melt(as.matrix(z$beta.bray)),c("site1","site2", "beta_div"))
  testo2<-setNames(melt(as.matrix(z$beta.bray.bal)),c("site1","site2", "turnover"))
  testo3<-setNames(melt(as.matrix(z$beta.bray.gra)),c("site1","site2", "nestedness"))
  testo<-left_join(testo,testo2)
  testo<-left_join(testo,testo3)
  testo$Distance<-test$Distance
  testo$metric<-"occurance"
  
  test<-rbind(test,testo)
}

mydat<-data.frame(site1=factor(), site2=factor(), beta_div=factor(), turnover=numeric(), nestedness=numeric(),
                  distance=numeric() ,metric=character(),status=character(),status=character,Ecoregion=character())

###loop 1
ecoregionz<-unique(d.env$NA_L3NAME) ##one has no numeric data
for (z in c(1:length(ecoregionz))){
  
  tlands<-dplyr::filter(d.env,US_L3NAME==ecoregionz[z]) ##
  
  dat<-filter(d, Plot %in% c(unique(tlands$Plot)))
  comms<-matricize(dat)
  occs<-ifelse(comms==0,comms,1)
  comms.out.abun<-beta.pair.abund(comms)
  comms.out.occs<-beta.pair.abund(occs)
  disty<-distancize(dat,comms)
  b<-geodist(disty,paired = TRUE,measure = "haversine")
  daty<-datacize(comms.out.abun,b,comms.out.occs)
  daty$status<-"introduced species included"
  daty$Ecoregion=ecoregionz[z]
  mydat<-rbind(daty,mydat)
  print(paste("done",which(ecoregionz==ecoregionz[z])))    
}

###now only natives
d.native<-filter(d,NativeStatus!="I")

mydat2<-data.frame(site1=factor(), site2=factor(), beta_div=factor(), turnover=numeric(), nestedness=numeric(),
                  distance=numeric() ,metric=character(),status=character(),status=character,Ecoregion=character())

###loop 2

for (z in c(1:length(ecoregionz))){
  
  tlands<-dplyr::filter(d.env,US_L3NAME==ecoregionz[z]) ##
  
  dat<-filter(d.native, Plot %in% c(unique(tlands$Plot)))
  comms<-matricize(dat)
  occs<-ifelse(comms==0,comms,1)
  comms.out.abun<-beta.pair.abund(comms)
  comms.out.occs<-beta.pair.abund(occs)
  disty<-distancize(dat,comms)
  b<-geodist(disty,paired = TRUE,measure = "haversine")
  daty<-datacize(comms.out.abun,b,comms.out.occs)
  daty$status<-"native species only"
  daty$Ecoregion=ecoregionz[z]
  mydat2<-rbind(daty,mydat2)
  print(paste("done",which(ecoregionz==ecoregionz[z])))    
}

mydatPlot<-rbind(mydat,mydat2)
mydatPlot$Distance<-mydatPlot$Distance/100

write.csv(mydatPlot,"compair300.csv")
save.image("compair300.Rda")

onetoone<-select(mydat,site1,site2,beta_div,status,Ecoregion)
onetoone<-tidyr::spread(onetoone,beta_div,status, -site1,-site2)


ggplot(mydatPlot,aes(Distance,beta_div))+geom_smooth(aes(color=status,fill=status,linetype=metric),method="glm",method.args = list(family = "quasibinomial"),size=1.2)+
  facet_wrap(~Ecoregion,scales="free")+ggthemes::theme_few(base_size = 11)+ylab("Bray-Curtis Dissimilarity")+xlab("Distance (km)")+
  scale_color_viridis_d(begin = .1,end=.7)+scale_fill_viridis_d(begin = .1,end=.7)


ggplot(mydatPlot,aes(Distance,turnover))+geom_smooth(aes(color=status,fill=status,linetype=metric),method="glm",method.args = list(family = "quasibinomial"),size=1.2)+
  facet_wrap(~Ecoregion,scales="free")+ggthemes::theme_few(base_size = 11)+ylab("Bray-Curtis Dissimilarity")+xlab("Distance (km)")+
  scale_color_viridis_d(begin = .1,end=.7)+scale_fill_viridis_d(begin = .1,end=.7)

ggplot(mydatPlot,aes(Distance,nestedness))+geom_smooth(aes(color=status,fill=status,linetype=metric),method="glm",method.args = list(family = "quasibinomial"),size=1.2)+
  facet_wrap(~Ecoregion,scales="free")+ggthemes::theme_few(base_size = 11)+ylab("Bray-Curtis Dissimilarity")+xlab("Distance (km)")+
  scale_color_viridis_d(begin = .1,end=.7)+scale_fill_viridis_d(begin = .1,end=.7)


BM<-filter(mydatPlot,Ecoregion=="High Plains")
ggplot(BM,aes(Distance,beta_div))+geom_smooth(aes(color=status,fill=status,linetype=metric),method="glm",method.args = list(family = "quasibinomial"),size=1.2)+
ggthemes::theme_few(base_size = 11)+ylab("Bray-Curtis Dissimilarity")+xlab("Distance (km)")+
  scale_color_viridis_d(begin = .1,end=.7)+scale_fill_viridis_d(begin = .1,end=.7)





BM$Dist.cent<-BM$Distance-mean(BM$Distance)
BM.occ<-filter(BM,metric=="occurance")
BM.occ$beta_div2<-ifelse(BM.occ$beta_div==1,.999999999999,BM.occ$beta_div)
library(brms)
hightest<-brm(beta_div2~Dist.cent+status,data=BM.occ,family=zero_inflated_beta(),chains=2,init=0)

stop()
##make input matrixes
abun<-matricize(d) 

occur<-ifelse(abun==0,abun,1)

#### Calculate dissimilarity using betapart package
abun.out<-betapart::beta.pair.abund(abun)
occur.out<-betapart::beta.pair.abund(occur)

disty<-distancize(d,abun) ###make spatial matrix
b<-geodist(disty,paired = TRUE,measure = "haversine") # make ditance matrix

### were going to put these into a data frame, so make sure they are the same dimensions
dim(b)
dim(as.matrix(abun.out$beta.bray))
dim(as.matrix(occur.out$beta.bray))

#########nice, Now build your datasett


daters<-datacize(abun.out,b,occur.out)


ggplot(daters,aes(Distance,beta_div))+geom_smooth(aes(color=metric,fill=metric),method="glm",method.args = list(family = "quasibinomial"),size=1.2)

#######
###now for plots without the invaders
d.native<-filter(d,NativeStatus!="I")

##make matrices
abun.native<-matricize(d.native)
occur.native<-ifelse(abun.native==0,abun.native,1)

##compute beta diversity
abun.native.out<-betapart::beta.pair.abund(abun.native)
occur.native.out<-betapart::beta.pair.abund(occur.native)

####check if new matrix matches previous created distance one
### this should occur of no plots are lost (e.g. no plots were 100% invaded)
dim(as.matrix(abun.native.out$beta.bray))
dim(b)

distyn<-distancize(d.native,abun.native)
bn<-geodist(distyn,paired = TRUE,measure = "haversine")
dim(bn)
####make you r data
daters.n<-datacize(abun.native.out,bn,occur.native.out)

daters$status<-"invaders and natives"
daters.n$status<-"natives only"

daters.plot<-rbind(daters,daters.n)
daters.plot$Distance<-daters.plot$Distance/100

jpeg("~/Documents/git/bioticHogs/plots/HIGHPLAINShomog.jpeg")
one<-ggplot(daters.plot,aes(Distance,beta_div))+geom_smooth(aes(color=status,fill=status),method="glm",method.args = list(family = "quasibinomial"),size=1.2)+
  facet_wrap(~metric)+ggthemes::theme_few(base_size = 11)+ylab("Bray-Curtis Dissimilarity")+xlab("Distance (km)")+
  scale_color_viridis_d(begin = .1,end=.7)+scale_fill_viridis_d(begin = .1,end=.7)
dev.off()

two<-ggplot(daters.plot,aes(Distance,nestedness))+geom_smooth(aes(color=status,fill=status),method="glm",method.args = list(family = "quasibinomial"),size=1.2)+
  facet_wrap(~metric)+ggthemes::theme_few(base_size = 11)+ylab("Bray-Curtis Dissimilarity")+xlab("Distance (km)")+
  scale_color_viridis_d(begin = .1,end=.7)+scale_fill_viridis_d(begin = .1,end=.7)


three<-ggplot(daters.plot,aes(Distance,turnover))+geom_smooth(aes(color=status,fill=status),method="glm",method.args = list(family = "quasibinomial"),size=1.2)+
  facet_wrap(~metric)+ggthemes::theme_few(base_size = 11)+ylab("Bray-Curtis Dissimilarity")+xlab("Distance (km)")+
  scale_color_viridis_d(begin = .1,end=.7)+scale_fill_viridis_d(begin = .1,end=.7)


jpeg("~/Documents/git/bioticHogs/plots/HIGHPLAINSparti.jpeg")
ggpubr::ggarrange(one,three,two,ncol=1,nrow=3, common.legend = TRUE, labels = c("Beta diversity","Turnover","Nestedness"))
dev.off()

###### part II Exploring functional Diversity
traits<-read.csv("Data/FULLDatabaseSppListTRY.csv")






stop()

mydat2<-data.frame(Distance=numeric(),beta_div_abun=numeric(),
                  turnover_abun=numeric(),nestedness_abun=numeric(),
                  beta_div_occ=numeric(),turnover_occ=numeric(),nestedness_occ=numeric(),
                  RowID=character(),ColumnID=character)

for(i in 1:nrow(b)) {
  for(j in 1:nrow(b)) {
    if(i<j) { #only include the lower diagonal of the distance matrix
      mydat2<-rbind(mydat2,data.frame(
        Distance=b[i,j],# distance
        beta_div_abun=as.matrix(abun.native.out$beta.bray)[i,j],
        turnover_abun=as.matrix(abun.native.out$beta.bray.bal)[i,j],
        nestedness_abun=as.matrix(abun.native.out$beta.bray.gra)[i,j],
        beta_div_occ=as.matrix(occur.native.out$beta.bray)[i,j],
        turnover_occ=as.matrix(occur.native.out$beta.bray.bal)[i,j],
        nestedness_occ=as.matrix(occur.native.out$beta.bray.gra)[i,j],
        RowID=as.character(i), # row ID
        ColumnID=as.character(j))) # column ID
    }
  }
}

mydat$status<-"Invasives"
mydat2$status<-"Natives only"

mydat3<-rbind(mydat,mydat2)
mydat3<-tidyr::gather()
?gather()
ggplot(mydat3,aes(Distance,beta_div_abun))+geom_smooth(aes(linetype=status),method="glm",method.args = list(family = "quasibinomial"),size=1.2)+
  geom_smooth(aes(Distance,beta_div_occ,linetype=status),method="glm",method.args = list(family = "quasibinomial"),size=1.2,color="red")+
  ggthemes::theme_few()

ggplot(mydat3,aes(Distance,nestedness_abun))+geom_smooth(aes(linetype=status),method="glm",method.args = list(family = "quasibinomial"),size=1.2)+
  geom_smooth(aes(Distance,nestedness_occ,linetype=status),method="glm",method.args = list(family = "quasibinomial"),size=1.2,color="red")+
  ggthemes::theme_few()

ggplot(mydat3,aes(Distance,turnover_abun))+geom_smooth(aes(linetype=status),method="glm",method.args = list(family = "quasibinomial"),size=1.2)+
  geom_smooth(aes(Distance,turnover_occ,linetype=status),method="glm",method.args = list(family = "quasibinomial"),size=1.2,color="red")+
  ggthemes::theme_few()


d.invasive<-filter(d,NativeStatus=="I")
#######make abundace maxtrix
dater.abun.invasive<-dplyr::select(d.invasive,Plot,AcceptedTaxonName,PctCov_100) ## rid of complilcating columns
dater.abun.invasive<-dplyr::filter(dater.abun.invasive,!is.na(dater.abun.invasive$AcceptedTaxonName)) ## spreading the data to long for (matrix) won't work if there are NA's in taxon name so get rid of them here
range(dater.abun.invasive$PctCov_100) #abun 122
abun.invasive<- tidyr::spread(dater.abun.invasive,AcceptedTaxonName,PctCov_100)  ## spread to long       
abun.invasive[is.na(abun.invasive)] <- 0 ## make na's zeros
abun.invasive<-as.data.frame(abun.invasive)
rownames(abun.invasive)<-abun.invasive$Plot # convert plot names to row names
abun.invasive<-dplyr::select(abun.invasive,-Plot) # remove plot column

####two data sets
abun.invasive<-as.matrix.data.frame(abun.invasive) # make it a matrix
occur.invasive<-ifelse(abun.invasive==0,abun,1)

abun.invasive.out<-betapart::beta.pair.abund(abun.invasive)
occur.invasive.out<-betapart::beta.pair.abund(occur.invasive)
dim(b)
dim(as.matrix(abun.invasive.out$beta.bray))


dater2.i<-dplyr::select(d.invasive,Plot,Original.Long,Original.Lat) ## subset to lat longs
dater2.i<-distinct(dater2.i) ## remove dplicates casue each species was a row
#remove<-setdiff(dater2$Plot,rownames(high))
#dater2<-filter(dater2,!Plot %in% c(remove))
colnames(dater2.i)<-c("Plot","Long","Lat")
b.i<-geodist(dater2.i,paired = TRUE,measure = "haversine")



mydat4<-data.frame(Distance=numeric(),beta_div_abun=numeric(),
                   turnover_abun=numeric(),nestedness_abun=numeric(),
                   beta_div_occ=numeric(),turnover_occ=numeric(),nestedness_occ=numeric(),
                   RowID=character(),ColumnID=character)

for(i in 1:nrow(b.i)) {
  for(j in 1:nrow(b.i)) {
    if(i<j) { #only include the lower diagonal of the distance matrix
      mydat4<-rbind(mydat4,data.frame(
        Distance=b.i[i,j],# distance
        beta_div_abun=as.matrix(abun.invasive.out$beta.bray)[i,j],
        turnover_abun=as.matrix(abun.invasive.out$beta.bray.bal)[i,j],
        nestedness_abun=as.matrix(abun.invasive.out$beta.bray.gra)[i,j],
        beta_div_occ=as.matrix(occur.invasive.out$beta.bray)[i,j],
        turnover_occ=as.matrix(occur.invasive.out$beta.bray.bal)[i,j],
        nestedness_occ=as.matrix(occur.invasive.out$beta.bray.gra)[i,j],
        RowID=as.character(i), # row ID
        ColumnID=as.character(j))) # column ID
    }
  }
}

mydat4$status<-"invasive only"

mydat5<-rbind(mydat4,mydat3)

ggplot(mydat5,aes(Distance,beta_div_abun))+geom_smooth(aes(color=status,fill=status),method="glm",method.args = list(family = "quasibinomial"),size=1.2)
+
  geom_smooth(aes(Distance,beta_div_occ,linetype=status),method="glm",method.args = list(family = "quasibinomial"),size=1.2,color="red")+
  ggthemes::theme_few()

clean4plot<-select(mydat5,Distance,beta_div_abun,beta_div_occ,status)
gimbo<-tidyr::gather(clean4plot,metric,measure,-Distance,-status)
gimbo2<-filter(gimbo,status!="invasive only")
ggplot(gimbo2,aes(Distance,measure))+geom_smooth(aes(color=metric,fill=metric,linetype=status),method="glm",method.args = list(family = "quasibinomial"),size=.5,alpha=0.6)
save.image("homoginizationbatholith.Rda")
