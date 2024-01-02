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


setwd("~/PowellCenter")
set.seed(3)

d<-read.csv("Data/FULLDatabase_05272022.csv")
d.env<-read.csv("Data/FULLDatabase_05272022_plotsenv4Jul2022.csv")
d.ref<-read.csv("Data/FullDatabase_diversity_July4.csv")


counters<- d.env %>% group_by(US_L3NAME) %>% count() %>% arrange(n)
#twohundo<-filter(counters,n<=200)
#twohundo<-filter(twohundo,n>=100)


d.env<-filter(d.env, US_L3NAME %in% c("Mojave Basin and Range"))
d<-filter(d, Plot %in% unique(d.env$Plot))

table(d$Resampled)

d<- d%>% group_by(Plot) %>% filter(Year==max(Year)) %>% ungroup()

d.ref2<-filter(d.ref,Plot %in% c(unique(d$Plot)))

d.ref2<-d.ref2 %>% group_by(Plot) %>% filter(Year==max(Year)) %>% ungroup()
d.ref2<-select(d.ref2,Plot,RelCov_I)

table(d$Dataset)
d<-filter(d,Dataset=="AIM")
d.env<-filter(d.env,Dataset=="AIM")
ecoregionz<-unique(d.env$US_L3NAME)

d.native<-filter(d,NativeStatus!="I")
mydat<-data.frame(local_sim=numeric(),Beta=numeric(),regional_sim=numeric(), Distance=numeric(),Invasion=numeric(),Ecoregion=character(),RowID=character(),ColumnID=character())


###loop 1
for (z in c(1:length(ecoregionz))){
  
  tlands<-dplyr::filter(d.env,US_L3NAME==ecoregionz[z]) ##
  
  dat<-filter(d, Plot %in% c(unique(tlands$Plot)))
  
  dater.high<-dplyr::select(dat,Plot,AcceptedTaxonName,PctCov_100) ## rid of complilcating columns
  dater.high<-dplyr::filter(dater.high,!is.na(dater.high$AcceptedTaxonName)) ## spreading the data to long for (matrix) won't work if there are NA's in taxon name so get rid of them here
  range(dater.high$PctCov_100) #high 122
  high<- tidyr::spread(dater.high,AcceptedTaxonName,PctCov_100)  ## spread to long       
  high[is.na(high)] <- 0 ## make na's zeros
  high<-as.data.frame(high)
  rownames(high)<-high$Plot # convert plot names to row names
  high<-dplyr::select(high,-Plot) # remove plot column
  high<-as.matrix.data.frame(high) # make it a matrix
  
  high.hill<-hill_taxa_parti_pairwise(high,q=2,output = "matrix")
  
  print(dim(high.hill$local_similarity))
  
   ######p2

  dater2<-dplyr::select(dat,Plot,Original.Long,Original.Lat) ## subset to lat longs
  dater2<-distinct(dater2) ## remove dplicates casue each species was a row
  remove<-setdiff(dater2$Plot,rownames(high))
  dater2<-filter(dater2,!Plot %in% c(remove))
  colnames(dater2)<-c("Plot","Long","Lat")
  b<-geodist(dater2,paired = TRUE,measure = "haversine")
  print(dim(b))
  

  
  ###p3
  #inv.data<-data.frame(Plot=rep(rownames(high),nrow(high)))
  
  #dim(inv.data)
  
  
  #inv.data<-left_join(inv.data,d.ref2)
  #dim(inv.data)
  #inv.data2<-data.frame(Plot=rep(rownames(high),each=nrow(high)))
  #inv.data2<-left_join(inv.data2,d.ref2)
  #dim(inv.data2)
  
#  inv.data<-cbind(inv.data,inv.data2)
 # colnames(inv.data)<-c("site1", "cov1","site2","cov2")
  
  #inv.data$cov1<-inv.data$cov1/100+.1
  #inv.data$cov2<-inv.data$cov2/100+.1
  #range(inv.data$cov1)
  
  #inv.data$reldiff<-abs((inv.data$cov1-inv.data$cov2)/.5*((inv.data$cov1+inv.data$cov2)))

  
  #goo<-select(inv.data,site1,site2,reldiff)
  #goober<-group_by(goo, site1,site2) %>% tally %>% arrange(n)
  
  #goo<-acast(goo,site1~site2,value.var="reldiff")
  
  #print(dim(goo))
   
  for(i in 1:nrow(b)) {
    for(j in 1:ncol(b)) {
      if(i<j) { #only include the lower diagonal of the distance matrix
        mydat<-rbind(mydat,data.frame(
          Distance=b[i,j],# distance
          local_sim=high.hill$local_similarity[i,j],
          regional_sim=high.hill$region_similarity[i,j],
          Beta=high.hill$TD_beta[i,j],
          invasion=NA,
          Ecoregion=unique(tlands$US_L3NAME),
          RowID=as.character(i), # row ID
          ColumnID=as.character(j))) # column ID
      }
    }
  }
}


mydat$Inv.z<-(mydat$invasion-mean(mydat$invasion))/sd(mydat$invasion)
mydat$Dist.z<-(mydat$Distance-mean(mydat$Distance))/sd(mydat$Distance)
#mydat$Sim<-ifelse(mydat$local_sim==1,.999999,mydat$local_sim)
range(mydat$invasion)
2.4/3
mydat$inv.bin<-ifelse(mydat$invasion<=.8,"low","intermediate")
mydat$inv.bin<-ifelse(mydat$invasion>=1.6,"high",mydat$inv.bin)
write.csv(mydat,"..//git/bioticHogs/twohundo_daters.csv")
colnames(mydat)


#mydat$Distance<-mydat$Distance*10

ggplot(mydat,aes(Distance,local_sim))+geom_smooth(method="glm",method.args = list(family = "quasibinomial"),aes(color=inv.bin))+
  ggthemes::theme_few()+facet_wrap(~Ecoregion,scales="free_x")+scale_fill_viridis_d()

test<-brm(Sim~Dist.z*Inv.z+(1|Ecoregion),data=mydat,family=zero_inflated_beta())
test2<-brm(Sim~Dist.z*Inv.z+(1|Ecoregion)+(1|mm(RowID,ColumnID)),data=mydat,family=zero_inflated_beta())

test3<-brm(Sim~Dist.z*Inv.z+(1|Ecoregion/RowID)+(1|Ecoregion/ColumnID),data=mydat,family=zero_inflated_beta())
library(lme4)

test3<-glmer(Sim~Dist.z*Inv.z+(1|Ecoregion/RowID)+(1|Ecoregion/ColumnID),data=mydat,family =binomial(link = "logit") )

??glm
