##play with indicies to chatacterize 
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

setwd("~Desktop/Powell/")
set.seed(3)

d<-read.csv("Data/FULLDatabase_05272022.csv")
d.env<-read.csv("Data/FULLDatabase_05272022_plotsenv4Jul2022.csv")
d.ref<-read.csv("Data/FullDatabase_diversity_July4.csv")

bos<-dplyr::filter(d.env,US_L3NAME== "Driftless Area")
dat.bos<-filter(d, Plot %in% c(unique(bos$Plot)))
dat.bos<-filter(dat.bos, Resampled=="N")
driftless.ref<-filter(d.ref,Plot %in% c(unique(dat.bos$Plot)))
#rownames(driftless.ref)<-driftless.ref$Plot
driftless.ref<-select(driftless.ref,Plot,RelCov_I)

#reldiff<-abs(apply(combn(driftless.ref$RelCov_I,2), 2, diff))/abs(apply(combn(driftless.ref$RelCov_I,2), 2, mean))




dater.high<-dplyr::select(dat.bos,Plot,AcceptedTaxonName,PctCov) ## rid of complilcating columns
dater.high<-dplyr::filter(dater.high,!is.na(dater.high$AcceptedTaxonName)) ## spreading the data to long for (matrix) won't work if there are NA's in taxon name so get rid of them here
range(dater.high$PctCov) #high 122
high<- tidyr::spread(dater.high,AcceptedTaxonName,PctCov)  ## spread to long       
high[is.na(high)] <- 0 ## make na's zeros
rownames(high)<-high$Plot # convert plot names to row names
high<-dplyr::select(high,-Plot) # remove plot column
high<-as.matrix.data.frame(high) # make it a matrix

high.hill<-hill_taxa_parti_pairwise(high,q=2,output = "matrix")

dim(high.hill$local_similarity)

dater2<-dplyr::select(dat.bos,Plot,Long,Lat) ## subset to lat longs
dater2<-distinct(dater2) ## remove dplicates casue each species was a row

#rownames(dater2.high)<-dater2.highPlot # convert plot names to row names
#dater2<-dplyr::select(dater2,-Plot)
setdiff(dater2$Plot,rownames(high))

dater2<-filter(dater2,Plot!="MISS_111")

b<-geodist(dater2,paired = TRUE,measure = "haversine")
dim(b)

inv.data<-data.frame(Plot=rep(rownames(high),nrow(high)))
inv.data<-left_join(inv.data,driftless.ref)
inv.data2<-data.frame(Plot=rep(rownames(high),each=nrow(high)))
inv.data2<-left_join(inv.data2,driftless.ref)

inv.data<-cbind(inv.data,inv.data2)
colnames(inv.data)<-c("site1", "cov1","site2","cov2")
inv.data$reldiff<-(abs(inv.data$cov1-inv.data$cov2))/(abs(inv.data$cov1+inv.data$cov2)/2)

goo<-inv.data %>% mutate_at(vars(reldiff), ~replace(., is.nan(.), 0))

goo<-select(goo,site1,site2,reldiff)


goo<-acast(goo,site1~site2,value.var="reldiff")

dim(goo)


mydat <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(mydat) <- c('Distance','Sim','Reldiff','RowID', 'ColumnID')
for(i in 1:nrow(b)) {
  for(j in 1:ncol(b)) {
    if(i<j) { #only include the lower diagonal of the distance matrix
      mydat<-rbind(mydat,data.frame(
        Distance=b[i,j],# distance
        Sim=high.hill$local_similarity[i,j],
        Reldiff=goo[i,j],
        RowID=as.character(i), # row ID
        ColumnID=as.character(j))) # column ID
    }
  }
}


ggplot(mydat,aes(Simround))+geom_histogram()
ggplot(mydat,aes(Distance))+geom_histogram()

mod1<-glm(Sim~Distance,data=mydat,family=gaussian(link=log),start=c(0,0))

summary(mod1)
summary(mod1)



library(brms)
get_prior(Sim~Distance,
          family=zero_inflated_beta_binomial(),
          data=mydat)


priorz<-c(prior("student_t(3, 0, 1.5", "b"),
          prior("uniform(.5,1)","Intercept"))
          

mydat$Sim.adj<-ifelse(mydat$Sim==0,0.000000000000001,mydat$Sim)
mydat$Simround<-as.integer(round(mydat$Sim*100))

ggplot(mydat,aes(Distance,Simround))+geom_point()
table(mydat$Simround)

mod1a<-brm(Simround~Distance,
              family=zero_inflated_negbinomial(),init=0,
              data=mydat,warmup=3000,iter=4000,chains=2)
summary(mod1a)
plot(mod1a)
pp_check(mod1a,ndraws = 100)

+(1|mm(RowID,ColumnID))







