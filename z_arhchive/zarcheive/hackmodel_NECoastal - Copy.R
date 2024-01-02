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

setwd("~PowellCenter")
set.seed(3)

d<-read.csv("Data/FULLDatabase_05272022.csv")
d.env<-read.csv("Data/FULLDatabase_05272022_plotsenv4Jul2022.csv")
d.ref<-read.csv("Data/FullDatabase_diversity_July4.csv")

bos<-dplyr::filter(d.env,US_L3NAME== "Northeastern Coastal Zone")
dat.bos<-filter(d, Plot %in% c(unique(bos$Plot)))
dat.bos<-filter(dat.bos, Resampled=="N")
driftless.ref<-filter(d.ref,Plot %in% c(unique(dat.bos$Plot)))
#rownames(driftless.ref)<-driftless.ref$Plot
driftless.ref<-select(driftless.ref,Plot,RelCov_I)






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
inv.data$cov1<-inv.data$cov1/100
inv.data$cov2<-inv.data$cov2/100

inv.data$reldiff<-ifelse(inv.data$cov1==0,inv.data$cov2*2,(abs(inv.data$cov1-inv.data$cov2))/(abs(inv.data$cov1+inv.data$cov2)/2))

fixcov<-which(inv.data$cov2==0)
inv.data$reldiff[fixcov] <-inv.data$cov1[fixcov]*2



goo<-select(inv.data,site1,site2,reldiff)


goo<-acast(goo,site1~site2,value.var="reldiff")

goo %>% 
  as.data.frame() %>%
  rownames_to_column("cov1") %>%
  pivot_longer(-c(cov1), names_to = "cov2", values_to = "reldiff") %>%
  ggplot(aes(x=cov2, y=cov1, fill=reldiff)) + 
  geom_raster()+
  scale_fill_viridis_c()

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

mydat$Diff.z<-(mydat$Reldiff-mean(mydat$Reldiff))/sd(mydat$Reldiff)
mydat$Dist.z<-(mydat$Distance-mean(mydat$Distance))/sd(mydat$Distance)



mod1<-glm(Sim~Distance+Reldiff,data=mydat,family=quasibinomial(link=logit))
summary(mod1)

mod2<-glm(Sim~Distance*Reldiff,data=mydat,family=gaussian(link=log),start=c(0,0,0,0))

quantile(mydat$Distance)

new.data<-data.frame(Distance=rep(c(0.00,  83922.93, 182868.45, 240369.10, 267699.66),each=5), Reldiff=rep(c(0,.5,1,1.5,2),5))

predy<-predict(mod1, new.data, type="response")
plotty<-cbind(new.data,predy)

predy2<-predict(mod2, new.data, type="response")
plotty2<-cbind(new.data,predy2)

cc<-ggplot()+
  geom_line(data=plotty,aes(Distance,predy,color=as.factor(Reldiff)),size=1.5)+
  ggthemes::theme_few()+ggtitle("quasibinomial")+scale_color_viridis_d("Invasion Difference")


mod1a<-glm(Sim~Distance*Reldiff,data=mydat,family=quasibinomial(link=logit))
summary(mod1a)

predya<-predict(mod1a, new.data, type="response")
plottya<-cbind(new.data,predya)



cc<-ggplot()+
  geom_line(data=plotty,aes(Distance,predy,color=as.factor(Reldiff)),size=1.5)+
  ggthemes::theme_few()+ggtitle("quasibinomial")+scale_color_viridis_d("Invasion Difference")



jpeg("..//git/bioticHogs/plots/modnopooling_NEcoast.jpeg")
cc
dev.off()
range(mydat$Sim)
mydat$Sim2<-ifelse(mydat$Sim==1,.9999999,mydat$Sim)
test<-brm(Sim2~Dist.z*Diff.z+(1|mm(RowID,ColumnID)),data=mydat,family=zero_inflated_beta())
summary(test)
834*834

save.image("../git/bioticHogs/NE.coastalincomplee.Rda")
