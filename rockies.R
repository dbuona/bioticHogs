###The beginning of this code is a loop to try and calculate all the pairwise dissimilarites
## for all level three eco regions, but it takes forever and needs a super computer.

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
library("brms")

setwd("~/PowellCenter/")
set.seed(3)

d<-read.csv("Data/FULLDatabase_05272022.csv")
d.env<-read.csv("Data/SPCIS_plots_env_17thJune2022.csv")
d.ref<-read.csv("Data/FullDatabase_diversity.csv")

#####to start with lets just do native species

d<-filter(d,Resampled=="N") ### we'll need to figure out how to deal with multiple timepoints later
d.nat<-filter(d,NativeStatus=="N")
NE<-unique(d.env$US_L3NAME)

df2<-data.frame(beta_bal=numeric(),beta_gra=numeric(),beta_div=numeric(), geo_dist=numeric(),eco_region3=character())

for (i in c(1:length(NE))){

tlands<-dplyr::filter(d.env,US_L3NAME==NE[i]) ##start with 2

dat<-filter(d.nat, Plot %in% c(unique(tlands$Plot)))

dater<-dplyr::select(dat,Plot,AcceptedTaxonName,PctCov) ## rid of complilcating columns
dater<-dplyr::filter(dater,!is.na(dater$AcceptedTaxonName)) ## spreading the data to long for (matrix) won't work if there are NA's in taxon name so get rid of them here
goo<- tidyr::spread(dater,AcceptedTaxonName,PctCov)  ## spread to long       
goo[is.na(goo)] <- 0 ## make na's zeros
rownames(goo)<-goo$Plot # convert plot names to row names
goo<-dplyr::select(goo,-Plot) # remove plot column
goo<-as.matrix.data.frame(goo) # make it a matrix

###I think you need "count data" for this packages to work so multply everything by 100
#goo<-round(goo*100)


### now try geogrpahic distance of all plots #####
dater2<-dplyr::select(dat,Plot,Long,Lat) ## subset to lat longs
dater2<-distinct(dater2) ## remove dplicates casue each species was a row
rownames(dater2)<-dater2$Plot # convert plot names to row names
dater2<-dplyr::select(dater2,-Plot)


###compute the matrices
d1<-beta.pair.abund(goo)$beta.bray
d1a<-beta.pair.abund(goo)$beta.bray.bal
d1b<-beta.pair.abund(goo)$beta.bray.gra
d2<-dist(dater2,method="euclidean")

dfhere2 <- data.frame(beta_bal=as.vector(d1a), beta_gra=as.vector(d1b),beta_div=as.vector(d1), geo_dist=as.vector(d2),eco_region3=NE[i])  ##generate fake data
df2 <- rbind(df2, dfhere2) 
}
unique(df2$eco_region3) # 14 eco regions in 30 hours

write.csv(df2,"bioticHogs/westernbetas.csv")
stop()
colnames(df2)
ggplot(df2,aes(geo_dist,beta_div))+
  geom_smooth(method="gam", formula= (y ~ s(log(x))),aes(color=eco_region3))

#nam<-paste("ecoregion",i,sep="")
#assign(nam,  decay.model(d1, d2, y.type="dissimilarities", model.type="exp", perm=100))
jpeg("~/git/bioticHogs/plots/Northeast_native_beta_occ.jpeg")
plot.decay(decay1, col="purple",ylim=c(0.8,1),xlim=c(0,7), remove.dots=TRUE)
plot.decay(decay2, col="blue", remove.dots=TRUE, add=TRUE)
plot.decay(decay3, col="darkgreen", remove.dots=TRUE, add=TRUE)
plot.decay(decay4, col="orange", remove.dots=TRUE, add=TRUE)
plot.decay(decay5, col="lightblue", remove.dots=TRUE, add=TRUE)
plot.decay(decay6, col="brown", remove.dots=TRUE, add=TRUE)
plot.decay(decay7, col="black", remove.dots=TRUE, add=TRUE)
legend("bottomright", legend=c("59","58","83","82","60","62","84"),title="Level 3 Ecoregion",
       col=c("purple", "blue","darkgreen","orange","lightblue","brown","black"), lty=c(1,1,1,1,1,1,1), cex=0.8)
dev.off()
