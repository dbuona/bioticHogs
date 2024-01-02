###for first meeting
### 1) compare how taxonomic diversity changes with levels of invasion.
### 2) Compare taxonomic and functional diversity
### 3) Compare R packages

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

setwd("~/PowellCenter/")
set.seed(3)

d<-read.csv("Data/FULLDatabase_05272022.csv")
d.env<-read.csv("Data/SPCIS_plots_env_17thJune2022.csv")
d.ref<-read.csv("Data/FullDatabase_diversity.csv")


### non in this are resalmpled

bos<-dplyr::filter(d.env,US_L3NAME== "Coast Range")
dat.bos<-filter(d, Plot %in% c(unique(bos$Plot)))
nrow(dat.bos)
table(dat.bos$Resampled)
dat.bos<-filter(dat.bos, Resampled=="N")

##start with 2


lowinv.bos<-dplyr::filter(d.ref,RelCov_I<= 15)
dat.bos.low<-filter(dat.bos, Plot %in% c(unique(lowinv.bos$Plot)))

dat.bos.pristine<-filter(dat.bos.low,NativeStatus=="N")

### make a community matrix
dater.pristine<-dplyr::select(dat.bos.pristine,Plot,AcceptedTaxonName,PctCov) ## rid of complilcating columns
dater.pristine<-dplyr::filter(dater.pristine,!is.na(dater.pristine$AcceptedTaxonName)) ## spreading the data to long for (matrix) won't work if there are NA's in taxon name so get rid of them here
range(dater.pristine$PctCov) #high 122
pristine<- tidyr::spread(dater.pristine,AcceptedTaxonName,PctCov)  ## spread to long       
pristine[is.na(pristine)] <- 0 ## make na's zeros
rownames(pristine)<-pristine$Plot # convert plot names to row names
pristine<-dplyr::select(pristine,-Plot) # remove plot column
pristine<-as.matrix.data.frame(pristine) # make it a matrix

### calculate betadiversity
d1<-beta.pair.abund(pristine,index.family = "bray")$beta.bray
d1a<-beta.pair.abund(pristine, index.family = "bray")$beta.bray.bal
d1b<-beta.pair.abund(pristine, index.family = "bray")$beta.bray.gra

#### calculate distance matrix
dater2.pri<-dplyr::select(dat.bos.pristine,Plot,Long,Lat) ## subset to lat longs
dater2.pri<-distinct(dater2.pri) ## remove dplicates casue each species was a row
rownames(dater2.pri)<-dater2.pri$Plot # convert plot names to row names
dater2.pri<-dplyr::select(dater2.pri,-Plot)

d2<-dist(dater2.pri,method="euclidean") ### distance matrix

beta.decay.pristine<-decay.model(d1, d2, y.type="dissimilarities", model.type="exp", perm=100)

plot.decay(beta.decay.pristine, col="purple",ylim=c(0.8,1), remove.dots=TRUE)

###### compare to high invasion plots
dater.high<-dplyr::select(dat.bos,Plot,AcceptedTaxonName,PctCov) ## rid of complilcating columns
dater.high<-dplyr::filter(dater.high,!is.na(dater.high$AcceptedTaxonName)) ## spreading the data to long for (matrix) won't work if there are NA's in taxon name so get rid of them here
range(dater.high$PctCov) #high 122
high<- tidyr::spread(dater.high,AcceptedTaxonName,PctCov)  ## spread to long       
high[is.na(high)] <- 0 ## make na's zeros
rownames(high)<-high$Plot # convert plot names to row names
high<-dplyr::select(high,-Plot) # remove plot column
high<-as.matrix.data.frame(high) # make it a matrix


### calculate betadiversity
d1.h<-beta.pair.abund(high,index.family = "bray")$beta.bray
d1a.h<-beta.pair.abund(high, index.family = "bray")$beta.bray.bal
d1b.h<-beta.pair.abund(high, index.family = "bray")$beta.bray.gra

nrow(high)
#### calculate distance matrix
dater2.high<-dplyr::select(dat.bos,Plot,Long,Lat) ## subset to lat longs
dater2.high<-distinct(dater2.high) ## remove dplicates casue each species was a row
dater2.high<-filter(dater2.high,Plot!="15052116152265152015-09-01")
rownames(dater2.high)<-dater2.high$Plot # convert plot names to row names
dater2.high<-dplyr::select(dater2.high,-Plot)
d2.h<-dist(dater2.high,method="euclidean",diag = TRUE) ### distance matrix


setdiff( rownames(dater2.high),rownames(high))



length(as.vector(d2.h))-length(as.vector(d1.h))

beta.decay.high<-decay.model(d1.h, d2.h, y.type="dissimilarities", model.type="exp", perm=100)

#### low invasion
dater.low<-dplyr::select(dat.bos.low,Plot,AcceptedTaxonName,PctCov) ## rid of complilcating columns
dater.low<-dplyr::filter(dater.low,!is.na(dater.low$AcceptedTaxonName)) ## spreading the data to long for (matrix) won't work if there are NA's in taxon name so get rid of them here
range(dater.low$PctCov) #high 122
low<- tidyr::spread(dater.low,AcceptedTaxonName,PctCov)  ## spread to long       
low[is.na(low)] <- 0 ## make na's zeros
rownames(low)<-low$Plot # convert plot names to row names
low<-dplyr::select(low,-Plot) # remove plot column
low<-as.matrix.data.frame(low) # make it a matrix

### calculate betadiversity
d1.l<-beta.pair.abund(low,index.family = "bray")$beta.bray
d1a.l<-beta.pair.abund(low, index.family = "bray")$beta.bray.bal
d1b.l<-beta.pair.abund(low, index.family = "bray")$beta.bray.gra

#### calculate distance matrix
dater2.low<-dplyr::select(dat.bos.low,Plot,Long,Lat) ## subset to lat longs
dater2.low<-distinct(dater2.low) ## remove dplicates casue each species was a row
rownames(dater2.low)<-dater2.low$Plot # convert plot names to row names
dater2.low<-dplyr::select(dater2.low,-Plot)

d2.l<-dist(dater2.low,method="euclidean") ### distance matrix

beta.decay.low<-decay.model(d1.l, d2.l, y.type="dissimilarities", model.type="exp", perm=100)

###invaders only
dat.bos.invonly<-filter(dat.bos.low,NativeStatus!="N")

dater.invonly<-dplyr::select(dat.bos.invonly,Plot,AcceptedTaxonName,PctCov) ## rid of complilcating columns
dater.invonly<-dplyr::filter(dater.invonly,!is.na(dater.invonly$AcceptedTaxonName)) ## spreading the data to long for (matrix) won't work if there are NA's in taxon name so get rid of them here
range(dater.invonly$PctCov) #high 30
invonly<- tidyr::spread(dater.invonly,AcceptedTaxonName,PctCov)  ## spread to long       
invonly[is.na(invonly)] <- 0 ## make na's zeros
rownames(invonly)<-invonly$Plot # convert plot names to row names
invonly<-dplyr::select(invonly,-Plot) # remove plot column
invonly<-as.matrix.data.frame(invonly) # make it a matrix

d1.i<-beta.pair.abund(invonly,index.family = "bray")$beta.bray
d1a.i<-beta.pair.abund(invonly, index.family = "bray")$beta.bray.bal
d1b.i<-beta.pair.abund(invonly, index.family = "bray")$beta.bray.gra

dater2.invonly<-dplyr::select(dat.bos.invonly,Plot,Long,Lat) ## subset to lat longs
dater2.invonly<-distinct(dater2.invonly) ## remove dplicates casue each species was a row
rownames(dater2.invonly)<-dater2.invonly$Plot # convert plot names to row names
dater2.invonly<-dplyr::select(dater2.invonly,-Plot)

d2.i<-dist(dater2.invonly,method="euclidean") ### distance matrix

beta.decay.invonly<-decay.model(d1.i, d2.i, y.type="dissimilarities", model.type="exp", perm=100)


###compare
jpeg("../git/bioticHogs/plots/coastrange_1_beta.jpeg")
plot.decay(beta.decay.pristine, col="dark green",ylim=c(0.8,1), remove.dots=TRUE)
plot.decay(beta.decay.low, col="orange",ylim=c(0.8,1), remove.dots=TRUE,add=TRUE)
plot.decay(beta.decay.high, col="red",ylim=c(0.8,1), remove.dots=TRUE,add=TRUE)
plot.decay(beta.decay.invonly, col="black",lty = 2,ylim=c(0.8,1), remove.dots=TRUE,add=TRUE)
legend("bottomright", legend=c("pristine","low","full","introduced only"),title="Invasion Level",
       col=c("darkgreen","orange","red","black"), lty=c(1,1,1,2), cex=.8)
title("Coast Range- Beta diversity")
dev.off()





###### balance

bal.decay.pristine<-decay.model(d1a, d2, y.type="dissimilarities", model.type="exp", perm=100)

bal.decay.low<-decay.model(d1a.l, d2.l, y.type="dissimilarities", model.type="exp", perm=100)

bal.decay.high<-decay.model(d1a.h, d2.h, y.type="dissimilarities", model.type="exp", perm=100)

bal.decay.invonly<-decay.model(d1a.i, d2.i, y.type="dissimilarities", model.type="exp", perm=100)


jpeg("../git/bioticHogs/plots/coastrange_1_bal.jpeg")
plot.decay(bal.decay.pristine, col="dark green",ylim=c(0.7,1), remove.dots=TRUE)
plot.decay(bal.decay.low, col="orange",ylim=c(0.7,1), remove.dots=TRUE,add=TRUE)
plot.decay(bal.decay.high, col="red",ylim=c(0.7,1), remove.dots=TRUE,add=TRUE)
plot.decay(bal.decay.invonly, col="black",lty = 2,ylim=c(0.7,1), remove.dots=TRUE,add=TRUE)
legend("bottomright", legend=c("pristine","low","full", "introduced only"),title="Invasion Level",
       col=c("darkgreen","orange","red","black"), lty=c(1,1,1,2), cex=.8)
title("Coastal Range- bal. diversity")
dev.off()

####
gra.decay.pristine<-decay.model(d1b, d2, y.type="dissimilarities", model.type="exp", perm=100)

gra.decay.low<-decay.model(d1b.l, d2.l, y.type="dissimilarities", model.type="exp", perm=100)

gra.decay.high<-decay.model(d1b.h, d2.h, y.type="dissimilarities", model.type="exp", perm=100)

gra.decay.invonly<-decay.model(d1b.i, d2.i, y.type="dissimilarities", model.type="exp", perm=100)

jpeg("../git/bioticHogs/plots/coastrange_1_gra.jpeg")
plot.decay(gra.decay.pristine, col="dark green",ylim=c(0,.11), remove.dots=TRUE)
plot.decay(gra.decay.low, col="orange",ylim=c(0,.11), remove.dots=TRUE,add=TRUE)
plot.decay(gra.decay.high, col="red",ylim=c(0,.11), remove.dots=TRUE,add=TRUE)
plot.decay(gra.decay.invonly, col="black",lty=2,ylim=c(0,.11), remove.dots=TRUE,add=TRUE)
legend("topright", legend=c("pristine","low","full","introduced only"),title="Invasion Level",
       col=c("darkgreen","orange","red","black"), lty=c(1,1,1,2), cex=.8)
title("coastrange- gra. diversity")
dev.off()

####try it with HillR

#pristin.hill<-hill_taxa_parti_pairwise(pristine,q=2)
pristin.hill.trix<-hill_taxa_parti_pairwise(pristine,q=2,output = "matrix")

#low.hill<-hill_taxa_parti_pairwise(low,q=2)
low.hill.trix<-hill_taxa_parti_pairwise(low,q=2,output = "matrix")

#high.hill<-hill_taxa_parti_pairwise(high,q=2)
high.hill.trix<-hill_taxa_parti_pairwise(high,q=2,output = "matrix")

#invonly.hill<-hill_taxa_parti_pairwise(invonly,q=2)
invonly.hill.trix<-hill_taxa_parti_pairwise(invonly,q=2,output="matrix")

head(d2.l)
colnames(dater2.pri)<-c("lon","lat")
colnames(dater2.low)<-c("lon","lat")
colnames(dater2.high)<-c("lon","lat")
colnames(dater2.invonly)<-c("lon","lat")
b<-geodist(dater2.pri,paired = TRUE,measure = "haversine")
b.l<-geodist(dater2.low,paired = TRUE,measure = "haversine")
b.h<-geodist(dater2.high,paired = TRUE,measure = "haversine")
b.i<-geodist(dater2.invonly,paired = TRUE,measure = "haversine")

length(b.h)
length(high.hill.trix$TD_beta)


mydat<-data.frame(dist=as.vector(b), beta=as.vector(pristin.hill.trix$TD_beta-1),invasion="pristine")
mydat2<-data.frame(dist=as.vector(b.l), beta=as.vector(low.hill.trix$TD_beta-1),invasion="low")
mydat3<-data.frame(dist=as.vector(b.h), beta=as.vector(high.hill.trix$TD_beta-1),invasion="high")
mydat4<-data.frame(dist=as.vector(b.i), beta=as.vector(invonly.hill.trix$TD_beta-1),invasion="introduced only")
mydatty<-rbind(mydat,mydat2,mydat3,mydat4)

jpeg("..//git/bioticHogs/plots/hillcoastalrange_1.jpeg")
ggplot(mydatty,aes(dist,beta))+ geom_smooth(method="glm",method.args = list(family = "quasibinomial"),aes(color=invasion,fill=invasion))+
  ggthemes::theme_few()+scale_color_manual(values = c("firebrick","black","orange","forestgreen"))+
  scale_fill_manual(values = c("firebrick","black","orange","forestgreen"))
dev.off()

mydat<-data.frame(dist=as.vector(b), beta=as.vector(pristin.hill.trix$local_similarity),invasion="pristine")
mydat2<-data.frame(dist=as.vector(b.l), beta=as.vector(low.hill.trix$local_similarity),invasion="low")
mydat3<-data.frame(dist=as.vector(b.h), beta=as.vector(high.hill.trix$local_similarity),invasion="high")
mydat4<-data.frame(dist=as.vector(b.i), beta=as.vector(invonly.hill.trix$local_similarity),invasion="introduced only")
mydatty2<-rbind(mydat,mydat2,mydat3,mydat4)

jpeg("..//git/bioticHogs/plots/hillmorsitahorn_1_coastalrange.jpeg")
ggplot(mydatty2,aes(dist,beta))+ geom_smooth(method="glm",method.args = list(family = "quasibinomial"),aes(color=invasion,fill=invasion))+
  ggthemes::theme_few()+scale_color_manual(values = c("firebrick","black","orange","forestgreen"))+
  scale_fill_manual(values = c("firebrick","black","orange","forestgreen"))+ylab("Morista-Horn")

dev.off()
##mantel tests
mantel(d1.h,d2.h,method = "spearman")
mantel(d1.i,d2.i,method = "spearman")
