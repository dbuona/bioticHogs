####started June 16, 2022 by Dan B.
#### some compuationsal things to discuss
### HillR packages converts to relative cover, so what happen is if you give it relative cover already?
## How does plot size effect these measures?
# these a computationally energetic
# are there random effects in this model? if so what? 

### potential model:
#dissimilarity~distance*invasionmean*invasiondiff


rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
graphics.off()

library(dplyr)
library(vegan)
library(geodist)
library(ggplot2)
library(hillR)
library(betapart)
library(reshape2)
library(brms)

setwd("~/PowellCenter/")
set.seed(3)

d<-read.csv("Data/FULLDatabase_05272022.csv")
d.ref<-read.csv("Data/FullDatabase_diversity.csv")
d.eco<-read.csv("Data/MergedDatasets_16July2021_EB_BB_hex_state_eco_EVT.csv") ###not sure if this corresponds
d.eco<-filter(d.eco,NA_L3NAME=="Atlantic Coastal Pine Barrens")

dat<-dplyr::filter(d,Dataset=="VNHP") ## subsemt to just virginia data for now
dat<-filter(dat,Resampled=="N")

ref.prestine<-filter(d.ref,RelCov_I==0)
unique(dat$Plot)
range(dat$Lat)
range(dat$Long)
dat<-filter(dat,Lat>=38)




use<-sample_n(as.data.frame(unique(dat$Plot)),500)
 ## choose 100 random site for now for computational purpos

dat.p<-filter(dat,Plot %in% c(unique(ref.prestine$Plot)))
use_prestine<-sample_n(as.data.frame(unique(dat.p$Plot)),500)
## remove site that that have been resampled...will need to think about how to deal with this later


daty<-dplyr::filter(dat,Plot %in% c(use$`unique(dat$Plot)`))
dater<-dplyr::select(daty,Plot,AcceptedTaxonName,PctCov_100) ## rid of complilcating columns
dater<-dplyr::filter(dater,!is.na(dater$AcceptedTaxonName)) ## spreading the data to long for (matrix) won't work if there are NA's in taxon name so get rid of them here
goo<- tidyr::spread(dater,AcceptedTaxonName,PctCov_100)  ## spread to long       
goo[is.na(goo)] <- 0 ## make na's zeros
rownames(goo)<-goo$Plot # convert plot names to row names
goo<-dplyr::select(goo,-Plot) # remove plot column
goo<-as.matrix.data.frame(goo) # make it a matrix



###prestine data)
daty.p<-dplyr::filter(dat.p,Plot %in% c(use$`unique(dat$Plot)`))
dater.p<-dplyr::select(daty.p,Plot,AcceptedTaxonName,PctCov_100) ## rid of complilcating columns
dater.p<-dplyr::filter(dater.p,!is.na(dater.p$AcceptedTaxonName)) ## spreading the data to long for (matrix) won't work if there are NA's in taxon name so get rid of them here
goo.p<- tidyr::spread(dater.p,AcceptedTaxonName,PctCov_100)  ## spread to long       
goo.p[is.na(goo.p)] <- 0 ## make na's zeros
rownames(goo.p)<-goo.p$Plot # convert plot names to row names
goo.p<-dplyr::select(goo.p,-Plot) # remove plot column
goo.p<-as.matrix.data.frame(goo.p) # make it a matrix


#a<-vegdist(goo,method="jaccard",binary=TRUE)# calculate pairwise dissiminarity with horn index


### now try geogrpahic distance of all plots #####
dater2<-dplyr::select(daty,Plot,Long,Lat) ## subset to lat longs
dater2<-distinct(dater2) ## remove dplicates casue each species was a row
rownames(dater2)<-dater2$Plot # convert plot names to row names
dater2<-dplyr::select(dater2,-Plot) # remove plot column


dater2.p<-dplyr::select(daty.p,Plot,Long,Lat) ## subset to lat longs
dater2.p<-distinct(dater2.p) ## remove dplicates casue each species was a row
rownames(dater2.p)<-dater2.p$Plot # convert plot names to row names
dater2.p<-dplyr::select(dater2.p,-Plot) # remove plot column

#######betapartpackages#######################
#remember, this doesnt work right on % data

d1<-beta.pair.abund(goo)$beta.bray
d2<-dist(dater2)

beta.multi.abund(goo)

beta.multi.abund(goo.p)

d1.p<-beta.pair.abund(goo.p)$beta.bray
d2.p<-dist(dater2.p)



BCI.decay.exp<-decay.model(d1, d2, y.type="dissimilarities", model.type="exp", perm=100)
BCI.decay.exp.p<-decay.model(d1.p, d2.p, y.type="dissimilarities", model.type="exp", perm=100)

?decay.model()
plot.decay(BCI.decay.exp, col=rgb(0,0,0,0.5),ylim=c(0.3,1))
plot.decay(BCI.decay.exp, col="red", remove.dots=TRUE, add=TRUE)
plot.decay(BCI.decay.exp.p, col="blue", remove.dots=TRUE, add=TRUE)




######################################





b<-geodist(dater2,measure="haversine") ## haversime distance

colnames(b)<-rownames(dater2)
rownames(b)<-rownames(dater2)
#b[upper.tri(b,diag= TRUE)] <- NA

b <- data.frame(site1=rep(row.names(b),ncol(b)),
                   site2=rep(colnames(b),each=nrow(b)),
                   distance=as.vector(b))




a<-hill_taxa_parti_pairwise(goo, q = 2,pairs = "full") ### caclulate pairwise hill number based on simpsons (abudance weighted)

mydata<-left_join(a,b) ## data frame with distance and dissimilarity



##compute invasion level differences
d.ref.here<-dplyr::filter(d.ref, Plot %in% c(use$`unique(dat$Plot)`) ) ### 

d.ref.here<-select(d.ref.here,Plot,RelCov_I)
#d.ref.here<-distinct(d.ref.here)
rownames(d.ref.here)<-d.ref.here$Plot
#d.ref.here<-select(d.ref.here,RelCov_I)

d.invasion<-data.frame(site1=rep(d.ref.here$Plot,100),site2=rep(d.ref.here$Plot,each=100),inv1=rep(d.ref.here$RelCov_I,100),inv2=rep(d.ref.here$RelCov_I,each=100))
d.invasion$inv1<-as.numeric(d.invasion$inv1)
d.invasion$inv2<-as.numeric(d.invasion$inv2)
d.invasion$diffinv<-abs(d.invasion$inv1-d.invasion$inv2)
d.invasion$meaninv<-(d.invasion$inv1+d.invasion$inv2)/2

d.invasion<-select(d.invasion,-inv1,-inv2)
mydata<-left_join(mydata,d.invasion)

mydata$diffinv.bin<-ifelse(mydata$diffinv<.2,"low","high")
mydata$meaninv.bin<-ifelse(mydata$meaninv<.2,"low","high")
mydata$inv_composite<-mydata$diffinv*mydata$meaninv
mydata$inv.bin<-ifelse(mydata$meaninv<4.178227e-01 ,"low","high")

quantile(mydata$inv_composite,na.rm=TRUE)

mydata<-filter(mydata,site1!=site2)


ggplot(mydata,aes(distance,local_similarity))+geom_smooth(aes(color=meaninv.bin))

ggplot(mydata,aes(distance,region_similarity))+
  geom_smooth(aes(color=inv.bin),method="gam")

ggplot(mydata,aes(distance,TD_beta))+geom_point(aes(color=inv_composite))+
  geom_smooth(aes(linetype=inv.bin),method="gam")

ggplot(mydata,aes(distance,region_similarity))+geom_point(aes(color=inv_composite))+
  geom_smooth(aes(linetype=inv.bin),method="gam")


mod1<-brm(local_similarity~distance,data=mydata,inits = 0,
          zero_inflated_beta(link = "logit", link_phi = "log", link_zi = "logit"))
stop()
###############################################
###Below is other ways of calulating diversity with different data subsets
aa<-hill_taxa_parti_pairwise(goo, q = 2,output = "matrix",pairs = "full")$TD_beta

c<-hill_taxa_parti_pairwise(goo, q = 2,output = "matrix",pairs = "full")$local_similarity
e<-hill_taxa_parti_pairwise(goo, q = 2,output = "matrix",pairs = "full")$region_similarity

data <- data.frame(betaDiv = as.vector(a),
                   Distance = as.vector(b))



data2 <- data.frame(betaDiv = as.vector(aa),
                   Distance = as.vector(b))

dat3<- data.frame(betaDiv = as.vector(c),
                  Distance = as.vector(b))

dat4<- data.frame(betaDiv = as.vector(e),
                  Distance = as.vector(b))





plota<-ggplot(data=data, aes(x = Distance,y = betaDiv ))+geom_point()+ geom_smooth(method="gam")+
  labs(x="Distance", y= "richness-based") 

plotb<-ggplot(data=data2, aes(x = Distance,y = betaDiv ))+geom_point()+ geom_smooth(method="gam")+
  labs(x="Distance", y= "abundance based") 


plotc<-ggplot(data=dat3, aes(x = Distance,y = betaDiv ))+geom_point()+ geom_smooth(method="gam")+
  labs(x="Distance", y= "local_sim") 

plote<-ggplot(data=dat4, aes(x = Distance,y = betaDiv ))+geom_point()+ geom_smooth(method="gam")+
  labs(x="Distance", y= "regional_sim") 

ggpubr::ggarrange(plota,plotb)
ggpubr::ggarrange(plotc,plote)

### how to interpret Hill numbers?
### the range is 1-2 
dummy = FD::dummy
goober<-hill_taxa_parti_pairwise(comm = dummy$abun, q = 2)
## from paper: local is sorenson and regional is jaccard for q=0
## local is morissita-Hornan d classic :"regional overlap" is regional

#How does sinvasion interaction with this?

date<-left_join(dater,d.ref)
date<-dplyr::filter(date,TotalPctCover_I==0)
datea<-dplyr::select(date,Plot,AcceptedTaxonName,PctCov_100) ## rid of complilcating columns

datea<-dplyr::filter(datea,!is.na(date$AcceptedTaxonName)) ## spreading the data to long for (matrix) won't work if there are NA's in taxon name so get rid of them here


goo2<- tidyr::spread(datea,AcceptedTaxonName,PctCov_100)  ## spread to long       
goo2[is.na(goo2)] <- 0 ## make na's zeros


rownames(goo2)<-goo2$Plot # convert plot names to row names
goo2<-dplyr::select(goo2,-Plot) # remove plot column
goo2<-as.matrix.data.frame(goo2) 

googoo2<-dplyr::select(date,Plot,Long,Lat) ## subset to lat longs

googoo2<-distinct(googoo2) ## remove dplicates casue each species was a row

rownames(googoo2)<-googoo2$Plot # convert plot names to row names
googoo2<-dplyr::select(googoo2,-Plot) # remove plot column

bee<-geodist(googoo2,measure="haversine")
aeh<-hill_taxa_parti_pairwise(goo2, q = 0,output = "matrix",pairs = "full")$TD_beta

cee<-hill_taxa_parti_pairwise(goo2, q = 2,output = "matrix",pairs = "full")$local_similarity
eeeh<-hill_taxa_parti_pairwise(goo2, q = 2,output = "matrix",pairs = "full")$region_similarity


dataa<- data.frame(betaDiv = as.vector(aeh),
                  Distance = as.vector(bee))

datathree<- data.frame(betaDiv = as.vector(cee),
                   Distance = as.vector(bee))

datafour<- data.frame(betaDiv = as.vector(eeeh),
                       Distance = as.vector(bee))


ggplot()+geom_point(data=data,aes(x = Distance,y = betaDiv),alpha=0.5,color="green")+
  geom_point(data=dataa,aes(x = Distance,y = betaDiv),alpha=0.5,color="blue")+
  geom_smooth(data=data,aes(x = Distance,y = betaDiv),color="forestgreen")+
  geom_smooth(data=dataa,aes(x = Distance,y = betaDiv),color="royalblue")



woo<-ggplot()+geom_point(data=dat3,aes(x = Distance,y = betaDiv),alpha=0.5,color="green")+
  
  geom_point(data=datathree,aes(x = Distance,y = betaDiv),alpha=0.5,color="lightblue")+
  geom_smooth(data=dat3,aes(x = Distance,y = betaDiv),color="forestgreen")+
  geom_smooth(data=datathree,aes(x = Distance,y = betaDiv),color="royalblue")+ylab("local similarity")

woot<-ggplot()+geom_point(data=dat4,aes(x = Distance,y = betaDiv),alpha=0.5,color="green")+
  
  geom_point(data=datafour,aes(x = Distance,y = betaDiv),alpha=0.5,color="lightblue")+
  geom_smooth(data=dat4,aes(x = Distance,y = betaDiv),color="forestgreen")+
  geom_smooth(data=datafour,aes(x = Distance,y = betaDiv),color="royalblue")+ylab("regional similarity")

jpeg("virginia_test.jpeg")
ggpubr::ggarrange(woo,woot,labels = c("blue=native only plot", "green: native +invaded") ) 
dev.off()
                  

woo2<-ggplot()+geom_point(data=dat3,aes(x = Distance,y = betaDiv),alpha=0.5,color="green")+
  
  geom_point(data=datathree,aes(x = Distance,y = betaDiv),alpha=0.5,color="lightblue")+
  geom_smooth(data=dat3,method="lm",aes(x = Distance,y = betaDiv),color="forestgreen")+
  geom_smooth(data=datathree,method="lm",aes(x = Distance,y = betaDiv),color="royalblue")+ylab("local similarity")

woot2<-ggplot()+geom_point(data=dat4,aes(x = Distance,y = betaDiv),alpha=0.5,color="green")+
  
  geom_point(data=datafour,aes(x = Distance,y = betaDiv),alpha=0.5,color="lightblue")+
  geom_smooth(data=dat4,method="lm",aes(x = Distance,y = betaDiv),color="forestgreen")+
  geom_smooth(data=datafour,method="lm",aes(x = Distance,y = betaDiv),color="royalblue")+ylab("regional similarity")
       
jpeg("virginia_test_linear.jpeg")
ggpubr::ggarrange(woo2,woot2,labels = c("blue=native only plot", "green: native +invaded") ) 
dev.off()           
