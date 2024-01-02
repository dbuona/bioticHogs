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
d.env<-read.csv("Data/SPCIS_plots_env_17thJune2022.csv")
d.ref<-read.csv("Data/FullDatabase_diversity.csv")

#####to start with lets just do native species

unique(d.env$US_L3NAME)
d<-filter(d,Resampled=="N") ### we'll need to figure out how to deal with multiple timepoints later

d.nat<-filter(d,NativeStatus=="N")

tlands<-dplyr::filter(d.env,US_L3NAME=="North Central Appalachians") ## subsemt to just virginia data for now

dat<-filter(d.nat, Plot %in% c(unique(tlands$Plot)))

dater<-dplyr::select(dat,Plot,AcceptedTaxonName,PctCov) ## rid of complilcating columns
dater<-dplyr::filter(dater,!is.na(dater$AcceptedTaxonName)) ## spreading the data to long for (matrix) won't work if there are NA's in taxon name so get rid of them here
goo<- tidyr::spread(dater,AcceptedTaxonName,PctCov)  ## spread to long       
goo[is.na(goo)] <- 0 ## make na's zeros
rownames(goo)<-goo$Plot # convert plot names to row names
goo<-dplyr::select(goo,-Plot) # remove plot column
goo<-as.matrix.data.frame(goo) # make it a matrix




#a<-vegdist(goo,method="jaccard",binary=TRUE)# calculate pairwise dissiminarity with horn index


### now try geogrpahic distance of all plots #####
dater2<-dplyr::select(dat,Plot,Long,Lat) ## subset to lat longs
dater2<-distinct(dater2) ## remove dplicates casue each species was a row
rownames(dater2)<-dater2$Plot # convert plot names to row names
dater2<-dplyr::select(dater2,-Plot) # remove plot column

#######betapartpackages#######################
#remember, this doesnt work right on % data

d1<-beta.pair.abund(goo)$beta.bray
d1a<-beta.pair.abund(goo)$beta.bray.bal
d1b<-beta.pair.abund(goo)$beta.bray.gra


mydat<-data.frame(sim=d1,distance=d2)

d1<-as.vector(d1)
d2<-as.vector(d2)

ggplot(mydat,aes(d2,d1))+geom_point()+geom_smooth()

BCI.decay.exp<-decay.model(d1, d2, y.type="dissimilarities", model.type="exp", perm=100)
BCI.decay.exp.bal<-decay.model(d1a, d2, y.type="dissimilarities", model.type="exp", perm=100)
BCI.decay.exp.gra<-decay.model(d1b, d2, y.type="dissimilarities", model.type="exp", perm=100)

summary(BCI.decay.exp)

plot.decay(BCI.decay.exp, col=rgb(0,0,0,0.5),ylim=c(0,1))
plot.decay(BCI.decay.exp, col="red", remove.dots=TRUE, add=TRUE)

plot.decay(BCI.decay.exp.bal, col=rgb(0,0,0,0.5),ylim=c(0,1))
plot.decay(BCI.decay.exp.bal, col="blue", remove.dots=TRUE, add=TRUE)

plot.decay(BCI.decay.exp.gra, col=rgb(0,0,0,0.5),ylim=c(0,1))
plot.decay(BCI.decay.exp.gra, col="green", remove.dots=TRUE, add=TRUE)



###############################gr#######




b<-dist(dater2,ethm="haversine") ## haversime distance

colnames(b)<-rownames(dater2)
rownames(b)<-rownames(dater2)
#b[upper.tri(b,diag= TRUE)] <- NA

b <- data.frame(site1=rep(row.names(b),ncol(b)),
                   site2=rep(colnames(b),each=nrow(b)),
                   distance=as.vector(b))




hill_taxa_parti_pairwise(goo, q = 2,pairs = "full") ### caclulate pairwise hill number based on simpsons (abudance weighted)

mydata<-left_join(a,d2) ## data frame with distance and dissimilarity







ggplot(mydata,aes(distance,local_similarity))+geom_smooth()

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
