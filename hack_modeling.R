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

setwd("~/PowellCenter/")
set.seed(3)

d<-read.csv("Data/FULLDatabase_05272022.csv")
d.env<-read.csv("Data/SPCIS_plots_env_17thJune2022.csv")
d.ref<-read.csv("Data/FullDatabase_diversity.csv")

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

high.hill<-hill_taxa_parti_pairwise(high,q=2,pairs="unique")
nrow(high.hill)

inv.data<-data.frame(Plot=high.hill$site1)
inv.data<-left_join(inv.data,driftless.ref)
inv.data2<-data.frame(Plot=high.hill$site2)
inv.data2<-left_join(inv.data2,driftless.ref)

inv.data<-cbind(inv.data,inv.data2)
colnames(inv.data)<-c("site1", "cov1","site2","cov2")
inv.data$reldiff<-abs((inv.data$cov1-inv.data$cov2))/(abs(inv.data$cov1+inv.data$cov2)/2)


goo<-inv.data %>% mutate_at(vars(reldiff), ~replace(., is.nan(.), 0))

goo<-select(goo,site1,site2,reldiff)
hill<-left_join(high.hill,goo)

dater2<-dplyr::select(dat.bos,Plot,Long,Lat) ## subset to lat longs
dater2<-distinct(dater2) ## remove dplicates casue each species was a row

#rownames(dater2.high)<-dater2.high$Plot # convert plot names to row names
#dater2<-dplyr::select(dater2,-Plot)




b<-geodist(dater2,paired = TRUE,measure = "haversine")
colnames(b)<-dater2$Plot
rownames(b)<-dater2$Plot
b2<-setNames(melt(b), c('site1', 'site2', 'distance'))
goony<-left_join(hill,b2)

nrow(goony)

ggplot(goony,aes(distance,local_similarity))+geom_point() +geom_smooth(method="glm",method.args = list(family = "quasibinomial"))+
  geom_smooth()
quantile(goony$distance,na.rm = TRUE)
quantile(goony$reldiff,na.rm = TRUE)

mod1<-glm(local_similarity~distance*reldiff,data=goony,family=quasibinomial(link=logit))
summary(mod1)

mod2<-glm(local_similarity~distance*reldiff,data=goony,family=gaussian(link=log),start=c(0,0,0,0))
S
mod2a<-brm(local_similarity~distance*reldiff, family=exgaussian(),
           data=goony)

summary(mod2a)


new.data<-data.frame(distance=rep(c(0.00 , 16078.82,  74454.06 ,137064.15, 329217.91),each=5), reldiff=rep(c(0,.5,1,1.5,2),5))

predy<-predict(mod1, new.data, type="response")
plotty<-cbind(new.data,predy)

predy2<-predict(mod2, new.data, type="response")
plotty2<-cbind(new.data,predy2)



ggplot()+
  geom_point(data=goony,aes(x=distance,y=local_similarity))+
  geom_smooth(data=plotty,aes(distance,predy,color=as.factor(reldiff)),method="glm",method.args = list(family = "quasibinomial"))+
  ggthemes::theme_few()

ggplot()+
geom_smooth(data=plotty,aes(distance,predy,color=as.factor(reldiff)),method="glm",method.args = list(family = "quasibinomial"))+
  ggthemes::theme_few()+ggtitle("quasibinomial")

cc<-ggplot()+
  geom_smooth(data=plotty2,aes(distance,predy,color=as.factor(reldiff)),method="glm",method.args = list(family = "quasibinomial"))+
  ggthemes::theme_few()+ggtitle("exponetial")

jpeg("..//git/bioticHogs/plots/firstlook_modcomp.jpeg")
ggpubr::ggarrange(aa,cc,common.legend = TRUE,labels=c("a)","b)"),widths=c(.6,.6))
dev.off()

jpeg("..//git/bioticHogs/plots/firstlook_mod.jpeg")
ggpubr::ggarrange(aa,bb,common.legend = TRUE,labels=c("a)","b)"),widths=c(.6,.5))  
dev.off()
summary(mod1)


###try again with everything as a matrix
