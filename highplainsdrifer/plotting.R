####see highplainsdrifter.R for most uptodat version of these loops.

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

setwd("~/Documents/git/bioticHogs/highplainsdrifer/")

d<-read.csv("highplainsdrifterdata.csv")
d.ref<-read.csv("~/Desktop/Powell/Data/FullDatabase_diversity_July4.csv")


oneway<-d %>% group_by(site1,Ecoregion,status,metric) %>% summarise(mean_beta=mean(beta_div))

d.ref<-select(d.ref,Plot,RelCov_I)
colnames(d.ref)[1]<-"site1"
oneway<-left_join(oneway,d.ref)
oneway<-filter(oneway,status!="native species only")
oneway.d<-d%>% group_by(site1,Ecoregion,status,metric) %>% summarise(mean_distance=mean(Distance))
oneway.d<-filter(oneway.d,status!="native species only")
goomba<-left_join(oneway.d,oneway)

ggplot(goomba,aes(mean_distance, mean_beta))+geom_point(aes(color=RelCov_I))+
  facet_wrap(~Ecoregion,scale="free_x")
goomba<-filter(goomba,metric=="occurance")
goomba<-filter(goomba,=="occurance")


mod1<-glm(mean_beta~mean_distance*RelCov_I,data=goomba, family=quasibinomial())
summary(mod1)

new.dat<-dataframe()


jpeg(filename = "highplainsdrifteroneway.jpeg",res=100)
ggplot(oneway,aes(RelCov_I,mean_beta))+geom_point(aes(shape=metric))+
  geom_smooth(method="lm",aes(linetype=metric))+ylab("Mean dissimilarity from focal site")+
  xlab("Invasive Relative Cover")+
  facet_wrap(~Ecoregion)+ggthemes::theme_few()

dev.off()



native<-filter(d,status=="native species only")
all<-filter(d,status!="native species only")

#native<-filter(native,metric!="occurance")
#all<-filter(all,metric!="occurance")


native<-select(native,,-status,-Distance)
colnames(native)[3:5]<-paste(colnames(native)[3:5],"native",sep = ".")

all<-select(all,,-status,-Distance)

mydat<-left_join(native,all,by=c("site1","site2","metric","Ecoregion"))
dev.off()
a<-ggplot(mydat,aes(beta_div.native,beta_div))+geom_point()+
  facet_grid(metric~Ecoregion)+geom_abline(intercept = 0,slope=1,color="red")+ggthemes::theme_few()+
  ylab("Pairwise dissimilarity native and introduced")+xlab("Pairwise dissimilarity native only")

b<-ggplot(mydat,aes(turnover.native,turnover))+geom_point()+
  facet_grid(metric~Ecoregion)+geom_abline(intercept = 0,slope=1,color="red")+ggthemes::theme_few()+
  ylab("Pairwise turnover native and introduced")+xlab("Pairwise turnover native only")


c<-ggplot(mydat,aes(nestedness.native,nestedness))+geom_point()+
  facet_grid(metric~Ecoregion)+geom_abline(intercept = 0,slope=1,color="red")+ggthemes::theme_few()+
  ylab("Pairwise nestedness native and introduced")+xlab("Pairwise nestedness native only")


jpeg(filename = "highplainsdrifterline.jpeg",width=6,height=12,unit='in',res=200)
ggpubr::ggarrange(a,b,c,ncol=1,nrow=3,labels = c("a)","b)","c)"))
dev.off()

mydat$H<-log((mydat$beta_div+.001)/(mydat$beta_div.native+.001))
mydf<- mydat[order(mydat$H),]


mydf = mydf[seq(1, nrow(mydf), 2), ]
colnames(d.ref)<-c("site1","RelCov1")
mydf<-left_join(mydf,d.ref)

d.ref2<-d.ref
colnames(d.ref2)<-c("site2","RelCov2")

mydf<-left_join(mydf,d.ref2)

mydf<-filter(mydf,RelCov1!=RelCov2)
mydf.high<-filter(mydf,Ecoregion=="High Plains")


ggplot(mydf.high,aes(H))+geom_histogram(bins=200)+facet_wrap(~Ecoregion,scales="free")




library(brms)
library(lme4)
mod1<-lmer(beta_div~beta_div.native*metric+(beta_div.native*metric|Ecoregion), data=mydat)
coef(mod1)
