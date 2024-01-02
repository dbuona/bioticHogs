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

setwd("~/git/bioticHogs/")

d<-read.csv("twohundo_daters.csv")
d$Beta2<-d$Beta-1
jpeg("..//bioticHogs/plots/twohundorev.jpeg",,width = 10, height=8,units = "in",res=200)
ggplot(d,aes(Distance,local_sim))+geom_point(aes(color=inv.bin),alpha=0.1,size=0.01)+geom_smooth(method="glm",method.args = list(family = "quasibinomial"),aes(color=inv.bin),size=1.2)+
  ggthemes::theme_few()+facet_wrap(~Ecoregion,scales="free_x")+scale_color_viridis_d(option = "inferno",begin = .1, end=.7)+
  scale_fill_viridis_d(option = "inferno",begin = .1, end=.7)
dev.off()

jpeg("..//bioticHogs/plots/twohundonovasion.jpeg",,width = 10, height=8,units = "in",res=200)
ggplot(d,aes(Distance,local_sim))+geom_point(alpha=0.1,size=0.01)+geom_smooth(method="glm",method.args = list(family = "quasibinomial"),aes(),size=1.2)+
  ggthemes::theme_few()+facet_wrap(~Ecoregion,scales="free_x")+scale_color_viridis_d(option = "inferno",begin = .1, end=.7)+
  scale_fill_viridis_d(option = "inferno",begin = .1, end=.7)
dev.off()

jpeg("..//bioticHogs/plots/twohundobeta.jpeg",width = 10, height=8,units = "in",res=200)
ggplot(d,aes(Distance,Beta2))+geom_point(aes(color=inv.bin),alpha=0.1,size=0.01)+geom_smooth(method="glm",method.args = list(family = "quasibinomial"),aes(color=inv.bin),size=1.2)+
  ggthemes::theme_few()+facet_wrap(~Ecoregion,scales="free_x")+scale_color_viridis_d(option = "inferno",begin = .1, end=.7)+
  scale_fill_viridis_d(option = "inferno",begin = .1, end=.7)
dev.off()

library(brms)

test3<-brm(local_sim~Dist.z*Inv.z+(1|Ecoregion/RowID:ColumnID),data=d,family =zero_inflated_beta() ) ### indicating that its crossed
summary(test3)


