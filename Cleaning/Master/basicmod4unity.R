

rm(list=ls())
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())
#rstan_options(auto_write = TRUE)

library(brms)
library(tidybayes)
if(FALSE){
t1<-read.csv("Analyses/Master/Input/pairwise_invaded_taxo.csv")
t1$community<-"invaded plots"

t2<-read.csv("Analyses/Master/Input/pairwise_pristine_taxo.csv")
t2$community<-"pristine plots"

t3<-read.csv("Analyses/Master/Input/pairwise_remnant_taxo.csv")
t3$community<-"remnant natives"



dat<-rbind(t1,t2)

dat$Dist_jit<-ifelse(dat$Distance==0,.00001,dat$Distance)
dat$logDist<-log(dat$Dist_jit)
p<-read.csv("Data/FULLDatabase_10272022_diversity_PctCov_100_10272022_INHABIT.csv")
p<-dplyr::select(p,Dataset,L4_KEY,L3_KEY, L2_KEY, L1_KEY)
p<-dplyr::distinct(p)
p$group<-paste(p$Dataset,p$L4_KEY,sep="_")
dat<-dplyr::left_join(dat,p)

dat<-filter(dat,community!="remnant natives")


unique(dat$community)
write.csv(dat,"Analyses/Master/Input/taxdataforunity.csv")

## this is for checking things on
library(dplyr)
dat<-sample_frac(dat,.0001)
}

setwd("Documents/git/bioticHogs/")
dat<-read.csv("Analyses/Master/Input/taxdataforunity.csv")
saveRDS(dat,"Analyses/Master/Input/taxdat4unity.rds")
colnames(dat)
dater<-dplyr::select(dat,site1,site2,local_similarity,logDist,Dataset,L4_KEY,L1_KEY)
write.csv(dater,"Analyses/Master/Input/taxdataforunity2.csv")
saveRDS(dater,"Analyses/Master/Input/taxdat4unity2.rds")
tax1mod<-brm(
  bf(local_similarity~logDist*community+(logDist*community|L1_KEY)+(1|L4_KEY)+(1|Dataset)+(1|mm(site1,site2)),
     phi ~ logDist*community+(logDist*community|L1_KEY)+(1|L4_KEY)+(1|Dataset)+(1|mm(site1,site2)),
  zoi ~ logDist*community+(logDist*community|L1_KEY)+(1|L4_KEY)+(1|Dataset)+(1|mm(site1,site2)),
  coi ~ logDist*community+(logDist*community|L1_KEY)+(1|L4_KEY)+(1|Dataset)+(1|mm(site1,site2))),
  data = dat,
  family = zero_one_inflated_beta(),
  control=list(adapt_delta=.9),
  chains = 2, iter = 2000, warmup = 1000,
  cores = 4, seed = 1234, init=0)      



if(FALSE){

new.data<-data.frame(L1_KEY=rep(unique(dat2$L1_KEY),each=3),community=rep(unique(dat2$community),8),logDist=rep(c(-1,1,2,3,4,5,6,7,8),each=24))

modpred<- tax1mod %>% 
  epred_draws(newdata =new.data,ndraws = 100,re_formula = ~(logDist*community|L1_KEY))

modpred$grouper<-paste(modpred$L1_KEY,modpred$community,modpred$.draw)

library(ggplot2)
ggplot(modpred,aes(logDist,.epred))+geom_line(aes(color=community,group=grouper))+facet_wrap(~L1_KEY)
ggplot(modpred,aes(logDist,.epred))+geom_smooth(aes(color=community),method="lm")+facet_wrap(~L1_KEY)


ggplot(modpred,aes(community,.epred))+stat_summary()+facet_wrap(~L1_KEY)
}