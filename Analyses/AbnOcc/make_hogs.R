rm(list=ls())
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())
#rstan_options(auto_write = TRUE)
graphics.off()

library(dplyr)
library(tidyr)
library(hillR)
library(ggplot2)
#require(V.PhyloMaker2)
library(ape)

setwd("~/Documents/git/bioticHogs/")# mac
if(FALSE){
d1<-read.csv("Data/FULLDatabase_10272022.csv")

use.list<-read.csv("Analyses/AbnOcc/Input/env_matched_plots.csv")


d<-d1 %>%
  group_by(Plot) %>%
  filter(Year==max(Year)) %>%
  ungroup() 

use.list<-dplyr::select(use.list,-X)
colnames(use.list)[4]<-"pair"



d<-left_join(d,use.list)

d<-dplyr::filter(d,!is.na(d$AcceptedTaxonName))

matricize<-function(x){
  temp<-dplyr::select(x,pair,AcceptedTaxonName,PctCov_100) ### select useful columns
  temp<-dplyr::filter(temp,!is.na(temp$AcceptedTaxonName))
  temp<- tidyr::spread(temp,AcceptedTaxonName,PctCov_100)
  temp[is.na(temp)] <- 0 
  temp<-as.data.frame(temp)
  rownames(temp)<-temp$pair# convert plot names to row names
  temp<-dplyr::select(temp,-pair)
  temp<-as.matrix.data.frame(temp)#
  
} 

count.use<-use.list%>%group_by(L4_KEY,status) %>% count()
count.use<-filter(count.use,n>1)


list.i<-filter(use.list, status==1)### invaded
list.p<-filter(use.list, status==0) ### native



natyA<-data.frame()
natyO<-data.frame()
inyA<-data.frame()
inyO<-data.frame()

ecoregionz<-sort(unique(count.use$L4_KEY)) ### 396 ecoregions

for (z in c(1:length(ecoregionz))){
  
  tlands<-dplyr::filter(d,L4_KEY==ecoregionz[z])
  #pairzz<-unique(tlands$pair)
  #use.pairz<-pairzz#sample(pairzz,2, replace=FALSE)
  # d.use<-filter(d,pair %in% tlands$pair)
  
  dat.i<-filter(tlands,Plot %in% list.i$Plot)
  dat.p<-filter(tlands,Plot %in% list.p$Plot) 
  
  
  commsall<-matricize(dat.i)
  commsnative<-matricize(dat.p)
  
  datyall.occ<-hill_taxa_parti_pairwise(commsall,q = 0)
  datyall.occ$L4_KEY<-ecoregionz[z]
  #datyall.occ<-dplyr::select(datyall.occ,site1,site2,local_similarity)
  #colnames(datyall.occ)<-c("site1","site2","beta_I_occ")
  
  datyall.abn<-hill_taxa_parti_pairwise(commsall,q = 1)
  datyall.abn$L4_KEY<-ecoregionz[z]
  #datyall.abn<-dplyr::select(datyall.abn,site1,site2,local_similarity)
  #colnames(datyall.abn)<-c("site1","site2","beta_I_abn")
  
  daty.naty.occ<-hill_taxa_parti_pairwise(commsnative,q = 0)
  daty.naty.occ$L4_KEY<-ecoregionz[z]
  #daty.naty.occ<-dplyr::select(daty.naty.occ,site1,site2,local_similarity)
  #colnames(daty.naty.occ)<-c("site1","site2","beta_N_occ")
  
  daty.naty.abn<-hill_taxa_parti_pairwise(commsnative,q = 1)
  daty.naty.abn$L4_KEY<-ecoregionz[z]
  #daty.naty.abn<-dplyr::select(daty.naty.abn,site1,site2,local_similarity)
  #colnames(daty.naty.abn)<-c("site1","site2","beta_N_abn")
  
  natyA<-rbind(daty.naty.abn,natyA)
  natyO<-rbind(daty.naty.occ,natyO)
  
  inyA<-rbind(datyall.abn,inyA)
  inyO<-rbind(datyall.occ,inyO)
  
  print(paste("done",which(ecoregionz==ecoregionz[z])))   
}

InyA<-dplyr::select(inyA,site1,site2,local_similarity,L4_KEY)
NatyA<-dplyr::select(natyA,site1,site2,local_similarity,L4_KEY)
colnames(InyA)[3]<-"invA"
colnames(NatyA)[3]<-"natA"
Abn.env<-left_join(InyA,NatyA)

InyO<-dplyr::select(inyO,site1,site2,local_similarity,L4_KEY)
NatyO<-dplyr::select(natyO,site1,site2,local_similarity,L4_KEY)
colnames(InyO)[3]<-"invO"
colnames(NatyO)[3]<-"natO"
Occ.env<-left_join(InyO,NatyO)
ENV<-left_join(Occ.env,Abn.env)


ENV$H.A<-ENV$invA-ENV$natA
ENV$H.O<-ENV$invO-ENV$natO

write.csv(ENV,"Analyses/AbnOcc/Input/abn_v_occ_environmental.csv")
}
env<-read.csv("Analyses/AbnOcc/Input/abn_v_occ_environmental.csv")
#env<-filter(env2,invA!=0&natA!=0&invO!=0&invO!=0) #
####recover original plot names
use.list<-read.csv("Analyses/AbnOcc/Input/env_matched_plots.csv")

invaded<-filter(use.list,status==1)
uninvaded<-filter(use.list,status!=1)

invaded<-dplyr::select(invaded,Plot,subclass)
uninvaded<-dplyr::select(uninvaded,Plot,subclass)

uninvaded1<-uninvaded
uninvaded2<-uninvaded  
colnames(uninvaded1)<-c("uninvaded.plot.1","site1")
colnames(uninvaded2)<-c("uninvaded.plot.2","site2") 

invaded1<-invaded
invaded2<-invaded  
colnames(invaded1)<-c("invaded.plot.1","site1")
colnames(invaded2)<-c("invaded.plot.2","site2") 

env<-left_join(env,uninvaded1)
env<-left_join(env,uninvaded2)

env<-left_join(env,invaded1)
env<-left_join(env,invaded2)
###add dataset
d1<-read.csv("Data/FULLDatabase_10272022.csv")
dset<-dplyr::select(d1,Plot,Dataset)
dset<-distinct(dset)

dsetu1<-dset
dsetu2<-dset

dseti1<-dset
dseti2<-dset

colnames(dsetu1)<-c("uninvaded.plot.1","dataset.u.plot.1")
colnames(dsetu2)<-c("uninvaded.plot.2","dataset.u.plot.2")

colnames(dseti1)<-c("invaded.plot.1","dataset.i.plot.1")
colnames(dseti2)<-c("invaded.plot.2","dataset.i.plot.2")

env<-left_join(env,dsetu1)
env<-left_join(env,dsetu2)

env<-left_join(env,dseti1)
env<-left_join(env,dseti2)


##now get dominate invader
d<-d1 %>%
  group_by(Plot) %>%
  filter(Year==max(Year)) %>%
  ungroup() 


d.inv<-filter(d,NativeStatus=="I")

d.inv<-dplyr::select(d.inv,Plot,AcceptedTaxonName,PctCov_100)  
  
dom<-d.inv %>% group_by(Plot) %>% slice(which.max(PctCov_100))

dom1<-dom
dom2<-dom

colnames(dom1)<-c("invaded.plot.1","dom.invader.1","cov.dom.1")
colnames(dom2)<-c("invaded.plot.2","dom.invader.2","cov.dom.2")

env<-left_join(env,dom1)
env<-left_join(env,dom2)

env$quad<-NA
env$quad[which(env$H.A>0 &env$H.O>0)]<-1
env$quad[which(env$H.A<=0&env$H.O>0)]<-4  
env$quad[which(env$H.A<0&env$H.O<0)]<-3 
env$quad[which(env$H.A>0&env$H.O<=0)]<-2
env<-filter(env,!is.na(quad))

env$dif<-abs(env$H.O-env$H.A)
round(mean(env$dif,na.rm=TRUE),3)
round(sd(env$dif,na.rm=TRUE),3)
env$number<-ifelse(env$dif>=.5,1,0)
round(table(env$number)[2]/(table(env$number)[1]+table(env$number)[2]),2)


cor.test(env$H.A,env$H.O ,method="pearson")

countit<-env %>%group_by(quad) %>% count()    
countit$perc<-round(countit$n/ sum(countit$n),3)
sum(countit$n)
sum(countit$perc)

heat.time<-env
heat.time$H.A<-round(heat.time$H.A,1)
heat.time$H.O<-round(heat.time$H.O,1)


heat.time$catcher<-abs(heat.time$H.A)+abs(heat.time$H.O)
#heat.time<-filter(heat.time,catcher!=0)

heat.test<-heat.time %>% group_by(H.O,H.A) %>% count()
heat.test$frequency<-heat.test$n/sum(heat.test$n)
heat.test$frequency<-round(heat.test$frequency,2)

heat.test$frequency2<-ifelse(heat.test$frequency==0,NA,heat.test$frequency)
heat.test$frequency3<-ifelse(is.na(heat.test$frequency2),"",paste(heat.test$frequency2*100,"%",sep=""))

jpeg("Analyses/AbnOcc/heatmap.jpeg",width=8,height=8,units = 'in',res=300)
ggplot()+
  geom_point(data=env,aes(x=H.O,y=H.A),size=0.01)+
geom_tile(data=heat.test,aes(x=H.O,y=H.A,fill=frequency),alpha=.8)+
  tidybayes::theme_tidybayes()+
  scale_fill_distiller(palette = "PuRd",direction=1 )+
  geom_hline(yintercept=0,color="darkgray")+geom_vline(xintercept=0,color="darkgray")+
  ylab("abundance")+xlab("occurence")+
  #geom_text(data=heat.test,aes(x=H.O,y=H.A,label=frequency3),size=2,color="navy")+
#geom_smooth(method="lm",color="yellow")
geom_abline(slope = 1,intercept = 0,color="black",linetype="dotdash",size=1)
dev.off()



ggplot()+
  geom_point(data=env,aes(x=H.O,y=H.A),size=0.01)+
  geom_tile(data=heat.test,aes(x=H.O,y=H.A,fill=frequency),alpha=.8)+
  tidybayes::theme_tidybayes()+
  scale_fill_distiller(palette = "PuRd",direction=1 )+
  geom_hline(yintercept=0,color="darkgray")+geom_vline(xintercept=0,color="darkgray")+
  ylab("abundance")+xlab("occurence")+
  geom_text(data=heat.test,aes(x=H.O,y=H.A,label=frequency3),size=2,color="navy")+
  #geom_smooth(method="lm",color="yellow")
  geom_abline(slope = 1,intercept = 0,color="black",linetype="dotdash",size=1)


d.ref<-read.csv("Data/FULLDatabase_10272022_diversity_PctCov_100_10272022_INHABIT.csv")
#d.ref<-read.csv("Data/FULLDatabase_10272022_diversity_PctCov_100_10272022_INHABIT.csv")
d.ref<-d.ref %>%
  group_by(Plot) %>%
  filter(Year==max(Year)) %>%
  ungroup() 

d.key<-dplyr::select(d.ref,L4_KEY,L3_KEY,L2_KEY,L1_KEY)
d.key<-distinct(d.key)
#rm(d.ref)
env<-left_join(env,d.key)


####question: what drive differences

env$same.dom<-ifelse(env$dom.invader.1==env$dom.invader.2,1,0)
env$same.dat.i<-ifelse(env$dataset.i.plot.1==env$dataset.i.plot.2,1,0)
env$same.dat.u<-ifelse(env$dataset.u.plot.1==env$dataset.u.plot.2,1,0)
env$same.all<-ifelse(env$dataset.u.plot.1==env$dataset.u.plot.2& env$dataset.u.plot.1==env$dataset.i.plot.1&env$dataset.i.plot.1==env$dataset.i.plot.2,1,0)
table(env$same.all)
env$dif_nat<-abs(env$natO-env$natA)

table(env$same.dat.u)
table(env$same.dom)
zscore.2 <- function(x){(x-mean(x,na.rm=TRUE))/(2*sd(x,na.rm=TRUE))} # zscore for 01 and continuous variables

ggplot(env,aes(as.factor(same.dom),dif))+stat_summary()
ggplot(env,aes(as.factor(same.dat.u),dif))+stat_summary()

env$dom.z<-zscore.2(env$same.dom)
env$dat.i.z<-zscore.2(env$same.dat.i)
env$dat.u.z<-zscore.2(env$same.dat.u)
env$dat.z<-zscore.2(env$same.all)
env$ndif.z<-zscore.2(env$dif_nat)
env$disagree<-ifelse(env$quad==2 |env$quad==4,1,0)

summary(lm(dif~same.dom+same.all+dif_nat,data=env))
summary(lm(dif~dat.z+dom.z+ndif.z,data=env))
summary(glm(disagree~dat.z+dom.z+ndif.z,data=env,family = binomial(link="logit")))


ex1<-ggplot(env,aes(as.factor(same.all),dif))+stat_summary(size=0.001,alpha=0.2,aes(group=L4_KEY),position="jitter",shape=1)+
  stat_summary(
  fun.data= "mean_sdl", fun.args = list(mult = 1),size=.75)+coord_cartesian(ylim=c(0,.5))+ggthemes::theme_few(base_size = 9)+ylab("absolute difference abn vs. occ")+
  scale_x_discrete(name="",labels=c("multiple datasets","same dataset"))

ex2<-ggplot(env,aes(as.factor(same.dom),dif))+stat_summary(size=0.001,alpha=0.2,aes(group=L4_KEY),position="jitter",shape=1)+
  stat_summary(
  fun.data= "mean_sdl", fun.args = list(mult = 1),size=.75)+coord_cartesian(ylim=c(0,.5))+ggthemes::theme_few(base_size = 9)+ylab("")+
  scale_x_discrete(name="",labels=c("different invaders","same invader"))+theme(axis.text.y = element_blank(),axis.ticks.y = element_blank())

jpeg("Analyses/AbnOcc/explaination1.jpeg",width=6,height=3,units = 'in',res=200)  
ggpubr::ggarrange(ex1,ex2,widths=c(.35,.3))
dev.off()


ex4<-ggplot(env,aes(as.factor(quad)))+geom_bar(position="fill",aes(fill=as.factor(same.dom)),width = .4)+
  scale_fill_manual(name="",values=c("grey40","lightgray"),labels=c("different invaders","same invader"))+ggthemes::theme_few(base_size = 9)+
scale_x_discrete(name="",labels=c("H abn \nH occ","H abn \nD occ","D abn \nD occ","D abn \nH occ"))+theme(legend.position = "top")+
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank())+ylab("")

ex3<-ggplot(env,aes(as.factor(quad)))+geom_bar(position="fill",aes(fill=as.factor(same.all)),width = .4,color="black")+
  scale_fill_manual(name="",values=c("black","white"),labels=c("multiple datasets","same dataset"))+ggthemes::theme_few(base_size = 9)+
  scale_x_discrete(name="",labels=c("H abn \nH occ","H abn \nD occ","D abn \nD occ","D abn \nH occ"))+theme(legend.position = "top")+ylab("proportion")

jpeg("Analyses/AbnOcc/explaination1.jpeg",width=7,height=7,units = 'in',res=200) 
ggpubr::ggarrange(ex1,ex2,ex3,ex4,widths=c(.35,.3),heights=c(.3,.35),labels=c("a)","b)","c)","d)"))
dev.off()  

ggplot(env,aes(as.factor(quad)))+geom_bar(position="fill",aes(fill=as.factor(same.all)))

help(geom_bar)
env.small<-sample_frac(env,0.1)

test.mod<-brm(dif~dom.z+dat.z+ndif.z+(1|mm(invaded.plot.1,invaded.plot.2,uninvaded.plot.1,uninvaded.plot.2)),data=env.small)
test.mod.2<-brm(disagree~dom.z+dat.z+ndif.z+(1|mm(invaded.plot.1,invaded.plot.2,uninvaded.plot.1,uninvaded.plot.2)),data=env.small,family="bernoulli")


mod.2<-brm(disagree~dom.z+dat.z+ndif.z+(1|mm(invaded.plot.1,invaded.plot.2,uninvaded.plot.1,uninvaded.plot.2)),data=env,family="bernoulli")
conditional_effects(test.mod)

conditional_effects(test.mod.2,points=TRUE)









resp<-env %>% group_by(L4_KEY)%>%summarise(mean_diff=mean(dif,na.rm=TRUE),mean_init=mean(dif_nat,na.rm=TRUE))

#invs.stats<-d.ref %>%group_by(L4_KEY) %>% summarise(meanT=mean(TMIN_Dec_Feb_1981_2018_mean.tif),meanEta=mean(ETa_Jun_Aug_2003_2018.tif,na.rm=TRUE),meanTree=mean( VCF_percent_treeCover_2000_2016_mean_integer.tif,na.rm=TRUE),mean_Rn=mean(Richness_N))
invs.stats<-d.ref %>%group_by(L4_KEY) %>% summarise(meanT=mean(TMIN_Dec_Feb_1981_2018_mean.tif),meanNMDI=mean(NDMI_median_Apr1_Sept30_1985_2020_90m.tif,na.rm=TRUE),
                                                      
                                                      SG_mean=mean(SG_n_tot_M_sl2_100m.tif,na.rm=TRUE),
                                               
                                                      mean_hum=mean(gHM.tif,na.rm=TRUE),meanTree=mean(VCF_percent_treeCover_2000_2016_mean_integer.tif,na.rm=TRUE),mean_Rn=mean(Richness_N))
dsets<-d.ref %>% group_by(L4_KEY,Dataset) %>% count()
dsets<-dsets %>% group_by(L4_KEY) %>% count()
samped<-d.ref %>% group_by(L4_KEY) %>% count()
colnames(dsets)[2]<-"n_dsets"

resp<-left_join(resp,invs.stats)
resp<-left_join(resp,samped)
resp<-left_join(resp,dsets)

resp<-filter(resp,!is.na(L4_KEY))

zscore <- function(x){(x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)}

colnames(resp)

resp$init.z<-zscore(resp$mean_init)
resp$n.z<-zscore(resp$n)
resp$dset.z<-zscore(resp$n_dsets)
resp$rich.z<-zscore(resp$mean_Rn)

resp$T.z<-zscore(resp$meanT)
resp$tree.z<-zscore(resp$meanTree)
resp$soil.z<-zscore(resp$SG_mean)
resp$NMDI.z<-zscore(resp$meanNMDI)
resp$hum.z<-zscore(resp$mean_hum)

summary(lm(resp$mean_init~resp$mean_Rn))


summary(lm(mean_diff~init.z+n.z+dset.z+init.z:n.z+n.z:dset.z+n.z+init.z:dset.z,data=resp))
summary(lm(mean_diff~T.z+tree.z+soil.z+NMDI.z+hum.z,data=resp))
resp<-left_join(resp,d.key)
library(brms)
modprops<-brm(mean_diff~init.z+n.z+dset.z+init.z:n.z+n.z:dset.z+n.z+init.z:dset.z+(1|L1_KEY),warmup=3000,iter=4000,control=list(adapt_delta=.9),data=resp)
modprops2<-brm(mean_diff~init.z+n.z+dset.z+rich.z+(1|L1_KEY),warmup=3000,iter=4000,control=list(adapt_delta=.99),data=resp)
modprops3<-brm(mean_diff~init.z+n.z+dset.z+(1|L1_KEY),warmup=3000,iter=4000,control=list(adapt_delta=.99),data=resp)

bayes_R2(modprops)
bayes_R2(modprops2)
bayes_R2(modprops3)

fixef(modprops,probs = c(.05,.95))
coef(modprops,probs = c(.05,.95))

eco<-brm(mean_init~(1|L1_KEY),data=resp)
coef(eco)

cor.test(resp$tree.z,resp$T.z)
modenvs<-brm(mean_diff~T.z+tree.z+soil.z+NMDI.z+hum.z+(T.z+tree.z+soil.z+NMDI.z+hum.z|L1_KEY),warmup=3000,iter=4000,control=list(adapt_delta=.99),data=resp)

library(tidybayes)


envout<-as.data.frame(fixef(modenvs))
envout$predictor<-rownames(envout)
envout<-filter(envout,predictor!="Intercept")
ggplot(envout,aes(Estimate,reorder(predictor,-Estimate)))+
  geom_point(size=3)+geom_errorbar(aes(xmin=`Q2.5`,xmax=`Q97.5`),width=0.05)+
  geom_vline(xintercept=0)+ylab("")+ggthemes::theme_few()

coef(modprops3)

propout<-as.data.frame(fixef(modprops3))
propout$predictor<-rownames(propout)
propout<-filter(propout,predictor!="Intercept")
ggplot(propout,aes(Estimate,reorder(predictor,-Estimate)))+
  geom_point(size=3)+geom_errorbar(aes(xmin=`Q2.5`,xmax=`Q97.5`),width=0.05)+
  geom_vline(xintercept=0,linetype="dotted")+ylab("")+xlab("estimated effect")+xlim(-.05,.1)+ggthemes::theme_few()+scale_y_discrete(labels=c("mean inital abn. vs. occ. differences/region","sample size/region","number of datasets/region"))

bayes_R2(modenvs)
bayes_R2(modprops2)

initmod<-brm(mean_init~mean_Rn+(mean_Rn|L1_KEY),data=resp)
fixef(initmod)

#invs.stats<-filter(invs.stats,L1_KEY!="15  TROPICAL WET FORESTS")
invs.stats<-filter(invs.stats,!is.na(L4_KEY))
#invs.stats$meanRich<-round(invs.stats$meanRich,1)
#invs.stats$meanCov<-round(invs.stats$meanCov,1)
ggplot(invs.stats,aes(meanEta,disagree))+geom_point()

invs.stats<-left_join(fiter,invs.stats)

ggplot(invs.stats,aes(log(meanTree),Estimate))+geom_point()
ggplot(invs.stats,aes(log(meanRich),Estimate))+geom_point()

summary(lm(Estimate~meanTree+meanEta+meanT,invs.stats))
