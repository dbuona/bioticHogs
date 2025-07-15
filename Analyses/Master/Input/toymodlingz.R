####do the modeling

rm(list=ls())
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())



library(tidyverse)
library(broom)
library(dplyr)
library(brms)
library(tidybayes)
setwd("/Users/danielbuonaiuto/Documents/git/bioticHogs/Analyses/Master")

if(FALSE){
tax1<-read.csv("Input/L3/q1_tax.csv")
tax0<-read.csv("Input/L3/q0_tax.csv")
tax2<-read.csv("Input/L3/q2_tax.csv")

#phy1<-read.csv("Input/L3/q1_phy.csv")
#phy0<-read.csv("Input/L3/q0_phy.csv")

set.seed(111)
tax1.small<-sample_frac(tax1,.005) ##sample .05 percept of the data for modeling at home
set.seed(111)
tax0.small<-sample_frac(tax0,.005)
set.seed(111)
tax2.small<-sample_frac(tax2,.005)
rm(tax1)
rm(tax2)
rm(tax0)

tax1.small<-filter(tax1.small,local_similarity>0) ## get rid of 0 and 1s for beta dist
tax1.small<-filter(tax1.small,local_similarity<1)

tax2.small<-filter(tax2.small,local_similarity>0)
tax2.small<-filter(tax2.small,local_similarity<1)

tax0.small<-filter(tax0.small,local_similarity>0)
tax0.small<-filter(tax0.small,local_similarity<1)

write.csv(tax0.small,"Input/L3/q0_tax_small.csv")
write.csv(tax1.small,"Input/L3/q1_tax_small.csv")
write.csv(tax2.small,"Input/L3/q2_tax_small.csv")
}

####next time just read in above

toy.mod.tax1.slope<- brm(
  bf(local_similarity ~ status*logDist+(status*logDist|NA_L3NAME)+(1|mm(site1,site2)),
     phi ~ status*logDist+(status*logDist|NA_L3NAME)),
  data = tax1.small,
  family = Beta(),
  control=list(adapt_delta=.9),
  chains = 4, iter = 2000, warmup = 1000,
  cores = 4, seed = 1234,backend = "cmdstanr")

toy.mod.tax2.slope<- brm(
  bf(local_similarity ~ status*logDist+(status*logDist|NA_L3NAME)+(1|mm(site1,site2)),
     phi ~ status*logDist+(status*logDist|NA_L3NAME)),
  data = tax2.small,
  family = Beta(),
  control=list(adapt_delta=.9),
  chains = 4, iter = 2000, warmup = 1000,
  cores = 4, seed = 1234,backend = "cmdstanr") 

toy.mod.tax0.slope<- brm(
  bf(local_similarity ~ status*logDist+(status*logDist|NA_L3NAME)+(1|mm(site1,site2)),
     phi ~ status*logDist+(status*logDist|NA_L3NAME)),
  data = tax0.small,
  family = Beta(),
  #control=list(adapt_delta=.9),
  chains = 4, iter = 2000, warmup = 1000,
  cores = 4, seed = 1234,backend = "cmdstanr") ### 4 divergent transitions, probably up the adapt_delta will handle

####try seperate models with diatance and 0
tax1.small.native<-filter(tax1.small,status=="native")
tax1.small.invaded<-filter(tax1.small,status!="native")
table(tax1.small.native$NA_L3NAME)

sorock<-filter(tax1.small.native,NA_L3NAME=="Middle Rockies")

sorock$binsim<-as.integer(sorock$local_similarity)
get_prior(local_similarity ~ s(Distance)+1,
              data = sorock)
          

sorock.native.lm<- lm(local_similarity ~ Distance+1,
                    data = sorock)


lm_shift_up <- glm(local_similarity ~ Distance +0 + 
                    offset(rep(1, nrow(sorock))),family="quasibinomial", 
                  data=sorock)

summary(lm_shift_up)

intercept <- 1.0
fit <- glm(I(local_similarity - intercept) ~ 0 + Distance, sorock,family="Gamma")
summary(fit)


     
  abline(intercept, coef(fit))

sorock.native<- brm(local_similarity ~ Distance+0,
  data = sorock,
  init = 1,
  chains = 4, iter = 2000, warmup = 1000,
  cores = 4, seed = 1234,backend = "cmdstanr") 

sorock %>%
  add_epred_draws(sorock.native, 
                  ndraws=10, 
                  allow_new_levels = TRUE) %>%
  ggplot(aes(x = Distance, 
             y = local_similarity)) +
  geom_point(data = sorock, color="black",size=0.01)



fit2<-brms::brm(bf(local_similarity ~ 
                     exp(-a1 * Distance) + b1,
                   a1 + b1~1 , nl=TRUE), 
                data = sorock,
                prior = c(
                  prior(normal(0, 1), nlpar = "a1"),
                  prior(normal(0, 1), nlpar = "b1")), 
                chains=4, 
                cores=4, 
                warmup = 1000,
                iter = 2000)


m2 <- brm(bf(local_similarity ~ s(Distance)),
          data = sorock, family = gaussian(), cores = 4, seed = 17,
          iter = 4000, warmup = 1000, thin = 10, refresh = 0)






save.image("Input/L3/L3toymodels_wslope.Rda") 

##posterior predictions
quantile(tax0.small$logDist,na.rm=TRUE)
quantile(tax1.small$logDist,na.rm=TRUE)
quantile(tax2.small$logDist,na.rm=TRUE)
length(unique(tax0.small$NA_L3NAME))
length(unique(tax1.small$NA_L3NAME))
length(unique(tax2.small$NA_L3NAME))

new.data2<-data.frame(logDist=rep(-2:9,each=2),status=rep(unique(tax1.small$status),12))

new.data3<-data.frame(NA_L3NAME=rep(unique(tax1.small$NA_L3NAME),each=24))
new.data5<-data.frame(NA_L3NAME=rep(unique(tax2.small$NA_L3NAME),each=24))
new.data4<-data.frame(NA_L3NAME=rep(unique(tax0.small$NA_L3NAME),each=24))

new.data6<-data.frame(NA_L3NAME=rep(unique(phy1.small$NA_L3NAME),each=24))
new.data7<-data.frame(NA_L3NAME=rep(unique(sla1.small$NA_L3NAME),each=24))

new.data3<-cbind(new.data3,new.data2)
new.data4<-cbind(new.data4,new.data2)
new.data5<-cbind(new.data5,new.data2)
new.data6<-cbind(new.data6,new.data2)
new.data7<-cbind(new.data7,new.data2)
#### predictions
toy_pred.tax1.slope<- toy.mod.tax1.slope %>% 
  epred_draws(newdata =new.data3,ndraws = 100,re_formula = ~(status*logDist|NA_L3NAME))

toy_pred.tax0.slope<- toy.mod.tax0.slope %>% 
  epred_draws(newdata =new.data4,ndraws = 100,re_formula = ~(status*logDist|NA_L3NAME))

toy_pred.tax2.slope<- toy.mod.tax2.slope %>% 
  epred_draws(newdata =new.data5,ndraws = 100,re_formula = ~(status*logDist|NA_L3NAME))

toy_pred.phy1.slope<- toy.mod.phy1.slope %>% 
  epred_draws(newdata =new.data6,ndraws = 100,re_formula = ~(status*logDist|NA_L3NAME))

toy_pred.sla1.slope<- toy.mod.sla1.slope %>% 
  epred_draws(newdata =new.data7,ndraws = 20,re_formula = ~(status*logDist|NA_L3NAME))




toy_pred.tax1<- toy.mod.tax1.slope %>% 
  epred_draws(newdata =new.data2,ndraws = 100,re_formula = NA)

toy_pred.tax0<- toy.mod.tax0.slope %>% 
  epred_draws(newdata =new.data2,ndraws = 100,re_formula =NA)

toy_pred.tax2<- toy.mod.tax2.slope %>% 
  epred_draws(newdata =new.data2,ndraws = 100,re_formula = NA)

toy_pred.phy1<- toy.mod.phy1.slope %>% 
  epred_draws(newdata =new.data2,ndraws = 100,re_formula = NA)

toy_pred.sla1<- toy.mod.sla1.slope %>% 
  epred_draws(newdata =new.data2,ndraws = 100,re_formula = NA)

toy_pred.height1<- toy.mod.height1.slope %>% 
  epred_draws(newdata =new.data2,ndraws = 100,re_formula = NA)

##group em
toy_pred.tax1.slope$grouper<-paste(toy_pred.tax1.slope$NA_L3NAME,toy_pred.tax1.slope$status,toy_pred.tax1.slope$.draw)
toy_pred.tax0.slope$grouper<-paste(toy_pred.tax0.slope$NA_L3NAME,toy_pred.tax0.slope$status,toy_pred.tax0.slope$.draw)
toy_pred.tax2.slope$grouper<-paste(toy_pred.tax2.slope$NA_L3NAME,toy_pred.tax2.slope$status,toy_pred.tax2.slope$.draw)


toy_pred.phy1.slope$grouper<-paste(toy_pred.phy1.slope$NA_L3NAME,toy_pred.phy1.slope$status,toy_pred.phy1.slope$.draw)
toy_pred.sla1.slope$grouper<-paste(toy_pred.sla1.slope$NA_L3NAME,toy_pred.sla1.slope$status,toy_pred.sla1.slope$.draw)


##group em
toy_pred.tax1$grouper<-paste(toy_pred.tax1$status,toy_pred.tax1$.draw)
toy_pred.tax0$grouper<-paste(toy_pred.tax0$status,toy_pred.tax0$.draw)
toy_pred.tax2$grouper<-paste(toy_pred.tax2$status,toy_pred.tax2$.draw)

toy_pred.phy1$grouper<-paste(toy_pred.phy1$status,toy_pred.phy1$.draw)
toy_pred.sla1$grouper<-paste(toy_pred.sla1$status,toy_pred.sla1$.draw)
toy_pred.height1$grouper<-paste(toy_pred.height1$status,toy_pred.height1$.draw)


ecoregionz<-read.csv("Input/ecoregions.csv")
ecoregionz<-dplyr::select(ecoregionz,-X,-Plot)
ecoregionz<-distinct(ecoregionz)

toy_pred.tax1.slope<-left_join(toy_pred.tax1.slope,ecoregionz)
toy_pred.tax2.slope<-left_join(toy_pred.tax2.slope,ecoregionz)
toy_pred.tax0.slope<-left_join(toy_pred.tax0.slope,ecoregionz)
toy_pred.phy1.slope<-left_join(toy_pred.phy1.slope,ecoregionz)

###plots
pa<-ggplot()+#geom_point(data=taxIN.small,aes(x=logDist,y=local_similarity,color=status),size=.01,alpha=0.2)+
  geom_line(data=toy_pred.tax1,aes(x=logDist,y=.epred,color=status,group=grouper),size=0.01)+
  geom_smooth(data=toy_pred.tax1,method="gam", formula = y ~ s(x, bs = "cs", k=9),aes(x=logDist,y=.epred,color=status))+ylab("similarity")+xlab("log(Distance)")+
  ggthemes::theme_few(base_size = 12)+scale_color_manual("invasion status",values=c("navyblue","forestgreen"),labels=c("invasive","native"))
#coord_cartesian(xlim=c(0.3,8),ylim=c(.1,.7))

pb<-ggplot()+#geom_point(data=taxIN.small,aes(x=logDist,y=local_similarity,color=status),size=.01,alpha=0.2)+
  geom_line(data=toy_pred.tax0,aes(x=logDist,y=.epred,color=status,group=grouper),size=0.01)+
  geom_smooth(data=toy_pred.tax0,method="gam", formula = y ~ s(x, bs = "cs", k=9),aes(x=logDist,y=.epred,color=status))+ylab("similarity")+xlab("log(Distance)")+
  ggthemes::theme_few(base_size = 12)+scale_color_manual("invasion status",values=c("navyblue","forestgreen"),labels=c("invasive","native"))
#coord_cartesian(xlim=c(0.3,8),ylim=c(.1,.7))

pc<-ggplot()+#geom_point(data=taxIN.small,aes(x=logDist,y=local_similarity,color=status),size=.01,alpha=0.2)+
  geom_line(data=toy_pred.tax2,aes(x=logDist,y=.epred,color=status,group=grouper),size=0.01)+
  geom_smooth(data=toy_pred.tax2,method="gam", formula = y ~ s(x, bs = "cs", k=9),aes(x=logDist,y=.epred,color=status))+ylab("similarity")+xlab("log(Distance)")+
  ggthemes::theme_few(base_size = 12)+scale_color_manual("invasion status",values=c("navyblue","forestgreen"),labels=c("invasive","native"))
#coord_cartesian(xlim=c(0.3,8),ylim=c(.1,.7))



phyp<-ggplot()+#geom_point(data=taxIN.small,aes(x=logDist,y=local_similarity,color=status),size=.01,alpha=0.2)+
  geom_line(data=toy_pred.phy1,aes(x=logDist,y=.epred,color=status,group=grouper),size=0.01)+
  geom_smooth(data=toy_pred.phy1,method="gam", formula = y ~ s(x, bs = "cs", k=9),aes(x=logDist,y=.epred,color=status))+ylab("similarity")+xlab("log(Distance)")+
  ggthemes::theme_few(base_size = 12)+scale_color_manual("invasion status",values=c("navyblue","forestgreen"),labels=c("invasive","native"))
#coord_cartesian(xlim=c(0.3,8),ylim=c(.1,.7))

heip<-ggplot()+#geom_point(data=taxIN.small,aes(x=logDist,y=local_similarity,color=status),size=.01,alpha=0.2)+
  geom_line(data=toy_pred.height1,aes(x=logDist,y=.epred,color=status,group=grouper),size=0.01)+
  geom_smooth(data=toy_pred.height1,method="gam", formula = y ~ s(x, bs = "cs", k=9),aes(x=logDist,y=.epred,color=status))+ylab("similarity")+xlab("log(Distance)")+
  ggthemes::theme_few(base_size = 12)+scale_color_manual("invasion status",values=c("navyblue","forestgreen"),labels=c("invasive","native"))
#coord_cartesian(xlim=c(0.3,8),ylim=c(.1,.7))

slap<-ggplot()+#geom_point(data=taxIN.small,aes(x=logDist,y=local_similarity,color=status),size=.01,alpha=0.2)+
  geom_line(data=toy_pred.sla1,aes(x=logDist,y=.epred,color=status,group=grouper),size=0.01)+
  geom_smooth(data=toy_pred.sla1,method="gam", formula = y ~ s(x, bs = "cs", k=9),aes(x=logDist,y=.epred,color=status))+ylab("similarity")+xlab("log(Distance)")+
  ggthemes::theme_few(base_size = 12)+scale_color_manual("invasion status",values=c("navyblue","forestgreen"),labels=c("invasive","native"))
#coord_cartesian(xlim=c(0.3,8),ylim=c(.1,.7))


jpeg("Plots/tax_maineffects_metrics.jpeg",width=8, height=6, units = 'in', res=200)
ggpubr::ggarrange(pa,phyp, heip, slap, labels=c("tax","phy","height","sla"))
dev.off()

jpeg("Plots/tax_maineffects.jpeg",width=8, height=4, units = 'in', res=200)
ggpubr::ggarrange(pb,pa,pc,common.legend = TRUE, labels = c("q=0","q=1","q=2"),ncol=3)
dev.off()


jpeg("Plots/tax_L1s.jpeg",width=8, height=7, units = 'in', res=200)
ggplot()+#geom_point(data=taxIN.small,aes(x=logDist,y=local_similarity,color=status),size=.01,alpha=0.2)+
  geom_line(data=toy_pred.tax1.slope,aes(x=logDist,y=.epred,color=status,group=grouper),size=0.01)+
  geom_smooth(data=toy_pred.tax1.slope,method="gam", formula = y ~ s(x, bs = "cs", k=9),aes(x=logDist,y=.epred,color=status))+ylab("similarity")+xlab("log(Distance)")+
  ggthemes::theme_few(base_size = 12)+scale_color_manual("invasion status",values=c("navyblue","forestgreen"),labels=c("invasive","native"))+
  facet_wrap(~NA_L1NAME,scale="free_y")
  #coord_cartesian(xlim=c(0.3,8),ylim=c(.1,.7))
dev.off()
ggplot()+#geom_point(data=taxIN.small,aes(x=logDist,y=local_similarity,color=status),size=.01,alpha=0.2)+
  geom_line(data=toy_pred.tax0.slope,aes(x=logDist,y=.epred,color=status,group=grouper),size=0.01)+
  geom_smooth(data=toy_pred.tax0.slope,method="gam", formula = y ~ s(x, bs = "cs", k=9),aes(x=logDist,y=.epred,color=status))+ylab("similarity")+xlab("log(Distance)")+
  ggthemes::theme_few(base_size = 12)+scale_color_manual("invasion status",values=c("navyblue","forestgreen"),labels=c("invasive","native"))+
  facet_wrap(~NA_L1NAME,scales="free_y")

ggplot()+#geom_point(data=taxIN.small,aes(x=logDist,y=local_similarity,color=status),size=.01,alpha=0.2)+
  geom_line(data=toy_pred.tax2.slope,aes(x=logDist,y=.epred,color=status,group=grouper),size=0.01)+
  geom_smooth(data=toy_pred.tax2.slope,method="gam", formula = y ~ s(x, bs = "cs", k=9),aes(x=logDist,y=.epred,color=status))+ylab("similarity")+xlab("log(Distance)")+
  ggthemes::theme_few(base_size = 12)+scale_color_manual("invasion status",values=c("navyblue","forestgreen"),labels=c("invasive","native"))+
  facet_wrap(~NA_L1NAME,scales="free_y")




ggplot()+#geom_point(data=taxIN.small,aes(x=logDist,y=local_similarity,color=status),size=.01,alpha=0.2)+
  geom_line(data=toy_pred.phy1.slope,aes(x=logDist,y=.epred,color=status,group=grouper),size=0.01)+
  geom_smooth(data=toy_pred.phy1.slope,method="gam", formula = y ~ s(x, bs = "cs", k=9),aes(x=logDist,y=.epred,color=status))+ylab("similarity")+xlab("log(Distance)")+
  ggthemes::theme_few(base_size = 12)+scale_color_manual("invasion status",values=c("navyblue","forestgreen"),labels=c("invasive","native"))+
  facet_wrap(~NA_L1NAME,scale="free_y")

ggplot()+#geom_point(data=taxIN.small,aes(x=logDist,y=local_similarity,color=status),size=.01,alpha=0.2)+
  geom_line(data=toy_pred.phy1.slope,aes(x=logDist,y=.epred,color=status,group=grouper),size=0.01)+
  geom_smooth(data=toy_pred.phy1.slope,method="gam", formula = y ~ s(x, bs = "cs", k=9),aes(x=logDist,y=.epred,color=status))+ylab("similarity")+xlab("log(Distance)")+
  ggthemes::theme_few(base_size = 12)+scale_color_manual("invasion status",values=c("navyblue","forestgreen"),labels=c("invasive","native"))+
  facet_wrap(~NA_L3NAME,scale="free_y")

ggplot()+#geom_point(data=taxIN.small,aes(x=logDist,y=local_similarity,color=status),size=.01,alpha=0.2)+
  geom_line(data=toy_pred.sla1.slope,aes(x=logDist,y=.epred,color=status,group=grouper),size=0.01)+
  geom_smooth(data=toy_pred.sla1.slope,method="gam", formula = y ~ s(x, bs = "cs", k=9),aes(x=logDist,y=.epred,color=status))+ylab("similarity")+xlab("log(Distance)")+
  ggthemes::theme_few(base_size = 12)+scale_color_manual("invasion status",values=c("navyblue","forestgreen"),labels=c("invasive","native"))+
  facet_wrap(~NA_L3NAME,scale="free_y")




jpeg("Plots/tax_L3s.jpeg",width=8, height=8, units = 'in', res=200)
ggplot()+#geom_point(data=taxIN.small,aes(x=logDist,y=local_similarity,color=status),size=.01,alpha=0.2)+
  geom_line(data=toy_pred.tax1.slope,aes(x=logDist,y=.epred,color=status,group=grouper),size=0.01)+
  geom_smooth(data=toy_pred.tax1.slope,method="gam", formula = y ~ s(x, bs = "cs", k=9),aes(x=logDist,y=.epred,color=status))+ylab("similarity")+xlab("log(Distance)")+
  ggthemes::theme_few(base_size = 7)+scale_color_manual("invasion status",values=c("navyblue","forestgreen"),labels=c("invasive","native"))+
  facet_wrap(~NA_L3NAME,scale="free_y")
dev.off()


ggplot()+#geom_point(data=taxIN.small,aes(x=logDist,y=local_similarity,color=status),size=.01,alpha=0.2)+
  geom_line(data=toy_pred.tax0.slope,aes(x=logDist,y=.epred,color=status,group=grouper),size=0.01)+
  geom_smooth(data=toy_pred.tax0.slope,method="gam", formula = y ~ s(x, bs = "cs", k=9),aes(x=logDist,y=.epred,color=status))+ylab("similarity")+xlab("log(Distance)")+
  ggthemes::theme_few(base_size = 12)+scale_color_manual("invasion status",values=c("navyblue","forestgreen"),labels=c("invasive","native"))+
  facet_wrap(~NA_L3NAME,scale="free_y")


ggplot()+#geom_point(data=taxIN.small,aes(x=logDist,y=local_similarity,color=status),size=.01,alpha=0.2)+
  geom_line(data=toy_pred.tax2.slope,aes(x=logDist,y=.epred,color=status,group=grouper),size=0.01)+
  geom_smooth(data=toy_pred.tax2.slope,method="gam", formula = y ~ s(x, bs = "cs", k=9),aes(x=logDist,y=.epred,color=status))+ylab("similarity")+xlab("log(Distance)")+
  ggthemes::theme_few(base_size = 9)+scale_color_manual("invasion status",values=c("navyblue","forestgreen"),labels=c("invasive","native"))+
  facet_wrap(~NA_L3NAME,scale="free_y")

pp_check(toy.mod.tax1.slope)

#####other metrics with q=1
if(FALSE){
phy1<-read.csv("Input/L3/q1_phy.csv")
set.seed(11)
phy1.small<-sample_frac(phy1,.005)
rm(phy1)

phy1.small<-filter(phy1.small,local_similarity>0)
phy1.small<-filter(phy1.small,local_similarity<1) 


height1<-read.csv("Input/L3/q1_height.csv")
set.seed(11)
height1.small<-sample_frac(height1,.005)
rm(height1)

height1.small<-filter(height1.small,local_similarity>0)
height1.small<-filter(height1.small,local_similarity<1) #8,008,228

sla1<-read.csv("Input/L3/q1_sla.csv")

set.seed(11)
sla1.small<-sample_frac(sla1,.005)
rm(sla1)
sla1.small<-filter(sla1.small,local_similarity>0)
sla1.small<-filter(sla1.small,local_similarity<1) #7,991,404



write.csv(phy1.small,"Input/L3/q1_phy_small.csv")
write.csv(height1.small,"Input/L3/q1_height_small.csv")
write.csv(sla1.small,"Input/L3/q1_sla_small.csv")
}
toy.mod.phy1.slope<- brm(
  bf(local_similarity ~ status*logDist+(status*logDist|NA_L3NAME)+(1|mm(site1,site2)),
     phi ~ status*logDist+(status*logDist|NA_L3NAME)),
  data = phy1.small,
  family = Beta(),
  control=list(adapt_delta=.9),
  chains = 4, iter = 2000, warmup = 1000,
  cores = 4, seed = 1234,backend = "cmdstanr")

toy.mod.height1.slope<- brm(
  bf(local_similarity ~ status*logDist+(status*logDist|NA_L3NAME)+(1|mm(site1,site2)),
     phi ~ status*logDist+(status*logDist|NA_L3NAME)),
  data = height1.small,
  family = Beta(),
  control=list(adapt_delta=.9),
  chains = 4, iter = 2000, warmup = 1000,
  cores = 4, seed = 1234,backend = "cmdstanr")

toy.mod.sla1.slope<- brm(
  bf(local_similarity ~ status*logDist+(status*logDist|NA_L3NAME)+(1|mm(site1,site2)),
     phi ~ status*logDist+(status*logDist|NA_L3NAME)),
  data = sla1.small,
  family = Beta(),
  control=list(adapt_delta=.9),
  chains = 4, iter = 2000, warmup = 1000,
  cores = 4, seed = 1234,backend = "cmdstanr")


save.image("Input/L3/L3toymodels_wslope.Rda") 


if(FALSE){
  toy.mod.tax1<- brm(
    bf(local_similarity ~ status*logDist+(1|NA_L3NAME)+(1|mm(site1,site2)),
       phi ~ status*logDist+(1|NA_L3NAME)),
    data = tax1.small,
    family = Beta(),
    chains = 4, iter = 2000, warmup = 1000,
    cores = 4, seed = 1234,backend = "cmdstanr") #618 seconds
  toy_pred.tax1<- toy.mod.tax1 %>% 
    epred_draws(newdata =new.data2,ndraws = 400,re_formula = NA )
  
  toy_pred.tax1$grouper<-paste(toy_pred.tax1$status,toy_pred.tax1$.draw)
  taxp<-ggplot()+geom_line(data=toy_pred.tax1,aes(x=logDist,y=.epred,color=status,group=grouper),size=0.01)+
    geom_smooth(data=toy_pred.tax1,method="gam", formula = y ~ s(x, bs = "cs", k=9),aes(x=logDist,y=.epred,color=status))+ylab("similarity")+xlab("log(Distance)")
  
  
  
  
  
  
  toy.mod.height1<- brm(
    bf(local_similarity ~ status*logDist+(1|NA_L3NAME)+(1|mm(site1,site2)),
       phi ~ status*logDist+(1|NA_L3NAME)),
    data = height1.small,
    family = Beta(),
    chains = 4, iter = 2000, warmup = 1000,
    cores = 4, seed = 1234,backend = "cmdstanr") #618 seconds
  
  toy_pred.height1<- toy.mod.height1 %>% 
    epred_draws(newdata =new.data2,ndraws = 400,re_formula = NA )
  
  toy_pred.height1$grouper<-paste(toy_pred.height1$status,toy_pred.height1$.draw)
  
  heip<-ggplot()+geom_line(data=toy_pred.height1,aes(x=logDist,y=.epred,color=status,group=grouper),size=0.01)+
    geom_smooth(data=toy_pred.height1,method="gam", formula = y ~ s(x, bs = "cs", k=9),aes(x=logDist,y=.epred,color=status))+ylab("similarity")+xlab("log(Distance)")
  
  toy.mod.sla1<- brm(
    bf(local_similarity ~ status*logDist+(1|NA_L3NAME)+(1|mm(site1,site2)),
       phi ~ status*logDist+(1|NA_L3NAME)),
    data = sla1.small,
    family = Beta(),
    chains = 4, iter = 2000, warmup = 1000,
    cores = 4, seed = 1234,backend = "cmdstanr") #618 seconds
  
  toy_pred.sla1<- toy.mod.sla1 %>% 
    epred_draws(newdata =new.data2,ndraws = 400,re_formula = NA )
  
  toy_pred.sla1$grouper<-paste(toy_pred.sla1$status,toy_pred.sla1$.draw)
  
  slap<-ggplot()+geom_line(data=toy_pred.sla1,aes(x=logDist,y=.epred,color=status,group=grouper),size=0.01)+
    geom_smooth(data=toy_pred.sla1,method="gam", formula = y ~ s(x, bs = "cs", k=9),aes(x=logDist,y=.epred,color=status))+ylab("similarity")+xlab("log(Distance)")
  
  toy.mod.phy1<- brm(
    bf(local_similarity ~ status*logDist+(1|NA_L3NAME)+(1|mm(site1,site2)),
       phi ~ status*logDist+(1|NA_L3NAME)),
    data = phy1.small,
    family = Beta(),
    chains = 4, iter = 2000, warmup = 1000,
    cores = 4, seed = 1234,backend = "cmdstanr") #618 seconds
  
  toy_pred.phy1<- toy.mod.phy1 %>% 
    epred_draws(newdata =new.data2,ndraws = 400,re_formula = NA )
  
  toy_pred.phy1$grouper<-paste(toy_pred.phy1$status,toy_pred.phy1$.draw)
  
  phyp<-ggplot()+geom_line(data=toy_pred.phy1,aes(x=logDist,y=.epred,color=status,group=grouper),size=0.01)+
    geom_smooth(data=toy_pred.phy1,method="gam", formula = y ~ s(x, bs = "cs", k=9),aes(x=logDist,y=.epred,color=status))+ylab("similarity")+xlab("log(Distance)")
  
  save.image("Input/L3/L3toymodels.Rda")  
}

