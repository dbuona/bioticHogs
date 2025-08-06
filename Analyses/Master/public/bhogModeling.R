rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())
graphics.off()

setwd("~/Documents/git/bioticHogs/")

library(brms)
library(dplyr)

d<-read.csv("Analyses/Master/public/biohogsData.csv")
#d$n_scale<-d$n/max(d$n)
taxdat2<-filter(d,metric=="taxonomic")
taxdat2<-filter(taxdat2,similarity>0)

phydat2<-filter(d,metric!="taxonomic")



tax.mod<- brm(
  bf(similarity|mi(se)~inv_prev*set+(inv_prev*set|NA_L1NAME),
     phi ~inv_prev*set+(inv_prev*set|NA_L1NAME)),
  data = taxdat2,
  family = Beta(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234)

round(fixef(tax.mod,probs = c(.25,.75,.05,.95)),2)

phy.mod<- brm(
  bf(similarity|mi(se)~inv_prev*set+(inv_prev*set|NA_L1NAME),
     phi ~inv_prev*set+(inv_prev*set|NA_L1NAME)),
  data = phydat2,
  family = Beta(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234)


round(fixef(phy.mod,probs = c(.25,.75,.05,.95)),2)

taxdat<-filter(d,set=="all species"& metric=="taxonomic")
phydat<-filter(d,set=="all species"& metric!="taxonomic")

nattaxdat<-filter(d,set=="natives"& metric=="taxonomic")
nattaxdat<-filter(nattaxdat,similarity>0)

natphydat<-filter(d,set=="natives"& metric!="taxonomic")



q1.tax<- brm(
  bf(similarity|mi(se)~inv_prev+(inv_prev|NA_L1NAME),
     phi ~inv_prev+(inv_prev|NA_L1NAME)),
  data = taxdat,
  family = Beta(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234)


fixef(q1.tax,probs = c(.25,.75))

q1.tax.gam<- brm(
  bf(similarity|mi(se)~scale(inv_gamma)*scale(inv_prev)+(1|US_L4NAME),
     phi ~scale(inv_gamma)*scale(inv_prev)+(1|US_L4NAME)),
  data = taxdat,
  family = Beta(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234)
coef(q1.tax.gam)


q1.tax.nat<- brm(
  bf(similarity|mi(se)~inv_prev+(inv_prev|NA_L1NAME),
     phi ~inv_prev+(inv_prev|NA_L1NAME)),
  data = nattaxdat,
  family = Beta(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 4321)

q1.phy.nat<- brm(
  bf(similarity|mi(se)~inv_prev+(inv_prev|NA_L1NAME),
     phi ~inv_prev+(inv_prev|NA_L1NAME)),
  data = natphydat,
  family = Beta(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 4321)


q1.phy<- brm(
  bf(similarity|mi(se)~inv_prev+(inv_prev|NA_L1NAME),
     phi ~inv_prev+(inv_prev|NA_L1NAME)),
  data = phydat,
  family = Beta(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234)


###predictions
new.tax<-data.frame(NA_L1NAME=rep(unique(taxdat2$NA_L1NAME),each=11),se=mean(taxdat2$se),inv_prev=rep(c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1),10),set=rep(c("all species","natives"),each=110))
new.phy<-data.frame(NA_L1NAME=rep(unique(phydat2$NA_L1NAME),each=11),se=mean(phydat2$se),inv_prev=rep(c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1),10),set=rep(c("all species","natives"),each=110))
#new.phy<-data.frame(NA_L1NAME=rep(unique(phydat$NA_L1NAME),each=11),se=mean(phydat$se),inv_prev=rep(c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1),10))

phypred<- q1.phy %>% 
  epred_draws(newdata =new.phy,ndraws = 1000)

library(tidybayes)
taxpred<- tax.mod %>% 
  epred_draws(newdata =new.tax,ndraws = 1000)

phypred<- phy.mod %>% 
  epred_draws(newdata =new.tax,ndraws = 1000)


#natphypred<- q1.phy.nat %>% 
 # epred_draws(newdata =new.phy,ndraws = 1000)

#nattaxpred<- q1.tax.nat %>% 
#  epred_draws(newdata =new.tax,ndraws = 1000)


#taxpred$class<-"taxonomic"
#phypred$class<-"phylogenetic"

#nattaxpred$class<-"taxonomic"
#natphypred$class<-"phylogenetic"



taxpred$grouper<-paste(taxpred$NA_L1NAME,taxpred$set,taxpred$.draw)
phypred$grouper<-paste(taxpred$NA_L1NAME,phypred$set,phypred$.draw)

#nattaxpred$grouper<-paste(nattaxpred$NA_L1NAME,nattaxpred$class,nattaxpred$.draw)
#natphypred$grouper<-paste(natphypred$NA_L1NAME,natphypred$class,natphypred$.draw)

library(ggplot2)
p1<-ggplot()+
  geom_point(data=taxdat2,aes(x=inv_prev*100,y=similarity,shape=set),size=.1)+
  geom_line(data=taxpred,aes(x=inv_prev*100,y=.epred, group=grouper, linetype=set),size=.005,alpha=.5)+
  coord_cartesian(ylim=c(0.15,.55))+ggthemes::theme_few()+theme(legend.position = "none")+
  geom_smooth(data=taxpred,aes(x=inv_prev*100,y=.epred,linetype=set),method="glm")+ylab("taxonomic similarity")+xlab("invasion prevelence (%)")



p2<-ggplot()+
  geom_point(data=phydat2,aes(x=inv_prev*100,y=similarity,shape=set),size=.1)+
  geom_line(data=phypred,aes(x=inv_prev*100,y=.epred, group=grouper, linetype=set),size=.005,alpha=.5)+
  coord_cartesian(ylim=c(0.65,1))+ggthemes::theme_few()+theme(legend.position = "none")+
  geom_smooth(data=phypred,aes(x=inv_prev*100,y=.epred,linetype=set),method="glm")+ylab("phylogenetic similarity")+xlab("invasion prevelence (%)")





#p3<-ggplot()+
 # geom_point(data=phydat,aes(x=inv_prev*100,y=similarity),size=.1)+
  #geom_line(data=phypred,aes(x=inv_prev*100,y=.epred, group=grouper),size=.01,alpha=.5)+
  #coord_cartesian(ylim=c(0,6))+ggthemes::theme_few()+
  #geom_smooth(data=phypred,aes(x=inv_prev*100,y=.epred),method="glm",color="grey60")+ylab("phylogenetic similarity")+xlab("invasion prevelence (%)")


order<-taxdat %>% group_by(NA_L1NAME) %>% 
  count() %>%
  arrange(-n)

order$NA_L1NAME
taxdat2$biome <- factor(taxdat2$NA_L1NAME, levels = c(order$NA_L1NAME))
taxpred$biome <- factor(taxpred$NA_L1NAME, levels = c(order$NA_L1NAME))

phydat2$biome <- factor(phydat2$NA_L1NAME, levels = c(order$NA_L1NAME))
phypred$biome <- factor(phypred$NA_L1NAME, levels = c(order$NA_L1NAME))



#nattaxpred$biome <- factor(nattaxpred$NA_L1NAME, levels = c(order$NA_L1NAME))
#natphypred$biome <- factor(natphypred$NA_L1NAME, levels = c(order$NA_L1NAME))



p3<-ggplot()+
  geom_point(data=taxdat2,aes(x=inv_prev*100,y=similarity,color=biome,shape=set),size=.5)+
  geom_line(data=taxpred,aes(x=inv_prev*100,y=.epred, group=grouper,color=biome,linetype=set),size=.005,alpha=0.6)+
  coord_cartesian(ylim=c(.15,.5))+ggthemes::theme_few()+facet_wrap(~biome,nrow=1)+theme(strip.text.x = element_blank(),legend.spacing.x = unit(0.1, "cm"))+
  geom_smooth(data=taxpred,aes(x=inv_prev*100,y=.epred,linetype=set,color=biome),method="glm")+ylab("taxonomic similarity")+xlab("")+
  scale_color_viridis_d(name="Biome",option = "C")+guides(linetype = "none")+guides(shape = "none")+
  scale_x_continuous(breaks=c(50,100))


p4<-ggplot()+
  geom_point(data=phydat2,aes(x=inv_prev*100,y=similarity,color=biome,shape=set),size=.5)+
  geom_line(data=phypred,aes(x=inv_prev*100,y=.epred, group=grouper,color=biome,linetype=set),size=.005,alpha=0.6)+
  coord_cartesian(ylim=c(0.65,1))+ggthemes::theme_few()+facet_wrap(~biome,nrow=1)+theme(strip.text.x = element_blank(),legend.spacing.x = unit(0.1, "cm"))+
  geom_smooth(data=phypred,aes(x=inv_prev*100,y=.epred,linetype=set,color=biome),method="glm")+ylab("phylogenetic similarity")+xlab("invasion prevelence (%)")+
  scale_color_viridis_d(name="Biome",option = "C")+guides(linetype = "none")+guides(shape = "none")+
  scale_x_continuous(breaks=c(50,100))

p1a<-ggpubr::ggarrange(p1,p2,common.legend = TRUE,ncol=2,legend = "top")
p1b<-ggpubr::ggarrange(p3,p4,common.legend = TRUE,ncol=1,legend = "right")

jpeg("Analyses/Master/Plots/full_species_linear_pred.jpg",width = 11,height=8,unit='in',res=200)
ggpubr::ggarrange(p1a,p1b,ncol=1)
dev.off()

library(tidybayes)
library(tidyverse)
library(bayesplot)



draws <- tax.mod %>%
  spread_draws(
    b_inv_prev,
   `b_inv_prev:setnatives`,
    r_NA_L1NAME[NA_L1NAME, term]
  )%>%
  filter(term %in% c("inv_prev", "inv_prev:setnatives")) %>%
  pivot_wider(names_from = term, values_from = r_NA_L1NAME, names_prefix = "r_") %>%
  mutate(
    slope_set0 = b_inv_prev + r_inv_prev,
    slope_set1 = b_inv_prev + `b_inv_prev:setnatives` + r_inv_prev + `r_inv_prev:setnatives`
  )

plot_data <- draws %>%
  pivot_longer(cols = c(slope_set0, slope_set1),
               names_to = "set",
               names_prefix = "slope_",
               values_to = "slope") %>%
  group_by(NA_L1NAME, set) %>%
  median_qi(slope, .width = c(0.5))


plot_data$biome<-NA
plot_data$NA_L1NAME <- gsub("\\.", " ", plot_data$NA_L1NAME)

plot_data$biome <- factor(plot_data$NA_L1NAME, levels = c(rev(order$NA_L1NAME)))
plot_data$community<-ifelse(plot_data$set=="set0","all species","natives only")

p3<-ggplot(plot_data, aes(x = slope, y = biome, color = community)) +
  geom_point(size=3,position = position_dodge(width = 0.5),aes(shape=community))+
  geom_errorbarh(aes(xmin = .lower, xmax = .upper), height = 0.3,
                 position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    x = "estimated effects of\n invasion prevenlence on \ntaxonomic similarity",
    y = ""
  ) +
  ggthemes::theme_few()+scale_color_viridis_d(name="community",begin = .3,end=.6)+theme(legend.position = "none")





draws2 <- phy.mod %>%
  spread_draws(
    b_inv_prev,
    `b_inv_prev:setnatives`,
    r_NA_L1NAME[NA_L1NAME, term]
  )%>%
  filter(term %in% c("inv_prev", "inv_prev:setnatives")) %>%
  pivot_wider(names_from = term, values_from = r_NA_L1NAME, names_prefix = "r_") %>%
  mutate(
    slope_set0 = b_inv_prev + r_inv_prev,
    slope_set1 = b_inv_prev + `b_inv_prev:setnatives` + r_inv_prev + `r_inv_prev:setnatives`
  )

plot_data2 <- draws2 %>%
  pivot_longer(cols = c(slope_set0, slope_set1),
               names_to = "set",
               names_prefix = "slope_",
               values_to = "slope") %>%
  group_by(NA_L1NAME, set) %>%
  median_qi(slope, .width = c(0.5))


plot_data2$biome<-NA
plot_data2$NA_L1NAME <- gsub("\\.", " ", plot_data2$NA_L1NAME)

plot_data2$biome <- factor(plot_data2$NA_L1NAME, levels = c(rev(order$NA_L1NAME)))
plot_data2$community<-ifelse(plot_data2$set=="set0","all species","natives only")

p4<-ggplot(plot_data2, aes(x = slope, y = biome, color = community)) +
  geom_point(size=3,position = position_dodge(width = 0.5),aes(shape=community))+
  geom_errorbarh(aes(xmin = .lower, xmax = .upper), height = 0.3,
                 position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    x = "estimated effects of\n invasion prevenlence on \nphylogenetic similarity",
    y = ""
  ) +
  ggthemes::theme_few()+scale_color_viridis_d(name="community",begin = .3,end=.6)+
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),legend.position = "none")



ggpubr::ggarrange(p3,p4,common.legend = TRUE,widths=c(3,2))

fixed_draws <- tax.mod %>%
  spread_draws(
    b_inv_prev,
    `b_inv_prev:setnatives`
  )

fixed_slopes <- fixed_draws %>%
  mutate(
    slope_set0 = b_inv_prev,
    slope_set1 = b_inv_prev + `b_inv_prev:setnatives`
  ) %>%
  pivot_longer(
    cols = starts_with("slope_set"),
    names_to = "setnatives",
    names_prefix = "slope_",
    values_to = "slope"
  )


plot_data3 <- fixed_slopes %>%
  group_by(setnatives) %>%
  median_qi(slope, .width = c(0.5))

plot_data3$community<-ifelse(plot_data3$setnatives=="set0","all species","natives only")

p1<-ggplot(plot_data3, aes(x = slope, y = community,color=community,shape=community)) +
  #stat_pointinterval() +
  geom_point(size=5,position = position_dodge(width = 0.5))+
  geom_errorbarh(aes(xmin = .lower, xmax = .upper), height = 0.1)+
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    x = "fixed effects of invasion prevenlence on \ntaxonomic similarity",
    y = "",
  ) +
  ggthemes::theme_few()+scale_color_viridis_d(,begin = .3,end=.6)




fixed_draws <- phy.mod %>%
  spread_draws(
    b_inv_prev,
    `b_inv_prev:setnatives`
  )

fixed_slopes <- fixed_draws %>%
  mutate(
    slope_set0 = b_inv_prev,
    slope_set1 = b_inv_prev + `b_inv_prev:setnatives`
  ) %>%
  pivot_longer(
    cols = starts_with("slope_set"),
    names_to = "setnatives",
    names_prefix = "slope_",
    values_to = "slope"
  )


plot_data4 <- fixed_slopes %>%
  group_by(setnatives) %>%
  median_qi(slope, .width = c(0.5))

plot_data4$community<-ifelse(plot_data4$setnatives=="set0","all species","natives only")

p2<-ggplot(plot_data4, aes(x = slope, y = community,color=community,shape=community)) +
  #stat_pointinterval() +
  geom_point(size=5,position = position_dodge(width = 0.5))+
  geom_errorbarh(aes(xmin = .lower, xmax = .upper), height = 0.1)+
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    x = "fixed effects of invasion prevenlence on \nphynogenetic similarity",
    y = "",
  ) +
  ggthemes::theme_few()+scale_color_viridis_d(,begin = .3,end=.6)+
  theme(axis.text.y=element_blank(),axis.ticks.y = element_blank())

fixed_draws <- tax.mod %>%
  spread_draws(
    b_inv_prev,
    `b_inv_prev:setnatives`
  )

fixed_slopes <- fixed_draws %>%
  mutate(
    slope_set0 = b_inv_prev,
    slope_set1 = b_inv_prev + `b_inv_prev:setnatives`
  ) %>%
  pivot_longer(
    cols = starts_with("slope_set"),
    names_to = "setnatives",
    names_prefix = "slope_",
    values_to = "slope"
  )


plot_data3 <- fixed_slopes %>%
  group_by(setnatives) %>%
  median_qi(slope, .width = c(0.5))

plot_data3$community<-ifelse(plot_data3$setnatives=="set0","all species","natives only")

ggplot(plot_data3, aes(x = slope, y = community,color=community,shape=community)) +
  #stat_pointinterval() +
  geom_point(size=5,position = position_dodge(width = 0.5))+
  geom_errorbarh(aes(xmin = .lower, xmax = .upper), height = 0.1)+
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    x = "fixed effects of invasion prevenlence on \ntaxonomic similarity",
    y = "",
  ) +
  ggthemes::theme_few()+scale_color_viridis_d(,begin = .3,end=.6)


pfir<-ggpubr::ggarrange(p1,p2,widths=c(1.25,1),common.legend = TRUE, labels=c("a)","b)"))
psec<-ggpubr::ggarrange(p3,p4,widths=c(2.2,1), labels=c("c)","d)"))


jpeg("Analyses/Master/Plots/updated_important.jpg",width = 8,height = 6, unit='in',res=300)
ggpubr::ggarrange(pfir,psec,ncol=1)
dev.off()

natphypred$set<-"natives"
nattaxpred$set<-"natives"

phypred$set<-"all species"
taxpred$set<-"all species"


ptop<-ggpubr::ggarrange(p1,p3,labels = c("a)","b)"))

pbot<-ggpubr::ggarrange(p2s,p4s,nrow=2, common.legend = TRUE,legend = "right",labels = c("c)","d)"))
jpeg("Analyses/Master/Plots/imprtantfigure.jpeg",width = 10,height = 9, unit="in",res=200)
ggpubr::ggarrange(ptop,pbot,nrow=2,heights = c(1,1.5))
dev.off()

combophypred<-rbind(phypred,natphypred)
combotaxpred<-rbind(taxpred,nattaxpred)

combotaxpred$grouper2<-paste(combotaxpred$NA_L1NAME,combotaxpred$class,combotaxpred$set,combotaxpred$.draw)

combophypred$grouper2<-paste(combophypred$NA_L1NAME,combophypred$class,combophypred$set,combophypred$.draw)

  pcomtax<-ggplot()+
  geom_line(data=combotaxpred,aes(x=inv_prev*100,y=.epred, group=grouper2,color=biome,fill=biome,linetype=set),size=.01, alpha=0.3)+
  coord_cartesian(ylim=c(0,1))+ggthemes::theme_few()+facet_wrap(~biome,nrow=2)+theme(strip.text.x = element_blank(),legend.spacing.x = unit(0.1, "cm"))+
  geom_smooth(data=combotaxpred,aes(x=inv_prev*100,y=.epred,linetype = set,color=biome),method="glm",se=FALSE)+ylab("taxonomic similarity")+xlab("invader prevelance (%)")+
  scale_color_viridis_d(name="Biome",option = "H")+scale_linetype_manual(values=c("solid","dotdash"))

pcomphy<-ggplot()+
  geom_line(data=combophypred,aes(x=inv_prev*100,y=.epred, group=grouper2,color=biome,fill=biome,linetype=set),size=.01, alpha=0.3)+
  coord_cartesian(ylim=c(0,1))+ggthemes::theme_few()+facet_wrap(~biome,nrow=2)+theme(strip.text.x = element_blank(),legend.spacing.x = unit(0.1, "cm"))+
  geom_smooth(data=combophypred,aes(x=inv_prev*100,y=.epred,linetype = set,color=biome),method="glm",se=FALSE)+ylab("phylogenetic similarity")+xlab("invader prevelance (%)")+
  scale_color_viridis_d(name="Biome",option = "H")+scale_linetype_manual(values=c("solid","dotdash"))

jpeg("Analyses/Master/Plots/allVnat.jpeg",width = 10,height = 9, unit="in",res=200)
ggpubr::ggarrange(pcomtax,pcomphy,nrow=2, common.legend = TRUE,legend = "right",labels = c("a)","b)"))
dev.off()


####explainers
gammamod<-brm(
  bf(inv_gamma ~ (1|NA_L1NAME)),
  data = taxdat,
  family =  gaussian(),
  control=list(adapt_delta=.95), chains = 4, iter = 4000, warmup = 3000, seed = 1234)


newgam<-data.frame(NA_L1NAME=unique(taxdat$NA_L1NAME))
library(tidybayes)
gampred<-epred_draws(gammamod,newdata=newgam)

ggplot(gampred,aes(.epred,NA_L1NAME))+stat_halfeye()

coef(gammamod)



taxdat$seI_alpha<-ifelse(taxdat$seI_alpha==0,0.0000000001,taxdat$seI_alpha)


alphamod<-brm(
  bf(meanI_alpha|mi(seI_alpha)~(1|NA_L1NAME)),
  data = taxdat,
  family =  gaussian(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 4000, warmup = 3000,cores = 4, seed = 1234)
