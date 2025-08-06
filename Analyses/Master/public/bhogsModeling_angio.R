rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())
graphics.off()

setwd("~/Documents/git/bioticHogs/")

library(brms)
library(dplyr)

d<-read.csv("Analyses/Master/public/bothspeciesBdata_angio.csv")

tr<-read.tree("Analyses/Master/Input/angiosperms.tre")
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
round(fixef(tax.mod,probs = c(.25,.75,.05,.95)),2)

coef(phy.mod,probs = c(.25,.75,.05,.95))
coef(tax.mod,probs = c(.25,.75,.05,.95))


order<-taxdat2 %>% group_by(NA_L1NAME) %>% 
  count() %>%
  arrange(-n)


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


draws<-tidyr::gather(draws,'seter','estimate',9:10)
draws$biome<-NA
draws$NA_L1NAME <- gsub("\\.", " ", draws$NA_L1NAME)
draws$biome <- factor(draws$NA_L1NAME, levels = c(rev(order$NA_L1NAME)))
draws$community<-ifelse(draws$seter=="slope_set0","all species","natives only")


p3<-ggplot(draws,aes(estimate,biome))+stat_eye(aes(fill=community,color=community,shape=community),.width=c(.005,.5),position = position_dodge(width = 0.5))+
  coord_cartesian(xlim=c(-1.5,1.5))+
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    x = "estimated effects of\n invasion prevenlence on \ntaxonomic similarity",
    y = ""
  ) +
  ggthemes::theme_few()+
  scale_fill_viridis_d(name="community",begin = .3,end=.6,alpha = 0.2)+scale_shape_manual(values=c(15,19))+
  scale_color_viridis_d(name="community",begin = .3,end=.6)+theme(legend.position = "none")




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


draws2<-tidyr::gather(draws2,'seter','estimate',9:10)
draws2$biome<-NA
draws2$NA_L1NAME <- gsub("\\.", " ", draws2$NA_L1NAME)

draws2$biome <- factor(draws2$NA_L1NAME, levels = c(rev(order$NA_L1NAME)))
draws2$community<-ifelse(draws2$seter=="slope_set0","all species","natives only")




p4<-ggplot(draws2,aes(estimate,biome))+stat_eye(aes(fill=community,color=community,shape=community),.width=c(0.005,.5),position = position_dodge(width = 0.5))+
  coord_cartesian(xlim=c(-1.5,1.5))+
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    x = "estimated effects of\n invasion prevenlence on \nphylogenetic similarity",
    y = ""
  ) +
  ggthemes::theme_few()+scale_color_viridis_d(name="community",begin = .3,end=.6)+scale_shape_manual(values=c(15,19))+
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),legend.position = "none")+
  scale_fill_viridis_d(name="community",begin = .3,end=.6,alpha = 0.2)



ggpubr::ggarrange(p3,p4,common.legend = TRUE,widths=c(3,2))

fixed_draws <- tax.mod %>%
  spread_draws(
    b_inv_prev,
    `b_inv_prev:setnatives`
  )


fixed_draws<-tidyr::gather(fixed_draws,'seter','estimate',4:5)
fixed_draws$community<-ifelse(fixed_draws$seter=="b_inv_prev","all species","natives only")




p1<-ggplot(fixed_draws,aes(estimate,community))+stat_eye(aes(fill=community,color=community,shape=community),.width=c(0.005,.5),position = position_dodge(width = 0.3))+
  coord_cartesian(xlim=c(-1,1))+

  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    x = "fixed effects of invasion prevenlence on \ntaxonomic similarity",
    y = "",
  ) +
  ggthemes::theme_few()+scale_color_viridis_d(,begin = .3,end=.6)+scale_shape_manual(values=c(15,19))+
  scale_fill_viridis_d(name="community",begin = .3,end=.6,alpha = 0.2)




fixed_draws2 <- phy.mod %>%
  spread_draws(
    b_inv_prev,
    `b_inv_prev:setnatives`
  )


fixed_draws2<-tidyr::gather(fixed_draws2,'seter','estimate',4:5)
fixed_draws2$community<-ifelse(fixed_draws2$seter=="b_inv_prev","all species","natives only")


p2<-ggplot(fixed_draws2,aes(estimate,community))+stat_eye(aes(fill=community,color=community,shape=community),.width=c(.05,.5),position = position_dodge(width = 0.3))+
  coord_cartesian(xlim=c(-1,1))+
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    x = "fixed effects of invasion prevenlence on \nphynogenetic similarity",
    y = "",
  ) +
  ggthemes::theme_few()+scale_color_viridis_d(begin = .3,end=.6)+
  theme(axis.text.y=element_blank(),axis.ticks.y = element_blank())+scale_shape_manual(values=c(15,19))+
  scale_fill_viridis_d(name="community",begin = .3,end=.6,alpha = 0.2)


pfir<-ggpubr::ggarrange(p1,p2,widths=c(1.25,1),common.legend = TRUE, labels=c("a)","b)"))
psec<-ggpubr::ggarrange(p3,p4,widths=c(2.2,1), labels=c("c)","d)"))


jpeg("Analyses/Master/Plots/angios_important.jpg",width = 8,height = 6, unit='in',res=300)
ggpubr::ggarrange(pfir,psec,ncol=1,heights=c(1,2))
dev.off()


inv<-read.csv("Analyses/Master/public/invaderBdata_angio.csv")
inv$se_sim<-ifelse(inv$se_sim==0,0.0000000001,inv$se_sim)
inv<-filter(inv,!is.na(se_sim))
inv<-filter(inv,similarity>0)
inv<-filter(inv,!is.na(NA_L1NAME))

inv.tax<-filter(inv,metric=="taxonomic")
inv.phy<-filter(inv,metric!="taxonomic")


gammamod<-brm(
  bf(gamma ~ (1|NA_L1NAME)),
  data = inv.tax,
  family =  gaussian(),
  control=list(adapt_delta=.95), chains = 4, iter = 4000, warmup = 3000, seed = 1234)



sim.mod<- brm(
  bf(similarity|mi(se_sim)~+(1|NA_L1NAME),
     phi ~+(1|NA_L1NAME)),
  data = inv.tax,
  family = Beta(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234)


sim.mod.phy<- brm(
  bf(similarity|mi(se_sim)~+(1|NA_L1NAME),
     phi ~+(1|NA_L1NAME)),
  data = inv.phy,
  family = Beta(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234)


sim.mod2<- brm(
  bf(similarity|mi(se_sim)~inv_prev+(inv_prev|NA_L1NAME),
     phi ~inv_prev+(inv_prev|NA_L1NAME)),
  data = inv.tax,
  family = Beta(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234)

new.inv<-data.frame(NA_L1NAME=rep(unique(inv.tax$NA_L1NAME),each=11),se_sim=mean(inv.tax$se_sim),inv_prev=rep(c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1),10))
new.phyr<-data.frame(NA_L1NAME=rep(unique(inv.phy$NA_L1NAME),each=11),se_sim=mean(inv.phy$se_sim),inv_prev=rep(c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1),10))

invypred<- sim.mod2 %>% 
  epred_draws(newdata =new.inv,ndraws = 100)

invypred2<- sim.modphy2 %>% 
  epred_draws(newdata =new.phyr,ndraws = 100)

ggplot(invypred,aes(inv_prev,.epred))+stat_smooth(aes(fill=NA_L1NAME),linewidth=0.01,alpha=.9)+scale_fill_viridis_d(option = "H")+ylim(0,.75)+
  ggthemes::theme_few()
ggplot(invypred2,aes(inv_prev,.epred))+stat_smooth(aes(fill=NA_L1NAME),linewidth=0.01,alpha=.9)+scale_fill_viridis_d(option = "H")+ylim(.25,1)+
  ggthemes::theme_few()


sim.modphy2<- brm(
  bf(similarity|mi(se_sim)~inv_prev+(inv_prev|NA_L1NAME),
     phi ~inv_prev+(inv_prev|NA_L1NAME)),
  data = inv.phy,
  family = Beta(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234)

coef(sim.mod2)
coef(sim.modphy2)


draws.sim <- sim.mod %>%
  spread_draws(
    r_NA_L1NAME[NA_L1NAME, term])%>%

  pivot_wider(names_from = term, values_from = r_NA_L1NAME, names_prefix = "r_") 



draws.sim$biome<-NA
draws.sim$NA_L1NAME <- gsub("\\.", " ", draws.sim$NA_L1NAME)
draws.sim$biome<-draws.sim$NA_L1NAME


poo<-ggplot(draws.sim, aes(x = r_Intercept, y = reorder(biome,r_Intercept))) + stat_eye(.width=c(0.05,.5))+
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    x = "non-native \ntaxonomic similarity",
    y = ""
  ) +
  ggthemes::theme_few()+scale_color_viridis_d(name="community",begin = .3,end=.6)+theme(legend.position = "none")




draws.sim2 <- sim.mod.phy %>%
  spread_draws(
    r_NA_L1NAME[NA_L1NAME, term])%>%
  
  pivot_wider(names_from = term, values_from = r_NA_L1NAME, names_prefix = "r_") 



draws.sim2$biome<-NA
draws.sim2$NA_L1NAME <- gsub("\\.", " ", draws.sim2$NA_L1NAME)
draws.sim2$biome<-draws.sim2$NA_L1NAME


poot<-ggplot(draws.sim2, aes(x = r_Intercept, y = reorder(biome,r_Intercept))) + stat_eye(.width=c(0.05,.5))+
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    x = "non-native \nphylogenetic similarity",
    y = ""
  ) +
  ggthemes::theme_few()+scale_color_viridis_d(name="community",begin = .3,end=.6)+theme(legend.position = "none")


jpeg("Analyses/Master/Plots/angios_invers.jpg",width = 8,height = 8, unit='in',res=300)
ggpubr::ggarrange(poo,poot,ncol=1,labels=c("a)","b)"))
dev.off()


library(ggtree)
dat<-read.csv("Data/SPCIS_plant_taxa.csv")
dat<-dplyr::select(dat,SpCode,NativeStatus)
dat<-filter(dat,SpCode %in% c(tr$tip.label))



p<-ggtree(tr) %<+% dat

dat$status<-ifelse(dat$NativeStatus=="I","non-native","native")
p2 <- p + geom_tippoint(aes(shape = status, color = status),size=0.1)



ggplot(invypred)
