rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())
graphics.off()

setwd("~/Documents/git/bioticHogs/")

library(brms)
library(tidyverse)
library(tidybayes)

d<-read.csv("Analyses/Master/public/bothspeciesBdata_angio.csv")


#co.tax<-filter(d,metric=="taxonomic")
#co.tax<-filter(co.tax,similarity>0)

#co.tax.mod<- brm(
# bf(similarity|mi(se)~inv_prev*set+(inv_prev*set|NA_L1NAME)+(1|regionID),
 #  phi ~inv_prev*set+(inv_prev*set|NA_L1NAME)+(1|regionID)),
#data = co.tax,
#family = Beta(),
##control=list(adapt_delta=.999),
#chains = 4, iter = 6000, warmup = 5000,
#cores = 4, seed = 1234)


d2<-filter(d,set=="all species")
d2<-filter(d2,similarity>0)

#all.mod<- brm(
 # bf(similarity|mi(se)~inv_prev*metric+(inv_prev*metric|NA_L1NAME)+(1|regionID),
  #   phi ~inv_prev*metric+(1|NA_L1NAME)+(1|regionID)),
  #d#ata = d2,
  #family = Beta(),
  #control=list(adapt_delta=.999),
  #chains = 4, iter = 6000, warmup = 5000,
  #cores = 4, seed = 1234)


taxdat2<-filter(d2,metric=="taxonomic")
taxdat2<-filter(taxdat2,similarity>0)

phydat2<-filter(d2,metric!="taxonomic")


###all species only
tax.mod<- brm(
  bf(similarity|mi(se)~inv_prev+(inv_prev|NA_L1NAME),
     phi ~inv_prev+(inv_prev|NA_L1NAME)),
  data = taxdat2,
  family = Beta(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234)



phy.mod<- brm(
  bf(similarity|mi(se)~inv_prev+(inv_prev|NA_L1NAME),
     phi ~inv_prev+(inv_prev|NA_L1NAME)),
  data = phydat2,
  family = Beta(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234)

round(fixef(tax.mod,probs = c(.25,.75,.05,.95)),2)
round(fixef(phy.mod,probs = c(.25,.75,.05,.95)),2)

phycof<-as.data.frame(coef(phy.mod,probs = c(.25,.75,.05,.95)))
taxcof<-as.data.frame(coef(tax.mod,probs = c(.25,.75,.05,.95)))

write.csv(phycof,"Analyses/Master/public/phylomodcoefs.csv")
write.csv(taxcof,"Analyses/Master/public/taxmodcoefs.csv")

order<-taxdat2 %>% group_by(NA_L1NAME) %>% 
  count() %>%
  arrange(-n)


####plots:
draws.tax <- tax.mod %>%
  spread_draws(
    b_inv_prev,
    r_NA_L1NAME[NA_L1NAME, term]
  )%>%
  filter(term %in% c("inv_prev")) %>%
  pivot_wider(names_from = term, values_from = r_NA_L1NAME, names_prefix = "r_") %>%
  mutate(
    effect_size = b_inv_prev + r_inv_prev,
  )


draws.tax$biome<-NA
draws.tax$NA_L1NAME <- gsub("\\.", " ", draws.tax$NA_L1NAME)
draws.tax$biome <- factor(draws.tax$NA_L1NAME, levels = c(rev(order$NA_L1NAME)))



draws.phy <- phy.mod %>%
  spread_draws(
    b_inv_prev,
    r_NA_L1NAME[NA_L1NAME, term]
  )%>%
  filter(term %in% c("inv_prev")) %>%
  pivot_wider(names_from = term, values_from = r_NA_L1NAME, names_prefix = "r_") %>%
  mutate(
    effect_size = b_inv_prev + r_inv_prev,
  )


draws.phy$NA_L1NAME <- gsub("\\.", " ", draws.phy$NA_L1NAME)
draws.phy$biome <- factor(draws.phy$NA_L1NAME, levels = c(rev(order$NA_L1NAME)))



draws.tax$metric<-"taxonomic"
draws.phy$metric<-"phylogenetic"
draws.biome<-rbind(draws.tax,draws.phy)



p2<-ggplot(draws.biome,aes(effect_size,biome))+stat_halfeye(aes(color=metric,shape=metric,fill=metric),.width=c(.05,.5),position = position_dodge(width = 0.5))+
  coord_cartesian(xlim=c(-1.5,1.5))+
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    x = "estimated effect (μ)",
    y = ""
  ) +
  ggthemes::theme_few()+
  scale_fill_viridis_d(name="community",begin = .3,end=.6,alpha = 0.4)+scale_shape_manual(values=c(15,19))+
  scale_color_viridis_d(name="community",begin = .3,end=.6)+theme(legend.position = "none")



fixed_draws <- tax.mod %>%
  spread_draws(
    b_inv_prev
  )

fixed_draws2 <- phy.mod %>%
  spread_draws(
    b_inv_prev
  )

fixed_draws$metric<-"taxonomic"
fixed_draws2$metric<-"phylogenetic"
fixed_draws3<-rbind(fixed_draws,fixed_draws2)



p1<-ggplot(fixed_draws3,aes(b_inv_prev,0))+stat_halfeye(aes(fill=metric,color=metric,shape=metric),.width=c(0.005,.5),position = position_dodge(width = 0.2))+
  coord_cartesian(xlim=c(-1,1.5),ylim=c(-0.1,0.4))+
  
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    x = "estimated effect (μ)",
    y = "",
  ) +
  ggthemes::theme_few()+scale_color_viridis_d(,begin = .3,end=.6)+scale_shape_manual(values=c(15,19))+
  scale_fill_viridis_d(name="metric",begin = .3,end=.6,alpha = 0.2)+theme(axis.text.y = element_blank(),axis.ticks.y=element_blank())


jpeg("Analyses/Master/Plots/angios_fig1.jpg",width = 8,height = 8, unit='in',res=300)
ggpubr::ggarrange(p1,p2,ncol=1,common.legend = TRUE,heights=c(1,2),labels=c("a)","b)"))
dev.off()



inv<-read.csv("Analyses/Master/public/invaderBdata_angio.csv")
inv$se_sim<-ifelse(inv$se_sim==0,0.0000000001,inv$se_sim)
inv$se_alph<-ifelse(inv$se_alph==0,0.0000000001,inv$se_alph)
inv<-filter(inv,!is.na(se_sim))
inv<-filter(inv,similarity>0)
inv<-filter(inv,!is.na(NA_L1NAME))

inv.tax<-filter(inv,metric=="taxonomic")
inv.phy<-filter(inv,metric!="taxonomic")




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


round(fixef(sim.mod2,probs = c(.25,.75,.05,.95)),2)
round(fixef(sim.modphy2,probs = c(.25,.75,.05,.95)),2)


invtaxcof<-as.data.frame(coef(sim.mod,probs = c(.25,.75,.05,.95)))
invphycof<-as.data.frame(coef(sim.mod.phy,probs = c(.25,.75,.05,.95)))

write.csv(invphycof,"Analyses/Master/public/inv_phylomodcoefs.csv")
write.csv(invtaxcof,"Analyses/Master/public/inv_taxmodcoefs.csv")


homog.tax<-draws.tax %>% group_by(biome) %>%
  summarise(homog=round(mean(effect_size),2)) %>% ungroup()

homog.phy<-draws.phy %>% group_by(biome) %>% 
  summarise(homog=round(mean(effect_size),2)) %>% ungroup()



new.inv<-data.frame(NA_L1NAME=unique(inv.tax$NA_L1NAME),se_sim=mean(inv.tax$se_sim))
new.phyr<-data.frame(NA_L1NAME=unique(inv.phy$NA_L1NAME),se_sim=mean(inv.phy$se_sim))

draws.sim<-epred_draws(sim.mod,newdata=new.inv,ndraws=1000)

draws.sim$biome<-NA
draws.sim$NA_L1NAME <- gsub("\\.", " ", draws.sim$NA_L1NAME)
draws.sim$biome<-draws.sim$NA_L1NAME


draws.sim.phy<-epred_draws(sim.mod.phy,newdata=new.phyr,ndraws=1000)

draws.sim.phy$biome<-NA
draws.sim.phy$NA_L1NAME <- gsub("\\.", " ", draws.sim.phy$NA_L1NAME)
draws.sim.phy$biome<-draws.sim.phy$NA_L1NAME





draws.sim$metric<-"taxonomic"
draw.sim<-left_join(draws.sim,homog.tax)


draws.sim.phy$metric<-"phylogenetic"
draw.sim.phy<-left_join(draws.sim.phy,homog.phy)


draws.invy<-bind_rows(draw.sim.phy,draw.sim)

save.image("biohug.Rda")

draws.invy<-left_join(draws.invy,htax)

jpeg("Analyses/Master/Plots/angios_fig2.jpg",width = 11,height = 6, unit='in',res=400)
ggplot(draws.invy,aes(.epred,reorder(biome,.epred)))+
  stat_halfeye(aes(fill=homog,shape=mean),.width=c(.005,.5),position = position_dodge(width = 0.5))+

  labs(
    x = "biotic similarity of\n non-native species",
    y = ""
  ) +
  ggthemes::theme_few()+
  scale_fill_distiller(name="similarity ~ invasion prevalence", type='div')+
 facet_wrap(~metric,scales="free_x")+theme(legend.position = "top")+scale_shape_binned(name="potential water deficit",n.breaks=4)
dev.off()


c<-read.csv("Data/FULLDatabase_10272022_diversity_PctCov_100_10272022_INHABIT.csv")


clims<-c %>% group_by(NA_L1NAME) %>% 
  summarise(mean=mean(PPT_ETo_Apr_Oct_1981_2018_mean.tif,na.rm=TRUE),sd=sd(PPT_ETo_Apr_Oct_1981_2018_mean.tif,na.rm=TRUE))
clims<-filter(clims,!is.na(NA_L1NAME))

colnames(clims)[1]<-'biome'
left_join(clims,htax)

htax <- homog.tax[order(homog.tax$biome), ]

htax<-left_join(clims,htax)
htax$metric<-"taxonomic"

ptax <- homog.phy[order(homog.tax$biome), ]
ptax<-left_join(clims,ptax)
ptax$metric<-"phylogenetic"

htax<-rbind(ptax,htax)

cor.test(htax$mean,htax$homog)


htax.tax<-filter(htax,metric=="taxonomic")
jpeg("Analyses/Master/Plots/pwd.jpg",width = 10,height = 5, unit='in',res=200)
ggplot(htax.tax,aes(mean,reorder(biome,-mean)))+geom_bar(stat='identity',color='black',aes(fill=homog))+
  geom_errorbarh(aes(xmin=mean-sd,xmax=mean+sd),height=.1)+
  ggthemes::theme_few()+geom_vline(xintercept=0,color='black')+ylab("")+xlab("potential water deficit \nApril-October")+
  scale_fill_distiller(name="taxonomic similarity ~ invasion prevalence", type='div')+theme(legend.position = 'top')
dev.off()



sim.mod2<- brm(
  bf(similarity|mi(se_sim)~inv_prev+(inv_prev|NA_L1NAME),
     phi ~inv_prev+(inv_prev|NA_L1NAME)),
  data = inv.tax,
  family = Beta(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234)

sim.modphy2<- brm(
  bf(similarity|mi(se_sim)~inv_prev+(inv_prev|NA_L1NAME),
     phi ~inv_prev+(inv_prev|NA_L1NAME)),
  data = inv.phy,
  family = Beta(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234)


new.inv2<-data.frame(NA_L1NAME=rep(unique(inv.tax$NA_L1NAME),each=11),se_sim=mean(inv.tax$se_sim),inv_prev=rep(c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1),10))
new.phyr2<-data.frame(NA_L1NAME=rep(unique(inv.phy$NA_L1NAME),each=11),se_sim=mean(inv.phy$se_sim),inv_prev=rep(c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1),10))

invypred<- sim.mod2 %>% 
  epred_draws(newdata =new.inv2,ndraws = 1000)

invypred2<- sim.modphy2 %>% 
  epred_draws(newdata =new.phyr2,ndraws = 1000)


invypred$biome <- factor(invypred$NA_L1NAME, levels = c(order$NA_L1NAME))
invypred2$biome <- factor(invypred2$NA_L1NAME, levels = c(order$NA_L1NAME))


pooa<-ggplot(invypred,aes(inv_prev,.epred))+stat_smooth(linewidth=0.1,alpha=1,color="black",se=TRUE,aes(fill=biome))+scale_fill_viridis_d(option = "H")+ylim(0,.65)+
  ggthemes::theme_few()+ylab("taxonomic similarity")+xlab("invasion prevalence")+guides(fill = guide_legend(ncol = 5))

prevcomp<-select(phydat2,regionID,NA_L1NAME,inv_prev)
prevcomp<-distinct(prevcomp)       

pcom<-prevcomp %>%group_by(NA_L1NAME) %>% summarise(mean=mean(inv_prev),sd=sd(inv_prev),
                                                    n = n(),
                                                    se = sd / sqrt(n))




poota<-ggplot(invypred2,aes(inv_prev,.epred))+stat_smooth(aes(fill=biome),linewidth=.1,alpha=.9,color="black")+scale_fill_viridis_d(option = "H")+ylim(.25,.9)+
  ggthemes::theme_few()+ylab("phylogenetic similarity")+xlab("invasion prevalence")


jpeg("Analyses/Master/Plots/angios_fig3.jpg",width =11,height = 5, unit='in',res=200)
ggpubr::ggarrange(poota,pooa,common.legend = TRUE,legend = "right")
dev.off()

































fixed_draws2 <- phy.mod %>%
  spread_draws(
    b_inv_prev
  )



p2<-ggplot(fixed_draws2,aes(b_inv_prev,0))+stat_halfeye(.width=c(.05,.5),position = position_dodge(width = 0.3))+
  coord_cartesian(xlim=c(-1,1.5))+scale_y_discrete()+
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    x = "fixed effects of invasion prevenlence on \nphynogenetic similarity",
    y = "",
  ) +
  ggthemes::theme_few()+scale_color_viridis_d(begin = .3,end=.6)+
  theme(axis.text.y=element_blank(),axis.ticks.y = element_blank())+scale_shape_manual(values=c(15,19))+
  scale_fill_viridis_d(name="community",begin = .3,end=.6,alpha = 0.2)


pfir<-ggpubr::ggarrange(p1,p2,widths=c(1,1),common.legend = TRUE, labels=c("a)","b)"))
psec<-ggpubr::ggarrange(p3,p4,widths=c(2.2,1), labels=c("c)","d)"))


jpeg("Analyses/Master/Plots/angios_important.jpg",width = 8,height = 6, unit='in',res=300)
ggpubr::ggarrange(pfir,psec,ncol=1,heights=c(1,2))
dev.off()






















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
inv$se_alph<-ifelse(inv$se_alph==0,0.0000000001,inv$se_alph)
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

alpa.mod<- brm(
  bf(alpha|mi(se_alph)~+(1|NA_L1NAME)),

  data = inv.tax,
  family = gaussian(),
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

invypred$biome <- factor(invypred$NA_L1NAME, levels = c(order$NA_L1NAME))
invypred2$biome <- factor(invypred2$NA_L1NAME, levels = c(order$NA_L1NAME))

invypred2<- sim.modphy2 %>% 
  epred_draws(newdata =new.phyr,ndraws = 100)

pooa<-ggplot(invypred,aes(inv_prev,.epred))+stat_smooth(aes(fill=biome),linewidth=0.1,alpha=.9,color="black")+scale_fill_viridis_d(option = "H")+ylim(0,.65)+
  ggthemes::theme_few()+ylab("taxonomic similarity")+xlab("invasion prevalence")+guides(fill = guide_legend(ncol = 3))

prevcomp<-select(phydat2,regionID,NA_L1NAME,inv_prev)
prevcomp<-distinct(prevcomp)       

pcom<-prevcomp %>%group_by(NA_L1NAME) %>% summarise(mean=mean(inv_prev),sd=sd(inv_prev),
                                                    n = n(),
                                                    se = sd / sqrt(n))




poota<-ggplot(invypred2,aes(inv_prev,.epred))+stat_smooth(aes(fill=biome),linewidth=.1,alpha=.9,color="black")+scale_fill_viridis_d(option = "H")+ylim(.25,.9)+
  ggthemes::theme_few()+ylab("phylogenetic similarity")+xlab("invasion prevalence")





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

gooby<-ggpubr::ggarrange(poo,poot,ncol=1,labels=c("a)","b)"))
gooby2<-ggpubr::ggarrange(pooa,poota,common.legend = TRUE,labels=c("c)","d)"))

jpeg("Analyses/Master/Plots/angios_invers.jpg",width = 9,height = 9, unit='in',res=300)
ggpubr::ggarrange(gooby,gooby2,ncol=1,heights=c(1.3,1))
dev.off()

jpeg("Analyses/Master/Plots/inv_levels.jpg",width = 6,height = 4, unit='in',res=200)
ggplot(pcom,aes(mean,reorder(NA_L1NAME,mean)))+geom_bar(stat='identity',fill='black',color="black")+
  geom_errorbarh(aes(xmin=mean,xmax=mean+se),height=0)+xlab("invasion prevalence")+
  ylab("")+ggthemes::theme_few()
dev.off()


p4<-ggplot(draws.phy,aes(effect_size,biome))+stat_halfeye(,.width=c(0.05,.5),position = position_dodge(width = 0.5))+
  coord_cartesian(xlim=c(-1,1.5))+
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    x = "estimated effects of\n invasion prevenlence on \nphylogenetic similarity",
    y = ""
  ) +
  ggthemes::theme_few()+scale_color_viridis_d(name="community",begin = .3,end=.6)+scale_shape_manual(values=c(15,19))+
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),legend.position = "none")+
  scale_fill_viridis_d(name="community",begin = .3,end=.6,alpha = 0.2)


ggpubr::ggarrange(p3,p4,widths=c(2,1))