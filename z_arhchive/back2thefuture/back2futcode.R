###temporal

####see highplainsdrifter.R for most uptodat version of these loops.

rm(list=ls())
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#graphics.off()

library(dplyr)
library(vegan)
library("geodist")
library(ggplot2)
library(hillR)
library(betapart)
library(reshape2)

load("~/Documents/git/bioticHogs/back2thefuture/back2thefuturemods.Rda")
setwd("~/Desktop/Powell/")# for mac
set.seed(3)

###read in data
d1<-read.csv("Data/FULLDatabase_05272022.csv")
d.env1<-read.csv("Data/FULLDatabase_05272022_plotsenv4Jul2022.csv")
d.ref<-read.csv("Data/FullDatabase_diversity_July4.csv")

d.env.time<-filter(d.env1,Resampled=="Y")
counters<- d.env.time%>% group_by(US_L3NAME) %>% count() %>% arrange(n)
hundo<-filter(counters,n>50)

d.env.time<-filter(d.env.time, US_L3NAME %in% c(hundo$US_L3NAME))

d<-filter(d1, Plot %in% unique(d.env.time$Plot)) 

d.min<- d%>% group_by(Plot) %>% 
  filter(Year==min(Year)) %>%
  ungroup() 



d.max<- d%>% group_by(Plot) %>% 
  filter(Year==max(Year)) %>%
  ungroup() 




d.mino<-select(d.min,Plot, Year)
d.mino<-distinct(d.mino)

d.maxo<-select(d.max,Plot, Year)
d.maxo<-distinct(d.maxo)

ggplot()+geom_histogram(data=d.maxo,bins=12,aes(Year),fill="firebrick1",alpha=0.5)+
  geom_histogram(data=d.mino,bins=12,aes(Year),fill="royalblue1",alpha=0.5)


ggplot()+geom_density(data=d.mino,aes(Year))+
  geom_density(data=d.maxo,aes(Year),color="firebrick")
colnames(d.maxo)[2]<-"Year.max"
maxomino<-left_join(d.maxo,d.mino)
maxomino$interval<-maxomino$Year.max-maxomino$Year
table(maxomino$interval)
a<-ggplot()+geom_histogram(data=maxomino,bins=20,aes(interval),fill="firebrick")+
  ggthemes::theme_few()+xlab("Return Interval")+ylab("frequency")+
  annotate("text",x = 15,y=380,label="n= 1460")




matricize<-function(x){
  temp<-dplyr::select(x,Plot,AcceptedTaxonName,PctCov_100) ### select useful columns
  temp<-dplyr::filter(temp,!is.na(temp$AcceptedTaxonName))
  temp<- tidyr::spread(temp,AcceptedTaxonName,PctCov_100)
  temp[is.na(temp)] <- 0 
  temp<-as.data.frame(temp)
  rownames(temp)<-temp$Plot # convert plot names to row names
  temp<-dplyr::select(temp,-Plot)
  temp<-as.matrix.data.frame(temp)# remove plot column
}

mydat<-data.frame(site1=character(),site2=character(),local_similarity=numeric(),Ecoregion=character())
ecoregionz<-c(unique(d.env.time$US_L3NAME))

for (z in c(1:length(ecoregionz))){
  
  tlands<-dplyr::filter(d.env.time,US_L3NAME==ecoregionz[z]) ##
  dat<-filter(d.min, Plot %in% c(unique(tlands$Plot)))
  comms<-matricize(dat)
 
daty<-hill_taxa_parti_pairwise(comms,q = 2)
daty<-dplyr::select(daty,site1,site2,local_similarity)
daty$Ecoregion<-ecoregionz[z]
mydat<-rbind(daty,mydat)
print(paste("done",which(ecoregionz==ecoregionz[z]))) 
}


mydat2<-data.frame(site1=character(),site2=character(),local_similarity=numeric(),Ecoregion=character())
for (z in c(1:length(ecoregionz))){
  
  tlands<-dplyr::filter(d.env.time,US_L3NAME==ecoregionz[z]) ##
  dat<-filter(d.max, Plot %in% c(unique(tlands$Plot)))
  comms<-matricize(dat)
  
  daty<-hill_taxa_parti_pairwise(comms,q = 2)
  daty<-dplyr::select(daty,site1,site2,local_similarity)
  daty$Ecoregion<-ecoregionz[z]
  mydat2<-rbind(daty,mydat2)
  print(paste("done",which(ecoregionz==ecoregionz[z]))) 
}

colnames(mydat)[3]<-"early"
colnames(mydat2)[3]<-"later"

mydat3<-left_join(mydat,mydat2)

ggplot(mydat3,aes(early,later))+geom_point(size=.1)+geom_abline(intercept=0,slope=1,color="royalblue")+facet_wrap(~Ecoregion)


mydat3$H<-log((mydat3$later+.001)/(mydat3$early+.001)) ##(positive log resp rations indicate homgenization)
mydat4<-filter(mydat3,H!=0)



###color change in invasion
d.ref.time<-filter(d.ref,Plot %in% c(unique(d.env.time$Plot)))

d.ref.time<-select(d.ref.time,Plot,Year,RelCov_I)


d.ref.min<- d.ref.time%>% group_by(Plot) %>% 
  filter(Year==min(Year)) %>%
  ungroup() 



d.ref.max<- d.ref.time%>% group_by(Plot) %>% 
  filter(Year==max(Year)) %>%
  ungroup() 

d.ref.min<-select(d.ref.min,Plot,RelCov_I)
d.ref.max<-select(d.ref.max,Plot,RelCov_I)

colnames(d.ref.max)[2]<-"RelCov_I_late"
I_maxomino<-left_join(d.ref.max,d.ref.min)

I_maxomino$I_change<-I_maxomino$RelCov_I_late-I_maxomino$RelCov_I
I_maxomino$direction<-NA
I_maxomino$direction<-ifelse(I_maxomino$I_change>0,"increase","no change")
I_maxomino$direction<-ifelse(I_maxomino$I_change<0,"decrease",I_maxomino$direction)


table(I_maxomino$direction)

b<-ggplot()+geom_histogram(data=I_maxomino,aes(I_change),bins=50)+
ggthemes::theme_few()+ylab("frequency")+xlab("change in invasive cover")+annotate("text",x=50,y=400,label=c("Decrease=552 \nIncrease= 726 \nno change= 182"))
 
?annotate() 


###plot
jpeg("~/Documents/git/bioticHogs/back2thefuture/data_params.jpeg",width=6, height=4, unit= "in", res=200)
ggpubr::ggarrange(a,b)
dev.off()

jpeg("~/Documents/git/bioticHogs/back2thefuture/onetoones.jpeg",width=6, height=6, unit= "in", res=200)
ggplot(mydat3,aes(early,later))+geom_point(size=0.2)+geom_abline(intercept=0,slope=1,color="royalblue")+facet_wrap(~Ecoregion)+
  ggthemes::theme_few()+ylab("Latest Survey")+xlab("Earlier Survey")
dev.off()


ggplot(mydat3,aes(H))+geom_histogram(bins=60)+facet_wrap(~Ecoregion,scales="free")+geom_vline(xintercept = 0,color="royalblue")+
  ggthemes::theme_few()+xlab("Homoginization Index")


#####statistics
inv<-I_maxomino
inv<-select(inv,Plot,I_change)
colnames(inv)<-c("site1","I_change1")
inv2<-inv
colnames(inv2)<-c("site2","I_change2")
mydat3<-left_join(mydat3,inv)
mydat3<-left_join(mydat3,inv2)



save.image("~/Documents/git/bioticHogs/back2thefuture/back2thefuturemods.Rda")


fixef(mod1)



coef(mod1)

coords<-select(d,Plot,Lat,Long)
coords<-distinct(coords)
colnames(coords)<-c("site1", "Lat.1", "Lon.1")
coords2<-coords
colnames(coords2)<-c("site2", "Lat.2", "Lon.2")

mydat3<-left_join(mydat3,coords)
mydat3<-left_join(mydat3,coords2)

cord1<-data.frame(Lon=mydat3$Lon.1,Lat=mydat3$Lat.1)
cord2<-data.frame(Lon=mydat3$Lon.2,Lat=mydat3$Lat.2)

disty1<-geodist(x=cord1,y=cord2,measure="haversine",paired = TRUE)
mydat3$Distance<-disty1/100

jpeg("~/Documents/git/bioticHogs/back2thefuture/distancedecay.jpeg",width=6, height=8, unit= "in", res=200)
ggplot()+
  geom_smooth(data=mydat3,aes(Distance,early),method="glm",method.args = list(family = "quasibinomial"),size=1,color="navyblue",fill="navyblue",alpha=0.2)+
  geom_smooth(data=mydat3,aes(Distance,later),method="glm",method.args = list(family = "quasibinomial"),size=1,color="firebrick",fill="firebrick",alpha=0.2)+
  facet_wrap(~Ecoregion,scales="free")+ylab("Taxonomic Similarity")+ggthemes::theme_few()
dev.off()




jpeg("~/Documents/git/bioticHogs/back2thefuture/relcov_H.jpeg",width=8, height=6, unit= "in", res=200)
ggplot(mydat3,aes(I_change1,H))+geom_point(size=.1)+facet_wrap(~Ecoregion,scales="free")+geom_smooth(method="lm")+
  ggthemes::theme_few()+ylab("Homoginization Index")+xlab("Change in relative invasive cover of site 1")
dev.off()

incre<-filter(mydat3,I_change1>=0)

ggplot(incre,aes(I_change1,H))+geom_point(size=.1)+facet_wrap(~Ecoregion,scales="free")+geom_smooth(method="lm")+
  ggthemes::theme_few()+ylab("Homoginization Index")+xlab("Change in relative invasive cover of site 1")


library(brms)
mod1<-brm(H~I_change1+I_change1:I_change2,data=incre)
fixef(mod1)
mod1a<-brm(H~(1|Ecoregion),data=mydat3)
summary(mod1a)
coef(mod1a)

library(tidybayes)
get_variables(mod1a)

yaya<-mod1a%>%
  spread_draws(r_Ecoregion[Ecoregion,H])

yaya2<-mod1a%>%
  spread_draws(b_Intercept)

yaya$r_Ecoregion2<-yaya2$b_Intercept+yaya$r_Ecoregion

jpeg("~/Documents/git/bioticHogs/back2thefuture/ecoregion_mus.jpeg",width=8, height=6, unit= "in", res=200)
ggplot()+ stat_pointinterval(data=yaya,aes(r_Ecoregion2,Ecoregion),fill="skyblue1",.width = c(.5,.9))+
  geom_vline(xintercept=0,linetype="dotted")+xlab("Homoginization Index")+ggthemes::theme_few()
dev.off()

library(maps)
jpeg("~/Documents/git/bioticHogs/back2thefuture/map.jpeg",width=8, height=8, unit= "in", res=200)
map("state", interior = FALSE)
points(cord1$Lon, cord1$Lat, pch=19, cex=0.1,col="purple")
dev.off()
