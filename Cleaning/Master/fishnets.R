rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())

library("hkdatasets")
library(sf)
library(dplyr)
library("mapview")
library("tmap")
library(maps)
library(usmap)
library(phytools)


setwd("~/Documents/git/bioticHogs/")
#load("fishnets.Rda")
d<-read.csv("Data/SPCIS_plant_taxa.csv") ### this is the public data for reproducibility
colnames(d)
p<-read.csv("Data/SPCIS_plots.csv")
regionz<-read.csv("Data/FULLDatabase_10272022_diversity_PctCov_100_10272022_INHABIT.csv")
regionz<-dplyr::select(regionz,Plot,NA_L1NAME,NA_L2NAME,NA_L3NAME,US_L4NAME)
regionz<-distinct(regionz)


d<-left_join(d,regionz)
d<-d %>% #select only the most recent survey for resampled plots
  group_by(Plot) %>%
  filter(Year==max(Year)) %>%
  ungroup() 



p<-left_join(p,regionz)

p<-p %>% #select only the most recent survey for resampled plots
  group_by(Plot) %>%
  filter(Year==max(Year)) %>%
  ungroup() 










p<-filter(p,!is.na(Long))
p<-filter(p,!is.na(Lat))


if(FALSE){

pts<-data.frame(lon=p$Long,lat=p$Lat)

eq_transformed <- usmap_transform(pts)
colnames(eq_transformed)[1:2]<-c("Long","Lat")
class(eq_transformed)

usmap_crs()
my.sf.point <- st_as_sf(x = eq_transformed, 
                        coords = c("x", "y"),
                        crs = "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +ellps=sphere +units=m +no_defs")



area_fishnet_grid_smallest = st_make_grid(my.sf.point,n = c(80,80),square = FALSE)

fishnet_grid_smallest = st_sf(area_fishnet_grid_smallest) %>%
  # add grid ID
  mutate(grid_id = 1:length(lengths(area_fishnet_grid_smallest)))

fishnet_grid_smallest$nplots<-lengths(st_intersects(fishnet_grid_smallest, my.sf.point))
fishnet_grid_smallest = filter(fishnet_grid_smallest, nplots > 9)


area_fishnet_grid_small = st_make_grid(my.sf.point,n = c(40,40),square=FALSE)

fishnet_grid_small = st_sf(area_fishnet_grid_small) %>%
  # add grid ID
  mutate(grid_id = 1:length(lengths(area_fishnet_grid_small)))

fishnet_grid_small$nplots<-lengths(st_intersects(fishnet_grid_small, my.sf.point))
fishnet_grid_small = filter(fishnet_grid_small, nplots > 19)

area_fishnet_grid_big = st_make_grid(my.sf.point,n = c(20,20),square=FALSE)

fishnet_grid_big = st_sf(area_fishnet_grid_big) %>%
  # add grid ID
  mutate(grid_id = 1:length(lengths(area_fishnet_grid_big)))

fishnet_grid_big$nplots<-lengths(st_intersects(fishnet_grid_big, my.sf.point))
fishnet_grid_big = filter(fishnet_grid_big, nplots > 39)


my.sf.point$smallest_grid<-st_intersects(my.sf.point,fishnet_grid_smallest)
my.sf.point$small_grid<-st_intersects(my.sf.point,fishnet_grid_small)
my.sf.point$big_grid<-st_intersects(my.sf.point,fishnet_grid_big)



p<-left_join(p,my.sf.point)













tmap_mode("view")

map_fishnet_smallest = tm_shape(fishnet_grid_smallest) +
  tm_fill(
    col = "nplots",
    palette = "Reds",
    style = "cont",
    title = "Plots",
    id = "grid_id",
    showNA = FALSE,
    alpha = 0.5,
    popup.vars = c(
      "Plots: " = "nplots"
    ),
    popup.format = list(
      nplots = list(format = "f", digits = 0)
    )
  ) +
  tm_borders(col = "grey40", lwd = 0.7)

map_fishnet_smallest



map_fishnet_small = tm_shape(fishnet_grid_small) +
  tm_fill(
    col = "nplots",
    palette = "Reds",
    style = "cont",
    title = "Plots",
    id = "grid_id",
    showNA = FALSE,
    alpha = 0.5,
    popup.vars = c(
      "Plots: " = "nplots"
    ),
    popup.format = list(
      nplots = list(format = "f", digits = 0)
    )
  ) +
  tm_borders(col = "grey40", lwd = 0.7)

map_fishnet_big = tm_shape(fishnet_grid_big) +
  tm_fill(
    col = "nplots",
    palette = "Reds",
    style = "cont",
    title = "Plots",
    id = "grid_id",
    showNA = FALSE,
    alpha = 0.5,
    popup.vars = c(
      "Plots: " = "nplots"
    ),
    popup.format = list(
      nplots = list(format = "f", digits = 0)
    )
  ) +
  tm_borders(col = "grey40", lwd = 0.7)

map_fishnet_smallest
map_fishnet_small
map_fishnet_big

gridmerge<-dplyr::select(p, Plot, smallest_grid,small_grid,big_grid)


dmer<-left_join(gridmerge,d)

#### gids
grid80tots<-p%>%group_by(smallest_grid) %>% count()

###invaded
grid80inv<-dmer %>% group_by(Plot,NativeStatus,smallest_grid)%>% summarize(RelCov=sum(PctCov_100,na.rm=TRUE))
grid80inv<-filter(grid80inv,NativeStatus=="I")
grid80inv<-filter(grid80inv,RelCov>5)
grid80inv<- grid80inv %>% group_by(smallest_grid) %>% count() ## plot where relcov I is >5

colnames(grid80inv)[2]<-"n_inv"
grid80tots<-left_join(grid80tots,grid80inv)
grid80tots$n_inv<-ifelse(is.na(grid80tots$n_inv),0,grid80tots$n_inv)
grid80tots<-filter(grid80tots,smallest_grid!=`integer(0)`)
grid80tots$smallest_grid<-as.character(grid80tots$smallest_grid)
class(grid80tots$smallest_grid)
grid80tots<-filter(grid80tots,smallest_grid!="integer(0)")
#L4tots<-filter(L4tots,n>19) ###must have at least 20 total plots
matricize<-function(x){
  temp<-dplyr::select(x,Plot,SpCode,PctCov_100) ### select useful columns
  temp<-dplyr::filter(temp,!is.na(temp$SpCode))
  temp<- tidyr::spread(temp,SpCode,PctCov_100)
  temp[is.na(temp)] <- 0 
  temp<-as.data.frame(temp)
  rownames(temp)<-temp$Plot # convert plot names to row names
  temp<-dplyr::select(temp,-Plot)
  temp<-as.matrix.data.frame(temp)#

}

dmer$smallest_grid<-as.character(dmer$smallest_grid)
class(dmer$smallest_grid)

library(hillR)
ecoregionz<-unique(grid80tots$smallest_grid)


taxdat<-data.frame()
for (z in c(1:length(ecoregionz))){
  dat<-dplyr::filter(dmer,smallest_grid==ecoregionz[z])
  
  dat<-filter(dat,!is.na(SpCode))
  dat<-dplyr::select(dat,Plot,SpCode,PctCov_100,smallest_grid)
  dat<-distinct(dat)
  comms<-matricize(dat)
  
  taxy0<-hill_taxa_parti(comms,q = 0)
  taxy1<-hill_taxa_parti(comms,q = 1)
  taxy2<-hill_taxa_parti(comms,q = 2)
  taxy<-rbind(taxy0,taxy1,taxy2)
  taxy$grid80<-ecoregionz[z]
  taxdat<-rbind(taxy,taxdat)
print("done")
  }  

colnames(grid80tots)[1]<-"grid80"
taxdat<-left_join(taxdat,grid80tots) 
taxdat$perc_invaded<-taxdat$n_inv/taxdat$n

ggplot(taxdat,aes(perc_invaded,local_similarity))+geom_point()+geom_smooth(method="lm")+facet_wrap(~q)

#40

grid40tots<-p%>%group_by(small_grid) %>% count()

###invaded
grid40inv<-dmer %>% group_by(Plot,NativeStatus,small_grid)%>% summarize(RelCov=sum(PctCov_100,na.rm=TRUE))
grid40inv<-filter(grid40inv,NativeStatus=="I")
grid40inv<-filter(grid40inv,RelCov>5)
grid40inv<- grid40inv %>% group_by(small_grid) %>% count() ## plot where relcov I is >5

colnames(grid40inv)[2]<-"n_inv"
grid40tots<-left_join(grid40tots,grid40inv)
grid40tots$n_inv<-ifelse(is.na(grid40tots$n_inv),0,grid40tots$n_inv)


grid40tots$small_grid<-as.character(grid40tots$small_grid)
#grid40tots<-filter(grid40tots,small_grid!=`integer(0)`)

class(grid40tots$small_grid)
#grid40tots<-filter(grid40tots,smallest_grid!="integer(0)")
#L4tots<-filter(L4tots,n>19) ###must have at least 20 total plots


dmer$small_grid<-as.character(dmer$small_grid)
class(dmer$smallest_grid)

library(hillR)
ecoregionz40<-unique(grid40tots$small_grid)


taxdat40<-data.frame()
for (z in c(1:length(ecoregionz40))){
  dat<-dplyr::filter(dmer,small_grid==ecoregionz40[z])
  
  dat<-filter(dat,!is.na(SpCode))
  dat<-dplyr::select(dat,Plot,SpCode,PctCov_100,small_grid)
  dat<-distinct(dat)
  comms<-matricize(dat)
  
  taxy0<-hill_taxa_parti(comms,q = 0)
  taxy1<-hill_taxa_parti(comms,q = 1)
  taxy2<-hill_taxa_parti(comms,q = 2)
  taxy<-rbind(taxy0,taxy1,taxy2)
  taxy$grid40<-ecoregionz[z]
  taxdat40<-rbind(taxy,taxdat40)
  print("done")
}  

colnames(grid40tots)[1]<-"grid40"
taxdat40<-left_join(taxdat40,grid40tots) 
taxdat40$perc_invaded<-taxdat40$n_inv/taxdat40$n



grid20tots<-p%>%group_by(big_grid) %>% count()

###invaded
grid20inv<-dmer %>% group_by(Plot,NativeStatus,big_grid)%>% summarize(RelCov=sum(PctCov_100,na.rm=TRUE))
grid20inv<-filter(grid20inv,NativeStatus=="I")
grid20inv<-filter(grid20inv,RelCov>5)
grid20inv<- grid20inv %>% group_by(big_grid) %>% count() ## plot where relcov I is >5

colnames(grid20inv)[2]<-"n_inv"
grid20tots$big_grid<-as.character(grid20tots$big_grid)
grid20tots<-left_join(grid20tots,grid20inv)
grid20tots$n_inv<-ifelse(is.na(grid20tots$n_inv),0,grid20tots$n_inv)


#grid20tots$big_grid<-as.character(grid20tots$big_grid)
#grid40tots<-filter(grid40tots,small_grid!=`integer(0)`)

class(grid20tots$big_grid)
#grid40tots<-filter(grid40tots,smallest_grid!="integer(0)")
#L4tots<-filter(L4tots,n>19) ###must have at least 20 total plots


dmer$big_grid<-as.character(dmer$big_grid)
class(dmer$big_grid)


ecoregionz20<-unique(grid20tots$big_grid)


taxdat20<-data.frame()
for (z in c(1:length(ecoregionz20))){
  dat<-dplyr::filter(dmer,big_grid==ecoregionz20[z])
  
  dat<-filter(dat,!is.na(SpCode))
  dat<-dplyr::select(dat,Plot,SpCode,PctCov_100,big_grid)
  dat<-distinct(dat)
  comms<-matricize(dat)
  
  taxy0<-hill_taxa_parti(comms,q = 0)
  taxy1<-hill_taxa_parti(comms,q = 1)
  taxy2<-hill_taxa_parti(comms,q = 2)
  taxy<-rbind(taxy0,taxy1,taxy2)
  taxy$grid20<-ecoregionz[z]
  taxdat20<-rbind(taxy,taxdat20)
  print("done")
}  

colnames(grid20tots)[1]<-"grid20"
taxdat20<-left_join(taxdat20,grid20tots) 

taxdat20$perc_invaded<-taxdat20$n_inv/taxdat20$n


ggplot(taxdat,aes(perc_invaded,local_similarity))+geom_point()+geom_smooth(method="lm")+facet_wrap(~q)
ggplot(taxdat40,aes(perc_invaded,local_similarity))+geom_point()+geom_smooth(method="lm")+facet_wrap(~q)

ggplot(taxdat20,aes(perc_invaded,local_similarity))+geom_point()+geom_smooth(method="lm")+facet_wrap(~q)




#####

}

phy<-read.tree("Analyses/Master/Input/FINAL_PHYLOGENY_SpCodes_Apr2023.tr")
d<-filter(d,SpCode %in% c(phy$tip.label)) 
###ecoregion
##Total # of plots per ecoregion
L4tots<-p%>%group_by(US_L4NAME) %>% count()

###invaded
L4invs<-d %>% group_by(Plot,NativeStatus,US_L4NAME)%>% summarize(RelCov=sum(PctCov_100,na.rm=TRUE))
L4invs<-filter(L4invs,NativeStatus=="I")
L4invs<-filter(L4invs,RelCov>5)
L4invsmean<- L4invs %>% group_by(US_L4NAME) %>% summarise(meanCov=mean(RelCov,na.rm=TRUE)) ## plot where relcov I is >5
L4invs<- L4invs %>% group_by(US_L4NAME) %>% count() ## plot where relcov I is >5

colnames(L4invs)[2]<-"n_inv"
L4tots<-left_join(L4tots,L4invs)
L4tots$n_inv<-ifelse(is.na(L4tots$n_inv),0,L4tots$n_inv)
L4tots<-filter(L4tots,n>9) ###must have at least 20 total plots

matricize<-function(x){
  temp<-dplyr::select(x,Plot,SpCode,PctCov_100) ### select useful columns
  temp<-dplyr::filter(temp,!is.na(temp$SpCode))
  temp<- tidyr::spread(temp,SpCode,PctCov_100)
  temp[is.na(temp)] <- 0 
  temp<-as.data.frame(temp)
  rownames(temp)<-temp$Plot # convert plot names to row names
  temp<-dplyr::select(temp,-Plot)
  temp<-as.matrix.data.frame(temp)#
}


library(hillR)
ecoregionz<-unique(L4tots$US_L4NAME)
ecoregionz<-intersect(ecoregionz,unique(d$US_L4NAME))
ecoregionz<-ecoregionz[-554]

taxdat<-data.frame()
phydat<-data.frame()
for (z in c(1:length(ecoregionz))){
  dat<-dplyr::filter(d,US_L4NAME==ecoregionz[z])
  
  comms<-matricize(dat)
  
  taxy0<-hill_taxa_parti(comms,q = 0)
  taxy1<-hill_taxa_parti(comms,q = 1)
  taxy2<-hill_taxa_parti(comms,q = 2)
  
  taxy<-rbind(taxy0,taxy1,taxy2)
  taxy$US_L4NAME<-ecoregionz[z]
  taxyclass<-"taxonomic"
  
  phy0<-hill_phylo_parti(comms,phy,q=0)
  phy1<-hill_phylo_parti(comms,phy,q=1)
  phy2<-hill_phylo_parti(comms,phy,q=2)
  
  
  phyt<-rbind(phy0,phy1,phy2)
  phyt$US_L4NAME<-ecoregionz[z]
  phyt$class<-"phylogenetic"
  
  taxdat<-rbind(taxy,taxdat)
  phydat<-rbind(phyt,phydat)

  }  

taxdat<-left_join(taxdat,L4tots) 
taxdat$perc_invaded<-taxdat$n_inv/taxdat$n

phydat<-left_join(phydat,L4tots) 
phydat$perc_invaded<-phydat$n_inv/phydat$n





library(ggplot2)



if(FALSE){

##Total # of plots per ecoregion
L3tots<-p%>%group_by(NA_L3NAME) %>% count()

###invaded
L3invs<-d %>% group_by(Plot,NativeStatus,NA_L3NAME)%>% summarize(RelCov=sum(PctCov_100,na.rm=TRUE))
L3invs<-filter(L3invs,NativeStatus=="I")
L3invs<-filter(L3invs,RelCov>5)
L3invs<- L3invs %>% group_by(NA_L3NAME) %>% count() ## plot where relcov I is >5

colnames(L3invs)[2]<-"n_inv"
L3tots<-left_join(L3tots,L3invs)
L3tots$n_inv<-ifelse(is.na(L3tots$n_inv),0,L3tots$n_inv)
L3tots<-filter(L3tots,n>19) ###must have at least 20 total plots

ecoregionz3<-unique(L3tots$NA_L3NAME)
ecoregionz3<-intersect(ecoregionz3,unique(d$NA_L3NAME))
ecoregionz3<-ecoregionz3[-82]

taxdat3<-data.frame()
for (z in c(1:length(ecoregionz3))){
  dat<-dplyr::filter(d,NA_L3NAME==ecoregionz3[z])
  dat<-filter(dat,!is.na(SpCode))
  dat<-dplyr::select(dat,Plot,SpCode,PctCov_100,NA_L3NAME)
  dat<-distinct(dat)
  comms<-matricize(dat)
  
  #taxy0<-hill_taxa_parti(comms,q = 0)
  taxy1<-hill_taxa_parti(comms,q = 1)
  #taxy2<-hill_taxa_parti(comms,q = 2)
  #taxy<-rbind(taxy0,taxy1,taxy2)
  taxy1$NA_L3NAME<-ecoregionz3[z]
  taxdat3<-rbind(taxy1,taxdat3)
}  


taxdat3<-left_join(taxdat3,L3tots) 
taxdat3$perc_invaded<-taxdat3$n_inv/taxdat3$n



key2<-dplyr::select(key,-US_L4NAME)
key2<-distinct(key2)

taxdat<-left_join(taxdat,key)
taxdat3<-left_join(taxdat3,key)
ggplot(taxdat,aes(perc_invaded,local_similarity))+
  geom_point(aes())+
  geom_smooth(method="lm",aes(),se=FALSE)+
  facet_wrap(~NA_L1NAME)
ggplot(taxdat3,aes(perc_invaded,local_similarity))+
  geom_point(aes())+
  geom_smooth(method="lm",aes(),se=FALSE)+
  facet_wrap(~NA_L1NAME)
}
key<-dplyr::select(regionz,-Plot)
key<-distinct(key)

key2<-dplyr::select(key,-US_L4NAME)
key2<-distinct(key2)
taxdat<-left_join(taxdat,key)
phydat<-left_join(phydat,key)


###todo add n as a covariate
taxdat$n_scale<-taxdat$n/max(taxdat$n)
phydat$n_scale<-phydat$n/max(phydat$n)



taxdat<-left_join(taxdat,L4invsmean)
phydat<-left_join(phydat,L4invsmean)


#taxdat<-filter(taxdat,!NA_L1NAME%in%c("TROPICAL WET FORESTS","SOUTHERN SEMI-ARID HIGHLANDS"))
#phydat<-filter(phydat,!NA_L1NAME%in%c("TROPICAL WET FORESTS","SOUTHERN SEMI-ARID HIGHLANDS"))

taxdatq0<-filter(taxdat,q==0)
taxdatq1<-filter(taxdat,q==1)
taxdatq2<-filter(taxdat,q==2)

phydatq0<-filter(phydat,q==0)
phydatq1<-filter(phydat,q==1)
phydatq2<-filter(phydat,q==2)
library(brms)
q0.tax<- brm(
  bf(local_similarity ~n_scale+ perc_invaded+(perc_invaded|NA_L1NAME),
     phi ~n_scale+ perc_invaded+(perc_invaded|NA_L1NAME)),
  data = taxdatq0,
  family = Beta(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 6000, warmup = 5000,
  cores = 4, seed = 4321,backend = "cmdstanr")


q1.tax<- brm(
  bf(local_similarity ~n_scale+ perc_invaded+(perc_invaded|NA_L1NAME),
     phi ~n_scale+ perc_invaded+(perc_invaded|NA_L1NAME)),
  data = taxdatq1,
  family = Beta(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr")


coef(q1.tax,prob=c(.25,.75,.055,.945))

coef(q0.phy,prob=c(.25,.75,.055,.945))

round(fixef(q1.tax,prob=c(.25,.75,.055,.945)),2)
round(fixef(q1.phy,prob=c(.25,.75,.055,.945)),2)

round(fixef(q0.phy,prob=c(.25,.75,.055,.945)),2)
round(fixef(q0.tax,prob=c(.25,.75,.055,.945)),2)

round(fixef(q2.phy,prob=c(.25,.75,.055,.945)),2)
round(fixef(q2.tax,prob=c(.25,.75,.055,.945)),2)

coef(q1.phy,prob=c(.25,.75,.055,.945))
q2.tax<- brm(
  bf(local_similarity ~n_scale+ perc_invaded+(perc_invaded|NA_L1NAME),
     phi ~n_scale+ perc_invaded+(perc_invaded|NA_L1NAME)),
  data = taxdatq2,
  family = Beta(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr")



q0.phy<- brm(
  bf(local_similarity ~n_scale+ perc_invaded+(perc_invaded|NA_L1NAME),
     phi ~ n_scale+perc_invaded+(perc_invaded|NA_L1NAME)),
  data = phydatq0,
  family = Beta(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr")


q1.phy<- brm(
  bf(local_similarity ~n_scale+ perc_invaded+(perc_invaded|NA_L1NAME),
     phi ~ n_scale+perc_invaded+(perc_invaded|NA_L1NAME)),
  data = phydatq1,
  family = Beta(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr")

q2.phy<- brm(
  bf(local_similarity ~n_scale+ perc_invaded+(perc_invaded|NA_L1NAME),
     phi ~ n_scale+perc_invaded+(perc_invaded|NA_L1NAME)),
  data = phydatq2,
  family = Beta(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr")





fixef(cept.only.tax)
fixef(cept.only.phy)
coef(cept.only.tax,probs = c(.25,.75,.05,.95))

new.data<-data.frame(NA_L1NAME=rep(unique(taxdat$NA_L1NAME),each=11),perc_invaded=rep(c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1),10))
new.data$n_scale<-median(taxdat$n_scale)

new.data2<-data.frame(perc_invaded=c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1))
new.data2$n_scale<-median(taxdat$n_scale)



#### predictions
library(tidybayes)
toy_pred.l4.phy<- q1.phy %>% 
  epred_draws(newdata =new.data,ndraws = 1000)
toy_pred.l4.tax<- q1.tax %>% 
  epred_draws(newdata =new.data,ndraws = 1000)


toy_pred.l4.phy.noreg<- q1.phy %>% 
  epred_draws(newdata =new.data2,ndraws = 1000,re_formula = NA)
toy_pred.l4.tax.noreg<- q1.tax %>% 
  epred_draws(newdata =new.data2,ndraws = 1000,re_formula = NA)


toy_pred.l4.phy2<- q2.phy %>% 
  epred_draws(newdata =new.data,ndraws = 1000)
toy_pred.l4.tax2<- q2.tax %>% 
  epred_draws(newdata =new.data,ndraws = 1000)

toy_pred.l4.phy2.noreg<- q2.phy %>% 
  epred_draws(newdata =new.data2,ndraws = 1000,re_formula = NA)
toy_pred.l4.tax2.noreg<- q2.tax %>% 
  epred_draws(newdata =new.data2,ndraws = 1000,re_formula=NA)


toy_pred.l4.tax2$class<-"taxonomic"
toy_pred.l4.phy2$class<-"phylogenetic"

toy_pred.l4.tax2.noreg$class<-"taxonomic"
toy_pred.l4.phy2.noreg$class<-"phylogenetic"

toy_pred.l4.phy2$grouper<-paste(toy_pred.l4.phy2$NA_L1NAME,toy_pred.l4.phy2$class,toy_pred.l4.phy2$.draw)
toy_pred.l4.tax2$grouper<-paste(toy_pred.l4.tax2$NA_L1NAME,toy_pred.l4.tax2$class,toy_pred.l4.tax2$.draw)

toy_pred.l4.phy2.noreg$grouper<-paste(toy_pred.l4.phy2.noreg$class,toy_pred.l4.phy2.noreg$.draw)
toy_pred.l4.tax2.noreg$grouper<-paste(toy_pred.l4.tax2.noreg$class,toy_pred.l4.tax2.noreg$.draw)




toy_pred.l4.tax$class<-"taxonomic"
toy_pred.l4.phy$class<-"phylogenetic"

toy_pred.l4.tax.noreg$class<-"taxonomic"
toy_pred.l4.phy.noreg$class<-"phylogenetic"

toy_pred.l4.phy$grouper<-paste(toy_pred.l4.phy$NA_L1NAME,toy_pred.l4.phy$class,toy_pred.l4.phy$.draw)
toy_pred.l4.tax$grouper<-paste(toy_pred.l4.tax$NA_L1NAME,toy_pred.l4.tax$class,toy_pred.l4.tax$.draw)

toy_pred.l4.phy.noreg$grouper<-paste(toy_pred.l4.phy.noreg$class,toy_pred.l4.phy.noreg$.draw)
toy_pred.l4.tax.noreg$grouper<-paste(toy_pred.l4.tax.noreg$class,toy_pred.l4.tax.noreg$.draw)



toy_pred.l4.phy0<- q0.phy %>% 
  epred_draws(newdata =new.data,ndraws = 1000)
toy_pred.l4.tax0<- q0.tax %>% 
  epred_draws(newdata =new.data,ndraws = 1000)


toy_pred.l4.phy0.noreg<- q0.phy %>% 
  epred_draws(newdata =new.data2,ndraws = 1000,re_formula = NA)
toy_pred.l4.tax0.noreg<- q0.tax %>% 
  epred_draws(newdata =new.data2,ndraws = 1000,re_formula = NA)


toy_pred.l4.tax0$class<-"taxonomic"
toy_pred.l4.phy0$class<-"phylogenetic"

toy_pred.l4.tax0.noreg$class<-"taxonomic"
toy_pred.l4.phy0.noreg$class<-"phylogenetic"

toy_pred.l4.phy0$grouper<-paste(toy_pred.l4.phy0$NA_L1NAME,toy_pred.l4.phy0$class,toy_pred.l4.phy0$.draw)
toy_pred.l4.tax0$grouper<-paste(toy_pred.l4.tax0$NA_L1NAME,toy_pred.l4.tax0$class,toy_pred.l4.tax0$.draw)

toy_pred.l4.phy0.noreg$grouper<-paste(toy_pred.l4.phy0.noreg$class,toy_pred.l4.phy0.noreg$.draw)
toy_pred.l4.tax0.noreg$grouper<-paste(toy_pred.l4.tax0.noreg$class,toy_pred.l4.tax0.noreg$.draw)




library(ggplot2)


aaa<-ggpubr::ggarrange(ggplot()+
                    geom_point(data=taxdatq1,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
                    
                    geom_line(data=toy_pred.l4.tax.noreg,aes(x=perc_invaded*100,y=.epred, group=grouper),size=.01,alpha=.5)+
                    coord_cartesian(ylim=c(.2,.9))+ggthemes::theme_few(base_size = 10)+geom_smooth(method="glm",color="black")+ylab("taxonomic similarity")+xlab("% invaded"),
                  
                  
                  ggplot()+
                    geom_point(data=phydatq1,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
                    
                    geom_line(data=toy_pred.l4.phy.noreg,aes(x=perc_invaded*100,y=.epred, group=grouper),size=.01,alpha=.5)+
                    coord_cartesian(ylim=c(.75,1))+ggthemes::theme_few(base_size = 10)+geom_smooth(method="glm",color="black")+ylab("phylogenetic similarity")+xlab("% invaded"),ncol=2)



bbb<-ggpubr::ggarrange(ggplot()+
                    geom_point(data=taxdatq0,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
                    
                    geom_line(data=toy_pred.l4.tax0.noreg,aes(x=perc_invaded*100,y=.epred, group=grouper),size=.01,alpha=.5)+
                    coord_cartesian(ylim=c(.4,1))+ggthemes::theme_few(base_size = 10)+geom_smooth(method="glm",color="black")+ylab("taxonomic similarity")+xlab("% invaded"),
                  
                  
                  ggplot()+
                    geom_point(data=phydatq0,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
                    
                    geom_line(data=toy_pred.l4.phy0.noreg,aes(x=perc_invaded*100,y=.epred, group=grouper),size=.01,alpha=.5)+
                    coord_cartesian(ylim=c(.75,1))+ggthemes::theme_few(base_size = 10)+geom_smooth(method="glm",color="black")+ylab("phylogenetic similarity")+xlab("% invaded"),ncol=2)


ccc<-ggpubr::ggarrange(ggplot()+
                    geom_point(data=taxdatq2,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
                    
                    geom_line(data=toy_pred.l4.tax2.noreg,aes(x=perc_invaded*100,y=.epred, group=grouper),size=.01,alpha=.5)+
                    coord_cartesian(ylim=c(.0,.75))+ggthemes::theme_few(base_size = 10)+geom_smooth(method="glm",color="black")+ylab("taxonomic similarity")+xlab("% invaded"),
                  
                  
                  ggplot()+
                    geom_point(data=phydatq2,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
                    
                    geom_line(data=toy_pred.l4.phy2.noreg,aes(x=perc_invaded*100,y=.epred, group=grouper),size=.01,alpha=.5)+
                    coord_cartesian(ylim=c(.6,1))+ggthemes::theme_few(base_size = 10)+geom_smooth(method="glm",color="black")+ylab("phylogenetic similarity")+xlab("% invaded"),ncol=2)



jpeg("Analyses/Master/Plots/q1main.jpeg",width = 6,height=3,units='in',res=200)
aaa
dev.off()

jpeg("Analyses/Master/Plots/q0main.jpeg",width = 6,height=3,units='in',res=200)
bbb
dev.off()


jpeg("Analyses/Master/Plots/q2main.jpeg",width = 6,height=3,units='in',res=200)
ccc
dev.off()


mainq1tax<-filter(taxdatq1,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))
mainq1phy<-filter(phydatq1,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))

toy_mainq1tax<-filter(toy_pred.l4.tax,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))
toy_mainq1phy<-filter(toy_pred.l4.phy,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))


mainq0tax<-filter(taxdatq0,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))
mainq0phy<-filter(phydatq0,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))

toy_mainq0tax<-filter(toy_pred.l4.tax0,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))
toy_mainq0phy<-filter(toy_pred.l4.phy0,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))


mainq2tax<-filter(taxdatq2,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))
mainq2phy<-filter(phydatq2,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))

toy_mainq2tax<-filter(toy_pred.l4.tax2,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))
toy_mainq2phy<-filter(toy_pred.l4.phy2,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))


jpeg("Analyses/Master/Plots/q1expansion.jpeg",width = 11,height=8,units='in',res=200)
ggpubr::ggarrange(ggplot()+
                    geom_point(data=taxdatq1,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
                    
                    geom_line(data=toy_pred.l4.tax,aes(x=perc_invaded*100,y=.epred, group=grouper),size=.01,alpha=.5)+facet_wrap(~NA_L1NAME,nrow=2)+
                    coord_cartesian(ylim=c(.25,.8))+ggthemes::theme_few(base_size = 9)+geom_smooth(method="glm",color="black")+ylab("taxonomic similarity")+xlab("% invaded"),
                  
                  
                  ggplot()+
                    geom_point(data=phydatq1,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
                    
                    geom_line(data=toy_pred.l4.phy,aes(x=perc_invaded*100,y=.epred, group=grouper),size=.01,alpha=.5)+facet_wrap(~NA_L1NAME,nrow=2)+
                    coord_cartesian(ylim=c(.8,1))+ggthemes::theme_few(base_size = 9)+geom_smooth(method="glm",color="black")+ylab("phylogenetic similarity")+xlab("% invaded"),ncol=1)

dev.off()


aa<-ggpubr::ggarrange(ggplot()+
                    geom_point(data=mainq1tax,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
                    
                    geom_line(data=toy_mainq1tax,aes(x=perc_invaded*100,y=.epred, group=grouper),size=.01,alpha=.5)+facet_wrap(~NA_L1NAME,nrow=5)+
                    coord_cartesian(ylim=c(.25,.8))+ggthemes::theme_few(base_size = 9)+geom_smooth(method="glm",color="black")+ylab("taxonomic similarity")+xlab("% invaded"),
                  
                  
                  ggplot()+
                    geom_point(data=mainq1phy,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
                    
                    geom_line(data=toy_mainq1phy,aes(x=perc_invaded*100,y=.epred, group=grouper),size=.01,alpha=.5)+facet_wrap(~NA_L1NAME,nrow=5)+
                    coord_cartesian(ylim=c(.8,1))+ggthemes::theme_few(base_size = 9)+geom_smooth(method="glm",color="black")+ylab("phylogenetic similarity")+xlab("% invaded"),ncol=2)


bb<-ggpubr::ggarrange(ggplot()+
                    geom_point(data=mainq0tax,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
                    
                    geom_line(data=toy_mainq0tax,aes(x=perc_invaded*100,y=.epred, group=grouper),size=.01,alpha=.5)+facet_wrap(~NA_L1NAME,nrow=5)+
                    coord_cartesian(ylim=c(.6,1))+ggthemes::theme_few(base_size = 9)+geom_smooth(method="glm",color="black")+ylab("taxonomic similarity")+xlab("% invaded"),
                  
                  
                  ggplot()+
                    geom_point(data=mainq0phy,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
                    
                    geom_line(data=toy_mainq0phy,aes(x=perc_invaded*100,y=.epred, group=grouper),size=.01,alpha=.5)+facet_wrap(~NA_L1NAME,nrow=5)+
                    coord_cartesian(ylim=c(.8,1))+ggthemes::theme_few(base_size = 9)+geom_smooth(method="glm",color="black")+ylab("phylogenetic similarity")+xlab("% invaded"),ncol=2)


cc<-ggpubr::ggarrange(ggplot()+
                    geom_point(data=mainq2tax,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
                    
                    geom_line(data=toy_mainq2tax,aes(x=perc_invaded*100,y=.epred, group=grouper),size=.01,alpha=.5)+facet_wrap(~NA_L1NAME,nrow=5)+
                    coord_cartesian(ylim=c(0,.6))+ggthemes::theme_few(base_size = 9)+geom_smooth(method="glm",color="black")+ylab("taxonomic similarity")+xlab("% invaded"),
                  
                  
                  ggplot()+
                    geom_point(data=mainq2phy,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
                    
                    geom_line(data=toy_mainq2phy,aes(x=perc_invaded*100,y=.epred, group=grouper),size=.01,alpha=.5)+facet_wrap(~NA_L1NAME,nrow=5)+
                    coord_cartesian(ylim=c(.5,1))+ggthemes::theme_few(base_size = 9)+geom_smooth(method="glm",color="black")+ylab("phylogenetic similarity")+xlab("% invaded"),ncol=2)



jpeg("Analyses/Master/Plots/q1big5jpeg",width = 6,height=6,units='in',res=200)
aa
dev.off()

jpeg("Analyses/Master/Plots/q0big5jpeg",width = 6,height=6,units='in',res=200)
bb
dev.off()

jpeg("Analyses/Master/Plots/q2big5jpeg",width = 6,height=6,units='in',res=200)
cc
dev.off()

jpeg("Analyses/Master/Plots/q2expansion.jpeg",width = 11,height=5,units='in',res=200)
ggpubr::ggarrange(ggplot()+
  geom_point(data=taxdatq2,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
  
  geom_line(data=toy_pred.l4.tax2,aes(x=perc_invaded*100,y=.epred, group=grouper),size=.01,alpha=.5)+facet_wrap(~NA_L1NAME,nrow=2)+
  ggthemes::theme_few(base_size = 10)+geom_smooth(method="glm",color="black")+ylab("taxonomic similarity")+xlab("% invaded")+coord_cartesian(ylim=c(.0,.5)),


ggplot()+
  geom_point(data=phydatq2,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
  
  geom_line(data=toy_pred.l4.phy2,aes(x=perc_invaded*100,y=.epred, group=grouper),size=.01,alpha=.5)+facet_wrap(~NA_L1NAME,nrow=2)+
ggthemes::theme_few(base_size = 10)+geom_smooth(method="glm",color="black")+ylab("phylogenetic similarity")+xlab("% invaded")+coord_cartesian(ylim=c(.5,1)),ncol=1)

dev.off()



jpeg("Analyses/Master/Plots/q0expansion.jpeg",width = 11,height=5,units='in',res=200)
ggpubr::ggarrange(ggplot()+
                    geom_point(data=taxdatq0,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
                    
                    geom_line(data=toy_pred.l4.tax0,aes(x=perc_invaded*100,y=.epred, group=grouper),size=.01,alpha=.5)+facet_wrap(~NA_L1NAME,nrow=1)+
                    ggthemes::theme_few(base_size = 10)+geom_smooth(method="glm",color="black")+ylab("taxonomic similarity")+xlab("% invaded"),
                  
                  
                  ggplot()+
                    geom_point(data=phydatq0,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
                    
                    geom_line(data=toy_pred.l4.phy0,aes(x=perc_invaded*100,y=.epred, group=grouper),size=.01,alpha=.5)+facet_wrap(~NA_L1NAME,nrow=1)+
                    ggthemes::theme_few(base_size = 10)+geom_smooth(method="glm",color="black")+ylab("phylogenetic similarity")+xlab("% invaded")+coord_cartesian(ylim=c(.8,1)),ncol=1)


dev.off()
d.invy<-filter(d,NativeStatus=="I")
dsumz<- d.invy %>% group_by(Plot,US_L4NAME,NA_L1NAME) %>% count()
dsumz<- dsumz %>% group_by(US_L4NAME,NA_L1NAME) %>% summarise(I_rich=mean(n))

dsumz.gamma<- d.invy %>% group_by(AcceptedTaxonName,US_L4NAME,NA_L1NAME) %>% count()
dsumz.gamma<-dplyr::select(dsumz.gamma,-n)
dsumz.gamma<-distinct(dsumz.gamma)
dsumz.gamma<- dsumz.gamma %>% group_by(US_L4NAME,NA_L1NAME) %>% count()
dsumz.gamma<-filter(dsumz.gamma,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))

dsumz<-filter(dsumz,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))

ggpubr::ggarrange(
ggplot(mainq1tax,aes(NA_L1NAME,perc_invaded))+stat_summary(),
ggplot(mainq1tax,aes(NA_L1NAME,meanCov))+stat_summary(),
ggplot(dsumz,aes(NA_L1NAME,I_rich))+stat_summary())
ggplot(dsumz.gamma,aes(NA_L1NAME,n))+stat_summary()

taxdatq1$Iabn<-taxdatq1$meanCov/100
taxdatq1$Iabn<-ifelse(taxdatq1$Iabn>1,1,taxdatq1$Iabn)
 ### if we like it redo these with just the two regions of interest
richmod<-brm(
  bf(I_rich ~(1|NA_L1NAME)),
  data = dsumz,
  family =  gaussian(),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr")

gammamod<-brm(
  bf(n ~(1|NA_L1NAME)),
  data = dsumz.gamma,
  family =  gaussian(),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr")

abnmod<-brm(
  bf(Iabn ~(1|NA_L1NAME),
     phi ~ (1|NA_L1NAME),
  zoi~ (1|NA_L1NAME),
  coi~ (1|NA_L1NAME)),
  data = taxdatq1,
  family =  zero_one_inflated_beta(),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr")



covmod<-brm(
  bf(perc_invaded ~(1|NA_L1NAME),
     phi ~ (1|NA_L1NAME),
     zoi~ (1|NA_L1NAME),
     coi~ (1|NA_L1NAME)),
  data = taxdatq1,
  family =  zero_one_inflated_beta(),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr")

coef(covmod)
covdat<-data.frame(NA_L1NAME=c('EASTERN TEMPERATE FORESTS','NORTH AMERICAN DESERTS'))
predcov<-covmod %>% 
  epred_draws(newdata =covdat,ndraws = 1000)

predabn<-abnmod %>% 
  epred_draws(newdata =covdat,ndraws = 1000)

predgamma<-gammamod %>% 
  epred_draws(newdata =covdat,ndraws = 1000)

predrich<-richmod %>% 
  epred_draws(newdata =covdat,ndraws = 1000)

jpeg("Analyses/Master/Plots/ETFvNAD.jpeg",width = 7,height=4,units='in',res=200)
ggpubr::ggarrange(ggplot(predcov,aes(.epred,NA_L1NAME))+stat_halfeye(.width = c(.5,.9))+ggthemes::theme_few(base_size = 9)+xlab("mean % invaded")+ylab(""),
ggplot(predabn,aes(.epred,NA_L1NAME))+stat_halfeye(.width = c(.5,.9))+ggthemes::theme_few(base_size = 9)+xlab("mean invasive cover/plot")+ylab(""),

ggplot(predrich,aes(.epred,NA_L1NAME))+stat_halfeye(.width = c(.5,.9))+ggthemes::theme_few(base_size = 9)+xlab("mean invasive richness/plot")+ylab(""),
ggplot(predgamma,aes(.epred,NA_L1NAME))+stat_halfeye(.width = c(.5,.9))+ggthemes::theme_few(base_size = 9)+xlab("invasive gamma diversity")+ylab(""))
dev.off()



inyonly<-filter(L4tots,n_inv>1)
ecoregionz<-unique(inyonly$US_L4NAME)
ecoregionz<-intersect(ecoregionz,unique(d.invy$US_L4NAME))
ecoregionz<-ecoregionz[-461]

taxdat.invy<-data.frame()
phydat.invy<-data.frame()
for (z in c(1:length(ecoregionz))){
  dat<-dplyr::filter(d.invy,US_L4NAME==ecoregionz[z])
  
  comms<-matricize(dat)
  
  taxy0<-hill_taxa_parti(comms,q = 0)
  taxy1<-hill_taxa_parti(comms,q = 1)
  taxy2<-hill_taxa_parti(comms,q = 2)
  
  taxy<-rbind(taxy0,taxy1,taxy2)
  taxy$US_L4NAME<-ecoregionz[z]
  taxyclass<-"taxonomic"
  
  phy0<-hill_phylo_parti(comms,phy,q=0)
  phy1<-hill_phylo_parti(comms,phy,q=1)
  phy2<-hill_phylo_parti(comms,phy,q=2)
  
  
  phyt<-rbind(phy0,phy1,phy2)
  phyt$US_L4NAME<-ecoregionz[z]
  phyt$class<-"phylogenetic"
  
  taxdat.invy<-rbind(taxy,taxdat.invy)
  phydat.invy<-rbind(phyt,phydat.invy)
  
}  

taxinvyq1<-filter(taxdat.invy,q==1)
taxinvyq1<-left_join(taxinvyq1,key)
taxinvyq1<-filter(taxinvyq1,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))

ggplot(taxinvyq1,aes(NA_L1NAME,TD_beta))+stat_summary()
ggplot(taxinvyq1,aes(NA_L1NAME,TD_alpha))+stat_summary()
ggplot(taxinvyq1,aes(NA_L1NAME,TD_gamma))+stat_summary()

taxyinvq1<-left_join(taxinvyq1,L4tots)
taxyinvq1$perc_invaded<-taxyinvq1$n_inv/taxyinvq1$n

phyinvq1<-left_join(phyinvyq1,L4tots)
phyinvq1$perc_invaded<-phyinvq1$n_inv/phyinvq1$n

ggplot(taxyinvq1,aes(perc_invaded,TD_gamma))+geom_smooth(method="lm")+facet_wrap(~NA_L1NAME,scales="free_y")+geom_point()
ggplot(phyinvq1,aes(perc_invaded,PD_gamma))+geom_smooth(method="lm")+facet_wrap(~NA_L1NAME)+geom_point()

ggplot(taxyinvq1,aes(perc_invaded,TD_alpha))+geom_smooth(method="lm")+facet_wrap(~NA_L1NAME,scales="free_y")+geom_point()
ggplot(taxyinvq1,aes(perc_invaded,TD_beta))+geom_smooth(method="lm")+facet_wrap(~NA_L1NAME,scales="free_y")+geom_point()

ggplot(taxyinvq1,aes(perc_invaded,local_similarity))+geom_smooth(method="lm")+facet_wrap(~NA_L1NAME,scales="free_y")+geom_point()
ggplot(phyinvq1,aes(perc_invaded,local_similarity))+geom_smooth(method="lm")+facet_wrap(~NA_L1NAME,scales="free_y")+geom_point()

taxyinvq1$n_scale<-taxyinvq1$n/max(taxyinvq1$n)
q1.tax.invy<- brm(
  bf(local_similarity ~n_scale+ perc_invaded+(perc_invaded|NA_L1NAME),
     phi ~n_scale+ perc_invaded+(perc_invaded|NA_L1NAME)),
  data = taxyinvq1,
  family = Beta(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr")
coef(q1.tax.invy)

toy_pred.l4.taxinvy<- q1.tax.invy %>% 
  epred_draws(newdata =new.data,ndraws = 1000,allow_new_levels=TRUE)

toy_pred.l4.taxinvy$grouper<-paste(toy_pred.l4.taxinvy$NA_L1NAME,toy_pred.l4.taxinvy$.draw)




library(ggplot2)


ggplot()+#geom_point(data=taxinvyq1,aes(x=perc_invaded*100,y=local_similarity),size=.1)+
                         
                         geom_line(data=toy_pred.l4.taxinvy,aes(x=perc_invaded*100,y=.epred, group=grouper),size=.01,alpha=.5)+facet_wrap(~NA_L1NAME)+
                         coord_cartesian(ylim=c(.2,.9))+ggthemes::theme_few(base_size = 10)+geom_smooth(method="glm",color="black")+ylab("taxonomic similarity")+xlab("% invaded")
                       



q1.gamma.invy<- brm(
  bf(TD_gamma ~n_scale+ perc_invaded+(perc_invaded|NA_L1NAME)),
  data = taxyinvq1,
  family = gaussian(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr")
coef(q1.gamma.invy,probs = c(.25,.75))


q1.alpha.invy<- brm(
  bf(TD_alpha ~n_scale+ perc_invaded+(perc_invaded|NA_L1NAME)),
  data = taxyinvq1,
  family = gaussian(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr")
coef(q1.alpha.invy)

toy_pred.l4.phy<- q1.phy %>% 
  epred_draws(newdata =new.data,ndraws = 1000)






phyinvyq1<-filter(phydat.invy,q==1)
phyinvyq1<-left_join(phyinvyq1,key)
phyinvyq1<-filter(phyinvyq1,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))


gammamod<-brm(
  bf(TD_gamma ~n_scale+(1|NA_L1NAME)),
  data = taxyinvq1,
  family =  gaussian(),
  control=list(adapt_delta=.95),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr")

gammamodphy<-brm(
  bf(PD_gamma ~(1|NA_L1NAME)),
  data = phyinvyq1,
  family =  gaussian(),
  control=list(adapt_delta=.95),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr")


alphamod<-brm(
  bf(TD_alpha ~n_scale+(1|NA_L1NAME)),
  data = taxyinvq1,
  family =  gaussian(),
  control=list(adapt_delta=.95),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr")

alphamodphy<-brm(
  bf(PD_alpha ~(1|NA_L1NAME)),
  data = phyinvyq1,
  family =  gaussian(),
  control=list(adapt_delta=.95),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr")

betamod<-brm(
  bf(TD_beta ~(1|NA_L1NAME)),
  data = taxinvyq1,
  family =  gaussian(),
  control=list(adapt_delta=.95),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr")

betamodphy<-brm(
  bf(PD_beta ~(1|NA_L1NAME)),
  data = phyinvyq1,
  family =  gaussian(),
  control=list(adapt_delta=.95),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr")

simmod<-brm(
  bf(local_similarity ~n_scale+(1|NA_L1NAME),
     phi ~ n_scale+(1|NA_L1NAME)),
  data = taxyinvq1,
  control=list(adapt_delta=.99),
  family =  Beta(),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr")

simmodr<-brm(
  bf(region_similarity ~(1|NA_L1NAME),
     phi ~ (1|NA_L1NAME)),
  data = taxinvyq1,
  control=list(adapt_delta=.99),
  family =  Beta(),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr")

simmodphy<-brm(
  bf(local_similarity ~(1|NA_L1NAME),
     phi ~ (1|NA_L1NAME)),
  data = phyinvyq1,
  control=list(adapt_delta=.99),
  family =  Beta(),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr")


expdat<-data.frame(NA_L1NAME=c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))
expdat$n_scale<-median(taxyinvq1$n_scale)
predgamma<-gammamod %>% 
  epred_draws(newdata =expdat,ndraws = 1000)

predgammaphy<-gammamodphy %>% 
  epred_draws(newdata =expdat,ndraws = 1000)

predalpha<-alphamod %>% 
  epred_draws(newdata =expdat,ndraws = 1000)

predalphaphy<-alphamodphy %>% 
  epred_draws(newdata =expdat,ndraws = 1000)

predbeta<-betamod %>% 
  epred_draws(newdata =expdat,ndraws = 1000)

predbetaphy<-betamodphy %>% 
  epred_draws(newdata =expdat,ndraws = 1000)

predsim<-simmod %>% 
  epred_draws(newdata =expdat,ndraws = 1000)

predsimr<-simmodr %>% 
  epred_draws(newdata =expdat,ndraws = 1000)

predsimphy<-simmodphy %>% 
  epred_draws(newdata =expdat,ndraws = 1000)

pa<-ggplot(predgamma,aes(.epred,NA_L1NAME))+stat_halfeye(.width = c(.5,.9))+ggthemes::theme_few(base_size = 9)+xlab("non-native gamma diversity")+ylab("")
paa<-ggplot(predgammaphy,aes(.epred,NA_L1NAME))+stat_halfeye(.width = c(.5,.9))+ggthemes::theme_few(base_size = 9)+xlab("invasive gamma diversity")+ylab("")

pb<-ggplot(predbeta,aes(.epred,NA_L1NAME))+stat_halfeye(.width = c(.5,.9))+ggthemes::theme_few(base_size = 9)+xlab("invasive beta diversity")+ylab("")
ggplot(predbetaphy,aes(.epred,NA_L1NAME))+stat_halfeye(.width = c(.5,.9))+ggthemes::theme_few(base_size = 9)+xlab("invasive beta diversity")+ylab("")


pc<-ggplot(predalpha,aes(.epred,NA_L1NAME))+stat_halfeye(.width = c(.5,.9))+ggthemes::theme_few(base_size = 9)+xlab("non-native alpha diversity")+ylab("")+theme(axis.text.y = element_blank(),axis.ticks.y = element_blank())
pcc<-ggplot(predalphaphy,aes(.epred,NA_L1NAME))+stat_halfeye(.width = c(.5,.9))+ggthemes::theme_few(base_size = 9)+xlab("invasive alpha diversity")+ylab("")

pd<-ggplot(predsim,aes(.epred,reorder(NA_L1NAME,-.epred)))+stat_halfeye(.width = c(.5,.9))+ggthemes::theme_few(base_size = 9)+xlab("non-native taxonomic similarity")+ylab("")


pdd<-ggplot(predsimphy,aes(.epred,reorder(NA_L1NAME,-.epred)))+stat_halfeye(.width = c(.5))+ggthemes::theme_few(base_size = 9)+xlab("invasive similarity diversity")+ylab("")

ggplot(predsimr,aes(.epred,reorder(NA_L1NAME,-.epred)))+stat_halfeye(.width = c(.5))+ggthemes::theme_few(base_size = 9)+xlab("invasive similarity diversity")+ylab("")
pe<-ggpubr::ggarrange(pa,pc,widths=c(.6,.3))

jpeg("Analyses/Master/Plots/invdiversity.jpeg",width = 6,height=5,units='in',res=200)
ggpubr::ggarrange(pd,pe,ncol=1)
dev.off()

library(sf)
#demo(nc, ask = FALSE, echo = FALSE)
#plot(nc["AREA"])
#View(nc)
regie<-st_read("Data/us_eco_l4 (1)/us_eco_l4_no_st.shx")
regie2<-left_join(regie,taxdatq1)

regiea <- regie2 %>% select(perc_invaded,local_similarity, geometry)
colnames(regiea)[2]<-"tax_similarity"

regie3<-left_join(regie,phydatq1)
regieb <- regie3 %>% select(perc_invaded,local_similarity, geometry)
colnames(regieb)[2]<-"phy_similarity"

regiec<-st_join(regiea, regieb,largest=TRUE)
regied<-regiec %>% gather(var,scale , -geometry)
regied<-filter(regied,var!="perc_invaded.y")
jpeg("Analyses/Master/Plots/maps1.jpeg",width = 10,height=9,unit='in',res=200)
ggplot() + 
  geom_sf(data = regied, aes(fill = scale),lwd = 0.01)+scale_fill_viridis_c(option="plasma")+
  ggthemes::theme_few()+theme(legend.position = "top")+
  facet_wrap(~var,ncol=3,labeller = labeller(var = c(`perc_invaded.x` = "% invaded",`phy_similarity`=" phylogenetic\nsimilarity",`tax_similarity`="taxonomic\nsimilarity")))
dev.off()

pp1<-ggplot() + 
  geom_sf(data = regiea, aes(fill = perc_invaded),lwd = 0.01)+scale_fill_viridis_c()+ggthemes::theme_few()+theme(legend.position = "top")

pp2<-    ggplot() + 
  geom_sf(data = regiea, aes(fill = tax_similarity),lwd = 0.01)+scale_fill_viridis_c(option = "A")+ggthemes::theme_few()+theme(legend.position = "top")

pp3<-    ggplot() + 
  geom_sf(data = regieb, aes(fill = phy_similarity),lwd = 0.01)+scale_fill_viridis_c(option = "C")+ggthemes::theme_few()+theme(legend.position = "top")

jpeg("Analyses/Master/Plots/maps2.jpeg",width = 10,height=9,unit='in',res=200)
ggpubr::ggarrange(pp1,pp2,pp3,nrow=1,ncol=3)
dev.off()

#### Naitive only
inyonly<-filter(L4tots,n_inv>1)
ecoregionz<-unique(inyonly$US_L4NAME)
ecoregionz<-intersect(ecoregionz,unique(d.invy$US_L4NAME))
ecoregionz<-ecoregionz[-461]

d.natty<-filter(d,NativeStatus=="N")
taxdat.nats<-data.frame()
phydat.nats<-data.frame()
for (z in c(1:length(ecoregionz))){
  dat<-dplyr::filter(d.natty,US_L4NAME==ecoregionz[z])
  
  comms<-matricize(dat)
  
  taxy0<-hill_taxa_parti(comms,q = 0)
  taxy1<-hill_taxa_parti(comms,q = 1)
  taxy2<-hill_taxa_parti(comms,q = 2)
  
  taxy<-rbind(taxy0,taxy1,taxy2)
  taxy$US_L4NAME<-ecoregionz[z]
  taxyclass<-"taxonomic"
  
  phy0<-hill_phylo_parti(comms,phy,q=0)
  phy1<-hill_phylo_parti(comms,phy,q=1)
  phy2<-hill_phylo_parti(comms,phy,q=2)
  
  
  phyt<-rbind(phy0,phy1,phy2)
  phyt$US_L4NAME<-ecoregionz[z]
  phyt$class<-"phylogenetic"
  
  taxdat.nats<-rbind(taxy,taxdat.nats)
  phydat.nats<-rbind(phyt,phydat.nats)
  
}  

taxnatq1<-filter(taxdat.nats,q==1)
head(taxnatq1)
taxnatq1<-left_join(taxnatq1,key)
taxnatq1<-filter(taxnatq1,NA_L1NAME %in% c('EASTERN TEMPERATE FORESTS','GREAT PLAINS','NORTH AMERICAN DESERTS','NORTHERN FORESTS','NORTHWESTERN FORESTED MOUNTAINS'))
taxnatq1<-left_join(taxnatq1,L4tots)
taxnatq1$n_scale<-taxnatq1$n/max(taxnatq1$n)
taxnatq1$perc_invaded<-taxnatq1$n_inv/taxnatq1$n
ggplot(taxnatq1,aes(perc_invaded,TD_alpha))+geom_smooth(method='lm')+facet_wrap(~NA_L1NAME)
