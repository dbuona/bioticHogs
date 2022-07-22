library(tidyverse)


d<-data.frame(cov1=rep(seq(0,1,by=0.1),each=11),cov2=rep(seq(0,1,by=0.1),11))

d$reldiffx<-abs(d$cov1-d$cov2)/(abs(d$cov1+d$cov2)/2)
d$reldiffw<-abs(d$cov1-d$cov2)/(abs(d$cov1+d$cov2)/2)

d$reldiff<-ifelse(d$cov1==0,d$cov2*2,(abs(d$cov1-d$cov2))/(abs(d$cov1+d$cov2)/2))

fixcov<-which(d$cov2==0)
d$reldiff[fixcov] <-d$cov1[fixcov]*2




dx<-acast(d,cov1~cov2,value.var="reldiffx")

d2<-acast(d,cov1~cov2,value.var="reldiff")
dw<-acast(d,cov1~cov2,value.var="reldiffw")
8

jpeg("..//git/bioticHogs/plots/xreldiff_heat_x.jpeg")
  dx %>% 
  as.data.frame() %>%
  rownames_to_column("cov1") %>%
  pivot_longer(-c(cov1), names_to = "cov2", values_to = "reldiff") %>%
  ggplot(aes(x=cov2, y=cov1, fill=reldiff)) + 
  geom_raster()+
  scale_fill_viridis_c()
dev.off()

jpeg("..//git/bioticHogs/plots/reldiff_heat_good.jpeg")
d2 %>% 
  as.data.frame() %>%
  rownames_to_column("cov1") %>%
  pivot_longer(-c(cov1), names_to = "cov2", values_to = "reldiff") %>%
  ggplot(aes(x=cov2, y=cov1, fill=reldiff)) + 
  geom_raster()+
  scale_fill_viridis_c()
dev.off()


dw %>% 
  as.data.frame() %>%
  rownames_to_column("cov1") %>%
  pivot_longer(-c(cov1), names_to = "cov2", values_to = "reldiff") %>%
  ggplot(aes(x=cov2, y=cov1, fill=reldiff)) + 
  geom_raster()+
  scale_fill_viridis_c()
dev.off()
