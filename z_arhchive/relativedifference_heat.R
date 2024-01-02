##play with indicies to chatacterize 
rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
graphics.off()
library(tidyverse)


d<-data.frame(cov1=rep(seq(.1,1.1,by=0.1),each=11),cov2=rep(seq(.1,1.1,by=0.1),11))


d$reldiff<-abs((d$cov1-d$cov2)/.5*((d$cov1+d$cov2)))




#d$reldiff<-ifelse(d$cov1==0,d$cov2*1.6,abs((d$cov1-d$cov2))/((d$cov1+d$cov2)))

#fixcov<-which(d$cov2==0)
#d$reldiff[fixcov] <-d$cov1[fixcov]*1.6

#d$reldiffw<-ifelse(d$reldiffx==0,-.047*d$cov1,d$reldiffx)




dx<-acast(d,cov1~cov2,value.var="reldiff")



#dx[1,]<-dx[2,]-(dx[2,4]-dx[2,5])
#dx[,1]<-dx[,2]-(dx[2,3]-dx[2,4])


jpeg("..//bioticHogs/plots/reldiff_heat.jpeg")
  dx %>% 
  as.data.frame() %>%
  rownames_to_column("cov1") %>%
  pivot_longer(-c(cov1), names_to = "cov2", values_to = "reldiff") %>%
  ggplot(aes(x=cov2, y=cov1, fill=reldiff)) + 
  geom_raster()+
  scale_fill_viridis_c()
dev.off()

dy<-dx
dy[1,1]<--0
dy[2,2]<--0.02
dy[3,3]<--0.04
dy[4,4]<--0.06
dy[5,5]<--0.08
dy[6,6]<--0.1
dy[7,7]<--0.12
dy[8,8]<--0.14
dy[9,9]<--0.16
dy[10,10]<--.18
dy[11,11]<--.2



0.047/10
jpeg("..//git/bioticHogs/plots/reldiff_heat_zeroadj.jpeg")
dy %>% 
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
