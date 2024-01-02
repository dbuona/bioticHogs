##play with indicies to chatacterize 
rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

graphics.off()

library(dplyr)
library(vegan)
library("geodist")
library(ggplot2)
library(hillR)
library(betapart)
library(reshape2)

setwd("~PowellCenter")
set.seed(3)

d<-read.csv("Data/FULLDatabase_05272022.csv")
d.env<-read.csv("Data/FULLDatabase_05272022_plotsenv4Jul2022.csv")
d.ref<-read.csv("Data/FullDatabase_diversity_July4.csv")

bos<-dplyr::filter(d.env,US_L3NAME== "Driftless Area")
dat.bos<-filter(d, Plot %in% c(unique(bos$Plot)))
dat.bos<-filter(dat.bos, Resampled=="N")
driftless.ref<-filter(d.ref,Plot %in% c(unique(dat.bos$Plot)))
#rownames(driftless.ref)<-driftless.ref$Plot
driftless.ref<-select(driftless.ref,Plot,RelCov_I)






dater.high<-dplyr::select(dat.bos,Plot,AcceptedTaxonName,PctCov) ## rid of complilcating columns
dater.high<-dplyr::filter(dater.high,!is.na(dater.high$AcceptedTaxonName)) ## spreading the data to long for (matrix) won't work if there are NA's in taxon name so get rid of them here
range(dater.high$PctCov) #high 122
high<- tidyr::spread(dater.high,AcceptedTaxonName,PctCov)  ## spread to long       
high[is.na(high)] <- 0 ## make na's zeros
rownames(high)<-high$Plot # convert plot names to row names
high<-dplyr::select(high,-Plot) # remove plot column
high<-as.matrix.data.frame(high) # make it a matrix

high.hill<-hill_taxa_parti_pairwise(high,q=2,output = "matrix")

dim(high.hill$local_similarity)

dater2<-dplyr::select(dat.bos,Plot,Long,Lat) ## subset to lat longs
dater2<-distinct(dater2) ## remove dplicates casue each species was a row

#rownames(dater2.high)<-dater2.highPlot # convert plot names to row names
#dater2<-dplyr::select(dater2,-Plot)
setdiff(dater2$Plot,rownames(high))

dater2<-filter(dater2,Plot!="MISS_111")

b<-geodist(dater2,paired = TRUE,measure = "haversine")
dim(b)

inv.data<-data.frame(Plot=rep(rownames(high),nrow(high)))
dim(inv.data)
inv.data<-left_join(inv.data,driftless.ref)
inv.data2<-data.frame(Plot=rep(rownames(high),each=nrow(high)))
inv.data2<-left_join(inv.data2,driftless.ref)
dim(inv.data2)

inv.data<-cbind(inv.data,inv.data2)
colnames(inv.data)<-c("site1", "cov1","site2","cov2")

inv.data$cov1<-inv.data$cov1/100+.1
inv.data$cov2<-inv.data$cov2/100+.1
range(inv.data$cov1)

inv.data$reldiff<-abs((inv.data$cov1-inv.data$cov2)/((inv.data$cov1+inv.data$cov2)))






goo<-select(inv.data,site1,site2,reldiff)


goo<-acast(goo,site1~site2,value.var="reldiff")

dim(goo)


mydat <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(mydat) <- c('Distance','Sim','Reldiff','RowID', 'ColumnID')
for(i in 1:nrow(b)) {
  for(j in 1:ncol(b)) {
    if(i<j) { #only include the lower diagonal of the distance matrix
      mydat<-rbind(mydat,data.frame(
        Distance=b[i,j],# distance
        Sim=high.hill$local_similarity[i,j],
        Reldiff=goo[i,j],
        RowID=as.character(i), # row ID
        ColumnID=as.character(j))) # column ID
    }
  }
}

mod1<-glm(Sim~Distance*Reldiff,data=mydat,family=quasibinomial(link=logit))
summary(mod1)

mydat$Diff.z<-(mydat$Reldiff-mean(mydat$Reldiff))/sd(mydat$Reldiff)
mydat$Dist.z<-(mydat$Distance-mean(mydat$Distance))/sd(mydat$Distance)

mod1z<-glm(Sim~Dist.z*Diff.z,data=mydat,family=quasibinomial(link=logit))
summary(mod1z)



new.data<-data.frame(Distance=rep(c(0.00 , 16078.82,  74454.06 ,137064.15, 329217.91),each=5), Reldiff=rep(c(0,.5,1,1.5,2),5))
quantile(mydat$Dist.z)
quantile(mydat$Diff.z)
new.data.z<- data.frame(Dist.z=rep(c(-1.08259238, -0.88218901, -0.07548773,  0.59924142,  2.96791370 ),each=5), Diff.z=rep(c(-0.7760023, -0.7646318, -0.5801083,  0.6356005,  2.0804269 ),5),RowID=rep(2,25),ColumnID=rep(1,25))



predy<-predict(mod1, new.data, type="response")
plotty<-cbind(new.data,predy)

predy.z<-predict(mod1z, new.data.z, type="response")
plottyz<-cbind(new.data.z,predy.z)



predy2<-predict(mod2, new.data, type="response")
plotty2<-cbind(new.data,predy2)



cc<-ggplot()+
  geom_line(data=plotty,aes(Distance,predy,color=as.factor(Reldiff)),size=1.5)+
  ggthemes::theme_few()+ggtitle("quasibinomial")+scale_color_viridis_d("Invasion Difference")

ggplot()+
  geom_line(data=plottyz,aes(Dist.z,predy,color=as.factor(Diff.z)),size=1.5)+
  ggthemes::theme_few()+scale_color_viridis_d("Invasion Difference")



aa<-ggplot()+
  geom_line(data=plotty2,aes(Distance,predy2,color=as.factor(Reldiff)),size=1.5)+
  ggthemes::theme_few()+ggtitle("exponetial")+scale_color_viridis_d("Invasion Difference")

bb<-ggplot()+geom_point(data=mydat,aes(Distance,Sim))+  geom_line(data=plotty2,aes(Distance,predy2,color=as.factor(Reldiff)),size=1.5)+
  geom_smooth(data=plotty2,aes(Distance,predy,color=as.factor(Reldiff)),method="glm",method.args = list(family = "quasibinomial"))+
  ggthemes::theme_few()+ggtitle("exponetial")+scale_color_viridis_d("Invasion Difference")


jpeg("..//git/bioticHogs/plots/firstlook_modcomp.jpeg")
ggpubr::ggarrange(cc,aa,bb, common.legend = TRUE,labels=c("a)","b)","c)"),widths=c(.6,.6,.4))
dev.off()



library(brms)
test<-brm(Sim~Dist.z*Diff.z+(1|mm(RowID,ColumnID)),data=mydat,family=zero_inflated_beta())
range(mydat$Dist.z)

newdaters<-data.frame(Diff.z=rep(c(-0.8465717, 2.4832547),5),Dist.z=rep(c(-1,0,1,2,2.5),each=2))

predyz<-fitted(test,newdata = newdaters,re_formula = NA)
plotz<-cbind(newdaters,predyz)
ggplot(plotz,aes(Dist.z,Estimate))+geom_line(aes(color=as.factor(Diff.z)))

output<-cbind(mydat,fitted(test))
ggplot(output,aes(Dist.z,Estimate))+geom_smooth()

summary(test)







jpeg("../git/bioticHogs/plots/ppcheck_zibeta.jpeg")
pp_check(test,ndraws=100)
dev.off()
jpeg("../git/bioticHogs/plots/mm_output.jpeg")
plot()

dev.off()


jpeg("../git/bioticHogs/plots/conditionaleffect_zibeta.jpeg")
pp<-conditional_effects(test,effect="Dist.z:Diff.z",prob = .5, int_conditions = list(Diff.z = quantile))


plot(pp,plot=FALSE)[[1]]+
  scale_color_viridis_d(direction=-1)+scale_fill_viridis_d(direction = -1)+ggthemes::theme_few()
dev.off()


save.image("../git/bioticHogs/drifless.Rda")







