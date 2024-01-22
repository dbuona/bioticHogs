##lit up by Dan Jan 4 2024
##Goal recalcuate plot richness, and relcover
##balance invasive and univaded plots based on richness

###get ready
rm(list=ls())
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())
setwd("/Users/danielbuonaiuto/Documents/git/bioticHogs")

library(tidyverse)

d<-read.csv("Analyses/Master/Input/specis_wtraitphylo.csv")
colnames(d)
##calculate the richness of each plot
richy<-d%>% group_by(Plot,NativeStatus,NA_L1NAME,NA_L2NAME,NA_L3NAME,US_L4NAME)%>%count() %>%spread(NativeStatus,n)
richy<-select(richy,-`<NA>`)
colnames(richy)[6:8]<-paste("Richness",colnames(richy)[6:8],sep="_")
#richy[is.na(richy)] <- 0
richy[, 6:8][is.na(richy[, 6:8])] <- 0

richy$Richness_total<-richy$Richness_I+richy$Richness_N+richy$Richness_NI
richy$Status<-ifelse(richy$Richness_I==0,0,1)

##remove plots that have only 1 species at it
richy<-filter(richy,Richness_total>1) 
richy<- richy[complete.cases(richy), ]#69,873 plots remainin
table(richy$Status)

##now balance invaded and uninvaded plots within ecoregions
##startwith L4
library(MatchIt)
set.seed(111)
msimp<-matchit(Status~US_L4NAME+Richness_total,data=richy,method = "nearest",distance="glm",exact=~US_L4NAME+Richness_total)
dat.msim<-match.data(msimp) #25,350 plots
table(dat.msim$NA_L1NAME)

set.seed(111)
msimpL3<-matchit(Status~NA_L3NAME+Richness_total,data=richy,method = "nearest",distance="glm",exact=~NA_L3NAME+Richness_total)
dat.msimL3<-match.data(msimpL3) #32764
table(dat.msimL3$NA_L1NAME)

counter<-dat.msimL3 %>%group_by(NA_L3NAME)%>% count()
counter2<-dat.msim %>%group_by(US_L4NAME)%>% count()

write.csv(dat.msimL3,"Analyses/Master/Input/balancedL3_guide.csv")
write.csv(dat.msim,"Analyses/Master/Input/balancedL4_guide.csv")
