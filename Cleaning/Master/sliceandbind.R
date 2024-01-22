##lit up by Dan Jan 18 2024
##Goal break up the big data sets into 1
##balance invasive and univaded plots based on richness

###get ready
rm(list=ls())
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())


library(data.table)

library(tidyverse)
library(broom)
library(dplyr)
library(brms)
setwd("/Users/danielbuonaiuto/Documents/git/bioticHogs/Analyses/Master/Input/L3/invaded")

temp.invaded<- list.files(pattern="*.csv")
inv.dat<- rbindlist(sapply(temp.invaded, fread,simplify = FALSE))

setwd("/Users/danielbuonaiuto/Documents/git/bioticHogs/Analyses/Master/Input/L3/uninvaded")

temp.uninvaded<- list.files(pattern="*.csv")
nat.dat<- rbindlist(sapply(temp.uninvaded, fread,simplify = FALSE))

### first seperarte the qs
inv.dat$status<-"invaded"
nat.dat$status<-"native"



q0_invaded<-filter(inv.dat,q=="0")
q1_invaded<-filter(inv.dat,q=="1")
q2_invaded<-filter(inv.dat,q=="2")


q0_invaded_tax<-filter(q0_invaded,class=="taxonomic")
q0_invaded_phy<-filter(q0_invaded,class=="phylogenetic")
q0_invaded_height<-filter(q0_invaded,class=="height")
q0_invaded_sla<-filter(q0_invaded,class=="SLA")

q1_invaded_tax<-filter(q1_invaded,class=="taxonomic")

q1_invaded_phy<-filter(q1_invaded,class=="phylogenetic")
q1_invaded_height<-filter(q1_invaded,class=="height")
q1_invaded_sla<-filter(q1_invaded,class=="SLA")

q2_invaded_tax<-filter(q2_invaded,class=="taxonomic")
q2_invaded_phy<-filter(q2_invaded,class=="phylogenetic")
q2_invaded_height<-filter(q2_invaded,class=="height")
q2_invaded_sla<-filter(q2_invaded,class=="SLA")

head(q2_invaded_sla)


q0_native<-filter(nat.dat,q=="0")
q1_native<-filter(nat.dat,q=="1")
q2_native<-filter(nat.dat,q=="2")


q0_native_tax<-filter(q0_native,class=="taxonomic")
q0_native_phy<-filter(q0_native,class=="phylogenetic")
q0_native_height<-filter(q0_native,class=="height")
q0_native_sla<-filter(q0_native,class=="SLA")

q1_native_tax<-filter(q1_native,class=="taxonomic")

q1_native_phy<-filter(q1_native,class=="phylogenetic")
q1_native_height<-filter(q1_native,class=="height")
q1_native_sla<-filter(q1_native,class=="SLA")

q2_native_tax<-filter(q2_native,class=="taxonomic")
q2_native_phy<-filter(q2_native,class=="phylogenetic")
q2_native_height<-filter(q2_native,class=="height")
q2_native_sla<-filter(q2_native,class=="SLA")



rm(inv.dat)
rm(nat.dat)
###make the data frame
q0_tax<-rbind(q0_native_tax,q0_invaded_tax)
q1_tax<-rbind(q1_native_tax,q1_invaded_tax)
q2_tax<-rbind(q2_native_tax,q2_invaded_tax)

q0_phy<-rbind(q0_native_phy,q0_invaded_phy)
q1_phy<-rbind(q1_native_phy,q1_invaded_phy)
q2_phy<-rbind(q2_native_phy,q2_invaded_phy)

q0_height<-rbind(q0_native_height,q0_invaded_height)
q1_height<-rbind(q1_native_height,q1_invaded_height)
q2_height<-rbind(q2_native_height,q2_invaded_height)


q0_sla<-rbind(q0_native_sla,q0_invaded_sla)
q1_sla<-rbind(q1_native_sla,q1_invaded_sla)
q2_sla<-rbind(q2_native_sla,q2_invaded_sla)

write.csv(q0_tax,"..//q0_tax.csv")
write.csv(q1_tax,"..//q1_tax.csv")
write.csv(q2_tax,"..//q2_tax.csv")

write.csv(q0_phy,"..//q0_phy.csv")
write.csv(q1_phy,"..//q1_phy.csv")
write.csv(q2_phy,"..//q2_phy.csv")


write.csv(q0_height,"..//q0_height.csv")
write.csv(q1_height,"..//q1_height.csv")
write.csv(q2_height,"..//q2_height.csv")

write.csv(q0_sla,"..//q0_sla.csv")
write.csv(q1_sla,"..//q1_sla.csv")
write.csv(q2_sla,"..//q2_sla.csv")

