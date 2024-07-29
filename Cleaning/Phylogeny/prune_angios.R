rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())
graphics.off()

library(sf)
library(dplyr)
library("mapview")
library("tmap")
library(maps)
library(usmap)
library(phytools)
library(hillR)
library(brms)
library(tidybayes)
library(ggplot2)

setwd("~/Documents/git/bioticHogs/")


phy<-read.tree("Analyses/Master/Input/FINAL_PHYLOGENY_SpCodes_Apr2023.tr") #read in tree
library(ggtree)

ggtree(phy,layout="circular")+geom_nodelab()
phy$node.label
library("RRphylo")
t2<-extract.clade(phy, "Angiospermae")

ggtree(t2)+geom_nodelab()

write.tree(t2,"Analyses/Master/Input/angiosperms.tre")

