### ignited Jan 4 2024 by Dan
##goal: subsample SPCIS to include only plots with acceptable trait and phylogeny coverage

###get ready
rm(list=ls())
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())
setwd("/Users/danielbuonaiuto/Documents/git/bioticHogs")

library(tidyverse)
library(ape)

#read in data neededx
d1<-read.csv("Data/FULLDatabase_10272022.csv") ## SPCIS
trait.cov<-read.csv("Analyses/Master/Input/SPCIS_TRY_coverage_2.csv")###this is trait coverage
phy<-read.tree("Analyses/Master/Input/FINAL_PHYLOGENY_SpCodes_Apr2023.tr") ### phylogeny

d1<-d1 %>% #select only the most recent survey for resampled plots
  group_by(Plot) %>%
  filter(Year==max(Year)) %>%
  ungroup() 


unique(trait.cov$Trait)
trait.cov<-trait.cov %>% #select only the most recent survey for resampled plots
 group_by(Plot) %>%
  filter(Year==max(Year)) %>%
 ungroup() 

trait.cov<-dplyr::filter(trait.cov,Trait %in% c("SLA_mm2/mg","heightveg_m")) ##choose are the traits we want
trait.cov<-filter(trait.cov,trait_cov>=.8) 

##now filter species not in phylogeny
d<-filter(d1,SpCode %in% c(phy$tip.label)) ## filter out species not in phylogeny, mostly NAs

d<-filter(d,Plot %in% trait.cov$Plot) ##now filter out plots that dont have 80% trait coverage in both traits

##how many plots are we left with currently?
length(unique(d$Plot)) #78,412

##merge with trait data
d1a<-read.csv("Data/SPCIS_traits.csv")###this is trait data
d.trait<-dplyr::select(d1a,sps_try_match,SLA_mm2.mg,heightveg_m)
d.trait<-filter(d.trait,!is.na(SLA_mm2.mg))
d.trait<-filter(d.trait,!is.na(heightveg_m))
d.trait<-distinct(d.trait)


nameref<-select(d,SpCode,AcceptedTaxonName)
nameref<-distinct(nameref)
library(tibble)
colnames(d.trait)[1]<-"AcceptedTaxonName"
d.trait<-left_join(d.trait,nameref)
d.trait<-select(d.trait,-AcceptedTaxonName)

d.trait<-arrange(d.trait,SpCode)
d.trait<-filter(d.trait,!is.na(SpCode)) ## remove NAs

##remove old forms
rm(d1)
rm(d1a)
rm(nameref)
rm(trait.cov)

d<-left_join(d,d.trait)

##write a table of ecoregion equiventns
d.ref<-read.csv("Data/FULLDatabase_10272022_diversity_PctCov_100_10272022_INHABIT.csv")
d.ref<-select(d.ref,Plot, NA_L1NAME,NA_L2NAME,NA_L3NAME,US_L4NAME)
d.ref<-distinct(d.ref)
#write.csv(d.ref,"Analyses/Master/Input/ecoregions.csv")

colnames(d)[9]<-"NA_L1NAME" ### give the columns the same name


intersect(d.ref$NA_L1NAME,d$NA_L1NAME)
d<-left_join(d,d.ref)
d<-filter(d,!is.na(NA_L1NAME))### filter unatrtributed ecoregions
d<-filter(d,!is.na(SLA_mm2.mg))
d<-filter(d,!is.na(heightveg_m))

length(unique(d$Plot)) ##75,475 plots

##so below should output a version of SPECIS that has
#only species with SLA and height info that are in the phylogeny
#only plots that have 80% trait coverage
## all ecoregions

write.csv(d,"Analyses/Master/Input/specis_wtraitphylo.csv")
