### Getting to know to R packages for beta diversity
##hillR uses hill numbers and has built ins for phylogenetic and functional diversity
##betapart allows for partitioning beta diversity into nestedness and turnover
###updated to include phylogenies and traits and distance decay model see:
#https://search.r-project.org/CRAN/refmans/betapart/html/betapart-package.html
#also had a good reference list
library(hillR)

dummy<-FD::dummy
dummy$trait
dummy$abun
### comes in as count data not percent
######################################################################
###Does it matter if you relativize then pool or pool the relativize?###
#####################################################################
hill_taxa_parti_pairwise(dummy$abun,q=2,rel_then_pool = TRUE)
hill_taxa_parti_pairwise(dummy$abun,q=2,rel_then_pool = FALSE)
### doesnt matter for this


#################################################################
### Q2: what if you have percent cover data rather than count data###
#############################################################
abun<-as.data.frame(dummy$abun)
abun$total<-rowSums(abun)
abun.rel<-abun[,1:8]/abun$total

hill_taxa_parti_pairwise(abun.rel,q=2,rel_then_pool = TRUE)
hill_taxa_parti_pairwise(dummy$abun,q=2,rel_then_pool = TRUE)
### you appear to get the same answer..yay

####Try it non-pairwise
hill_taxa_parti(dummy$abun,q = 2)
hill_taxa_parti(abun.rel, q=2)
#A data frame with one row (across all sites) and six columns:q, gamma diversity, alpha diversity,
#beta diversity, local similarity (similar to Sorensen), and region similarity (similar to Jaccard).

####we probably need to go through the paper in detail to actually understand what we are reporting

#Q how do you interpret 2:4: Effective species numbers?
###fundamentally what is q 

library("betapart")
#Computes 3 distance matrices accounting for the (i) balanced variation in abundances, (ii) abundance gradients, and (iii) total dissimilarity (i.e. the sum of both components)
###paritioning beta diversity between nestedness and turnover
#https://besjournals.onlinelibrary.wiley.com/doi/10.1111/j.2041-210X.2012.00224.x

####below it seems like this packages can only really deal with count data
#as beta.bray.gra zeros out for percent data

beta.pair.abund(abun.rel, index.family="bray")
beta.pair.abund(dummy$abun, index.family="bray")

### you can see when you look at all the sites
beta.multi.abund(dummy$abun, index.family="bray")
beta.multi.abund(abun.rel, index.family="bray")
abun.rel.count<-round(abun.rel*100)
beta.pair.abund(abun.rel.count, index.family="bray")

beta.pair(dummy$abun, index.family="sorensen") #would need to binarize doent run
