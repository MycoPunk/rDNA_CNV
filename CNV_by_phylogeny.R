#set directory
#built from tutorial http://www.francoiskeck.fr/phylosignal/demo_plots.html
setwd("~/Desktop/Project_ITS_CNV/CNV_Tree_NEW")

#set libraries
library(phylosignal)
library(adephylo)
library(ape)
library(phylobase)

##read in file
CNV_df<- read.csv("NEXUS_for_phylosignal_with_branchlengths_no_outlier.csv")

#attach species names to train numbers (ITS CN)
CNV<- CNV_df$rDNA_CN
names(CNV) = CNV_df$X

#read in the tree
nexus<- as.character(CNV_df[1,3])

#make the .phylo file
tre.w.names<- read.tree(text = c(nexus, CNV))

#create a dataframe of three traits for the 34 species ITS, Random values, Simulated values under a Brownian Motion model along the tree
dat.w.names <- data.frame(CNV)
dat.w.names$random <- rnorm(dim(dat.w.names)[1], sd = 10)
dat.w.names$bm <- rTraitCont(tre.w.names)

#create .p4d object
p4d_named <- phylo4d(tre.w.names, dat.w.names)

#plot phylo4d object
barplot(p4d_named, 
        trait.cex = .8, 
        tip.cex = .3, 
        tree.ladderize = TRUE)
dotplot(p4d_named, 
        trait.cex = .8, 
        tip.cex = .3, 
        tree.ladderize = TRUE)
gridplot(p4d_named, 
         trait.cex = .8, 
         tip.cex = .3, 
         tree.ladderize = TRUE)

#phyloweights 
weights.df<- phyloWeights(p4d_named, dist.phylo = "patristic", 
                          method= "inverse", 
                          #alpha = 3, 
                          dmax = 3)
weights.df


#plot which traits to plot and their order 
#take off negative scale
barplot(p4d_named, center = FALSE, trait = c("CNV"))

#stats time
phyloSignal(p4d_named)

#compare random ...
corr.random<- phyloCorrelogram(p4d_named, trait = "random")
plot(corr.random, main = "random-no outlier")

#to bm (should be positive and significant)
corr.bm<- phyloCorrelogram(p4d_named, trait = "bm")
plot(corr.bm, main = "BM-no outlier")

#..to the data
corr.data<- phyloCorrelogram(p4d_named, trait = "CNV")
plot(corr.data, main = "CNV-no outlier")

#get phylogenetic distance matrix
distance<-adephylo::distTips(p4d_named)
dist_matrix<- as.matrix(distance)
View(dist_matrix)




#plot with CNV trait and flipped nodes to make pretty
dotplot(p4d_named, tree.ladderize = TRUE, center = FALSE, trait = "CNV")

#change axes on trait 
max.x<- max(CNV_df$rDNA_CN)
dotplot.phylo4d(p4d_named, 
                tree.ladderize = TRUE, 
                center = FALSE, trait = "CNV", 
                data.xlim=c(0, max.x), 
                scale = FALSE, 
                trait.cex = .8, 
                tip.cex = .3, 
                dot.cex = 1.5)




#compute lipaMoran I
local.i<- lipaMoran(p4d_named, 
                    trait = "CNV", 
                    prox.phylo = "nNodes",
                    as.p4d = TRUE)


points.col<- lipaMoran(p4d_named, trait = "CNV", 
                       prox.phylo = "nNodes")$p.value


points.col<-ifelse(points.col<0.05, "red", "black")
dotplot.phylo4d(local.i, dot.col = points.col, main = "no_outlier", 
                tree.ladderize = TRUE, 
                scale = FALSE, 
                trait.cex = .8, 
                tip.cex = .3, 
                dot.cex = 1.5)
