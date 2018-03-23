#set directory
setwd("~/Desktop/CNV_Tree")

#set libraries
library(phylosignal)
library(adephylo)
library(ape)
library(phylobase)

##read in file
CNV_df<- read.csv("CNV_for_tree.csv")

#attach species names to train numbers (ITS CN)
CNV<- CNV_df$ITS_CN_AVE_ITS_LSU
names(CNV) = CNV_df$X

#read in the tree
nexus<- as.character(CNV_df[1,6])

#make the .phylo file
tre.w.names<- read.tree(text = c(nexus, CNV))

#create a dataframe of three traits for the 34 species ITS, Random values, Simulated values under a Brownian Motion model along the tree
dat.w.names <- data.frame(CNV)
dat.w.names$random <- rnorm(dim(dat.w.names)[1], sd = 10)
dat.w.names$bm <- rTraitCont(tre.w.names)

#create .p4d object
p4d_named <- phylo4d(tre.w.names, dat.w.names)

#plot phylo4d object
barplot(p4d_named)
dotplot(p4d_named)
gridplot(p4d_named)

#to control which traits to plot and their order 
barplot(p4d_named, trait = c("CNV"))

#take off negative scale
barplot(p4d_named, center = FALSE, trait = c("CNV"))

#stats time
phyloSignal(p4d_named)

#compare random ...
corr.random<- phyloCorrelogram(p4d_named, trait = "random")
plot(corr.random)

#to bm (shoule be positive and significant)
corr.bm<- phyloCorrelogram(p4d_named, trait = "bm")
plot(corr.bm)

#..to the data
corr.data<- phyloCorrelogram(p4d_named, trait = "CNV")
plot(corr.data)

#plot with CNV trait and flipped nodes to make pretty
dotplot(p4d_named, tree.ladderize = TRUE, center = FALSE, trait = "CNV")

#change graph peram. 
max.x<- max(CNV_df$ITS_CN_AVE_ITS_LSU)
dotplot.phylo4d(p4d_named, 
                tree.ladderize = TRUE, 
                center = FALSE, trait = "CNV", 
                data.xlim=c(0, max.x), 
                scale = FALSE, 
                trait.cex = .8, 
                tip.cex = .6, 
                dot.cex = 1.5)


