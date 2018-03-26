#set directory
#built from tutorial http://www.francoiskeck.fr/phylosignal/demo_plots.html
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

#change axes on trait 
max.x<- max(CNV_df$ITS_CN_AVE_ITS_LSU)
dotplot.phylo4d(p4d_named, 
                tree.ladderize = TRUE, 
                center = FALSE, trait = "CNV", 
                data.xlim=c(0, max.x), 
                scale = FALSE, 
                trait.cex = .8, 
                tip.cex = .6, 
                dot.cex = 1.5)




#make scatter plot with pairwise distance matrix
#make distance matrix (phylo.distance)
dist_phylo<- cophenetic(tre.w.names)

#make difference matrix (CN)
dist_CNV<- as.matrix(dist(CNV))

#put the matrices in the same order
dist_phylo_in_row_order<- dist_phylo[order(rownames(dist_phylo)),] 
dist_CNV_in_row_order<- dist_CNV[order(rownames(dist_CNV)),] 
dist_phylo_in_row_col_order<- dist_phylo_in_row_order[,order(colnames(dist_phylo_in_row_order))]
dist_CNV_in_row_col_order<- dist_CNV_in_row_order[,order(colnames(dist_CNV_in_row_order))]


#make plot 
plot(dist_phylo_in_row_col_order ~ dist_CNV_in_row_col_order)


#first, reshape matrices. note, this will also remove 0 values (self-comparisons), 
#also removes duplicates that are present as a consequence of the original distance matrix
dist_phylo_in_row_col_order_m<- data.frame(dist_phylo_in_row_col_order)
reshaped_phylo<- data.frame( t(combn(names(dist_phylo_in_row_col_order_m),2)), dist=t(dist_phylo_in_row_col_order_m)[lower.tri(dist_phylo_in_row_col_order_m)] )
dist_CNV_in_row_col_order_m<- data.frame(dist_CNV_in_row_col_order)
reshaped_CNV<- data.frame( t(combn(names(dist_CNV_in_row_col_order_m),2)), dist=t(dist_CNV_in_row_col_order_m)[lower.tri(dist_CNV_in_row_col_order_m)] )

#print to a file and code in color catagories by phylogeny 
write.table(reshaped_phylo, "reshaped_phylo.csv", append = TRUE, na = "NA", quote =FALSE, col.names = FALSE, row.names = FALSE, sep = ",")
write.table(reshaped_CNV, "reshaped_CNV.csv", append = TRUE, na = "NA", quote =FALSE, col.names = FALSE, row.names = FALSE, sep = ",")

#reload data 
color_catagories<- read.csv(file = "color_catagories.csv", sep = ",")

#plot with colors
#color code by belonging in one of four lists 

#Interspecific = catagory 1 = #CCCC66
#Intergeneric = catagory 2 = #FF9933
#Interphylum = catagory 3 = #006666 
#Interdomain = catagory 4 = #FF6633

#plot without jitter
#plot(color_catagories$diff_CNV, color_catagories$diff_phylo, 
#     pch =19, 
#     col = ifelse(color_catagories$color_cat== 1, '#CCCC66', 
#                  ifelse(color_catagories$color_cat== 2, '#FF9933', 
#                         ifelse(color_catagories$color_cat== 3, '#006666','#FF6633'))))

#plot with jitter
plot(color_catagories$diff_CNV, jitter(color_catagories$diff_phylo, 2), 
     pch =19, 
     cex = .8,
     xlab = "difference in copy number",
     ylab = "phylogenetic distance",
     col = ifelse(color_catagories$color_cat== 1, '#CCCC66', 
                  ifelse(color_catagories$color_cat== 2, '#FF9933', 
                         ifelse(color_catagories$color_cat== 3, '#006666','#FF6633'))))


#add regression lines for each catagory 
#Total
abline((lm(color_catagories$diff_phylo ~ color_catagories$diff_CNV)), lwd=5, col = '#666666')
all_stat<- lm(color_catagories$diff_phylo ~ color_catagories$diff_CNV)
all_stat_sum<- summary(all_stat)
r2_all = all_stat_sum$adj.r.squared
p_val_all = all_stat_sum$coefficients[2,4] 

rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2_all,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(p_val_all, digits = 2)))[2]

#interspecific
abline(lm(color_catagories$diff_phylo[color_catagories$color_cat == 1] ~ color_catagories$diff_CNV[color_catagories$color_cat == 1]), col = '#CCCC66', lwd=2, lty="dotted")
IS_stat<- lm(color_catagories$diff_phylo[color_catagories$color_cat == 1] ~ color_catagories$diff_CNV[color_catagories$color_cat == 1])
IS_stat_sum<- summary(IS_stat)
r2_IS = IS_stat_sum$adj.r.squared
p_val_IS = IS_stat_sum$coefficients[2,4] 

rp2 = vector('expression',2)
rp2[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2_IS,dig=3)))[2]
rp2[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(p_val_IS, digits = 2)))[2]

#intergeneric
abline(lm(color_catagories$diff_phylo[color_catagories$color_cat == 2] ~ color_catagories$diff_CNV[color_catagories$color_cat == 2]), col = '#FF9933', lwd=2, lty="dotted")
IG_stat<- lm(color_catagories$diff_phylo[color_catagories$color_cat == 2] ~ color_catagories$diff_CNV[color_catagories$color_cat == 2])
IG_stat_sum<- summary(IG_stat)
r3_IG = IG_stat_sum$adj.r.squared
p_val_IG = IG_stat_sum$coefficients[2,4] 

rp3 = vector('expression',2)
rp3[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                    list(MYVALUE = format(r3_IG,dig=3)))[2]
rp3[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                    list(MYOTHERVALUE = format(p_val_IG, digits = 2)))[2]


#interphilum
abline(lm(color_catagories$diff_phylo[color_catagories$color_cat == 3] ~ color_catagories$diff_CNV[color_catagories$color_cat == 3]), col = '#006666', lwd=2, lty="dotted")
IP_stat<- lm(color_catagories$diff_phylo[color_catagories$color_cat == 3] ~ color_catagories$diff_CNV[color_catagories$color_cat == 3])
IP_stat_sum<- summary(IP_stat)
r4_IP = IP_stat_sum$adj.r.squared
p_val_IP = IP_stat_sum$coefficients[2,4] 

rp4 = vector('expression',2)
rp4[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                    list(MYVALUE = format(r4_IP,dig=3)))[2]
rp4[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                    list(MYOTHERVALUE = format(p_val_IP, digits = 2)))[2]

#interdomain
abline(lm(color_catagories$diff_phylo[color_catagories$color_cat == 4] ~ color_catagories$diff_CNV[color_catagories$color_cat == 4]), col = '#FF6633', lwd=2, lty="dotted")
ID_stat<- lm(color_catagories$diff_phylo[color_catagories$color_cat == 4] ~ color_catagories$diff_CNV[color_catagories$color_cat == 4])
ID_stat_sum<- summary(ID_stat)
r5_ID = ID_stat_sum$adj.r.squared
p_val_ID = ID_stat_sum$coefficients[2,4] 

rp5 = vector('expression',2)
rp5[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                    list(MYVALUE = format(r5_ID,dig=3)))[2]
rp5[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                    list(MYOTHERVALUE = format(p_val_ID, digits = 2)))[2]

legend('topright', legend = c(rp5, rp4, rp, rp3, rp2), bty = 'n',
       text.col = c('#FF6633', '#FF6633', '#006666', '#006666', '#666666','#666666', '#FF9933', '#FF9933', '#CCCC66', '#CCCC66'))

