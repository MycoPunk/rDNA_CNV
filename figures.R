setwd("")
library(data.table)
library(outliers)
#hash second line work with full data, unhash to work with data minus the outlier 
rDNA_by_taxa<- read.csv("genome_size.csv", sep = ",", header = TRUE)
rDNA_by_taxa<- rDNA_by_taxa[ rDNA_by_taxa$Genus!="Basidiobolus",]

#mean CN
mean(rDNA_by_taxa$CN)
median(rDNA_by_taxa$CN)

#min and max
range(rDNA_by_taxa$CN)
#n by phyla
sum(rDNA_by_taxa$group == "Ascomycota")
sum(rDNA_by_taxa$group == "Basidiomycota")
sum(rDNA_by_taxa$group == "Lower")


#suillus only
suillus<- rDNA_by_taxa[rDNA_by_taxa$Genus == "Suillus",]
#no brev.
suillus.wo.brev<- suillus[suillus$SE != "brevipes",]
range(suillus.wo.brev$CN)


#figure 1.d (Phylum)
boxplot(rDNA_by_taxa$CN ~rDNA_by_taxa$group, main = "no_outlier") 
cochran.test(CN~group, rDNA_by_taxa, inlying=FALSE) 
aov<- aov(rDNA_by_taxa$CN~rDNA_by_taxa$group, rDNA_by_taxa)
summary(aov(rDNA_by_taxa$CN~rDNA_by_taxa$group, rDNA_by_taxa))
hsd<- TukeyHSD(aov)
hsd

###Fig2
#figure 2.a by trophic mode 
#order goups for plotting
rDNA_by_taxa$TrophicMode=factor(rDNA_by_taxa$TrophicMode,c("Pathotroph","Saprotroph","Symbiotroph","Multi"))

#get n
b <- boxplot(CN ~TrophicMode, data=rDNA_by_taxa, plot=0) 
#set pars
par(mar=c(7,5,1,1))
par(cex.axis=0.7)
boxplot(CN ~TrophicMode, data=rDNA_by_taxa, 
        range = 3, 
        outline=TRUE, 
        names=paste(b$names, "(n=", b$n, ")"), 
        las = 2, 
        main = "Trophic mode", 
        col = "#5D2D8B", 
        ylim = c(0,350))

#Fig. 2.b Guild for Ascos, without outlier 
#make subset dataframe of just ascos 
#see how many catagories have 5 or more 
Asco_df<- rDNA_by_taxa[rDNA_by_taxa$Phylum == "Ascomycota",]
b <- boxplot(CN ~ Lifestyle, data=Asco_df, plot=0)
par(mar=c(7,5,1,1))
par(cex.axis=0.7)
boxplot(CN ~Lifestyle, data=Asco_df, 
        range = 3, 
        outline=TRUE, 
        names=paste(b$names, "(n=", b$n, ")"), 
        las = 2, 
        main = "CN by trophic mode - no outlier")

#only Pathogen, SAP S/L/O and PAthogen / SAP S/L/O do 
#make new DF
#Asco_df_with_enough_reps<- Asco_df[Asco_df$Lifestyle == "SAP S/L/O" | Asco_df$Lifestyle == "Pathogen" | Asco_df$Lifestyle == "SAP S/L/O / pathogen",]
Asco_df_with_enough_reps<- droplevels(subset(Asco_df[(Asco_df$Lifestyle == "SAP S/L/O") | (Asco_df$Lifestyle == "Pathogen") | (Asco_df$Lifestyle == "SAP S/L/O / pathogen"),]))

c <- boxplot(Asco_df_with_enough_reps$CN ~ Asco_df_with_enough_reps$Lifestyle, data=Asco_df_with_enough_reps, plot=0)
par(mar=c(7,5,1,1))
par(cex.axis=0.7)
boxplot(Asco_df_with_enough_reps$CN ~ Asco_df_with_enough_reps$Lifestyle, data=Asco_df_with_enough_reps, 
        range = 3, 
        outline=TRUE, 
        names=paste(c$names, "(n=", c$n, ")"), 
        las = 2, 
        main = "Ascomycota - no outlier")

cochran.test(CN ~ Lifestyle, data=Asco_df_with_enough_reps, inlying=FALSE)

#stripchart of above with no sap SLO/pathogen 
how_many<- boxplot(Asco_df_with_enough_reps$CN ~ Asco_df_with_enough_reps$Lifestyle, data=Asco_df_with_enough_reps, plot=0)

Asco_df_with_enough_reps<- droplevels(subset(Asco_df[(Asco_df$Lifestyle == "SAP S/L/O") | (Asco_df$Lifestyle == "Pathogen"),]))
stripchart(Asco_df_with_enough_reps$CN ~ Asco_df_with_enough_reps$Lifestyle, 
           main = "Guild - Ascomycota", 
           method="jitter", 
           pch = 16, 
           cex = 1.8, 
           vertical=TRUE, 
           col = "#006666", 
           at=c(1.1,1.9), 
           ylim = c(0, 350))

how_many$n

#Fig. 2.c Guild for Basids, without outlier 
#see how many catagories have 5 or more BASIDS 
Basid_df<- rDNA_by_taxa[rDNA_by_taxa$Phylum == "Basidiomycota",]
b <- boxplot(CN ~ Lifestyle, data=Basid_df, plot=0)
par(mar=c(7,5,1,1))
par(cex.axis=0.7)
boxplot(CN ~Lifestyle, data=Basid_df, 
        range = 3, 
        outline=TRUE, 
        names=paste(b$names, "(n=", b$n, ")"), 
        las = 2, 
        main = "CN by trophic mode - no outlier")




#only ECM, SAP BR and SAP WR 
#make new DF
Basid_df_with_enough_reps<- droplevels(subset(Basid_df[(Basid_df$Lifestyle == "ECM") | (Basid_df$Lifestyle == "SAP BR") | (Basid_df$Lifestyle == "SAP WR"),]))

#strip chart 
Basid_df_with_enough_reps$Lifestyle <- as.character(Basid_df_with_enough_reps$Lifestyle)
Basid_df_with_enough_reps.2<- replace(Basid_df_with_enough_reps$Lifestyle, Basid_df_with_enough_reps$Lifestyle=="SAP BR", "Wood Rot")
Basid_df_with_enough_reps.3<- replace(Basid_df_with_enough_reps.2, Basid_df_with_enough_reps.2=="SAP WR", "Wood Rot")
Basid_df_with_enough_reps.4<- cbind(Basid_df_with_enough_reps, Basid_df_with_enough_reps.3)

stripchart(Basid_df_with_enough_reps.4$CN ~ Basid_df_with_enough_reps.4$Basid_df_with_enough_reps.3, 
           main = "Guild - Basidiomycota", 
           method="jitter", 
           pch = 16, 
           cex = 1.8, 
           vertical=TRUE, 
           col = "#CCCB64", 
           at=c(1.1,1.9), 
           ylim = c(0,350))

#get n
sum(Basid_df_with_enough_reps.4$Basid_df_with_enough_reps.3 == "ECM")
sum(Basid_df_with_enough_reps.4$Basid_df_with_enough_reps.3 == "Wood Rot")


#run stats
#2.a by Trophic mode 
cochran.test(CN ~ TrophicMode, data=rDNA_by_taxa, inlying=FALSE)
#cochrans is significant. have to run KW test rather than anova. 
#summary(aov(CN ~TrophicMode, data=rDNA_by_taxa))
kruskal.test(CN ~TrophicMode, data=rDNA_by_taxa) # kw not significant. 


#2.b by Trophic mode ASCOS
cochran.test(CN ~ Lifestyle, data=Asco_df_with_enough_reps, inlying=FALSE)
#not sig, run anova
summary(aov(CN ~ Lifestyle, data=Asco_df_with_enough_reps))
#not sig dif. 

#2.c by Trophic mode BASIDS
cochran.test(CN ~ Lifestyle, data=Basid_df_with_enough_reps, inlying=FALSE)
#not sig, run anova
summary(aov(CN ~ Lifestyle, data=Basid_df_with_enough_reps))
#not sig dif.



###For fig S1
setwd()
library(data.table)
library(outliers)
rDNA_by_taxa<- read.csv("genome_size.csv", sep = ",", header = TRUE)
rDNA_by_taxa<- rDNA_by_taxa[ rDNA_by_taxa$Genus!="Basidiobolus",] #remove outlier

#at the species level (collapse species with more than one representation into a random representation of that species)
get_HSD_species <- function() { 
  split.df.by.SE<- split(rDNA_by_taxa, rDNA_by_taxa$SE)
  randomRows = function(df,n){
    return(df[sample(nrow(df),n),])
  }
  single_species_representation.df <-lapply(names(split.df.by.SE), function(x) randomRows(df = split.df.by.SE[[x]], n = 1))
  single_species_representation.df.df <- data.frame(matrix(unlist(single_species_representation.df), nrow=length(single_species_representation.df), byrow=T))
  aov<- aov(single_species_representation.df.df$X12~single_species_representation.df.df$X13, single_species_representation.df.df)
  summary(aov)
  hsd<- TukeyHSD(aov)
  TK<-(hsd)
  TK_data<-as.data.frame(TK[1:1])
  results<-data.frame(TK_data[1:3, 4], row.names = row.names(TK_data))
  return(results)
}

#run 1000 times 
AOV_results_df_SE<- as.data.frame(replicate(10000, get_HSD_species(), simplify = "array"))

#get results of the table, numeric
Basidiomycota_Ascomycota_SE<- sum(AOV_results_df_SE[1,1:ncol(AOV_results_df_SE)] < 0.05)
Lower_Ascomycota_SE<- sum(AOV_results_df_SE[2,1:ncol(AOV_results_df_SE)] < 0.05)
Lower_Basidiomycota_SE <- sum(AOV_results_df_SE[3,1:ncol(AOV_results_df_SE)] < 0.05)

Basidiomycota_Ascomycota_SE
Lower_Ascomycota_SE
Lower_Basidiomycota_SE

#at the genus level
get_HSD_genus <- function() { 
  split.df.by.genus<- split(rDNA_by_taxa, rDNA_by_taxa$Genus)
  randomRows = function(df,n){
    return(df[sample(nrow(df),n),])
  }
  single_genus_representation.df <-lapply(names(split.df.by.genus), function(x) randomRows(df = split.df.by.genus[[x]], n = 1))
  single_genus_representation.df.df <- data.frame(matrix(unlist(single_genus_representation.df), nrow=length(single_genus_representation.df), byrow=T))
  aov<- aov(single_genus_representation.df.df$X12~single_genus_representation.df.df$X13, single_genus_representation.df.df)
  summary(aov)
  hsd<- TukeyHSD(aov)
  TK<-(hsd)
  TK_data<-as.data.frame(TK[1:1])
  results<-data.frame(TK_data[1:3, 4], row.names = row.names(TK_data))
  return(results)
}


AOV_results_df_genus<- as.data.frame(replicate(10000, get_HSD_genus(), simplify = "array"))

#get results of the table, numeric
Basidiomycota_Ascomycota_genus<- sum(AOV_results_df_genus[1,1:ncol(AOV_results_df_genus)] < 0.05)
Lower_Ascomycota_genus<- sum(AOV_results_df_genus[2,1:ncol(AOV_results_df_genus)] < 0.05)
Lower_Basidiomycota_genus <- sum(AOV_results_df_genus[3,1:ncol(AOV_results_df_genus)] < 0.05)

Basidiomycota_Ascomycota_genus
Lower_Ascomycota_genus
Lower_Basidiomycota_genus

#at the family level
get_HSD_family <- function() { 
  split.df.by.family<- split(rDNA_by_taxa, rDNA_by_taxa$Family)
  randomRows = function(df,n){
    return(df[sample(nrow(df),n),])
  }
  single_family_representation.df <-lapply(names(split.df.by.family), function(x) randomRows(df = split.df.by.family[[x]], n = 1))
  single_family_representation.df.df <- data.frame(matrix(unlist(single_family_representation.df), nrow=length(single_family_representation.df), byrow=T))
  aov<- aov(single_family_representation.df.df$X12~single_family_representation.df.df$X13, single_family_representation.df.df)
  summary(aov)
  hsd<- TukeyHSD(aov)
  TK<-(hsd)
  TK_data<-as.data.frame(TK[1:1])
  results<-data.frame(TK_data[1:3, 4], row.names = row.names(TK_data))
  return(results)
}

AOV_results_df_family<- as.data.frame(replicate(10000, get_HSD_family(), simplify = "array"))

#get results of the table, numeric
Basidiomycota_Ascomycota_family<- sum(AOV_results_df_family[1,1:ncol(AOV_results_df_family)] < 0.05)
Lower_Ascomycota_family<- sum(AOV_results_df_family[2,1:ncol(AOV_results_df_family)] < 0.05)
Lower_Basidiomycota_family <- sum(AOV_results_df_family[3,1:ncol(AOV_results_df_family)] < 0.05)

Basidiomycota_Ascomycota_family
Lower_Ascomycota_family
Lower_Basidiomycota_family


#at the order level
get_HSD_order <- function() { 
  split.df.by.order<- split(rDNA_by_taxa, rDNA_by_taxa$Order)
  randomRows = function(df,n){
    return(df[sample(nrow(df),n),])
  }
  single_order_representation.df <-lapply(names(split.df.by.order), function(x) randomRows(df = split.df.by.order[[x]], n = 1))
  single_order_representation.df.df <- data.frame(matrix(unlist(single_order_representation.df), nrow=length(single_order_representation.df), byrow=T))
  aov<- aov(single_order_representation.df.df$X12~single_order_representation.df.df$X13, single_order_representation.df.df)
  summary(aov)
  hsd<- TukeyHSD(aov)
  TK<-(hsd)
  TK_data<-as.data.frame(TK[1:1])
  results<-data.frame(TK_data[1:3, 4], row.names = row.names(TK_data))
  return(results)
}


AOV_results_df_order<- as.data.frame(replicate(10000, get_HSD_order(), simplify = "array"))

#get results of the table, numeric
Basidiomycota_Ascomycota_order<- sum(AOV_results_df_order[1,1:ncol(AOV_results_df_order)] < 0.05)
Lower_Ascomycota_order<- sum(AOV_results_df_order[2,1:ncol(AOV_results_df_order)] < 0.05)
Lower_Basidiomycota_order <- sum(AOV_results_df_order[3,1:ncol(AOV_results_df_order)] < 0.05)

Basidiomycota_Ascomycota_order
Lower_Ascomycota_order
Lower_Basidiomycota_order

#at the class level 
get_HSD_class <- function() { 
  split.df.by.class<- split(rDNA_by_taxa, rDNA_by_taxa$Class)
  randomRows = function(df,n){
    return(df[sample(nrow(df),n),])
  }
  single_class_representation.df <-lapply(names(split.df.by.class), function(x) randomRows(df = split.df.by.class[[x]], n = 1))
  single_class_representation.df.df <- data.frame(matrix(unlist(single_class_representation.df), nrow=length(single_class_representation.df), byrow=T))
  aov<- aov(single_class_representation.df.df$X12~single_class_representation.df.df$X13, single_class_representation.df.df)
  summary(aov)
  hsd<- TukeyHSD(aov)
  TK<-(hsd)
  TK_data<-as.data.frame(TK[1:1])
  results<-data.frame(TK_data[1:3, 4], row.names = row.names(TK_data))
  return(results)
}


AOV_results_df_class<- as.data.frame(replicate(10000, get_HSD_class(), simplify = "array"))

#get results of the table, numeric
Basidiomycota_Ascomycota_class<- sum(AOV_results_df_class[1,1:ncol(AOV_results_df_class)] < 0.05)
Lower_Ascomycota_class<- sum(AOV_results_df_class[2,1:ncol(AOV_results_df_class)] < 0.05)
Lower_Basidiomycota_class <- sum(AOV_results_df_class[3,1:ncol(AOV_results_df_class)] < 0.05)

Basidiomycota_Ascomycota_class
Lower_Ascomycota_class
Lower_Basidiomycota_class


#mosaic plots
par(mfrow=c(1,5))
#species
p.vad.dist.df.species<- rbind("A-B" = c(Basidiomycota_Ascomycota_SE, 10000-Basidiomycota_Ascomycota_SE))
p.vad.dist.df.species<- rbind(p.vad.dist.df.species, "A-L" = c(Lower_Ascomycota_SE, 10000-Lower_Ascomycota_SE))
p.vad.dist.df.species<- rbind(p.vad.dist.df.species,"L-B" = c(Lower_Basidiomycota_SE, 10000- Lower_Basidiomycota_SE))
rownames(p.vad.dist.df.species)<- c("A-B", "A-L", "L-B")
mosaicplot(p.vad.dist.df.species, ylab = "relative frequency", main = "Species", col=c("black", "dark grey"), border = "white", margin = NULL, cex.axis = 0.66, off= 0)

#genus
p.vad.dist.df.genus<- rbind("A-B" = c(Basidiomycota_Ascomycota_genus, 10000-Basidiomycota_Ascomycota_genus))
p.vad.dist.df.genus<- rbind(p.vad.dist.df.genus, "A-L" = c(Lower_Ascomycota_genus, 10000-Lower_Ascomycota_genus))
p.vad.dist.df.genus<- rbind(p.vad.dist.df.genus,"L-B" = c(Lower_Basidiomycota_genus, 10000- Lower_Basidiomycota_genus))
rownames(p.vad.dist.df.genus)<- c("A-B", "A-L", "L-B")
mosaicplot(p.vad.dist.df.genus, ylab = "relative frequency", main = "Genus", col=c("black", "dark grey"), border = "white", margin = NULL, cex.axis = 0.66, off= 0)

#family
p.vad.dist.df.family<- rbind("A-B" = c(Basidiomycota_Ascomycota_family, 10000-Basidiomycota_Ascomycota_family))
p.vad.dist.df.family<- rbind(p.vad.dist.df.family, "A-L" = c(Lower_Ascomycota_family, 10000-Lower_Ascomycota_family))
p.vad.dist.df.family<- rbind(p.vad.dist.df.family,"L-B" = c(Lower_Basidiomycota_family, 10000- Lower_Basidiomycota_family))
rownames(p.vad.dist.df.family)<- c("A-B", "A-L", "L-B")
mosaicplot(p.vad.dist.df.family, ylab = "relative frequency", main = "Family", col=c("black", "dark grey"), border = "white", margin = NULL, cex.axis = 0.66, off= 0)

#order
p.vad.dist.df.order<- rbind("A-B" = c(Basidiomycota_Ascomycota_order, 10000-Basidiomycota_Ascomycota_order))
p.vad.dist.df.order<- rbind(p.vad.dist.df.order, "A-L" = c(Lower_Ascomycota_order, 10000-Lower_Ascomycota_order))
p.vad.dist.df.order<- rbind(p.vad.dist.df.order,"L-B" = c(Lower_Basidiomycota_order, 10000- Lower_Basidiomycota_order))
rownames(p.vad.dist.df.order)<- c("A-B", "A-L", "L-B")
mosaicplot(p.vad.dist.df.order, ylab = "relative frequency", main = "Order", col=c("black", "dark grey"), border = "white", margin = NULL, cex.axis = 0.66, off= 0)

#class
p.vad.dist.df.class<- rbind("A-B" = c(Basidiomycota_Ascomycota_class, 10000-Basidiomycota_Ascomycota_class))
p.vad.dist.df.class<- rbind(p.vad.dist.df.class, "A-L" = c(Lower_Ascomycota_class, 10000-Lower_Ascomycota_class))
p.vad.dist.df.class<- rbind(p.vad.dist.df.class,"L-B" = c(Lower_Basidiomycota_class, 10000- Lower_Basidiomycota_class))
rownames(p.vad.dist.df.class)<- c("A-B", "A-L", "L-B")
mosaicplot(p.vad.dist.df.class, main = "Class", col=c("black", "dark grey"), border = "white", margin = NULL, cex.axis = 0.66, off= 0)


#get n for each category
#species

randomRows = function(df,n){
  return(df[sample(nrow(df),n),])
}
split.df.by.SE<- split(rDNA_by_taxa, rDNA_by_taxa$SE)
single_species_representation.df <-lapply(names(split.df.by.SE), function(x) randomRows(df = split.df.by.SE[[x]], n = 1))
single_species_representation.df.df <- data.frame(matrix(unlist(single_species_representation.df), nrow=length(single_species_representation.df), byrow=T))
n.speies<- nrow(single_species_representation.df.df)
n.speies
n.asco<-sum(single_species_representation.df.df$X13 == "Ascomycota")
n.asco
n.basid<- sum(single_species_representation.df.df$X13 == "Basidiomycota")
n.basid  
n.lower<- sum(single_species_representation.df.df$X13 == "Lower")
n.lower 


 
split.df.by.genus<- split(rDNA_by_taxa, rDNA_by_taxa$Genus)
single_genus_representation.df <-lapply(names(split.df.by.genus), function(x) randomRows(df = split.df.by.genus[[x]], n = 1))
single_genus_representation.df.df <- data.frame(matrix(unlist(single_genus_representation.df), nrow=length(single_genus_representation.df), byrow=T))
n.genus<- nrow(single_genus_representation.df.df)
n.genus
n.asco<-sum(single_genus_representation.df.df$X13 == "Ascomycota")
n.asco
n.basid<- sum(single_genus_representation.df.df$X13 == "Basidiomycota")
n.basid  
n.lower<- sum(single_genus_representation.df.df$X13 == "Lower")
n.lower 


split.df.by.family<- split(rDNA_by_taxa, rDNA_by_taxa$Family)
single_family_representation.df <-lapply(names(split.df.by.family), function(x) randomRows(df = split.df.by.family[[x]], n = 1))
single_family_representation.df.df <- data.frame(matrix(unlist(single_family_representation.df), nrow=length(single_family_representation.df), byrow=T))
n.family<- nrow(single_family_representation.df.df)
n.family
n.asco<-sum(single_family_representation.df.df$X13 == "Ascomycota")
n.asco
n.basid<- sum(single_family_representation.df.df$X13 == "Basidiomycota")
n.basid  
n.lower<- sum(single_family_representation.df.df$X13 == "Lower")
n.lower 

split.df.by.order<- split(rDNA_by_taxa, rDNA_by_taxa$Order)
single_order_representation.df <-lapply(names(split.df.by.order), function(x) randomRows(df = split.df.by.order[[x]], n = 1))
single_order_representation.df.df <- data.frame(matrix(unlist(single_order_representation.df), nrow=length(single_order_representation.df), byrow=T))
n.order<- nrow(single_order_representation.df.df)
n.order
n.asco<-sum(single_order_representation.df.df$X13 == "Ascomycota")
n.asco
n.basid<- sum(single_order_representation.df.df$X13 == "Basidiomycota")
n.basid  
n.lower<- sum(single_order_representation.df.df$X13 == "Lower")
n.lower 

split.df.by.class<- split(rDNA_by_taxa, rDNA_by_taxa$Class)
single_class_representation.df <-lapply(names(split.df.by.class), function(x) randomRows(df = split.df.by.class[[x]], n = 1))
single_class_representation.df.df <- data.frame(matrix(unlist(single_class_representation.df), nrow=length(single_class_representation.df), byrow=T))
n.class<- nrow(single_class_representation.df.df)
n.class
n.asco<-sum(single_class_representation.df.df$X13 == "Ascomycota")
n.asco
n.basid<- sum(single_class_representation.df.df$X13 == "Basidiomycota")
n.basid  
n.lower<- sum(single_class_representation.df.df$X13 == "Lower")
n.lower

