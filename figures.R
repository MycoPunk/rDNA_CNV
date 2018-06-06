setwd("~/Desktop/Project_ITS_CNV")
library(data.table)

#hash second line work with full data, unhash to work with data minus the outlier 
rDNA_by_taxa<- read.csv("rDNA_by_taxa.csv", sep = ",", header = TRUE)
rDNA_by_taxa<- rDNA_by_taxa[ rDNA_by_taxa$Genus!="Basidiobolus",] #check this. 


#mean CN
mean(rDNA_by_taxa$CN)
median(rDNA_by_taxa$CN)

#figure 1.d (Phylum)
boxplot(rDNA_by_taxa$CN ~rDNA_by_taxa$group, main = "no_outlier") 

cochran.test(CN~group, rDNA_by_taxa, inlying=FALSE) 
aov<- aov(rDNA_by_taxa$CN~rDNA_by_taxa$group, rDNA_by_taxa)
summary(aov(rDNA_by_taxa$CN~rDNA_by_taxa$group, rDNA_by_taxa))
hsd<- TukeyHSD(aov)

###Fig2
#figure 2.a by trophic mode 
#order goups to put milti last
rDNA_by_taxa$TrophicMode=factor(rDNA_by_taxa$TrophicMode, levels=levels(rDNA_by_taxa$TrophicMode)[c(4,3,2,1)])
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


