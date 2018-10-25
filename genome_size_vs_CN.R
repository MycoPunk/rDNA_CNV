setwd("~/Desktop/J.Mol.Ecol.Revised")

#Read in file
df<- read.csv("genome_size.csv")

#exclue Basidiobolus outlier 
df_noB<- df[ df$Genus!="Basidiobolus",] 
#exclue NA's
df_noB_noOut<- df_noB[ !is.na(df_noB$genome.size..Mbp.),] 


#normallize genome size by the size of rDNA cassette 
cassette_size_in_bp<- df_noB_noOut$CN*3900 #assume 3900bp for average single cassette size
genome_size_in_pb<- df_noB_noOut$genome.size..Mbp.*1000000 # Mbp convert to bp
genome_length_wo_cassette<- genome_size_in_pb - cassette_size_in_bp #subtract casette size from genome size
df_with_rDNA_size<-cbind(df_noB_noOut, genome_length_wo_cassette )

#seperate Basids, Ascos, Lower
asco.df<- df_with_rDNA_size[df_with_rDNA_size$Phylum == "Ascomycota",] 
basid.df<- df_with_rDNA_size[df_with_rDNA_size$Phylum == "Basidiomycota",] 
lower.df<- df_with_rDNA_size[df_with_rDNA_size$Phylum != "Basidiomycota", ]
lower.df<- lower.df[lower.df$Phylum != "Ascomycota",]

#check assumptions for Pearson's 
#linear? Ascos
plot(asco.df$genome_length_wo_cassette, asco.df$CN)
abline(lm(asco.df$CN ~ asco.df$genome_length_wo_cassette), col='red',lwd=2)
#linear? Basids
plot(basid.df$genome_length_wo_cassette, basid.df$CN)
abline(lm(basid.df$CN ~ basid.df$genome_length_wo_cassette), col='red',lwd=2)
#linear? lower
plot(lower.df$genome_length_wo_cassette, lower.df$CN)
abline(lm(lower.df$CN ~ lower.df$genome_length_wo_cassette), col='red',lwd=2)

#normal? Ascos
hist(asco.df$genome_length_wo_cassette, breaks = 30)
shapiro.test(asco.df$genome_length_wo_cassett) #not normal genome length (<.05)
hist(asco.df$CN, breaks = 30)
shapiro.test(asco.df$CN) #normal CN (>.05)

#normal? Basids
hist(basid.df$genome_length_wo_cassette, breaks = 30)
shapiro.test(basid.df$genome_length_wo_cassett) #not normal genome length (<.05)
hist(basid.df$CN, breaks = 30)
shapiro.test(basid.df$CN) #normal CN (>.05)

#normal? Lower
hist(lower.df$genome_length_wo_cassette)
shapiro.test(lower.df$genome_length_wo_cassett) #not normal genome length(<.05)
hist(lower.df$CN, breaks = 30)
shapiro.test(lower.df$CN) #normal CN (>.05)

#normal after transformation? Ascos
hist(log(asco.df$genome_length_wo_cassette), breaks = 30)
shapiro.test(log(asco.df$genome_length_wo_cassette)) #normal genome length (>.05)
qqplot(log(asco.df$genome_length_wo_cassette), log(asco.df$CN))
hist(log(asco.df$CN), breaks = 30)
shapiro.test(log(asco.df$CN)) #normal CN (>.05)

#normal after transformation? Basids
hist(log(basid.df$genome_length_wo_cassette), breaks = 30)
shapiro.test(log(basid.df$genome_length_wo_cassett)) #normal (>.05)
qqplot(log(basid.df$genome_length_wo_cassette), log(basid.df$CN))
hist(log(basid.df$CN), breaks = 30)
shapiro.test(log(basid.df$CN)) #not normal (<.05)- don't transform

#normal after transformation? Lower
hist(log(lower.df$genome_length_wo_cassette))
shapiro.test(log(lower.df$genome_length_wo_cassett)) #not normal (<.05)
qqplot(log(lower.df$genome_length_wo_cassette), log(lower.df$CN))
hist(log(lower.df$CN))
shapiro.test(log(lower.df$CN)) #normal (>.05)

##Conclusion:
#transform all except CN for ascos (normal before transformation, not normal after)
#Question: what to do about the genome length for Lower? transformation with log didn't help

#Pearsons's 
#Ascos, log length, but not CN
cor.test(log(asco.df$genome_length_wo_cassette), asco.df$CN, method = "pearson") #not significant
#Basids log length, and log CN
cor.test(log(basid.df$genome_length_wo_cassette), log(basid.df$CN), method = "pearson") #not significant
#Lower log length (?still not normal), log CN
cor.test(log(lower.df$genome_length_wo_cassette), lower.df$CN, method = "pearson") #not significant

#Spearman's
#Ascos, log length but not CN
cor.test(log(asco.df$genome_length_wo_cassette), asco.df$CN, method = "spearman") #not significant
#Basids log length, and CN
cor.test(log(basid.df$genome_length_wo_cassette), log(basid.df$CN), method = "spearman") #not significant
#Lower log length (?still not normal), log CN
cor.test(log(lower.df$genome_length_wo_cassette), lower.df$CN, method = "spearman") #not significant


#Kendall's
#Ascos, log length but not CN
cor.test(log(asco.df$genome_length_wo_cassette), asco.df$CN, method = "kendall") #not significant
#Basids log length, and CN
cor.test(log(basid.df$genome_length_wo_cassette), log(basid.df$CN), method = "kendall") #not significant
#Lower log length (?still not normal), log CN
cor.test(log(lower.df$genome_length_wo_cassette), lower.df$CN, method = "kendall") #not significant



###try without log transforming anyhing 
#Pearsons's 
#Ascos no logs
cor.test(asco.df$genome_length_wo_cassette, asco.df$CN, method = "pearson") #not significant
#Basids no logs
cor.test(basid.df$genome_length_wo_cassette, basid.df$CN, method = "pearson") #not significant
#Lower no logs
cor.test(lower.df$genome_length_wo_cassette, lower.df$CN, method = "pearson") #not significant

#Spearman's
#Ascos, no logs
cor.test(asco.df$genome_length_wo_cassette, asco.df$CN, method = "spearman") #not significant
#Basids no logs
cor.test(basid.df$genome_length_wo_cassette, basid.df$CN, method = "spearman") #not significant
#Lower no logs
cor.test(lower.df$genome_length_wo_cassette, lower.df$CN, method = "spearman") #not significant


#Kendall's
#Ascos, no logs
cor.test(asco.df$genome_length_wo_cassette, asco.df$CN, method = "kendall") #not significant
#Basids no logs
cor.test(basid.df$genome_length_wo_cassette, basid.df$CN, method = "kendall") #not significant
#Lower no logs
cor.test(lower.df$genome_length_wo_cassette, lower.df$CN, method = "kendall") #not significant

#Ok Cool. Genome size independent of CN and CN are NOT correlated in ant way.
#Lets plot this to show it. 

#get plot range
CN_range<- range(df_with_rDNA_size$CN)
genome_size_range<- range(df_with_rDNA_size$genome_length_wo_cassette)

plot(asco.df$genome_length_wo_cassette, asco.df$CN,col = "#006666", pch = 20, xlab = "Geome size", ylab = "CN", 
     main = "corelation between genome size and CN", xlim = genome_size_range, ylim = CN_range)
abline(lm(asco.df$CN ~ asco.df$genome_length_wo_cassette), col = '#006666')
points(basid.df$genome_length_wo_cassette, basid.df$CN, col = '#CCCB64', pch = 20)
abline(lm(basid.df$CN ~ basid.df$genome_length_wo_cassette), col = '#CCCB64')
points(lower.df$genome_length_wo_cassette, lower.df$CN, col = '#F16639', pch = 20)
abline(lm(lower.df$CN ~ lower.df$genome_length_wo_cassette), col = '#F16639')
abline(lm(df_with_rDNA_size$CN ~ df_with_rDNA_size$genome_length_wo_cassette), col ="#5D2D88")


#How does this look as the log?
plot(log(asco.df$genome_length_wo_cassette), log(asco.df$CN),col = "#006666", pch = 20, xlab = "log Geome size", ylab = "log CN", 
     main = "genome size vs. CN", xlim = log(genome_size_range), ylim = log(CN_range))
abline(lm(log(asco.df$CN) ~ log(asco.df$genome_length_wo_cassette)), col = '#006666')
points(log(basid.df$genome_length_wo_cassette), log(basid.df$CN), col = '#CCCB64', pch = 20)
abline(lm(log(basid.df$CN) ~ log(basid.df$genome_length_wo_cassette)), col = '#CCCB64')
points(log(lower.df$genome_length_wo_cassette), log(lower.df$CN), col = '#F16639', pch = 20)
abline(lm(log(lower.df$CN) ~ log(lower.df$genome_length_wo_cassette)), col = '#F16639')
abline(lm(log(df_with_rDNA_size$CN) ~ log(df_with_rDNA_size$genome_length_wo_cassette)), col ="#5D2D88")

       
