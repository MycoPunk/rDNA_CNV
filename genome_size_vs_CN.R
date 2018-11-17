setwd("")

#Read in file
df<- read.csv("genome_size.csv")

#exclue Basidiobolus outlier 
df_noB<- df[ df$Genus!="Basidiobolus",] 
#exclue NA's
df_noB_noOut<- df_noB[ !is.na(df_noB$genome.size..Mbp.),] 

#normallize genome size by the size of rDNA cassette 
cassette_size_in_bp<- df_noB_noOut$CN*9100 #assume 9.1kb for average single cassette size
genome_size_in_bp<- df_noB_noOut$genome.size..Mbp.*1000000 # Mbp convert to bp
genome_length_w_cassette<- genome_size_in_bp + cassette_size_in_bp #add the casette size to the genome size
df_with_rDNA_size1<-cbind(df_noB_noOut, genome_size_in_bp)
df_with_rDNA_size<-cbind(df_with_rDNA_size1, genome_length_w_cassette)

#seperate Basids, Ascos, Lower
asco.df<- df_with_rDNA_size[df_with_rDNA_size$Phylum == "Ascomycota",] 
basid.df<- df_with_rDNA_size[df_with_rDNA_size$Phylum == "Basidiomycota",] 
lower.df<- df_with_rDNA_size[df_with_rDNA_size$Phylum != "Basidiomycota", ]
lower.df<- lower.df[lower.df$Phylum != "Ascomycota",]

###first look at CN by genoms size w/o incuding length added by CN (MC regions not in JGI est of genome size)
#check assumptions for Pearson's 
#linear? Ascos
plot(asco.df$genome_size_in_bp, asco.df$CN)
abline(lm(asco.df$CN ~ asco.df$genome_size_in_bp), col='red',lwd=2)
#linear? Basids
plot(basid.df$genome_size_in_bp, basid.df$CN)
abline(lm(basid.df$CN ~ basid.df$genome_size_in_bp), col='red',lwd=2)
#linear? lower
plot(lower.df$genome_size_in_bp, lower.df$CN)
abline(lm(lower.df$CN ~ lower.df$genome_size_in_bp), col='red',lwd=2)

#normal? Ascos
hist(asco.df$genome_size_in_bp, breaks = 30)
shapiro.test(asco.df$genome_length_w_cassett) #NOT normal genome length (<.05)
hist(asco.df$CN, breaks = 30)
shapiro.test(asco.df$CN) #normal CN (>.05)

#normal? Basids
hist(basid.df$genome_size_in_bp, breaks = 30)
shapiro.test(basid.df$genome_length_w_cassett) #NOT normal genome length (<.05)
hist(basid.df$CN, breaks = 30)
shapiro.test(basid.df$CN) #NOT normal CN (<.05)

#normal? Lower
hist(lower.df$genome_size_in_bp)
shapiro.test(lower.df$genome_length_w_cassett) #NOT normal genome length(<.05)
hist(lower.df$CN, breaks = 30)
shapiro.test(lower.df$CN) #normal CN (>.05)

#normal after transformation? Ascos
hist(log(asco.df$genome_size_in_bp), breaks = 30)
shapiro.test(log(asco.df$genome_size_in_bp)) #normal genome length (>.05)
qqplot(log(asco.df$genome_size_in_bp), log(asco.df$CN))
hist(log(asco.df$CN), breaks = 30)
shapiro.test(log(asco.df$CN)) #NOT normal CN (<.05)

#normal after transformation? Basids
hist(log(basid.df$genome_size_in_bp), breaks = 30)
shapiro.test(log(basid.df$genome_length_w_cassett)) #normal (>.05)
qqplot(log(basid.df$genome_size_in_bp), log(basid.df$CN))
hist(log(basid.df$CN), breaks = 30)
shapiro.test(log(basid.df$CN)) #normal CN (>.05)

#normal after transformation? Lower
hist(log(lower.df$genome_size_in_bp))
shapiro.test(log(lower.df$genome_length_w_cassett)) #NOT normal (<.05)
qqplot(log(lower.df$genome_size_in_bp), log(lower.df$CN))
hist(log(lower.df$CN))
shapiro.test(log(lower.df$CN)) #normal (>.05)

##Conclusion:
#Can use Pearson's on log-transformed data (for GS in Ascos, SG and CN in Basids and GS in Lower)
#Pearson's
#Ascos, 
cor.test(log(asco.df$genome_size_in_bp), asco.df$CN, method = "pearson") #not significant
#Basids 
cor.test(log(basid.df$genome_size_in_bp), log(basid.df$CN), method = "pearson") #not significant
#Lower 
cor.test(log(lower.df$genome_size_in_bp), lower.df$CN, method = "pearson") #not significant


#Spearman's
#Ascos, 
cor.test(asco.df$genome_size_in_bp, asco.df$CN, method = "spearman") #not significant
#Basids 
cor.test(basid.df$genome_size_in_bp, basid.df$CN, method = "spearman") #not significant
#Lower 
cor.test(lower.df$genome_size_in_bp, lower.df$CN, method = "spearman") #not significant

#Kendall's
#Ascos
cor.test(asco.df$genome_size_in_bp, asco.df$CN, method = "kendall") #not significant
#Basids
cor.test(basid.df$genome_size_in_bp, basid.df$CN, method = "kendall") #not significant
#Lower 
cor.test(lower.df$genome_size_in_bp, lower.df$CN, method = "kendall") #not significant

###Conclusion: genome size (without considering length added from CN) and CN is not significant.

###Now look as CN by genome size while considering the influence of the cassette to genome length

#check assumptions for Pearson's 
#linear? Ascos
plot(asco.df$genome_length_w_cassette, asco.df$CN)
abline(lm(asco.df$CN ~ asco.df$genome_length_w_cassette), col='red',lwd=2)
#linear? Basids
plot(basid.df$genome_length_w_cassette, basid.df$CN)
abline(lm(basid.df$CN ~ basid.df$genome_length_w_cassette), col='red',lwd=2)
#linear? lower
plot(lower.df$genome_length_w_cassette, lower.df$CN)
abline(lm(lower.df$CN ~ lower.df$genome_length_w_cassette), col='red',lwd=2)

#normal? Ascos
hist(asco.df$genome_length_w_cassette, breaks = 30)
shapiro.test(asco.df$genome_length_w_cassett) #NOT normal genome length (<.05)
hist(asco.df$CN, breaks = 30)
shapiro.test(asco.df$CN) #normal CN (>.05)

#normal? Basids
hist(basid.df$genome_length_w_cassette, breaks = 30)
shapiro.test(basid.df$genome_length_w_cassett) #NOT normal genome length (<.05)
hist(basid.df$CN, breaks = 30)
shapiro.test(basid.df$CN) #NOT normal CN (<.05)

#normal? Lower
hist(lower.df$genome_length_w_cassette)
shapiro.test(lower.df$genome_length_w_cassett) #NOT normal genome length(<.05)
hist(lower.df$CN, breaks = 30)
shapiro.test(lower.df$CN) #normal CN (>.05)

#normal after transformation? Ascos
hist(log(asco.df$genome_length_w_cassette), breaks = 30)
shapiro.test(log(asco.df$genome_length_w_cassette)) #normal genome length (>.05)
qqplot(log(asco.df$genome_length_w_cassette), log(asco.df$CN))
hist(log(asco.df$CN), breaks = 30)
shapiro.test(log(asco.df$CN)) #normal CN (>.05)

#normal after transformation? Basids
hist(log(basid.df$genome_length_w_cassette), breaks = 30)
shapiro.test(log(basid.df$genome_length_w_cassett)) #normal (>.05)
qqplot(log(basid.df$genome_length_w_cassette), log(basid.df$CN))
hist(log(basid.df$CN), breaks = 30)
shapiro.test(log(basid.df$CN)) #normal (>.05)

#normal after transformation? Lower
hist(log(lower.df$genome_length_w_cassette))
shapiro.test(log(lower.df$genome_length_w_cassett)) #normal (>.05)
qqplot(log(lower.df$genome_length_w_cassette), log(lower.df$CN))
hist(log(lower.df$CN))
shapiro.test(log(lower.df$CN)) #normal (>.05)

##Conclusion:
#Can use non-parametric methods (Pearson's), on log transformed data (GS for Ascos, GS and CN for basids, and GS for Lower)

#Pearson's
#Ascos
cor.test(log(asco.df$genome_length_w_cassette), asco.df$CN, method = "pearson") #not significant
#Basids 
cor.test(log(basid.df$genome_length_w_cassette), log(basid.df$CN), method = "pearson") #not significant
#Lower 
cor.test(log(lower.df$genome_length_w_cassette), lower.df$CN, method = "pearson") #not significant

#Spearman's
#Ascos, 
cor.test(asco.df$genome_length_w_cassette, asco.df$CN, method = "spearman") #not significant
#Basids 
cor.test(basid.df$genome_length_w_cassette, basid.df$CN, method = "spearman") #not significant
#Lower 
cor.test(lower.df$genome_length_w_cassette, lower.df$CN, method = "spearman") #not significant

#Kendall's
#Ascos
cor.test(asco.df$genome_length_w_cassette, asco.df$CN, method = "kendall") #not significant
#Basids
cor.test(basid.df$genome_length_w_cassette, basid.df$CN, method = "kendall") #not significant
#Lower 
cor.test(lower.df$genome_length_w_cassette, lower.df$CN, method = "kendall") #not significant

#Conclusion:
#Genome size independent of CN and CN are NOT correlated in any way.

#Lets plot it (GSclength w/ cassette)
CN_range<- range(df_with_rDNA_size$CN)
genome_size_range<- range(df_with_rDNA_size$genome_length_w_cassette)

plot(asco.df$genome_length_w_cassette, asco.df$CN,col = "#006666", pch = 19, xlab = "Genome size", ylab = "CN", 
     main = "corelation between genome size and CN w/ cassette", xlim = genome_size_range, ylim = CN_range)
abline(lm(asco.df$CN ~ asco.df$genome_length_w_cassette), col = '#006666')
points(basid.df$genome_length_w_cassette, basid.df$CN, col = '#CCCB64', pch = 19)
abline(lm(basid.df$CN ~ basid.df$genome_length_w_cassette), col = '#CCCB64')
points(lower.df$genome_length_w_cassette, lower.df$CN, col = '#F16639', pch = 19)
abline(lm(lower.df$CN ~ lower.df$genome_length_w_cassette), col = '#F16639')
abline(lm(df_with_rDNA_size$CN ~ df_with_rDNA_size$genome_length_w_cassette), col ="#5D2D88")

#as a log? (GSclength w/ cassette)
plot(log(asco.df$genome_length_w_cassette), log(asco.df$CN),col = "#006666", pch = 19, xlab = "log Genome size", ylab = "log CN", 
     main = "corelation between genome size and CN w/ cassette")
abline(lm(log(asco.df$CN) ~ log(asco.df$genome_length_w_cassette)), col = '#006666')
points(log(basid.df$genome_length_w_cassette), log(basid.df$CN), col = '#CCCB64', pch = 19)
abline(lm(log(basid.df$CN) ~ log(basid.df$genome_length_w_cassette)), col = '#CCCB64')
points(log(lower.df$genome_length_w_cassette), log(lower.df$CN), col = '#F16639', pch = 19)
abline(lm(log(lower.df$CN) ~ log(lower.df$genome_length_w_cassette)), col = '#F16639')
abline(lm(log(df_with_rDNA_size$CN) ~ log(df_with_rDNA_size$genome_length_w_cassette)), col ="#5D2D88")

#Lets plot it (GSclength w/o cassette)
CN_range<- range(df_with_rDNA_size$CN)
genome_size_range<- range(df_with_rDNA_size$genome_size_in_bp)

plot(asco.df$genome_size_in_bp, asco.df$CN,col = "#006666", pch = 19, xlab = "Genome size", ylab = "CN", 
     main = "corelation between genome size and CN w/o cassette", xlim = genome_size_range, ylim = CN_range)
abline(lm(asco.df$CN ~ asco.df$genome_size_in_bp), col = '#006666')
points(basid.df$genome_size_in_bp, basid.df$CN, col = '#CCCB64', pch = 19)
abline(lm(basid.df$CN ~ basid.df$genome_size_in_bp), col = '#CCCB64')
points(lower.df$genome_size_in_bp, lower.df$CN, col = '#F16639', pch = 19)
abline(lm(lower.df$CN ~ lower.df$genome_size_in_bp), col = '#F16639')
abline(lm(df_with_rDNA_size$CN ~ df_with_rDNA_size$genome_size_in_bp), col ="#5D2D88")

#as a log? (GSclength w/o cassette)
plot(log(asco.df$genome_size_in_bp), log(asco.df$CN),col = "#006666", pch = 19, xlab = "log Genome size", ylab = "log CN", 
     main = "corelation between genome size and CN w/o cassette")
abline(lm(log(asco.df$CN) ~ log(asco.df$genome_size_in_bp)), col = '#006666')
points(log(basid.df$genome_size_in_bp), log(basid.df$CN), col = '#CCCB64', pch = 19)
abline(lm(log(basid.df$CN) ~ log(basid.df$genome_size_in_bp)), col = '#CCCB64')
points(log(lower.df$genome_size_in_bp), log(lower.df$CN), col = '#F16639', pch = 19)
abline(lm(log(lower.df$CN) ~ log(lower.df$genome_size_in_bp)), col = '#F16639')
abline(lm(log(df_with_rDNA_size$CN) ~ log(df_with_rDNA_size$genome_size_in_bp)), col ="#5D2D88")
       
#difference when accountign for CN length is almost neglabable, you can see it shift just slightly here 
plot(df_with_rDNA_size$genome_size_in_bp)
plot(df_with_rDNA_size$genome_length_w_cassette)
