setwd("~/Desktop/Project_ITS_CNV/Mock_Community")


#read in data
otu_table_rep1 <- read.csv("mock_community_analysis_rep_1.csv")
otu_table_rep2 <- read.csv("mock_community_analysis_rep_2.csv")
otu_table_ITS2<- read.csv("ITS2_mock_community_OTU_TABLE.csv")

#parse to only the species put into the mock + negative controls 
otu_table_relevent_1<- otu_table_rep1[c(1:12,105,1517), c(2:5,7:10)]
otu_table_relevent_2<- otu_table_rep2[c(1:12,18,166), c(4:17)]
otu_table_relevent_ITS2<- otu_table_ITS2[,2:6]

#assign taxonomy as rownames 
rownames1<- as.matrix(otu_table_relevent_1$taxonomy)
rownames(otu_table_relevent_1) <-rownames1
rownames2<- as.matrix(otu_table_relevent_2$taxonomy)
rownames(otu_table_relevent_2) <-rownames2
rownames3<- as.matrix(otu_table_relevent_ITS2$taxonomy)
rownames(otu_table_relevent_ITS2) <-rownames3

#remove taxonomy column 
otu_table_relevent_1<- otu_table_relevent_1[, 1:7]
otu_table_relevent_2<- otu_table_relevent_2[,1:14]
otu_table_relevent_ITS2<- otu_table_relevent_ITS2[,1:4]

View(otu_table_relevent_1)
View(otu_table_relevent_subtract_rep1.2)
#subtrat reads from negative controls 
otu_table_relevent_subtract_rep1<-  otu_table_relevent_1[,1:4] - otu_table_relevent_1[,6]
otu_table_relevent_subtract_rep1.2<-  otu_table_relevent_subtract_rep1[,1:4] - otu_table_relevent_1[,7]
otu_table_relevent_subtract_rep1.2<- as.matrix(otu_table_relevent_subtract_rep1.2)
#remove values < 0
otu_table_relevent_subtract_no_neg_rep1<- ifelse(otu_table_relevent_subtract_rep1.2 < 0.1, 0, otu_table_relevent_subtract_rep1.2)

otu_table_relevent_subtract_rep2<-  otu_table_relevent_2[,1:12] - otu_table_relevent_2[,13]
otu_table_relevent_subtract_rep2.2<- as.matrix(otu_table_relevent_subtract_rep2)
#remove values < 0
otu_table_relevent_subtract_no_neg_rep2<- ifelse(otu_table_relevent_subtract_rep2.2 < 0.1, 0, otu_table_relevent_subtract_rep2.2)

#note, no negative control for two-step ITS2, just for one step. 

#normalize data from 0-1
norm.reads_rep1<- t(t(otu_table_relevent_subtract_no_neg_rep1)/rowSums(t(otu_table_relevent_subtract_no_neg_rep1)))
as.data.frame(norm.reads_rep1_100<- norm.reads_rep1 * 100)
#check that they == 100
colSums(norm.reads_rep1_100)
#...cool 

norm.reads_rep2<- t(t(otu_table_relevent_subtract_no_neg_rep2)/rowSums(t(otu_table_relevent_subtract_no_neg_rep2)))
as.data.frame(norm.reads_rep2_100<- norm.reads_rep2 * 100)
#check that they == 100
colSums(norm.reads_rep2_100)

norm.reads_repITS2<- t(t(otu_table_relevent_ITS2)/rowSums(t(otu_table_relevent_ITS2)))
as.data.frame(norm.reads_repITS2_100<- norm.reads_repITS2 * 100)
#check that they == 100
colSums(norm.reads_repITS2_100)

##RUN STATS
#load metadata
summary_table_lg <- read.csv("summary_table_new.csv")

#shrink data to species of interest
summary_table <- summary_table_lg[c(1,3,5,7,9:18),] 


#reassign taxonomy as row names 
rownames4<- as.matrix(summary_table$Species)
rownames(summary_table) <-rownames4
summary_table_renamed<- summary_table[,2:9]

#sort tables by speces name
norm.reads_rep1_100.sorted <- norm.reads_rep1_100[order(row.names(norm.reads_rep1_100)),] 
norm.reads_rep2_100.sorted <- norm.reads_rep2_100[order(row.names(norm.reads_rep2_100)),] 
norm.reads_repITS2_100.sorted <- norm.reads_repITS2_100[order(row.names(norm.reads_repITS2_100)),] 
summary_table_renamed.sorted <- summary_table_renamed[order(row.names(summary_table_renamed)),] 

#convert to dataframes
norm.reads_rep1_100.sorted<- data.frame(norm.reads_rep1_100.sorted)
norm.reads_rep2_100.sorted<- data.frame(norm.reads_rep2_100.sorted)
norm.reads_repITS2_100.sorted<- data.frame(norm.reads_repITS2_100.sorted)

ACT_all_sp<- data.frame(cbind(norm.reads_rep1_100.sorted[,c(1,3)]), norm.reads_rep2_100.sorted[,c(1:2)])
CTA_all_sp<- data.frame(cbind(norm.reads_rep1_100.sorted[,c(2,4)]), norm.reads_rep2_100.sorted[,c(3:4)])
ATC_suillus_only<- data.frame(norm.reads_rep2_100.sorted[,c(5,6)])
CTA_suillus_only<- data.frame(norm.reads_rep2_100.sorted[,c(7,8)])
ATC_no_suillus<- data.frame(norm.reads_rep2_100.sorted[,c(9,10)])
CTA_no_suillus<- data.frame(norm.reads_rep2_100.sorted[,c(11,12)])
ATC_all_sp_ITS2<- data.frame(norm.reads_repITS2_100.sorted[,c(1,2)])
CTA_all_sp_ITS2<- data.frame(norm.reads_repITS2_100.sorted[,c(3,4)])


#get means for rep and rep 2- all species #ATC=amplify then combine, #CTA = combine then amplify 
mean_ATC<-data.frame(rowMeans(ACT_all_sp))
colnames(mean_ATC) <-"ave"
mean_CTA<-data.frame(rowMeans(CTA_all_sp))
colnames(mean_CTA) <-"ave"
mean_ATC_suillus_only<-data.frame(rowMeans(ATC_suillus_only))
colnames(mean_ATC_suillus_only) <-"ave"
mean_CTA_suillus_only<-data.frame(rowMeans(CTA_suillus_only))
colnames(mean_CTA_suillus_only) <-"ave"
mean_ATC_no_suillus<-data.frame(rowMeans(ATC_no_suillus))
colnames(mean_ATC_no_suillus) <-"ave"
mean_CTA_no_suillus<-data.frame(rowMeans(CTA_no_suillus))
colnames(mean_CTA_no_suillus) <-"ave"
#for ITS2
mean_ATC_all_sp_ITS2<-data.frame(rowMeans(ATC_all_sp_ITS2))
colnames(mean_ATC_all_sp_ITS2) <-"ave"
mean_CTA_all_sp_ITS2<-data.frame(rowMeans(CTA_all_sp_ITS2))
colnames(mean_CTA_all_sp_ITS2) <-"ave"


##build model of read abundance (ITS1 averaged between lanes, no LSU) explained by CNV, GC, amplicon length

#ITS1 all ATC (n = 4)
lm_all_ATC <- lm(mean_ATC$ave ~ summary_table_renamed.sorted$Ave_CN_by.Species_ITS_ONLY + 
                  summary_table_renamed.sorted$GC_ITS1 + 
                  summary_table_renamed.sorted$ITS1_length, 
                  mean_ATC)
summary(lm_all_ATC )
#Mean_ATC by CN
plot(mean_ATC$ave ~ summary_table_renamed.sorted$Ave_CN_by.Species_ITS_ONLY, mean_ATC)

#CN Alone
lm_all_ATC_CN <- lm(mean_ATC$ave ~ summary_table_renamed.sorted$Ave_CN_by.Species_ITS_ONLY, 
   mean_ATC)

summary(lm_all_ATC_CN)

#Mean_ATC by GC
plot(mean_ATC$ave ~ summary_table_renamed.sorted$GC_ITS1, mean_ATC)

#GC alone
lm_all_ATC_GC <- lm(mean_ATC$ave ~ summary_table_renamed.sorted$GC_ITS1, 
                    mean_ATC)

summary(lm_all_ATC_GC)

#Mean_ATC by Length
plot(mean_ATC$ave ~ summary_table_renamed.sorted$ITS1_length, mean_ATC)

#Length alone
lm_all_ATC_Length <- lm(mean_ATC$ave ~ summary_table_renamed.sorted$ITS1_length, 
                    mean_ATC)

summary(lm_all_ATC_Length)


#ITS1 all CTA (n = 4)
lm_all_CTA <- lm(mean_CTA$ave ~ summary_table_renamed.sorted$Ave_CN_by.Species_ITS_ONLY + 
                summary_table_renamed.sorted$GC_ITS1 + 
                summary_table_renamed.sorted$ITS1_length, 
                mean_CTA)
summary(lm_all_CTA)

plot(mean_CTA$ave ~ summary_table_renamed.sorted$GC_ITS1, mean_CTA)




#####
#ITS1 Suillus only ACT (n = 2)
lm_ATC_suillus_only <- lm(mean_ATC_suillus_only$ave ~ summary_table_renamed.sorted$Ave_CN_by.Species_ITS_ONLY + 
                   summary_table_renamed.sorted$GC_ITS1 + 
                   summary_table_renamed.sorted$ITS1_length, 
                   mean_ATC_suillus_only)
summary(lm_ATC_suillus_only)

#ITS1 Suillus only CTA (n = 2)
lm_CTA_suillus_only <- lm(mean_CTA_suillus_only$ave ~ summary_table_renamed.sorted$Ave_CN_by.Species_ITS_ONLY + 
                            summary_table_renamed.sorted$GC_ITS1 + 
                            summary_table_renamed.sorted$ITS1_length, 
                          mean_CTA_suillus_only)
summary(lm_CTA_suillus_only)
#####

View(summary_table_renamed.sorted)
#ITS1 no Suillus ATC (n = 2)
lm_ATC_no_suillus <- lm(mean_ATC_no_suillus$ave ~ summary_table_renamed.sorted$Ave_CN_by.Species_ITS_ONLY + 
                            summary_table_renamed.sorted$GC_ITS1 + 
                            summary_table_renamed.sorted$ITS1_length, 
                        mean_ATC_no_suillus)
summary(lm_ATC_no_suillus)

#ITS1 no Suillus CTA (n = 2)
lm_CTA_no_suillus <- lm(mean_CTA_no_suillus$ave ~ summary_table_renamed.sorted$Ave_CN_by.Species_ITS_ONLY + 
                          summary_table_renamed.sorted$GC_ITS1 + 
                          summary_table_renamed.sorted$ITS1_length, 
                        mean_CTA_no_suillus)
summary(lm_CTA_no_suillus)



#####

#ITS2, all sp, ATC
lm_ATC_all_sp_ITS2 <- lm(mean_ATC_all_sp_ITS2$ave ~ summary_table_renamed.sorted$Ave_CN_by.Species_ITS_ONLY + 
                          summary_table_renamed.sorted$GC_ITS2 + 
                          summary_table_renamed.sorted$ITS2_length, 
                         mean_ATC_all_sp_ITS2)
summary(lm_ATC_all_sp_ITS2)

#ITS2, all sp, CTA
lm_CTA_all_sp_ITS2 <- lm(mean_CTA_all_sp_ITS2$ave ~ summary_table_renamed.sorted$Ave_CN_by.Species_ITS_ONLY + 
                           summary_table_renamed.sorted$GC_ITS2 + 
                           summary_table_renamed.sorted$ITS2_length, 
                         mean_CTA_all_sp_ITS2)
summary(lm_CTA_all_sp_ITS2)



#figure for ITS2
barplot(as.matrix(ATC_all_sp_ITS2[,c(1,2)]), main="Mock Community Analysis #1",
        ylab="read abundance", col=c("#C02E1D", "#D93807", "#D94E1F", "#F16C20", "#EF8B2C", "#ECAA38", "#EBC844", "#CCC155", "#A2B86C", "#5CA793", "#1395BA", "#117899", "#0F5B78", "#0D3C55") 
        , beside=TRUE,cex.names=0.8, ylim=c(0, 15))
legend("topleft", legend = fungi_names, fill=c("#C02E1D", "#D93807", "#D94E1F", "#F16C20", "#EF8B2C", "#ECAA38", "#EBC844", "#CCC155", "#A2B86C", "#5CA793", "#1395BA", "#117899", "#0F5B78", "#0D3C55"),
       ncol = 1, cex = .4, lwd = 0, box.lwd = 0, box.lty =0, box.col =0, xjust = 1)

#figure for ITS2
barplot(as.matrix(CTA_all_sp_ITS2[,c(1,2)]), main="Mock Community Analysis #1",
        ylab="read abundance", col=c("#C02E1D", "#D93807", "#D94E1F", "#F16C20", "#EF8B2C", "#ECAA38", "#EBC844", "#CCC155", "#A2B86C", "#5CA793", "#1395BA", "#117899", "#0F5B78", "#0D3C55") 
        , beside=TRUE,cex.names=0.8, ylim=c(0, 35))
legend("topleft", legend = fungi_names, fill=c("#C02E1D", "#D93807", "#D94E1F", "#F16C20", "#EF8B2C", "#ECAA38", "#EBC844", "#CCC155", "#A2B86C", "#5CA793", "#1395BA", "#117899", "#0F5B78", "#0D3C55"),
       ncol = 1, cex = .4, lwd = 0, box.lwd = 0, box.lty =0, box.col =0, xjust = 1)

