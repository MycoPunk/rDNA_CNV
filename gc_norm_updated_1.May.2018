

#set directory
setwd("~/Desktop/Project_ITS_CNV/genes and depth files/Rots_and_saps/")

#set libraries
library(Biostrings)

##read in fasta files
#single-copy
ELF1_fasta<- readDNAStringSet("C_quercuum_ELF1_IR.fasta")[[1]]
G6PDH_fasta<- readDNAStringSet("C_quercuum_G6PDH_IR.fasta")[[1]] 
GAPDH_fasta<- readDNAStringSet("C_quercuum_GAPDH_IR.fasta")[[1]] 
GH63_fasta<- readDNAStringSet("C_quercuum_GH63_IR.fasta")[[1]] 
LYS2_fasta<- readDNAStringSet("C_quercuum_LYS2_IR.fasta")[[1]] 
MCM7_fasta<- readDNAStringSet("C_quercuum_MCM7_IR.fasta")[[1]] 
MLS_fasta<- readDNAStringSet("C_quercuum_MLS_IR.fasta")[[1]] 
RPB1_fasta<- readDNAStringSet("C_quercuum_RPB1_IR.fasta")[[1]] 
RPB2_fasta<- readDNAStringSet("C_quercuum_RPB2_IR.fasta")[[1]] 
TOP2_fasta<- readDNAStringSet("C_quercuum_TOP2_IR.fasta")[[1]] 
#multi-copy
ITS_fasta<- readDNAStringSet("C_quercuum_ITS_IR.fasta")[[1]] 
LSU_fasta<- readDNAStringSet("C_quercuum_LSU_IR.fasta")[[1]] 

##read in depth files
#single-copy
ELF1_table<- read.table("C_quercuum_ELF1_depth.txt")
G6PDH_table<- read.table("C_quercuum_G6PDH_depth.txt")
GAPDH_table<- read.table("C_quercuum_GAPDH_depth.txt")
GH63_table<- read.table("C_quercuum_GH63_depth.txt")
LYS2_table<- read.table("C_quercuum_LYS2_depth.txt")
MCM7_table<- read.table("C_quercuum_MCM7_depth.txt")
MLS_table<- read.table("C_quercuum_MLS_depth.txt")
RPB1_table<- read.table("C_quercuum_RPB1_depth.txt")
RPB2_table<- read.table("C_quercuum_RPB2_depth.txt")
TOP2_table<- read.table("C_quercuum_TOP2_depth.txt")
#multi-copy
ITS_table<- read.table("C_quercuum_ITS_depth.txt")
LSU_table<- read.table("C_quercuum_LSU_depth.txt")


##attach bp ident. to table
#single-copy
ELF1_table<-cbind(ELF1_table, ELF1_fasta)
G6PDH_table<-cbind(G6PDH_table, G6PDH_fasta)
GAPDH_table<-cbind(GAPDH_table, GAPDH_fasta)
GH63_table<-cbind(GH63_table, GH63_fasta)
LYS2_table<-cbind(LYS2_table, LYS2_fasta)
MCM7_table<-cbind(MCM7_table, MCM7_fasta)
MLS_table<-cbind(MLS_table, MLS_fasta)
RPB1_table<-cbind(RPB1_table, RPB1_fasta)
RPB2_table<-cbind(RPB2_table, RPB2_fasta)
TOP2_table<-cbind(TOP2_table, TOP2_fasta)
#multi-copy
ITS_table<-cbind(ITS_table, ITS_fasta)
LSU_table<-cbind(LSU_table, LSU_fasta)


##function to calculate GC within 100 bp of the target bp
sliding.bin.ave <- function(DNA_seq){
  # this function calculates the GC % of each target bin, with bin size 100, and target bp centered
  n      <- 50
  total  <- length(DNA_seq)
  bins   <- seq(from=1, to=(total))
  result <- rep(NA,times = length(DNA_seq))
  thisGC <- as.integer(DNA_seq %in% c("G", "C"))
  for(i in 51:length(bins)){
    result[i] <- sum(thisGC[(bins[i]-n):(bins[i]+n)])
  }
  return(result)
}


##run sliding.bin.ave on all genes and attach results to each table
#single-copy
ELF1_table<- cbind(ELF1_table, sliding.bin.ave(DNA_seq = ELF1_table[,4]))
G6PDH_table<- cbind(G6PDH_table, sliding.bin.ave(DNA_seq = G6PDH_table[,4]))
GAPDH_table<- cbind(GAPDH_table, sliding.bin.ave(DNA_seq = GAPDH_table[,4]))
GH63_table<- cbind(GH63_table, sliding.bin.ave(DNA_seq = GH63_table[,4]))
LYS2_table<- cbind(LYS2_table, sliding.bin.ave(DNA_seq = LYS2_table[,4]))
MCM7_table<- cbind(MCM7_table, sliding.bin.ave(DNA_seq = MCM7_table[,4]))
MLS_table<- cbind(MLS_table, sliding.bin.ave(DNA_seq = MLS_table[,4]))
RPB1_table<- cbind(RPB1_table, sliding.bin.ave(DNA_seq = RPB1_table[,4]))
RPB2_table<- cbind(RPB2_table, sliding.bin.ave(DNA_seq = RPB2_table[,4]))
TOP2_table<- cbind(TOP2_table, sliding.bin.ave(DNA_seq = TOP2_table[,4]))
#multi-copy
ITS_table<- cbind(ITS_table, sliding.bin.ave(DNA_seq = ITS_table[,4]))
LSU_table<- cbind(LSU_table, sliding.bin.ave(DNA_seq = LSU_table[,4]))

#clean it up #1
#delete row if there's an ambiguity code 
BP_list<- c("A", "C", "T", "G")
ELF1_table.1<- ELF1_table[ELF1_table$ELF1_fasta %in% BP_list, ]
G6PDH_table.1<- G6PDH_table[G6PDH_table$G6PDH_fasta %in% BP_list, ]
GAPDH_table.1<- GAPDH_table[GAPDH_table$GAPDH_fasta %in% BP_list, ]
GH63_table.1<- GH63_table[GH63_table$GH63_fasta %in% BP_list, ]
LYS2_table.1<- LYS2_table[LYS2_table$LYS2_fasta %in% BP_list, ]
MCM7_table.1<- MCM7_table[MCM7_table$MCM7_fasta %in% BP_list, ]
MLS_table.1<- MLS_table[MLS_table$MLS_fasta %in% BP_list, ]
RPB1_table.1<- RPB1_table[RPB1_table$RPB1_fasta %in% BP_list, ]
RPB2_table.1<- RPB2_table[RPB2_table$RPB2_fasta %in% BP_list, ]
TOP2_table.1<- TOP2_table[TOP2_table$TOP2_fasta %in% BP_list, ]
#multi copy
ITS_table.1<- ITS_table[ITS_table$ITS_fasta %in% BP_list, ]
LSU_table.1<- LSU_table[LSU_table$LSU_fasta %in% BP_list, ]

#clean it up #2
#delete row if there are depth values of zero, or NA's
ELF1_table.2<- na.omit(ELF1_table.1[!(ELF1_table.1$V3==0),])
G6PDH_table.2<- na.omit(G6PDH_table.1[!(G6PDH_table.1$V3==0),])
GAPDH_table.2<- na.omit(GAPDH_table.1[!(GAPDH_table.1$V3==0),])
GH63_table.2<- na.omit(GH63_table.1[!(GH63_table.1$V3==0),])
LYS2_table.2<- na.omit(LYS2_table.1[!(LYS2_table.1$V3==0),])
MCM7_table.2<- na.omit(MCM7_table.1[!(MCM7_table.1$V3==0),])
MLS_table.2<- na.omit(MLS_table.1[!(MLS_table.1$V3==0),])
RPB1_table.2<- na.omit(RPB1_table.1[!(RPB1_table.1$V3==0),])
RPB2_table.2<- na.omit(RPB2_table.1[!(RPB2_table.1$V3==0),])
TOP2_table.2<- na.omit(TOP2_table.1[!(TOP2_table.1$V3==0),])
#multi copy
ITS_table.2<- na.omit(ITS_table.1[!(ITS_table.1$V3==0),])
LSU_table.2<- na.omit(LSU_table.1[!(LSU_table.1$V3==0),])

##change the names so that they match and can be merged 
NEW.names <- c("sp.", "bp.pos", "long.depth","bp", "gc.bin")
#single-copy
colnames(ELF1_table.2) <- NEW.names
colnames(G6PDH_table.2) <- NEW.names
colnames(GAPDH_table.2) <- NEW.names
colnames(GH63_table.2) <- NEW.names
colnames(LYS2_table.2) <- NEW.names
colnames(MCM7_table.2) <- NEW.names
colnames(MLS_table.2) <- NEW.names
colnames(RPB1_table.2) <- NEW.names
colnames(RPB2_table.2) <- NEW.names
colnames(TOP2_table.2) <- NEW.names
##multi-copy
colnames(ITS_table.2) <- NEW.names
colnames(LSU_table.2) <- NEW.names

##chop first and last 50 pb to use for bin ave. calculation 
#single-copy
ELF1_table_chop<- ELF1_table.2[50:(nrow(ELF1_table.2) -50),]
G6PDH_table_chop<- G6PDH_table.2[50:(nrow(G6PDH_table.2) -50),]
GAPDH_table_chop<- GAPDH_table.2[50:(nrow(GAPDH_table.2) -50),]
GH63_table_chop<- GH63_table.2[50:(nrow(GH63_table.2) -50),]
LYS2_table_chop<- LYS2_table.2[50:(nrow(LYS2_table.2) -50),]
MCM7_table_chop<- MCM7_table.2[50:(nrow(MCM7_table.2) -50),]
MLS_table_chop<- MLS_table.2[50:(nrow(MLS_table.2) -50),]
RPB1_table_chop<- RPB1_table.2[50:(nrow(RPB1_table.2) -50),]
RPB2_table_chop<- RPB2_table.2[50:(nrow(RPB2_table.2) -50),]
TOP2_table_chop<- TOP2_table.2[50:(nrow(TOP2_table.2) -50),]
##multi-copy
ITS_table_chop<- ITS_table.2[50:(nrow(ITS_table.2) -50),]
LSU_table_chop<- LSU_table.2[50:(nrow(LSU_table.2) -50),]

#STOP HERE: check to see how good the data looks
#SCG's
barplot(ELF1_table_chop$long.depth, main = "ELF1, pre-normalization")
barplot(G6PDH_table_chop$long.depth, main = "G6PDH, pre-normalization")
barplot(GAPDH_table_chop$long.depth, main = "GAPDH, pre-normalization")
barplot(GH63_table_chop$long.depth, main = "GH63, pre-normalization")
barplot(LYS2_table_chop$long.depth, main = "LYS2, pre-normalization")
barplot(MCM7_table_chop$long.depth, main = "MCM7, pre-normalization")
barplot(MLS_table_chop$long.depth, main = "MLS, pre-normalization")
barplot(RPB1_table_chop$long.depth, main = "RPB1, pre-normalization")
barplot(RPB2_table_chop$long.depth, main = "RPB2, pre-normalization")
barplot(TOP2_table_chop$long.depth, main = "TOP2, pre-normalization")
#MCG's
barplot(LSU_table_chop$long.depth, main = "LSU, pre-normalization")
barplot(ITS_table_chop$long.depth, main = "ITS, pre-normalization")


##are any of the "SC" genes likely to be "MC"? If they're more than 1 sd outside the median, don't include downstream.
#get means
long.means.ELF1<- mean(ELF1_table_chop$long.depth)
long.means.G6PDH<- mean(G6PDH_table_chop$long.depth)
long.means.GAPDH<- mean(GAPDH_table_chop$long.depth)
long.means.GH63<- mean(GH63_table_chop$long.depth)
long.means.LYS2<- mean(LYS2_table_chop$long.depth)
long.means.MCM7<- mean(MCM7_table_chop$long.depth)
long.means.MLS<- mean(MLS_table_chop$long.depth)
long.means.RPB1<- mean(RPB1_table_chop$long.depth)
long.means.RPB2<- mean(RPB2_table_chop$long.depth)
long.means.TOP2<- mean(TOP2_table_chop$long.depth)

#make a list of means
long.means<- c(long.means.ELF1,
               long.means.G6PDH,
               long.means.GAPDH,
               long.means.GH63,
               long.means.LYS2,
               long.means.MCM7,
               long.means.MLS,
               long.means.RPB1,
               long.means.RPB2,
               long.means.TOP2)

#name the genes
names(long.means)<- c( "ELF1", "G6PDH", "GAPDH", "GH63", "LYS2", "MCM7","MLS", "RPB1","RPB2","TOP2")

#get cutoffs 
medianSCG<- median(long.means)
sdSCG<- sd(long.means)
upper<- medianSCG + sdSCG
lower<- medianSCG - sdSCG

#don't include genes in the next step that return 'FALSE', (above abline in graph). 
only.genes.in.range<- long.means < upper & long.means > lower
only.genes.in.range

barplot(long.means, las = 2, ylim =c(0,400))
abline(a= upper, b= 0)
abline(a= lower, b= 0)


#combine all sc genes, #hash out genes above or below ab-line in graph
all_SC_gene_tables<- rbind(ELF1_table_chop, 
                           G6PDH_table_chop,
                           GAPDH_table_chop,
                           GH63_table_chop,
                           LYS2_table_chop,
                           MCM7_table_chop,
                           MLS_table_chop,
                           RPB1_table_chop,
                           RPB2_table_chop,
                           TOP2_table_chop)


#RD median for all bins
all_SC_bins<- median(all_SC_gene_tables[,3])

#get depth medians for each unique bin 
depth.ea.bin.all<- aggregate(data = all_SC_gene_tables, all_SC_gene_tables[,3]~all_SC_gene_tables[,5], FUN = median) 

#Function to adjust depth by GC bins 
cor.depth <- function(depth.df){
  result  <- seq(from=1, to=(nrow(depth.df))) 
  x <- all_SC_bins #median RD over all single copy genes
  for(i in 1:nrow(depth.df)){ 
    for (j in 1:nrow(depth.ea.bin.all)){
      if (depth.ea.bin.all[j,1] == depth.df[i,5]) {
        result[i] <- depth.df[i,3]*(x /(depth.ea.bin.all[j,2]))
      }
    }
  }
  return(result)
}


#run function on each single copy region over the relevent bp's (100 chopped from either side), and get mean
#hash out genes above or below ab-line in graph
adj.ELF1<- mean(cor.depth(depth.df = ELF1_table_chop))
adj.G6PDH<- mean(cor.depth(depth.df = G6PDH_table_chop))  
adj.GAPDH<- mean(cor.depth(depth.df = GAPDH_table_chop))  
adj.GH63<- mean(cor.depth(depth.df = GH63_table_chop)) 
adj.LYS2<- mean(cor.depth(depth.df = LYS2_table_chop))  
adj.MCM7<- mean(cor.depth(depth.df = MCM7_table_chop))
adj.MLS<- mean(cor.depth(depth.df = MLS_table_chop))
adj.RPB1<- mean(cor.depth(depth.df = RPB1_table_chop))
adj.RPB2<- mean(cor.depth(depth.df = RPB2_table_chop))
adj.TOP2<- mean(cor.depth(depth.df = TOP2_table_chop))  



#make a list of adjusted means 
#hash out genes above or below ab-line in graph
adjusted.means<- c(adj.ELF1,
                   adj.G6PDH,
                   adj.GAPDH,
                   adj.GH63,
                   adj.LYS2,
                   adj.MCM7,
                   adj.MLS,
                   adj.RPB1,
                   adj.RPB2,
                   adj.TOP2)

names<-          c("ELF1=",
                   "G6PDH=",
                   "GAPDH=",
                   "GH63=",
                   "LYS2=",
                   "MCM7=",
                   "MLS=",
                   "RPB1=",
                   "RPB2=",
                   "TOP2=")

##to adjust multi-copy region of interest (ITS). 
cor.depth.multi_ITS <- function(depth.df){
  #this function takes a dataframe (depth.df) with read depth at each bp, and the average GC content (%) of the 100 bps surrounding the bp of interest
  result  <- seq(from=1, to=(nrow(depth.df)))
  #ave RD over ITS unadjusted bins
  x <- mean(ITS_table_chop$long.depth) 
  #median read depth for each unique GC bin (percent GC over 100 bp, centered on the bp of interest)
  depth.ea.bin.ITS<- aggregate(data = ITS_table_chop, ITS_table_chop[,3]~ITS_table_chop[,5], FUN = median)
  for(i in 1:nrow(depth.df)){
    for (j in 1:nrow(depth.ea.bin.ITS)){
      #for each unique GC % median
      if (depth.ea.bin.ITS[j,1] == depth.df[i,5]) {
        #for each pb position, multiply the depth by (the average over all positions devided by the median for that unique GC %)
        result[i] <- depth.df[i,3]*(x /(depth.ea.bin.ITS[j,2]))
      }
    }
  }
  return(result)
}

corrected.ITS.dep<- cor.depth.multi_ITS(depth.df = ITS_table_chop)
mean.ITS.corrected<- mean(corrected.ITS.dep)
mean.ITS.corrected


##to adjust multi-copy region of interest (LSU). 
cor.depth.multi_LSU <- function(depth.df){
  #this function takes a dataframe (depth.df) with read depth at each bp, and the average GC content (%) of the 100 bps surrounding the bp of interest
  result  <- seq(from=1, to=(nrow(depth.df)))
  #ave RD over LSU unadjusted bins
  x <- mean(LSU_table_chop$long.depth)
  #median read depth for each unique GC bin (percent GC over 100 bp, centered on the bp of interest)
  depth.ea.bin.LSU<- aggregate(data = LSU_table_chop, LSU_table_chop[,3]~LSU_table_chop[,5], FUN = median) 
  for(i in 1:nrow(depth.df)){
    for (j in 1:nrow(depth.ea.bin.LSU)){
      #for each unique GC % median
      if (depth.ea.bin.LSU[j,1] == depth.df[i,5]) {
        #for each pb position, multiply the depth by (the average over all positions devided by the median for that unique GC %)
        result[i] <- depth.df[i,3]*(x /(depth.ea.bin.LSU[j,2]))
      }
    }
  }
  return(result)
}

corrected.LSU.dep<- cor.depth.multi_LSU(depth.df = LSU_table_chop)
mean.LSU.corrected<- mean(corrected.LSU.dep)
mean.LSU.corrected

#check that you trust the data
barplot(ITS_table_chop$long.depth, main = "ITS before correction")
barplot(corrected.ITS.dep, main = "ITS after correction")

barplot(LSU_table_chop$long.depth, main = "LSU before correction")
barplot(corrected.LSU.dep, main = "LSU after correction")

#Average depth
SCG_ave<- mean(adjusted.means)
MC_ave<- (mean.LSU.corrected + mean.ITS.corrected) /2

#Estimated Copy Number of rDNA cassette 
CN_est_ITS<- round(mean.ITS.corrected / SCG_ave)
CN_est<- round(MC_ave / SCG_ave)

#percent difference between ITS and LSU depth 
per_diff<- round(abs(mean.LSU.corrected - mean.ITS.corrected) / ((mean.LSU.corrected + mean.ITS.corrected) / 2)*100, digits = 3)

#generate output
#individual totals
df<- data.frame(cbind(names, adjusted.means))
write.table(df, file = "C_quercuum_depth_totals.txt",
            append = TRUE, sep = "\t", row.names=FALSE, col.names=FALSE, quote = FALSE)
names.2 <- "ITS="
ITS.df<- data.frame(cbind(names.2, mean.ITS.corrected))
write.table(ITS.df, file = "C_quercuum_depth_totals.txt",
            append = TRUE, sep = "\t", row.names=FALSE, col.names=FALSE, quote = FALSE)
names.3 <- "LSU="
LSU.df<- data.frame(cbind(names.3, mean.LSU.corrected))
write.table(LSU.df, file = "C_quercuum_depth_totals.txt",
            append = TRUE, sep = "\t", row.names=FALSE, col.names=FALSE, quote = FALSE)
#summary totals 
names.4 <- "SCG depth average="
SCG.df<- data.frame(cbind(names.4, round(SCG_ave)))
write.table(SCG.df, file = "C_quercuum_depth_totals.txt",
            append = TRUE, sep = "\t", row.names=FALSE, col.names=FALSE, quote = FALSE)
names.5 <- "estimated rDNA CN (ITS only)="
CN.df.ITS<- data.frame(cbind(names.5, CN_est_ITS))
write.table(CN.df.ITS, file = "C_quercuum_depth_totals.txt",
            append = TRUE, sep = "\t", row.names=FALSE, col.names=FALSE, quote = FALSE)
names.6 <- "estimated rDNA CN (ITS / LSU)="
CN.df<- data.frame(cbind(names.6, CN_est))
write.table(CN.df, file = "C_quercuum_depth_totals.txt",
            append = TRUE, sep = "\t", row.names=FALSE, col.names=FALSE, quote = FALSE)
names.7 <- "% diff ITS / LSU="
CN.df<- data.frame(cbind(names.7, per_diff))
write.table(CN.df, file = "C_quercuum_depth_totals.txt",
            append = TRUE, sep = "\t", row.names=FALSE, col.names=FALSE, quote = FALSE)

