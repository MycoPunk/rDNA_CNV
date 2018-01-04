
#set directory
setwd("~/Desktop/Project_ITS_CNV/genes and depth files/A_muscaria_1904")

#set libraries
library(Biostrings)

##read in fasta files
#single-copy
ELF1_fasta<- readDNAStringSet("A_muscaria_ELF1_IR.fasta")[[1]]
G6PDH_fasta<- readDNAStringSet("A_muscaria_G6PDH_IR.fasta")[[1]] 
GAPDH_fasta<- readDNAStringSet("A_muscaria_GAPDH_IR.fasta")[[1]] 
GH63_fasta<- readDNAStringSet("A_muscaria_GH63_IR.fasta")[[1]] 
LYS2_fasta<- readDNAStringSet("A_muscaria_LYS2_IR.fasta")[[1]] 
MCM7_fasta<- readDNAStringSet("A_muscaria_MCM7_IR.fasta")[[1]] 
MLS_fasta<- readDNAStringSet("A_muscaria_MLS_IR.fasta")[[1]] 
RPB1_fasta<- readDNAStringSet("A_muscaria_RPB1_IR.fasta")[[1]] 
RPB2_fasta<- readDNAStringSet("A_muscaria_RPB2_IR.fasta")[[1]] 
TOP2_fasta<- readDNAStringSet("A_muscaria_TOP2_IR.fasta")[[1]] 
#multi-copy
ITS_fasta<- readDNAStringSet("A_muscaria_ITS_IR.fasta")[[1]] 
LSU_fasta<- readDNAStringSet("A_muscaria_LSU_IR.fasta")[[1]] 

##read in depth files
#single-copy
ELF1_table<- read.table("A_muscaria_1904_ELF1_depth_raw.txt")
G6PDH_table<- read.table("A_muscaria_1904_G6PDH_depth_raw.txt")
GAPDH_table<- read.table("A_muscaria_1904_GAPDH_depth_raw.txt")
GH63_table<- read.table("A_muscaria_1904_GH63_depth_raw.txt")
LYS2_table<- read.table("A_muscaria_1904_LYS2_depth_raw.txt")
MCM7_table<- read.table("A_muscaria_1904_MCM7_depth_raw.txt")
MLS_table<- read.table("A_muscaria_1904_MLS_depth_raw.txt")
RPB1_table<- read.table("A_muscaria_1904_RPB1_depth_raw.txt")
RPB2_table<- read.table("A_muscaria_1904_RPB2_depth_raw.txt")
TOP2_table<- read.table("A_muscaria_1904_TOP2_depth_raw.txt")
#multi-copy
ITS_table<- read.table("A_muscaria_1904_ITS_depth_raw.txt")
LSU_table<- read.table("A_muscaria_1904_LSU_depth_raw.txt")

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


##run sliding.bin.ave on all single copy genes 
#single-copy
ELF1.bins<- sliding.bin.ave(DNA_seq = ELF1_table[,4])
G6PDH.bins<- sliding.bin.ave(DNA_seq = G6PDH_table[,4])
GAPDH.bins<- sliding.bin.ave(DNA_seq = GAPDH_table[,4])
GH63.bins<- sliding.bin.ave(DNA_seq = GH63_table[,4])
LYS2.bins<- sliding.bin.ave(DNA_seq = LYS2_table[,4])
MCM7.bins<- sliding.bin.ave(DNA_seq = MCM7_table[,4])
MLS.bins<- sliding.bin.ave(DNA_seq = MLS_table[,4])
RPB1.bins<- sliding.bin.ave(DNA_seq = RPB1_table[,4])
RPB2.bins<- sliding.bin.ave(DNA_seq = RPB2_table[,4])
TOP2.bins<- sliding.bin.ave(DNA_seq = TOP2_table[,4])
#multi-copy
ITS.bins<- sliding.bin.ave(DNA_seq = ITS_table[,4])
LSU.bins<- sliding.bin.ave(DNA_seq = LSU_table[,4])

##add bin averages to depth tables and remove NA values 
#single-copy
ELF1_table<- na.omit(cbind(ELF1_table, ELF1.bins))
G6PDH_table<- na.omit(cbind(G6PDH_table, G6PDH.bins))
GAPDH_table<- na.omit(cbind(GAPDH_table, GAPDH.bins))
GH63_table<- na.omit(cbind(GH63_table, GH63.bins))
LYS2_table<- na.omit(cbind(LYS2_table, LYS2.bins))
MCM7_table<- na.omit(cbind(MCM7_table, MCM7.bins))
MLS_table<- na.omit(cbind(MLS_table, MLS.bins))
RPB1_table<- na.omit(cbind(RPB1_table, RPB1.bins))
RPB2_table<- na.omit(cbind(RPB2_table, RPB2.bins))
TOP2_table<- na.omit(cbind(TOP2_table, TOP2.bins))
#multi-copy
ITS_table<- na.omit(cbind(ITS_table, ITS.bins))
LSU_table<- na.omit(cbind(LSU_table, LSU.bins))

##changes the names so that they match and can be merged 
new.names <- c("sp.", "bp.pos", "long.depth","bp", "gc.bin")
#single-copy
colnames(ELF1_table) <- new.names
colnames(G6PDH_table) <- new.names
colnames(GAPDH_table) <- new.names
colnames(GH63_table) <- new.names
colnames(LYS2_table) <- new.names
colnames(MCM7_table) <- new.names
colnames(MLS_table) <- new.names
colnames(RPB1_table) <- new.names
colnames(RPB2_table) <- new.names
colnames(TOP2_table) <- new.names
##multi-copy
colnames(ITS_table) <- new.names
colnames(LSU_table) <- new.names

##chop first and last 50 pb to use for bin ave. calculation 
#single-copy
ELF1_table_chop<- ELF1_table[50:(nrow(ELF1_table) -50),]
G6PDH_table_chop<- G6PDH_table[50:(nrow(G6PDH_table) -50),]
GAPDH_table_chop<- GAPDH_table[50:(nrow(GAPDH_table) -50),]
GH63_table_chop<- GH63_table[50:(nrow(GH63_table) -50),]
LYS2_table_chop<- LYS2_table[50:(nrow(LYS2_table) -50),]
MCM7_table_chop<- MCM7_table[50:(nrow(MCM7_table) -50),]
MLS_table_chop<- MLS_table[50:(nrow(MLS_table) -50),]
RPB1_table_chop<- RPB1_table[50:(nrow(RPB1_table) -50),]
RPB2_table_chop<- RPB2_table[50:(nrow(RPB2_table) -50),]
TOP2_table_chop<- TOP2_table[50:(nrow(TOP2_table) -50),]
##multi-copy
ITS_table_chop<- ITS_table[50:(nrow(ITS_table) -50),]
LSU_table_chop<- LSU_table[50:(nrow(LSU_table) -50),]

##are any of the "SC" genes likely to be "MC"? If they're more than 1 sd outside the median, don't include downstream.
#get long means
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

#make a list of long means
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

barplot(long.means, las = 2)
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
                           #RPB1_table_chop,
                           RPB2_table_chop,
                           TOP2_table_chop)


#RD median for all bins
all_SC_bins<- median(all_SC_gene_tables[,3])

#get depth medians for each unique bin 
depth.ea.bin.all<- aggregate(data = all_SC_gene_tables, all_SC_gene_tables[,3]~all_SC_gene_tables[,5], FUN = median) 

#ugly function to adjust depth by GC bins 
cor.depth <- function(depth.df){
  #this function takes ...
  result  <- seq(from=1, to=(nrow(depth.df))) 
  x <- all_SC_bins #ave RD over all single copy genes
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
#adj.RPB1<- mean(cor.depth(depth.df = RPB1_table_chop))
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
                   #adj.RPB1,
                   adj.RPB2,
                   adj.TOP2)

names<-          c("ELF1=",
                   "G6PDH=",
                   "GAPDH=",
                   "GH63=",
                   "LYS2=",
                   "MCM7=",
                   "MLS=",
                   #"RPB1=",
                   "RPB2=",
                   "TOP2=")

##to adjust multi-copy region of interest (ITS). 
cor.depth.multi_ITS <- function(depth.df){
  #this function takes ...
  result  <- seq(from=1, to=(nrow(depth.df)))
  x <- mean(ITS_table_chop$long.depth) #ave RD over ITS unadjusted bins
  depth.ea.bin.ITS<- aggregate(data = ITS_table_chop, ITS_table_chop[,3]~ITS_table_chop[,5], FUN = median)
  for(i in 1:nrow(depth.df)){
    for (j in 1:nrow(depth.ea.bin.ITS)){
      if (depth.ea.bin.ITS[j,1] == depth.df[i,5]) {
        result[i] <- depth.df[i,3]*(x /(depth.ea.bin.ITS[j,2]))
      }
    }
  }
  return(result)
}

corrected.ITS.dep<- cor.depth.multi_ITS(depth.df = ITS_table_chop)
mean.ITS.corrected<- mean(corrected.ITS.dep)
mean.ITS.corrected

# #to adjust multi-copy region of interest (LSU). 
cor.depth.multi_LSU <- function(depth.df){
  #this function takes ...
  result  <- seq(from=1, to=(nrow(depth.df)))
  x <- mean(LSU_table_chop$long.depth) #ave RD over LSU unadjusted bins
  depth.ea.bin.LSU<- aggregate(data = LSU_table_chop, LSU_table_chop[,3]~LSU_table_chop[,5], FUN = median)
  for(i in 1:nrow(depth.df)){
    for (j in 1:nrow(depth.ea.bin.LSU)){
      if (depth.ea.bin.LSU[j,1] == depth.df[i,5]) {
        result[i] <- depth.df[i,3]*(x /(depth.ea.bin.LSU[j,2]))
      }
    }
  }
  return(result)
}

corrected.LSU.dep<- cor.depth.multi_LSU(depth.df = LSU_table_chop)
mean.LSU.corrected<- mean(corrected.LSU.dep)
mean.LSU.corrected




#generate output
df<- data.frame(cbind(names, adjusted.means))
write.table(df, file = "A_muscaria_1904_Corrected_totals_sdabv_sdbel.txt",
            append = TRUE, sep = "\t", row.names=FALSE, col.names=FALSE, quote = FALSE)
names.2 <- "ITS="
ITS.df<- data.frame(cbind(names.2, mean.ITS.corrected))
write.table(ITS.df, file = "A_muscaria_1904_Corrected_totals_sdabv_sdbel.txt",
            append = TRUE, sep = "\t", row.names=FALSE, col.names=FALSE, quote = FALSE)
names.3 <- "LSU="
LSU.df<- data.frame(cbind(names.3, mean.LSU.corrected))
write.table(LSU.df, file = "A_muscaria_1904_Corrected_totals_sdabv_sdbel.txt",
            append = TRUE, sep = "\t", row.names=FALSE, col.names=FALSE, quote = FALSE)

