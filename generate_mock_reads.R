#This script takes in the mock genome and takes random draws to simulate illumina-like 101bp reads
#which can then be used in read depth analysis. This script generates 1M interleaved paired-end (-1, -2) reads.
#the .pbs submission scrip (included below) was submitted using job chaining, in 100 itterations to produce a total 100M paired end reads. 


##################################chained submission 
#!/bin/bash
one=$(qsub generate_mock_reads.pbs)
echo $one 
for id in seq 2 100; do 
 two=$(qsub -W depend=afterok:$one generate_mock_reads.pbs)
 one=$two
done

##################################.pbs script
#!/bin/bash -l 
#PBS -l walltime=00:08:00,mem=62gb,nodes=1:ppn=20
#PBS -m abe 
#PBS -M llofgren@umn.edu

cd /home/kennedyp/llofgren/genomes/Synth_genome
module load R
R CMD BATCH generate_mock_reads.R
wait
#################################


#setwd("~/Desktop")
.libPaths(new='~/R/x86_64-pc-linux-gnu-library/3.2')
options(stringsAsFactors = FALSE)

library(doParallel)
require(doParallel)
library(stringi)
library(foreach)
library(compiler)
library(data.table) 
library(R.utils) 


# Create cluster with desired number of cores
cl <- makeCluster(20)
clusterEvalQ(cl, .libPaths(new='~/R/x86_64-pc-linux-gnu-library/3.2'))
# Register cluster
registerDoParallel(cl)


#read in the synth genome file
whole_synth_genome <- readChar("synth_genome.txt",nchars=52107286) 


#brake synth genome into chunks to speed the draws

#what size chunks?
#52107286 / 10 = 5210729
wsg_1<- substr(whole_synth_genome, start = 1, stop = 5210729)
#5210729 + 5210729
wsg_2<- substr(whole_synth_genome, start = 5210730, stop = 10421458)
#10421458 + 5210729
wsg_3<- substr(whole_synth_genome, start = 10421459, stop = 15632187)
#15632187+5210729
wsg_4<- substr(whole_synth_genome, start = 15632188, stop = 20842916)
#20842916 + 5210729
wsg_5<- substr(whole_synth_genome, start = 20842917, stop = 26053645)
#26053645 + 5210729
wsg_6<- substr(whole_synth_genome, start = 26053646, stop = 31264374)
#31264374 + 5210729
wsg_7<- substr(whole_synth_genome, start = 31264375, stop = 36475103)
#36475103 + 5210729
wsg_8<- substr(whole_synth_genome, start = 36475104, stop = 41685832)
#41685832 + 5210729
wsg_9<- substr(whole_synth_genome, start = 41685833, stop = 46896561)
#46896561 + 5210729
wsg_10<- substr(whole_synth_genome, start = 46896562, stop = 52107286)

#clean up garbage 
rm(whole_synth_genome)
gc()

#function to make unique id line
Header_ID = function(){
  header = c("@HISEQ15:100:AXXXXXXXX:1:1000:", sample(1000:2000, 1), ":", sample(1000:2000, 1), " ", "1:N:0:XXXXXX")
  paste(header, collapse="")
}


#function to take n 500 bp draws based on input n  
Fragment = function(string) {
  nStart <- sample(1:(nchar(string) -500), 1)
  samp <- substr(string, nStart, nStart + 499)
} 


#function to assemble fastq files using above functions
make_fastq_ill <- function(df) { 
  forward.read  <- NULL 
  reverse.read  <- NULL
  fragment <- seq(from=1, to=(nrow(df)))
  first.101<- seq(from=1, to=101)
  for(i in 1:nrow(df)){
    header <- Header_ID()
    fragment[i] <- df[i,1]
    first.101 <- substring(fragment,1,101)
    last.101 <- substring(fragment,400,501)
    q_line <- rep("~",101)
    q_line2 <- paste(q_line, collapse="")
    x <- ("+")
    x2 <- paste(x, collapse="")
    complement <- chartr("ATGC","TACG", last.101)
    reverse <- stri_reverse(complement)
    reverse_header <- gsub("1:N:0", "2:N:0", header)
    forward.read[i] <- cat(header, first.101[i], x2, q_line2, sep = "\n")
    reverse.read[i] <- cat(reverse_header, reverse[i], x2, q_line2, sep = "\n")
  }}


#make dataframe of n 500 bp draws (this is the 1:n, where n is what you change to get different # of fastq files)
ill_df_1<- foreach(n = 1:10000, .combine = rbind)  %dopar%  Fragment(string= wsg_1)
ill_df_1<- data.frame(ill_df_1)
ill_df_1_fasta<- data.table(captureOutput(make_fastq_ill(df = (ill_df_1))))
rm(ill_df_1, wsg_1)


#then print that shit #1
write.table(ill_df_1_fasta, "synthreads_1M_w_cleanup.fasta", append = TRUE, na = "NA", quote =FALSE, col.names = FALSE, row.names = FALSE)
rm(ill_df_1_fasta)
gc()

#2
ill_df_2<- foreach(n = 1:10000, .combine = rbind)  %dopar%  Fragment(string= wsg_2)
ill_df_2<- data.frame(ill_df_2)
ill_df_2_fasta<- data.table(captureOutput(make_fastq_ill(df = (ill_df_2))))
rm(ill_df_2, wsg_2)


#then print that shit #2
write.table(ill_df_2_fasta, "synthreads_1M_w_cleanup.fasta", append = TRUE, na = "NA", quote =FALSE, col.names = FALSE, row.names = FALSE)
rm(ill_df_2_fasta)
gc()

#3
ill_df_3<- foreach(n = 1:10000, .combine = rbind)  %dopar%  Fragment(string= wsg_3)
ill_df_3<- data.frame(ill_df_3)
ill_df_3_fasta<- data.table(captureOutput(make_fastq_ill(df = (ill_df_3))))
rm(ill_df_3, wsg_3)


#then print that shit #3
write.table(ill_df_3_fasta, "synthreads_1M_w_cleanup.fasta", append = TRUE, na = "NA", quote =FALSE, col.names = FALSE, row.names = FALSE)
rm(ill_df_3_fasta)
gc()

#4
ill_df_4<- foreach(n = 1:10000, .combine = rbind)  %dopar%  Fragment(string= wsg_4)
ill_df_4<- data.frame(ill_df_4)
ill_df_4_fasta<- data.table(captureOutput(make_fastq_ill(df = (ill_df_4))))
rm(ill_df_4, wsg_4)


#then print that shit #4
write.table(ill_df_4_fasta, "synthreads_1M_w_cleanup.fasta", append = TRUE, na = "NA", quote =FALSE, col.names = FALSE, row.names = FALSE)
rm(ill_df_4_fasta)
gc()

#5
ill_df_5<- foreach(n = 1:10000, .combine = rbind)  %dopar%  Fragment(string= wsg_5)
ill_df_5<- data.frame(ill_df_5)
ill_df_5_fasta<- data.table(captureOutput(make_fastq_ill(df = (ill_df_5))))
rm(ill_df_5, wsg_5)


#then print that shit #5
write.table(ill_df_5_fasta, "synthreads_1M_w_cleanup.fasta", append = TRUE, na = "NA", quote =FALSE, col.names = FALSE, row.names = FALSE)
rm(ill_df_5_fasta)
gc()

#6
ill_df_6<- foreach(n = 1:10000, .combine = rbind)  %dopar%  Fragment(string= wsg_6)
ill_df_6<- data.frame(ill_df_6)
ill_df_6_fasta<- data.table(captureOutput(make_fastq_ill(df = (ill_df_6))))
rm(ill_df_6, wsg_6)


#then print that shit #6
write.table(ill_df_6_fasta, "synthreads_1M_w_cleanup.fasta", append = TRUE, na = "NA", quote =FALSE, col.names = FALSE, row.names = FALSE)
rm(ill_df_6_fasta)
gc()

#7
ill_df_7<- foreach(n = 1:10000, .combine = rbind)  %dopar%  Fragment(string= wsg_7)
ill_df_7<- data.frame(ill_df_7)
ill_df_7_fasta<- data.table(captureOutput(make_fastq_ill(df = (ill_df_7))))
rm(ill_df_7, wsg_7)


#then print that shit #7
write.table(ill_df_7_fasta, "synthreads_1M_w_cleanup.fasta", append = TRUE, na = "NA", quote =FALSE, col.names = FALSE, row.names = FALSE)
rm(ill_df_7_fasta)
gc()

#8
ill_df_8<- foreach(n = 1:10000, .combine = rbind)  %dopar%  Fragment(string= wsg_8)
ill_df_8<- data.frame(ill_df_8)
ill_df_8_fasta<- data.table(captureOutput(make_fastq_ill(df = (ill_df_8))))
rm(ill_df_8, wsg_8)


#then print that shit #8
write.table(ill_df_8_fasta, "synthreads_1M_w_cleanup.fasta", append = TRUE, na = "NA", quote =FALSE, col.names = FALSE, row.names = FALSE)
rm(ill_df_8_fasta)
gc()

#9
ill_df_9<- foreach(n = 1:10000, .combine = rbind)  %dopar%  Fragment(string= wsg_9)
ill_df_9<- data.frame(ill_df_9)
ill_df_9_fasta<- data.table(captureOutput(make_fastq_ill(df = (ill_df_9))))
rm(ill_df_9, wsg_9)

#then print that shit #9
write.table(ill_df_9_fasta, "synthreads_1M_w_cleanup.fasta", append = TRUE, na = "NA", quote =FALSE, col.names = FALSE, row.names = FALSE)
rm(ill_df_9_fasta)
gc()

#10
ill_df_10<- foreach(n = 1:10000, .combine = rbind)  %dopar%  Fragment(string= wsg_10)
ill_df_10<- data.frame(ill_df_10)
ill_df_10_fasta<- data.table(captureOutput(make_fastq_ill(df = (ill_df_10))))
rm(ill_df_10, wsg_10)


#then print that shit #10
write.table(ill_df_10_fasta, "synthreads_1M_w_cleanup.fasta", append = TRUE, na = "NA", quote =FALSE, col.names = FALSE, row.names = FALSE)
rm(ill_df_10_fasta)
gc()

stopCluster(cl)
quit()
