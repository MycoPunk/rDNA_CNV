
#everything old is new again backup 

#This script takes in PacBio .fastq consessus files and chops them into illumina-like 101bp reads
#whcih can then be used in read depth analysis. This creates paired-end (-1, -2) reads with 2X coverage of each fragment. 


.libPaths(new='~/R/x86_64-pc-linux-gnu-library/3.2')
#setwd("~/Desktop")
options(stringsAsFactors = FALSE)
library(doParallel)
library(foreach)

library(stringi)
require(compiler)
library(compiler)
library(data.table) 
library(R.utils) 

enableJIT(3)
#enableJIT(0) #to trun off Jit

# Create cluster with desired number of cores
cl <- makeCluster(20)

clusterEvalQ(cl, .libPaths(new='~/R/x86_64-pc-linux-gnu-library/3.2'))

# Register cluster
registerDoParallel(cl)

#set path and read in the PacBio .fastq
P_bio_seqs <- data.table::fread("corrected_pbio-a.fastq", header = FALSE, quote="", sep="}")

#P_bio_seqs <- data.table::fread("mock_pbio_reads.fastq", header = FALSE, quote="", sep="}")

#transpose 
Seq <- P_bio_seqs[seq(2, nrow(P_bio_seqs), 4), ]
q_score <- P_bio_seqs[seq(4, nrow(P_bio_seqs), 4), ]
PBIO_df<- data.frame(Seq, q_score)

#1: get number of draws for ea. fragment
frag_length<- 500 
PBIO_df<- data.table(PBIO_df, (round((nchar(PBIO_df[,1]) / frag_length) * 3, digits = 0)))
colnames(PBIO_df) <- c("BP", "q-score", "num_draws")

#function to take n 500 bp draws based on length of fragment  
Fragment = function(string, n) {
  replicate(n= n,  {nStart <- sample(1:(nchar(string) -500), 1)
  samp <- substr(string, nStart, nStart + 499)
  })   
  
}


#2: draw a different number of 500 pb samples based on length of ea fragment 
draws<- apply(PBIO_df, 1, function(x) data.frame(lapply(x[1:2], Fragment, n = x[3])))
draws.df<- rbindlist(draws)


#3: create function to make unique id line
Header_ID = function(){
  header = c("@HISEQ15:100:AXXXXXXXX:1:1000:", sample(1000:2000, 1), ":", sample(1000:2000, 1), " ", "1:N:0:XXXXXX")
  paste(header, collapse="")
}



#4: assemble fastq file:
make_fastq <- function(df) { 
  forward.read  <- NULL 
  reverse.read  <- NULL
  fragment <- seq(from=1, to=(nrow(df)))
  q_line <- seq(from=1, to=(nrow(df)))
  first.101<- seq(from=1, to=101)
  for(i in 1:nrow(df)){
    header <- Header_ID()
    fragment[i] <- df[i,1]
    first.101 <- substring(fragment,1,101)
    last.101 <- substring(fragment,400,501)
    q_line[i] <- df[i,2]
    first.q <- substring(q_line,1,101)
    last.q <- substring(q_line,400,501)
    x <- ("+")
    x2 <- paste(x, collapse="")
    complement <- chartr("ATGC","TACG", last.101)
    reverse <- stri_reverse(complement)
    reverse.q <- stri_reverse(last.q)
    reverse_header <- gsub("1:N:0", "2:N:0", header)
    forward.read[i] <- cat(header, first.101[i], x2, first.q[i], sep = "\n")
    reverse.read[i] <- cat(reverse_header, reverse[i], x2, reverse.q[i], sep = "\n")
  }}


make_fastq_comp <- cmpfun(make_fastq)

#Print the new Illumina-style reads 
#foreach(x=(capture.output(apply(draws.df, 1, function(x) make_fastq_comp(df = draws.df)), file = "Illumina_type_Pbio_reads_capture.output.txt", append = TRUE))) %dopar% {x}


chunk <- 100
#chunk <- 100
n <- nrow(draws.df)
r  <- rep(1:ceiling(n/chunk),each=chunk)[1:n]
d <- split(draws.df,r)
#access each sub data.table like this:
#df1 <- as.data.frame(d[[1]])


#split print
for(i in 1:length(d)){
  captureOutput(lapply(d[i], function(x) make_fastq_comp(df = d[[i]])), 
                file = "PRINT_WITH_SPLIT_from_3_a.fastq", append = TRUE)
}


#print without split 
#captureOutput(make_fastq_comp(df = draws.df), 
#              file = "PRINT_WITHOUT_split_from_a.fastq", append = TRUE)

quit()
