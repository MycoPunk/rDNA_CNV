
#This script takes in PacBio .fastq consessus files and chops them into illumina-like 101bp reads
#whcih can then be used in read depth analysis. This creates paired-end (-1, -2) reads with XX coverage of each fragment. 

'''
#NOTE: before running this script, you need to remove any pbio reads < 500 bp. 

module load cutadapt

cutadapt --minimum-length 500 <input>.fastq -o <output>.fastq
'''


.libPaths(new='~/R/x86_64-pc-linux-gnu-library/3.2')

options(stringsAsFactors = FALSE)
library(doParallel)
library(foreach)

library(stringi)
require(compiler)
library(compiler)
library(data.table) 
library(R.utils) 

#turn on compiler
enableJIT(2)

# Create cluster with desired number of cores
cl <- makeCluster(20)

# Register cluster
registerDoParallel(cl)

#set path and read in the PacBio .fastq
P_bio_seqs <- data.table::fread("pbio-1487.13649.fastq", header = FALSE, quote="", sep="}")

#transpose into dataframe
Seq <- P_bio_seqs[seq(2, nrow(P_bio_seqs), 4), ]
q_score <- P_bio_seqs[seq(4, nrow(P_bio_seqs), 4), ]
PBIO_df<- data.table(Seq, q_score)

#get number of draws for ea. fragment
frag_length<- 500 
PBIO_df<- data.table(PBIO_df, (round((nchar(PBIO_df[,1]) / frag_length) * .5, digits = 0)))
colnames(PBIO_df) <- c("BP", "q-score", "num_draws")


#function to make unique id line
Header_ID = function(){
  header = c("@HISEQ15:100:AXXXXXXXX:1:1000:", sample(1000:2000, 1), ":", sample(1000:2000, 1), " ", "1:N:0:XXXXXX")
  paste(header, collapse="")
}

#function to take n = 500 bp draws based on length of PacBio fragment  
Fragment = function(string, n) {
  replicate(n= n,  {nStart <- sample(1:(nchar(string) -500), 1)
  samp <- substr(string, nStart, nStart + 499)
  })   
  
}


#draw a different number 500 bp draws based on length of PacBio fragment  
draws<- apply(PBIO_df, 1, function(x) as.data.table(lapply(x[1:2], Fragment, n = x[3])))
draws.df<- rbindlist(draws)

#function to assemble fastq short reads
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

#cmp for speed (byte code compiler)
make_fastq_comp <- cmpfun(make_fastq)

#print output to file - note captureOutput will add NULL lines, but these will be removed when splitting foward from reverse reads
#and it's (a lot) faster than other methods
for(i in 1:nrow(draws.df)){
  captureOutput(apply(draws.df[i], 1, function(x) make_fastq_comp(df = (draws.df[i]))), file = "Pbio_to_ill_no_redundancies_1.fastq", append = TRUE)
}


quit()

