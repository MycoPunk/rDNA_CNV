#This script takes in a mock genome and pulls random draws to simulate illumina-like 101bp reads
#whcih can then be used in read depth analysis. This script creates 10M, interlaced paired-end (-1, -2) reads.

.libPaths(new='~/R/x86_64-pc-linux-gnu-library/3.2')
options(stringsAsFactors = FALSE)
library(doParallel)
library(stringi)
library(foreach)
library(compiler)
library(data.table) 
library(R.utils) 

# Create cluster with desired number of cores
cl <- makeCluster(24)

# Register cluster
registerDoParallel(cl)

#read in the synth genome file
whole_synth_genome <- readChar("synth_genome.txt",nchars=52107286) 
#whole_synth_genome <- readChar("synth_genome.txt",nchars=15585)

#function to make unique id line
Header_ID = function(){
  header = c("@HISEQ15:100:AXXXXXXXX:1:1000:", sample(1000:2000, 1), ":", sample(1000:2000, 1), " ", "1:N:0:XXXXXX")
  paste(header, collapse="")
}

#function to take n 500 bp draws based on input n  
Fragment = function(string, n) {
  replicate(n= n,  {nStart <- sample(1:(nchar(string) -500), 1)
  samp <- substr(string, nStart, nStart + 499)
  }) 
  
}


#make dataframe of n 500 bp draws (this is the n you change to get different # of fastq files - n=1M, yields 10M reads)
ill_df<- data.frame(Fragment(string= whole_synth_genome, n=1000000))


#split the table into subtables to increase output speed
subtable <- 10
rows <- nrow(ill_df)
groups  <- rep(1:ceiling(rows/subtable),each=subtable)[1:rows]
split_table <- split(ill_df,groups)


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
  }
  return(forward.read)
  return(reverse.read)
}



#cmp for speed
make_fastq_comp <- cmpfun(make_fastq_ill)


#capture the output of each subtable 
for(i in 1:length(split_table)){
  captureOutput(apply(split_table[[i]], 1, function(x) make_fastq_comp(df = (split_table[[i]]))), file = "Mock_illumina_reads_1M.fasta", append = TRUE)
}


quit()
