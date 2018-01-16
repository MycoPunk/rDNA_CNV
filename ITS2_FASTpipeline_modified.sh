#!/bin/bash -l 
#PBS -l walltime=0:30:00,mem=16gb,nodes=2:ppn=8
#PBS -m abe 
#PBS -M llofgren@umn.edu

#modified version of FASTpipeline to handle Kennedy Lab ITS2 dual index amplicons 

#Move to specific working directory in Projects. Double check the "read1" and "read2" folders make sure they have the same number of files in them. 
cd /panfs/roc/groups/5/kennedyp/shared/Projects/LotusAmanda/ITS2

#Set short cut to FAST - this allows access to all the scripts in the FAST directory in the commands below.  Make sure this is up to date with the latest FAST release.  
FAST='/panfs/roc/groups/5/kennedyp/shared/FAST_ITS1'

#Load the relevant analysis modules
module load python-epd
module load cutadapt
module load pear
module load vsearch 

#Generate sample mapping files:
python $FAST/fast.py -generate_mapping -i read1 -o read1_map.txt
python $FAST/fast.py -generate_mapping -i read2 -o read2_map.txt

#Label sequences:
python $FAST/fast.py -add_labels -m read1_map.txt -i read1 -o read1_labeled -t 24
python $FAST/fast.py -add_labels -m read2_map.txt -i read2 -o read2_labeled -t 24

#Merge all files:
python $FAST/fast.py -merge_seqs -i read1_labeled -o read1.fastq
python $FAST/fast.py -merge_seqs -i read2_labeled -o read2.fastq

#Trim off OUR sequencing primers:
cutadapt -a TTGGTCATTTAGAGGAAGTAAAAGTCGTAACAAGGTTTCC -A CGTTCTTCATCGATGCVAGARCCAAGAGATC -o read1.cut.fastq -p read2.cut.fastq read1.fastq read2.fastq -m 50

#Merge Read1 and Read2 sequences:
pear -f read1.cut.fastq -r read2.cut.fastq -o merge.pear.its2 -k -j 24

#Remove SSU and 5.8S regions:
python $FAST/fast.py -nucl_freq -i merge.pear.its2.assembled.fastq -o merge.pear.head.txt
python $FAST/fast.py -nucl_freq -i merge.pear.its2.assembled.fastq -o merge.pear.tail.txt -tail

#remove first 55 bp's (cut off 5.8s)
#in the form: cutadapt -u 55 -o trimmed.fastq reads.fastq
cutadapt -u 55 -o merge.pear.its2.cut_fr1.fastq merge.pear.its2.assembled.fastq

#remove the last 60pb's (cut off LSU)
cutadapt -u -60 -o merge.pear.its2.cut_fr.fastq merge.pear.its2.cut_fr1.fastq

#Discard low quality sequences: need to check you are using the right file
vsearch --fastq_filter merge.pear.its2.cut_fr.fastq --fastq_maxee 1 --fastaout merge.pear.its2.maxee1.fasta --fasta_width 0 --thread 16

#Dereplicate sequences (for Mesabi, -t 1 is actually the fastest way).
python $FAST/fast.py -dereplicate -i merge.pear.its2.maxee1.fasta -o raw.qc.derep -t 1

#Discard singletons:
python $FAST/fast.py -filter_otu_map -i raw.qc.derep.txt -o raw.qc.derep.size2.txt -min_size 2
python $FAST/fast.py -pick_seqs -i raw.qc.derep.fasta -map raw.qc.derep.size2.txt -o raw.qc.derep.size2.fasta -sizeout

#Chimera checking with VSEARCH using UNITE as reference
vsearch --uchime_ref raw.qc.derep.size2.fasta --nonchimeras raw.qc.derep.size2.uchime.fasta --db /panfs/roc/groups/5/kennedyp/shared/Tools/database/sh_general_release_dynamic_22.08.2016.fasta --sizeout --fasta_width 0 --thread 16

#Cluster OTU at 97% similarity using the greedy algorithm:
vsearch --cluster_size raw.qc.derep.size2.uchime.fasta --centroids raw.qc.vsearch.fasta --fasta_width 0 -id 0.97 --sizein --uc raw.qc.uc.txt --threads 16

#Parse the UC output into OTU map:
python $FAST/fast.py -parse_uc_cluster -i raw.qc.uc.txt -o raw.qc.vsearch.txt

#Combine the dereplicate map and sequences:
python $FAST/fast.py -generate_fast_map -map raw.qc.derep.size2.txt -seq raw.qc.derep.size2.uchime.fasta -o fast.derep.txt -derep

#Combine the OTU map and sequences:
python $FAST/fast.py -generate_fast_map -map raw.qc.vsearch.txt -seq raw.qc.vsearch.fasta -o fast.otu.txt -otu

#Combine to FAST derep map and OTU map into a single hybrid:
python $FAST/fast.py -combine_fast_map -derep_map fast.derep.txt -otu_map fast.otu.txt -o fast.hybrid.txt

#Rename the OTUs, so them will start with OTU_:
python $FAST/fast.py -rename_otu_map -fast_map fast.hybrid.txt -o fast.hybrid.otu.txt

#Generate the OTU table from the FAST hybrid map, along with the representative sequences:
python $FAST/fast.py -make_otu_table -fast_map fast.hybrid.otu.txt -o otu_table.txt -rep rep_seq.fasta

#Using VSEARCH, search UNITE database (v7.1) for records with similarity >= 60%
vsearch --usearch_global rep_seq.fasta -db /panfs/roc/groups/5/kennedyp/shared/Tools/database/sh_general_release_dynamic_22.08.2016.fasta --userout taxa.vsearch.txt --userfields query+target+ql+pairs+id --id 0.6

#Assign taxonomy using VSEARCH and keep only OTUs likely to be a fungus: filter at match length >= 70% and similarity >= 75% for fungal kingdom (from Tedersoo et al. 2014).  Also divide taxonomy scores over and under 90% (very roughly family level confidence).
python $FAST/fast.py -assign_taxonomy -otu otu_table.txt -tax taxa.vsearch.txt -o otu_table.taxa.txt -min_match 0.75 -min_pident 70 -pident 90 -scores

#Get a summary of the OTU table:
python $FAST/fast.py -summary_otu_table -otu otu_table.taxa.txt -o otu_report.txt

