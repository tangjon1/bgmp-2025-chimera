#!/bin/bash

#activate the environment with necessary packages (see READme)
mamba activate chimera_test

#First step is bbmerge. This combines the reads from the R1 and R2 files together to get complete reads
#Takes R1 and R2 files from Illumina sequencing for inputs. 
#Specify merged file with out=. unmerged files can be specified for reads that are unmerged. Command returns stats file.
bbmerge.sh in1=R1_file.fastq.gz\
    in2=R2_file.fastq.gz\
    out=merged_read_file.fastq.gz outu1=R1_unmerged.fastq.gz outu2=R2_unmerged.fastq.gz > your_bbmerge_stats.txt

#This removes the primers. hts primers allows for searching for forward and reverse orientations along with the variable length inserts used for library preparation.
#Input merged read file from bbmerge along with forward and reverse primer fasta files.
#Output merged/flipped/trimmed file from hts-primer results. can also specify stats file
hts_Primers -U merged_read_file.fastq.gz\
    -f merged_flipped_trimmed_file.fastq.gz\
    -P forward_primers.fasta\
    -Q reverse_primers.fasta\
    -l 5 -x -e 6 -d 6 -L your_flip_stats.txt -F 

#next step trims the read to grab the barcode, getting rid of any extra bases that might be in the data at this point
#input merged/flipped/trimmed file from the hts_primers results.
#specify output as a txt file to hold all the cut barcode sequences found in the data.
gzip -cd merged_flipped_trimmed_file.fastq.gz\
    | awk 'NR%4 ==2' | cut -c1-24 | awk 'BEGIN{FIELDWIDTHS="24"} {print $1}' > merged_flipped_trim.cut.txt 

#now find the counts of each trimmed barcode in the data
#input all cut barcode sequences from the last step.
#output as a text file which now contains each barcode with its count separated by a tab
cat merged_flipped_trim.cut.txt\
    | sort -T ./ --parallel=40 | uniq -c | awk '{print $2 "\t" $1}' > your_barcodes.count.txt

#Finally, STARCODE collapses barcodes based on calculated similarity to account for sequencing errors in barcodes
#input the counted barcodes from previous step
#output will return the collapsed barcode counts in tab separated format
starcode -d1 --sphere -i "your_barcodes.count.txt"\
    --output "your_collapsed_counts.tsv" --print-clusters