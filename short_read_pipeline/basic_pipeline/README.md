# Short-Read Sequencing Pipeline as a Bash Script (No Nextflow)

This folder contains the short read sequencing pipeline as a Bash Script. Good for running tests or smaller datasets, but it is not parallelized for larger datasets. 

## Requirements

|Package|Version|Availability|
|---|---|---|
|bbtools| 37.62|agbiome|
|htstream| 1.4.1|bioconda|
|starcode|1.4|bioconda|

## Inputs
R1 and R2 fastq files in .fastq.gz format
Fasta files for forward and reverse primer sets

## Pipeline Steps

1. First step is bbmerge. This combines the reads from the R1 and R2 files together to get complete reads

2. Second is hts_Primers. Removes primers. hts primers allows for searching for forward and reverse orientations along with the variable length inserts used for library preparation.

3. Third step trims the reads to grab the barcode, getting rid of any extra bases that might be in the data at this point

4. Fourth, find the counts of each barcode in the data set.

5. Finally, Starcode collapses barcodes based on sequence similarity to help account for mutations or sequencing errors.
