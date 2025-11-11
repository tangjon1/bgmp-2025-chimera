# Chimera Notebook 2025 - Makayla Schacht

My work for this project can be found here: /projects/bgmp/shared/groups/2025/chimera/mlscha/bgmp-2025-chimera

## October 10/15/2025 - Make combine_reads.sh
Part 2 of the chimera project involves matching barcodes to their ligand binding conditions. To start we said we needed to combine all of the R1 and R2 files of the illumina sequencing data. 

I want to know that when I combine my files, I didn't just overwrite it so let's find out how many lines are in the R1 and R2 files of one experiment run (DcuS_Asp10_S36):

```bash
zcat DcuS_Asp10_S36_L004_R1_001.fastq.gz |wc -l
zcat DcuS_Asp10_S36_L004_R2_001.fastq.gz |wc -l
```
There are 280,334,248 lines in the first read file and the same goes for the 2nd file. So the total should be 560,668,496 lines. 

started by simplimply trying to combine these 2 files into a new file as seen in the code below:
```bash
zcat /projects/bgmp/shared/groups/2025/chimera/shared/calin/GC3F_Novaseq_run_8255/DcuS_Asp10_S36_L004_R1_001.fastq.gz /projects/bgmp/shared/groups/2025/chimera/shared/calin/GC3F_Novaseq_run_8255/DcuS_Asp10_S36_L004_R2_001.fastq.gz >> /projects/bgmp/shared/groups/2025/chimera/shared/NovaSeq_merged/Asp10_S36_combined.fastq.gz
```

check to see if the new file has been created while I wait... yes! the file is there.

How large are these 2 files? How large would I expect the file to be after concatenating them together?

Used the command below to visualize sizes:
```bash
ls -lah
```
The read 1 file is 4G and the read 2 file is 3.3G so we should expect a size of 7.3G after these have been combined, but  the combined file seems to be 51G instead of 7.3G. It seems like the file is no longer zipped. The line count is as expected (560,668,496). So all we need to do after combining these is to re-zip them.

Want to concatenate all corresponding R1 and R2 files in that directory in a loop:
1. Get a list of all unique file parts
```
ls | sed 's|DcuS_||g' | sed -r 's|_L.+||g'| uniq |head -n 12
```
output:
Asp10_S36
Asp5_S32
Asp6_S31
Asp7_S33
Asp8_S34
Asp9_S35
NoLig10_S30
NoLig5_S26
NoLig6_S25
NoLig7_S27
NoLig8_S28
NoLig9_S29

Need to save this as a file.

```bash
ls | sed 's|DcuS_||g' | sed -r 's|_L.+||g'| uniq | head -n 12 > /projects/bgmp/shared/groups/2025/chimera/shared/NovaSeq_merged/conditions.txt
```

 Then for each file that contains one of these and ends in .gz, concatenate them and output to a file with a name containing this 

```bash
while read line; do
    zcat *"$line"* > $line.fastq;
done < conditions.txt
```

use this to make conditions instead so that uf we use it in the future, we can use it for future novaseq data file name formatting (hopefully):
```bash
ls | grep "fastq.gz" | sed -r 's|(.+)_L.+|\1|' 
```
Used the line below lets you store this as a temporary file, which gets deleted sometime later after the run.
```bash
ls $NovaSeq| grep "fastq.gz" | sed -r 's|(.+)_L.+|\1|' | uniq > $tmp_conditions
```
After running the concatenation loop, this is is some of the resource information from that run:
```bash
(base) [mlscha@n0349 GC3F_Novaseq_run_8255]$ sacct -j 39150618 -o jobname,elapsed,MaxRSS,exitcode
   JobName    Elapsed     MaxRSS ExitCode 
---------- ---------- ---------- -------- 
combine_r+   00:43:07                 0:0 
     batch   00:43:07       832K      0:0 
    extern   00:43:07       100K      0:0 
```
You can't get CPU use from sacct though, but I can atleast recover this since I did not use usr/bin/time -v on this loop.

## 10/16/25 - zip all of the combined reads from yesterday
In the same shell script (combine_reads.sh found here: /projects/bgmp/shared/groups/2025/chimera/mlscha/bgmp-2025-chimera), I have a loop to zip all of the files that were unzipped from using zcat (see below).

```bash
#read each line of the conditions file again to find the correspinding fastq file in the new merged directory and zip them
#note that this will always use however many CPU for pigz that is allotted to this bash script already
while read line; do
    /usr/bin/time -v pigz -p $SLURM_CPUS_PER_TASK $Merged/*$line*
done < $tmp_conditions
```
It's taken about an hour to zip half of the files so it will likely take another hour or so. 

So I now understand that when Calin said "merge" the R1 and R2 files, he meant using something like bbmerge to overlapp the corresponding reads from R1 and R2 and keep them, the rest will get dumped in another output file. He did not literally mean to concatenate the R1 and R2 files.

Therefore, I renamed the directory where these are to NovaSeq_combined_dont_use/, but I don't expect we will need this data later. However, I may be able to alter the script I made earlier to help easily run some of the pipeline until we get it implemented with NextFlow. 

### Actual Pipeline
Here is the pipeline we should follow (from Calin in the Slack channel):

Here is an example of our standard bioinformatic pipeline for getting these types of Illumina barcode counts from sequencing data. I've attached two small files here with 10,000 reads each as a test. I've also attached a set of degenerate primers for these small test files which you can use with hts_Primers (note these are not the same as the ones use for your data, I'll post those later)
Software used:
https://jgi.doe.gov/data-and-tools/software-tools/bbtools/
https://s4hts.github.io/HTStream/
awk
https://github.com/gui11aume/starcode
Most of these can be installed with conda (https://anaconda.org/).
#Step 1:
#Merge reads with bbmerge (from bbtools)
```bash
bbmerge.sh in1=G01_test_R1.fastq.gz in2=G01_test_R2.fastq.gz out=G01_mrg.fastq.gz outu1=G01_R1_unmerged.fastq.gz outu2=G01_R2_unmerged.fastq.gz > G01_bbmerge_stats.txt
#in the output you will see a line that says Joined: and a percentage. The percentage should be ideally in high 90%s otherwise something may be wrong.
```
#Step 2:
#flip reads to correct orientation and trim with hts_Primers (from HTStream) (Note: this is the slowest step and is not multithreaded) (Note 2: it can be a pain to get this working the first time)
#note the & at the end of the command makes the command run in the background
./hts_Primers -U G01_mrg.fastq.gz -f G01_mrg_flip_trim -P pSG1stub_FWD.fasta -Q pSG1stub_REV.fasta -l 5 -x -e 6 -d 6 -L G01_flip_stats.txt -F &
#Step 3 (optional)
#make some new directories to store files
mkdir -p cut1
mkdir -p cut1/Final_BC
#Step 4
#decompress the data with gzip
#discard the quality data (this turns the fastq into fasta)
#take the first 20 bases using the cut command (BGMP: note in your case the barcode is 24 bp)
#save these to a new file
gzip -cd G01_mrg_flip_trim_SE.fastq.gz | awk 'NR%4 ==2' | cut -c1-20 | awk 'BEGIN{FIELDWIDTHS="20"} {print $1}' > cut1/G01_mrg_flip_trim_SE.cut.txt &
#Step 5
#count how many times each barcode appears
#set the number (40) after parallel= to be the number of cores you want to use based on your computer
#the & at the end of the command makes it run in the background
cat "cut1/G01_mrg_flip_trim_SE.cut.txt" | sort -T ./ --parallel=40 | uniq -c | awk '{print $2 "\t" $1}' > cut1/G01.count.txt &
#Step 6
#look at the first line of the file
#is the barcode empty?
#note this step has not been automated yet
head cut1/G01.count.txt
#Step 7
#if the answer to Step 6 is yes than run this command:
#this removes the first line from the file
tail -n +2 cut1/G01.count.txt > cut1/G01.fix.count.txt
#if the answer to Step 6 is No, then copy the data as is:
cp cut1/G01.count.txt cut1/G01.fix.count.txt
#Step 8
#starcode collapse on distance 1
starcode -d1 --sphere -i "cut1/G01.fix.count.txt" --output "cut1/Final_BC/G01_collapse_d1.tsv" --print-clusters
If you did everything right the first few lines of the file cut1/Final_BC/G01_collapse_d1.tsv should look like this:
CGGTACCTAAGTGTCGCTGC    21      CGGTACCTAAGTGTCGCTGC,AGGTACCTAAGTGTCGCTGC,CGGAACCTAAGTGTCGCTGC,CGGTACCTAAGTGTCGCTGT,GGGTACCTAAGTGTCGCTGC,TGGTACCTAAGTGTCGCTGC
TCGGTACCTAAGTGTCGCTG    9       TCGGTACCTAAGTGTCGCTG,TGGGTACCTAAGTGTCGCTG,TTGGTACCTAAGTGTCGCTG,ACGGTACCTAAGTGTCGCTG
AATGGTACCTAAGTGTCGCT    7       AATGGTACCTAAGTGTCGCT,ACTGGTACCTAAGTGTCGCT,ATTGGTACCTAAGTGTCGCT
GGGTACCTAAGTGTCGCTGT    6       GGGTACCTAAGTGTCGCTGT,AGGTACCTAAGTGTCGCTGT,TGGTACCTAAGTGTCGCTGT
AGGCGGTACCTAAGTGTCGC    5       AGGCGGTACCTAAGTGTCGC,AGGGGGTACCTAAGTGTCGC,AGCCGGTACCTAAGTGTCGC
ATCGGTACCTAAGTGTCGCT    5       ATCGGTACCTAAGTGTCGCT,TTCGGTACCTAAGTGTCGCT,CTCGGTACCTAAGTGTCGCT
AGCTAAGACAGAGATCATAT    4       AGCTAAGACAGAGATCATAT
That's it! Next step is usually importing this data into R.

So I should start with setting up bbmerge first...

### bbmerge

The below is the code I would use for bbmerge on one set of R1 and R2 reads:
```bash
/usr/bin/time -v bbmerge.sh in1=$NovaSeq/DcuS_Asp10_S36_L004_R1_001.fastq.gz in2=$NovaSeq/DcuS_Asp10_S36_L004_R2_001.fastq.gz out=G01_mrg.fastq.gz outu1=$Merged/DcuS_Asp10_S36_R1_unmerged.fastq.gz outu2=$Merged/DcuS_Asp10_S36_R2_unmerged.fastq.gz ihist=ihist.txt > DcuS_Asp10_S36_bbmerge_stats.txt
```


## 10/20/25 - try bbmerge

create new environment called "NovaSeq" and install bbtools
```bash
conda create -n NovaSeq
conda activate NovaSeq
conda install agbiome::bbtools
```
How much time it took to run:
```bash
 Command being timed: "bbmerge.sh in1=/projects/bgmp/shared/groups/2025/chimera/shared/calin/GC3F_Novaseq_run_8255/DcuS_Asp10_S36_L004_R1_001.fastq.gz in2=/projects/bgmp/shared/groups/2025/chimera/shared/calin/GC3F_Novaseq_run_8255/DcuS_Asp10_S36_L004_R2_001.fastq.gz out=G01_mrg.fastq.gz outu1=/projects/bgmp/shared/groups/2025/chimera/shared/NovaSeq_getBarcodes/DcuS_Asp10_S36_R1_unmerged.fastq.gz outu2=/projects/bgmp/shared/groups/2025/chimera/shared/NovaSeq_getBarcodes/DcuS_Asp10_S36_R2_unmerged.fastq.gz ihist=ihist.txt"
	User time (seconds): 2120.39
	System time (seconds): 79.37
	Percent of CPU this job got: 292%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 12:31.21
	Maximum resident set size (kbytes): 678824
	Exit status: 0
```

## 10/22/25 -- check bbmerge and make shared environment for Illumina pipeline

create a shared environment called "Illumina_env"
```bash
mamba create -p /projects/bgmp/shared/groups/2025/chimera/envs/Illumina_env agbiome::bbtool
mamba activate /gpfs/projects/bgmp/shared/groups/2025/chimera/envs/Illumina_env
/gpfs/projects/bgmp/shared/groups/2025/chimera/envs/Illumina_env/bbtools/lib/bbversion.sh
mamba install bioconda::htstream
bbversion.sh
hts_Primers -v
```
bbtools version = 37.62
hts_Primers Version: 1.4.1

bbmerge did not work. I did not specify the merged output file and I had the wrong syntax for the histogram. I still got the histogram/other stats though. Here they are:

Below are based on the identified insert sizes from bbmerge:
#Mean	74.726
#Median	75
#Mode	75
#STDev	2.179
#PercentOfPairs	95.952

We got the expected pairs of >90%! This is a good start. 

The starting input files have reads that are 159 bases long. These get merged with bbmerge so the lengths couls be longer than 159 bp unless there is a LOT of overlap. Given that our average length of insert is 74.726, this still does not account for the adapters so we still need to use a similar command to test the lengths after bbmerge.

Check for the lengths of the resulting merged file:


Check the lengths after hts_Primers:
```bash
zcat DcuS_Asp10_S36_flip_trim_SE.fastq.gz | grep -A 1 "^@" | grep -v "^@" | grep -v "-" | awk '{print length}' | head
```

## 10/23/25 nextflow for Illumina pipeline

install nextflow into our illumina environment
```bash
mamba activate /gpfs/projects/bgmp/shared/groups/2025/chimera/envs/Illumina_env
conda install bioconda::nextflow
nextflow -v
```
nextflow version 0.30.1.4844

Nexflow is throwing an error at line 2. It could be a version issue. Trying to install the latest version in a new environment to fix it...
```bash
conda create --name nf-env bioconda::nextflow
conda activate nf-env
```
the new environment seems to run nextflow! It does not see bbtools to run bbmerge command so install bbtools

conda install agbiome::bbtools

try re-running

ERROR!

/gpfs/home/mlscha/miniforge3/envs/nf-env/lib/jvm/bin/java: symbol lookup error: /gpfs/home/mlscha/miniforge3/envs/nf-env/lib/jvm/bin/java: undefined symbol: JLI_StringDup
NOTE: Nextflow needs a Java virtual machine to run. To this end:
 - make sure a `java` command can be found; or
 - manually define the variables JAVA_HOME to point to an existing installation; or
 - install a Java virtual machine, for instance through https://sdkman.io (read the docs);
 - for more details please refer to the Nextflow Get Started page at http://docs.nextflow.io.
NOTE: Nextflow is trying to use the Java VM defined by the following environment variables:
 JAVA_CMD: /gpfs/home/mlscha/miniforge3/envs/nf-env/lib/jvm/bin/java
 NXF_OPTS: 

 Jason passed a link to the syntax that will let me call an environment within nextflow so hopefully that will work! dropping the link here:
https://www.nextflow.io/docs/latest/conda.html


## 10/24/25 Make nextflow shared env and bbmerge and hts primers env

Make nextflow shared env
```bash
mamba create -p /projects/bgmp/shared/groups/2025/chimera/envs/Nextflow_env bioconda::nextflow
mamba activate /gpfs/projects/bgmp/shared/groups/2025/chimera/envs/Nextflow_env
```
Make bbmerge and hts primers env
```bash

```

## 10/26/25 try installing more recent bbtools (uses Java 17+ just like Nextflow!)
Try bbmap to be able to put bbmerge, hts primers, AND NEXTFLOW!! IN THE SAME ENVIRONMENT
OLD bbtools version = 37.62

1. update get barcodes bbtools
```bash
#cd into shared env directory for getBarcodes
cd /gpfs/projects/bgmp/shared/groups/2025/chimera/envs/getBarcodes

#remove old bbtools
conda remove bbtools

# Download the latest version of bbtools
wget https://sourceforge.net/projects/bbmap/files/latest/download -O BBTools.tar.gz
tar -xzf BBTools.tar.gz
cd bbmap

# Test installation
./bbversion.sh

#check java version
java -Xmx90m -version
```
NEW bbtools version = 39.38!
openjdk version "1.8.0_462"
OpenJDK Runtime Environment (build 1.8.0_462-b08)
OpenJDK 64-Bit Server VM (build 25.462-b08, mixed mode)

got nextflow up to the counting barcodes command

(Nextflow_env) [mlscha@n0349 test_outfiles]$ cat G01_test_merged_flip_trim_SE.cut.txt | sort -d -T ./ --parallel=40 | head

AAAAACGACGCCTGGCTCTT
AAAAACGACGCCTGGCTCTT
AAAAACGACGCCTGGCTCTT
AAAAACGACGCCTGGCTCTT
AAAAACGACGCCTGGCTCTT
AAAAACGACGCCTGGCTCTT
AAAAACGACGCCTGGCTCTT
AAAAAGGTACCTAAGTGTCG
AAAAAGGTACCTAAGTGTCG


## 10/28/25 finish the nextflow pipeline running on the practice data

install starcode
```bash
conda create -p /projects/bgmp/shared/groups/2025/chimera/envs/starcode_env bioconda::starcode
conda activate /projects/bgmp/shared/groups/2025/chimera/envs/starcode_env
starcode -v
```
Starcode version - starcode-v1.4
starcode github was last updated in 2020 and I think uses java, but it was hard to tell by just looking at the github. There also may be a newer version that runs on a newer version of java, but I believe it requires paying for a license.

Nextflow pipeline is complete, but the output does not line up. My bbmerge output file does not match Ryan's so I think it was because I updated the bbmerge to the new bbmerge version. I'm going to try re-running nextflow using the old bbmerge and see if the output seems closer.

Old bbversion (./bbversion.sh)
37.62

This matches the previous bbtools old version. 

This still does not match Ryan's bbmerge output. I'll continue looking into this next time.

ack, ripgrep, silver searcher, or find and xargs, but it means making a file or a temp file--try this out for "grep parallelizing to remove the empty line" 

## 10/30/25 update the pipeline for 1 set of data

1. Run only bbmerge on 8cpus on srun on DcuS_Asp10_S36 paired files
Completed at: 30-Oct-2025 23:20:07
Duration    : 6m 14s
CPU hours   : 0.1
Succeeded   : 1

2. sbatch run on ENTIRE pipeline with **16cpu** on only one pair (DcuS_Asp10_S36 paired files)
all results are going to "onePair_outfiles" directory
Job #: 39677424 -- getBarcodes_39677424


## 11/01/25 

AAAACCTAAGTGTCGCTGCC - seems to be a barcode on 179 of the illumina reads after the illumina pipeline

Where is this in the original read?


## 11/03/25
previous run from 10/30 did not take a barcode size of 24, rerun un 110325_nf_onepair


## 11/04/25 

check if the run from yesterday keeps the correct size barcodes:
```bash
cat DcuS_Asp10_S36_L004_R_merged_flip_trim_SE_collapse_d1.tsv | cut -f1 | awk '{print length}' | sort -n | uniq -c
      1 1
      2 3
      1 4
     19 5
      3 6
      8 7
     18 8
     19 9
     10 10
     19 11
     34 12
     40 13
     15 14
     46 15
    221 16
     90 17
    114 18
    269 19
    889 20
    712 21
   1481 22
   2183 23
 129206 24
 ```
 it does!

 Run on 2 files and have each set output into it's own results directory

## 11/06/25
can't backup data
didn't fix the "push to repo" issue

Try to run short data files from real data in "ShortIlluminaData" Folder. This is the latest error. 
Command error:
  ERROR: Error: Cannot write to /projects/bgmp/shared/groups/2025/chimera/mlscha/bgmp-2025-chimera/110625_shortdata/hts_primers_out/DcuS_Asp10_S36_L004_merged_flip_stats.txt: No such file or directory

  I thought it could write a new, non-existing directory....


## 11/07/2025 backing up original data, fixing my github repo, running nextflow on 2 small files

making a backup of the original data on the external hard drive (1TB)

```bash
cd /mnt/d/ChimeraBackup2025
rsync -av -e ssh mlscha@login.talapas.uoregon.edu:/projects/bgmp/shared/groups/2025/chimera/shared/calin .
```
this started around 10:40 pm and seems to take a while given the large data files.

fixing the github repo
1. made a copy of the current repo (with a new name)
2. re-cloned the repo
3. checked that I could push to the repo
4. moved the .git file from the working repo to the previous copy
5. move all other repo copies to a folder
6. rename the current working repo to "bgmp-2025-chimera-fixed"

running nextflow on 2 small files AND keeping each output file in a directory named for it's corresponding process
    ex. dcus.merged.fq would be saved in the bbmerge_out folder 
    and the dcus.merged.cut.trim.fq would be in the hts_out folder

