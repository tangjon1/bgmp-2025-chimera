# Short-Read Sequencing Pipeline using Nextflow

This folder contains the nextflow script to run the short read sequencing pipeline. 

## Requirements

This script requires the following packages and versions. Additionally, HTSPrimers cannot be loaded into the same environment as Nextflow as HTS Primers runs on a much older version of Java. In the Nextflow script, most of these tools are loaded into their own conda environments, then called throughout the nextflow script as each tool is needed. 

|Package|Version|Availability|
|---|---|---|
|Nextflow|25.10.0.10289|conda|
|BBTools|39.38|conda|
|HTSPrimers|1.4.1|conda|
|Starcode|1.4|conda|

## Pipeline Steps

1. Input
- The input is the directory and the naming convention of the read 1 and read 2 files in the direcory that are to be merged.

2. Merge the read files using BBmerge from BBTools
- The corresponding read 1 and read 2 files are merged in parallel during this first step.
- This outputs 1 merged file and two unmerged files for each read file to be used in the next step.

3. Flip and trim reads using HTS Primers
- The merged files are then 
