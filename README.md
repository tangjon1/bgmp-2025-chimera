# BGMP 2025 Chimera Project

This repository comprises the data-processing scripts and pipelines developed in conjunction with the Plesa Lab at the University of Oregon in order to investigate novel ligand binding specificity in mutagenized sensor proteins using a fluorescent model.

## Overview

There are 3 pipelines described in this repository. One is for parsing long-read PacBio sequencing data to extract gene sequences, barcodes, and perform consensus calling. The other is for the short-read Illumina sequencing data to merge paired-end sequencing data, extract barcodes, and cluster barcodes. The short-read sequencing pipeline has the option to run each of the steps within a Bash script OR run the workflow in a more parallelized manner within Nextflow. 

The long- and short-read data processing pipelines were created in the context of barcode-associated variants, resulting from mutagenesis. These pipelines work in tandem to demultiplex scalable experimental procedures on the part of the Plesa Lab, minimizing the amount of sequencing necessary while still providing comprehensive data. By sequencing the variants and their associated barcodes prior to running the experiment (the input of the long-read pipeline), this allows for the sensor variants expressed by the FACS-sorted cells (where fluorescence levels are a function of binding affinity) to be identified by the barcode alone (the input of the short-read pipeline).

After running these pipelines, for additional insight, there is a sample R pipeline that explains the calculations for the median fluorescence of the variants.

## How to Set Up and Run

Below is an outline of each pipeline in this repository; `README.md` files within each describe how to set up and run the corresponding section:

1. Long Read Sequencing Pipeline
    1. [Extract](./extract/) - Extract a sequence and its associated barcode based on flanking motifs
    2. [Consensus](./consensus/) - Call consensus

2. Short Read Sequencing Pipeline
    - OPTION 1: [Run in Bash](./short_read_pipeline/basic_pipeline/)
    - OPTION 2: [Run using Nextflow](./short_read_pipeline/nextflow_pipeline/)

3. Calculating Median Fluorescence Pipeline
    - [Median Fluorescence](./median_fluorescence/)

## Acknowledgements

We would especially like to thank Dr. Calin Plesa and Luca Lippert for their mentorship. We also thank the faculty and staff from the Knight Campus Graduate Internship Program - Bioinformatics track including Dr. Leslie Coonrod, Prof. Jason Sydes, Dr. Hope Healey, and Dr. Maxine Wren. In particular, we thank Prof. Jason Sydes for acting as our non-technical liaison.

This project benefited from computing resources provided by Research Advanced Computing Services at the University of Oregon.

[Plesa Lab](https://www.plesalab.org/)

[More about the team](https://linktr.ee/mlscha?utm_source=qr_code)
