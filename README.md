# BGMP 2025 Chimera Project


There are 2 pipelines outlined in this repository. One is for parsing long-read PacBio sequencing data to extract gene sequences, barcodes, and perform barcode consensus calling. The other is for the short-read Illumina sequencing data to merge paired-end sequencing data, extract barcodes, and cluster barcodes. The short-read sequencing pipeline has the option to run each of the steps within a Bash script OR run the workflow in a more parallelized manner within Nextflow.

## How to Set Up and Run
There are 3 pipelines outlined in this repository. One is for parsing long-read PacBio Sequencing data to extract gene sequences, barcodes, and perform barcode consensus calling. Another is for the short-read Illumina Sequencing data to merge paired end sequencing, extract barcodes, and cluster barcodes. The short read sequencing pipeline has the option to run each of the steps within a bash script OR run the workflow in a more parallelized manner within Nextflow. After running these pipelines, if you are interested, there is an example pipeline in R that explains the calculations for the median fluorescence of the mutant variants.

Below are links to the README files in this repository that describe how to set up and run each pipeline:

1. Long Read Sequencing Pipeline
    1. [Extract](./extract/README.md) - Extract a sequence and its associated barcode based on flanking motifs
    2. [Consensus](./consensus/README.md)

2. Short Read Sequencing Pipeline
    - [Run in Bash](./short_read_pipeline/basic_pipeline/README.md)
    - [Run using Nextflow](./short_read_pipeline/nextflow_pipeline/README.md)

## Methodology

The long- and short-read data processing pipelines were created in the context of barcode-associated variants, derived from mutagenesis. 

3. Calculating Median Fluorescence Pipeline
- [`Median Fluorescence`](./median_fluorescence/README.md)

## Acknowledgements

This project was developed in conjunction with the Plesa Lab at the University of Oregon. We would especially like to thank Dr. Calin Plesa and Luca Lippert for their mentorship. We also thank the faculty and staff from the Knight Campus Graduate Internship Program - Bioinformatics track including Dr. Leslie Coonrod, Prof. Jason Sydes, Dr. Hope Healey, and Dr. Maxine Wren. In particular, we thank Prof. Jason Sydes for acting as our non-technical liaison.

This project benefited from computing resources provided by Research Advanced Computing Services at the University of Oregon.

[Plesa Lab](https://www.plesalab.org/)

[More about the team](https://linktr.ee/mlscha?utm_source=qr_code)
