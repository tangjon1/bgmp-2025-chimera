# BGMP 2025 Chimera Project

This repository comprises the data-processing scripts and pipelines developed in conjunction with the Plesa Lab at the University of Oregon in order to investigate novel ligand binding specificity in mutagenized sensor proteins using a fluorescent model.

## Background

The Plesa Lab [Website Here](https://www.plesalab.org/research.html) works on understanding and engineering biological systems, and along with projects like this seeks to further understand sequence-function relationships in proteins. Specifically, this project focuses on the bacterial sensor protein family of Sensor Histidine Kinases (SHKs). While many of these proteins exist and allow bacteria to respond to their environment, not many of them have been fully characterized. This project looks to close that gap, allowing us to better understand and characterize existing proteins as well as enable the engineering of new ones. To do this, this project analyzes the response of mutagenized DcuS variants to binding a novel ligand to discover the links between the protein's sequence and function.

## Experimental Design

Here we analyze the response of the Dcus variants to their native ligand, fumarate, and a novel ligand, aspartate. By observing the binding affinity of different variants to both fumarate and aspartate, we hope to characterize crucial parts of the protein sequence which determine function. The Plesa Lab created a plasmid library containing mutant sequences of the original DcuS sequence. Each plasmid also contains a unique 24 bp barcode sequence which is used to identify the variants after tbey are sorted by the assay. The assay uses fluorescent response to sort the variants, allow us to assess binding function. To create this fluorescent response, the bacterium in the assay has been engineered to produce a fluorescent response to the EnVZ intermembrane protein. By fusing the DcuS variants to this EnvZ protein, we are therefore able to analyze the fluorescent response.

Here we perform sequencing analysis in two parts on the plasmid library. The first uses longer read sequencing to identify each variant with its corresponding barcode sequence. The second uses shorter read sequencing to identify the variants after they have been sorted according to their fluorescent response in the assay. By combining these two techniques, we can see how changes in variant sequence lead to changes in ligand binding.

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
