# BGMP 2025 Chimera Project
Full documentation coming soon!

## How to Set Up and Run
There are 2 pipelines outlined in this repository. One is for parsing long-read PacBio Sequencing data to extract gene sequences, barcodes, and perform barcode consensus calling. The other is for the short-read Illumina Sequencing data to merge paired end sequencing, extract barcodes, and cluster barcodes. The short read sequencing pipeline has the option to run each of the steps within a bash script OR run the workflow in a more parallelized manner within Nextflow.

 Below are links to the README files in this repository that describe how to set up and run each pipeline:

1. Long Read Sequencing Pipeline
- [`Extract`](./extract/README.md)
- [`Consensus`](./consensus/README.md)

2. Short Read Sequencing Pipeline
- [`Run in Bash`](./short_read_pipeline/basic_pipeline/README.md)
- [`Run using Nextflow`](./short_read_pipeline/nextflow_pipeline/README.md)

## Acknowledgements

This project was made in conjunction with the Plesa Lab at the University of Oregon. We would like to especially thank Dr. Calin Plesa and Luca Lippert for their mentorship. We also thank the faculty and staff from the Knight Campus Graduate Internship Program - Bioinformatics track including Dr. Leslie Coonrod, Mr. Jason Sydes, Dr. Hope Healy, and Dr. Maxine Wren. In particular, we thank Mr. Jason Sydes for acting as our non-technical leison.

 This project benefitted from computing resources provided by Research Advanced Computing Services at the University of Oregon.

[`Plesa Lab`](https://www.plesalab.org/)

[`More about the team`](https://linktr.ee/mlscha?utm_source=qr_code)
