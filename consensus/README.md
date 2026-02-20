# Consensus Calling

To address the possibility of multiple variant sequences being associated with the same barcode following [motif extraction](../extract/), a singular consensus sequence must be chosen among the variants. This is done to minimize the impact of false positives in downstream analysis, affirming that each barcode corresponds to only one variant sequence.

## Requirements

This program requires Python 3.14+, in addition to the following packages:

|Package|Version|Availability|
|---|---|---|
|[`biopython`](https://biopython.org)|1.86|`conda-forge`|[Biopython website]|
|[`mafft`](https://mafft.cbrc.jp/)|7.526|`conda-forge`|[MAFFT website]|
|[`editdistance`](https://pypi.org/project/editdistance/)|0.8.1|PyPI (`pip`)|

## Usage

For efficient consensus calling functionality en masse, the input must be a FASTA/Q file where records are sorted by header. Sequences with the same header will be considered for consensus calling together—a single sequence is selected from the group.

### Sorting

To sort such files, use the `sort` subcommand:

```bash
call_consensus.py sort \
    --input-file data.fastq \
    --output-file data.sorted.fastq \
    --cluster-file clusters.tsv # Optional clustering file
```

#### Clustering

To ease the process of renaming headers, a `--cluster-file` can be provided, which is a TSV with the following format:

 - **Column 1:** Collapsed header *(does not need to exist in the FASTA/Q file)* to be used as the header in the final output file
 - **Column 2:** Comma-delimited series of headers (i.e., the cluster) *that should exist in the FASTA/Q file* to be mapped to the header in column 1

**Example format:**
```
HEADER_1    HEADER_1,HEADER_2,HEADER_3
HEADER_5    HEADER_4,HEADER_5
```

In this example, `HEADER_1`, `HEADER_2`, and `HEADER_3` will all be renamed to `HEADER_1` and sorted accordingly. This will allow these records to be considered for consensus calling together downstream; the same applies for `HEADER_4` and `HEADER_5`, which are both renamed to `HEADER_5`. It is not required that the cluster (column 2) contain the collapsed header (column 1).

> [!NOTE]
> While this script was written in the context of barcode-associated variants (where the headers are themselves nucleotide barcode sequences), information in the header is not parsed beyond a comparative level—as such, the headers themselves can contain arbitrary information so long as the sequences to be consensus called together share identical headers.

Use `call_consensus.py sort -h` for additional information.

### Consensus Calling

With the sorted file, the `consensus` subcommand can be used:

```bash
call_consensus.py consensus \
    --input-file data.sorted.fastq \
    --output-file data.consensus.tsv
```

All records that share the same header in the `--input-file` will be evaluated together as a cluster for consensus calling; the highest-scoring sequence will be chosen to represent the cluster in the `--output-file`.
