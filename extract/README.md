# Motif Extraction

Following the construction of a plasmid mutagenesis library, the variant sequence and its associated barcode must be extracted from the sequencing data based on the engineered motifs surrounding the regions.

This functionality is provided by the script [`extract_hk_seqs_by_barcode.py`](./extract_hk_seqs_by_barcode.py), which accepts FASTA/Q sequencing formats:

```bash
extract_hk_seqs_by_barcode.py \
    --input-file reads.fastq \
    --output-dir extraction_results/ \
    --quality \ # Preserve seqs' quality scores
```

This will provide an output in FASTA/Q format where—for each record—the header is the barcode and the sequence is its associated variant sequence from the read. Additionally, a barcodes counts TSV file is also created which contains the number of times each barcode appeared throughout the file.