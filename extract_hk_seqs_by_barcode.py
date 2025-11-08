#!/usr/bin/env python

# CLI
import argparse

# Extracting barcodes by motif
import regex

# Parsing FASTQ records
import gzip

# For writing output
from csv import writer

# Multithreading
import concurrent.futures
import threading

# Type annotations (and NamedTuple)
from typing import cast, NamedTuple, Any, IO, Generator, TextIO
from io import TextIOWrapper

# Misc
from os import path

# Establish the motif region to extract tha barcode from
BARCODE_MOTIF = r'((GTCGCTGCCGAACAGC)(........................)(AGGAGAAGAGCGCACG)){e<4}'

# Complement table for reverse complementing
COMPLEMENT_TABLE = str.maketrans("ATGC", "TACG")

# NamedTuple whose instances contain a nucleotide sequence and its corresponding quality sequence
class FastqSeqs(NamedTuple):
    nt: str
    qual: str


# NamedTuple whose instances contain the data to be returned after processing a read
class ProcessedSeqs(NamedTuple):
    barcode: str | None
    nt_seq: str | None
    qual_seq: str | None
    is_reverse: bool
    has_startseq: bool
    has_endseq: bool


def reverse_complement(sequence: str) -> str:
    '''Returns the reverse complement of a provided sequence'''
    return sequence.translate(COMPLEMENT_TABLE)[::-1]


def open_read(file_path: str) -> IO[Any] | TextIOWrapper:
    '''Takes an input file path and returns the file handle using the appropriate open() function, based on whether or not the file is compressed (ends with '.gz') or not.'''
    if file_path.endswith(".gz"):
        return gzip.open(file_path, 'rt')
    return open(file_path, 'r')


def get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(prog=path.basename(__file__), description="Extract and count barcode sequences")
    parser.add_argument('-i', '--input-file', help="Input FASTQ file to extract barcodes and associated DcuS sequence variant from", required=True)
    parser.add_argument('-o', '--output-dir', help="Directory for output files", required=True)
    parser.add_argument('-q', '--quality', help="Export sequence data in FASTQ format (preserve quality scores)", action='store_true')
    parser.add_argument('-t', '--threads', help="The number of threads to utilize", type=int, default=1)
    return parser.parse_args()


def get_fastq_seqs(fastq_file_handle: IO[Any] | TextIO) -> Generator[FastqSeqs, None, None]:
    '''Iteratively retrieve the nucleotide and quality score sequences from a file'''
    # Store the nucleotide and quality sequences in a list
    fastq_seqs_list: list[str] = list()
    for i, line in enumerate(fastq_file_handle):
        # If this is the start of a new record, yield the record seqs
        if i > 0 and i % 4 == 0:
            yield FastqSeqs(fastq_seqs_list[0], fastq_seqs_list[1])
            fastq_seqs_list.clear()
        # Otherwise, we are only interested in the sequence and quality line
        elif i % 4 == 1 or i % 4 == 3:
            fastq_seqs_list.append(line.strip('\n'))
    # Yield the stored sequences from the final record
    yield FastqSeqs(fastq_seqs_list[0], fastq_seqs_list[1])


class BarcodeSeqExtractor:
    def __init__(self, input_file_path: str, output_dir_path: str, fastq_export: bool):
        '''Class to process and track barcode/variant sequences'''
        # Dictionary to store barcodes as keys and the nucleotide/quality sequences as values
        self.sequences_dict: dict[str, list[FastqSeqs]] = dict()
        
        # Dictionary to store various frequency statistics from processing
        self.stats_dict: dict[str, int] = {
            "barcode_match": 0,
            "is_reverse": 0,
            "no_startseq": 0,
            "no_endseq": 0
        }

        # Store file paths
        self.input_file_path = input_file_path
        self.output_dir_path = output_dir_path

        # Store the flag to determine output format (FASTQ or FASTA)
        self.fastq_export = fastq_export

        # Threading locks to prevent race conditions
        self.seq_dict_lock: threading.Lock = threading.Lock()
        self.stats_dict_lock: threading.Lock = threading.Lock()
        self.output_write_lock: threading.Lock = threading.Lock()


    def process_seqs(self, seqs: FastqSeqs) -> ProcessedSeqs:
        '''
        Extracts sequences from a record, tracking them in the class
        '''
        # Initialize the return sequences to None. If these are still None
        # when returning, the calling thread should treat this record as 'invalid';
        # e.g. Did not have all motifs or contained Ns
        barcode: str | None = None
        return_nt_seq: str | None = None
        return_qual_seq: str | None = None

        # Used to pass stats to the worker
        is_reverse = False
        has_startseq = False
        has_endseq = False

        # Search the nt sequence for the motif
        nt_seq: str = seqs.nt.upper()
        motif_match: regex.Match[str] | None = regex.search(BARCODE_MOTIF, nt_seq, regex.BESTMATCH)

        # Check the reverse strand if there was no match on the forward strand
        if motif_match is None:
            nt_seq = reverse_complement(nt_seq)
            motif_match = regex.search(BARCODE_MOTIF, nt_seq, regex.BESTMATCH)
            if motif_match is not None:
                is_reverse = True

        # Search for the variant region
        if motif_match is not None:
            barcode = motif_match.group(3)
            barcode_index: int = motif_match.start() + 16
            sequence = nt_seq[:barcode_index]

            startseq_index = sequence.find("CATATG")
            if startseq_index != -1:
                has_startseq = True
                endseq_index = sequence.find("GACGACCGCACGCTGCTG")
                if endseq_index != -1:
                    has_endseq = True
                    seq_between_ndei_kpni = sequence[(startseq_index + 6):endseq_index]
                    if seq_between_ndei_kpni != "" and 'N' not in seq_between_ndei_kpni:
                        return_nt_seq = seq_between_ndei_kpni
                        return_qual_seq = seqs.qual[:barcode_index][(startseq_index + 6):endseq_index]

        return ProcessedSeqs(barcode, return_nt_seq, return_qual_seq, is_reverse, has_startseq, has_endseq)
        

    def worker(self, seqs: FastqSeqs, output_file: TextIO) -> None:
        '''TODO'''
        # Process the passed "raw" sequence
        processed_seqs: ProcessedSeqs = self.process_seqs(seqs)
        
        # Track the stats in the stats dictionary
        with self.stats_dict_lock:
            self.stats_dict['barcode_match'] += processed_seqs.nt_seq is not None
            self.stats_dict['is_reverse'] += processed_seqs.is_reverse
            self.stats_dict['no_startseq'] += not processed_seqs.has_startseq
            self.stats_dict['no_endseq'] += not processed_seqs.has_endseq
        
        # If the sequences could be successfully processed, update the dictionary and write the record to the output file
        if processed_seqs.nt_seq is not None:
            final_seqs: FastqSeqs = FastqSeqs(processed_seqs.nt_seq, cast(str, processed_seqs.qual_seq))

            # Update that barcode in the dictionary of sequences
            with self.seq_dict_lock:
                self.sequences_dict.setdefault(cast(str, processed_seqs.barcode), list())
                self.sequences_dict[cast(str, processed_seqs.barcode)].append(final_seqs)

            # Render the sequences into the appropriate format and export
            write_record: str = ""

            # Use FASTQ format if flag is set, otherwise FASTA
            if self.fastq_export:
                write_record = f"@{processed_seqs.barcode}\n{processed_seqs.nt_seq}\n+\n{processed_seqs.qual_seq}\n"
            else:
                write_record = f">{processed_seqs.barcode}\n{processed_seqs.nt_seq}\n"

            with self.output_write_lock:
                output_file.write(write_record)


    def extract_seqs(self, num_threads: int, max_pending: int=3500):
        '''Main process'''
        # Establish output files
        input_file_basename = path.basename(self.input_file_path)
        output_filepath_sequences = path.join(self.output_dir_path, f"{input_file_basename}.fast{'q' if self.fastq_export else 'a'}")

        # Run 
        with (
            concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as executor,
            open_read(self.input_file_path) as input_file,
            open(output_filepath_sequences, 'w') as output_file_sequences
        ):
            # Track the futures
            futures = set()

            # Get the sequences of each read
            for i, fastq_seqs in enumerate(get_fastq_seqs(input_file)):
                # Keep no more than max_pending futures in memory
                if len(futures) >= max_pending:
                    # Wait for one future to complete before adding a new one
                    _, futures = concurrent.futures.wait(futures, return_when=concurrent.futures.FIRST_COMPLETED)

                # Process the record
                futures.add(executor.submit(self.worker, fastq_seqs, output_file_sequences))
                
                # Inform progress
                #if (i + 1) % 10 == 0:
                #    print(f"Submitted {i + 1} reads")

            # Wait for remaining tasks to finish
            for _ in concurrent.futures.as_completed(futures):
                pass

            self.total_reads: int = i
        

    def summarize(self) -> None:
        '''Print different statistics regarding the processed records'''
        print(f"total_reads\t{self.total_reads}")
        print(f"unique_barcode_count\t{len(self.sequences_dict)}")
        for stat in self.stats_dict:
            print(f"{stat}\t{self.stats_dict[stat]}")
        print(f"")


if __name__ == '__main__':
    args = get_args()
    bcsx = BarcodeSeqExtractor(args.input_file, args.output_dir, args.quality)
    bcsx.extract_seqs(args.threads)
