#!/usr/bin/env python

# CLI
import argparse

# Extracting barcodes by motif
import regex

# Parsing FASTQ records
import gzip

# Multithreading
import concurrent.futures
import multiprocessing
from multiprocessing.managers import DictProxy

# Type annotations (and NamedTuple)
from typing import NamedTuple, Any, IO, Generator, TextIO
from io import TextIOWrapper

# Misc
from os import path

# Establish the motif region to extract the barcode from
BARCODE_MOTIF: regex.Pattern = regex.compile(r'((GTCGCTGCCGAACAGC)(........................)(AGGAGAAGAGCGCACG)){e<4}', regex.BESTMATCH | regex.IGNORECASE)

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
    has_barcode_motif: bool
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
    '''Defines the command line arguments using argparse'''
    parser = argparse.ArgumentParser(prog=path.basename(__file__), description="Extract and count barcode sequences")
    parser.add_argument('-i', '--input-file', help="Input FASTQ file to extract barcodes and associated DcuS sequence variant from", required=True)
    parser.add_argument('-o', '--output-dir', help="Directory for output files", required=True)
    parser.add_argument('-q', '--quality', help="Export sequence data in FASTQ format (preserve quality scores)", action='store_true')
    parser.add_argument('-t', '--threads', help="The number of threads to utilize", type=int, required=True)
    parser.add_argument('-m', '--max-pending', help="A multiplier applied to the number of threads to determine the total number of records to ever load into memory at once, including those currently processing.", type=int, default=10)
    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    '''Validate the user's configuration of arguments; raises an error if the current configuration leads to an incompatibility'''
    if args.threads < 2:
        raise argparse.ArgumentError(args.threads, "Number of threads must be ≥2")
    if args.max_pending < 1:
        raise argparse.ArgumentError(args.max_pending, "Maximum pending processes multiplier must be ≥1")


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
        # Store file paths
        self.input_file_path = input_file_path
        self.output_dir_path = output_dir_path
        self.input_file_basename = path.basename(self.input_file_path)

        # Store the flag to determine output format (FASTQ or FASTA)
        self.fastq_export = fastq_export
        
        # Concurrency
        self.manager = multiprocessing.Manager()
        
        # Dictionary to store barcodes as keys and the nucleotide/quality sequences as values
        self.barcode_freq_dict: DictProxy[str, int] = self.manager.dict()
        
        # Dictionary to store various frequency statistics from processing
        self.stats_dict: DictProxy[str, int] = self.manager.dict({
            'has_barcode_motif': 0,
            'is_reverse': 0,
            'no_startseq': 0,
            'no_endseq': 0,
            'unambiguous_site': 0
        })


    @staticmethod
    def process_seqs(seqs: FastqSeqs) -> ProcessedSeqs:
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
        nt_seq: str = seqs.nt
        motif_match: regex.Match[str] | None = regex.search(BARCODE_MOTIF, nt_seq)

        # Check the reverse strand if there was no match on the forward strand
        if motif_match is None:
            nt_seq = reverse_complement(nt_seq)
            motif_match = regex.search(BARCODE_MOTIF, nt_seq)
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

        return ProcessedSeqs(barcode, return_nt_seq, return_qual_seq, (motif_match is not None), is_reverse, has_startseq, has_endseq)


    @staticmethod
    def writer_worker(result_queue: multiprocessing.Queue, output_path: str, fastq_export: bool, barcode_freq_dict: DictProxy[str, int], stats_dict: DictProxy[str, int]):
        '''Writer process that consumes results from the queue and writes to file'''
        with open(output_path, 'w') as output_file:
            while True:
                processed_seqs: ProcessedSeqs = result_queue.get()
                if processed_seqs is None:
                    break

                stats_dict['has_barcode_motif'] += processed_seqs.has_barcode_motif
                stats_dict['is_reverse'] += processed_seqs.is_reverse
                stats_dict['no_startseq'] += processed_seqs.has_barcode_motif and not processed_seqs.has_startseq
                stats_dict['no_endseq'] += processed_seqs.has_barcode_motif and not processed_seqs.has_endseq and processed_seqs.has_startseq
                stats_dict['unambiguous_site'] += processed_seqs.nt_seq is not None

                # Update the barcode frequency
                if processed_seqs.barcode is not None:
                    barcode_freq_dict[processed_seqs.barcode] = barcode_freq_dict.get(processed_seqs.barcode, 0) + 1

                if processed_seqs.nt_seq is not None:
                    if fastq_export:
                        output_file.write(f"@{processed_seqs.barcode}\n{processed_seqs.nt_seq}\n+\n{processed_seqs.qual_seq}\n")
                    else:
                        output_file.write(f">{processed_seqs.barcode}\n{processed_seqs.nt_seq}\n")


    def extract_seqs(self, num_threads: int, max_pending: int):
        '''Main process'''
        # Establish output file
        output_filepath_sequences = path.join(self.output_dir_path, f"{self.input_file_basename}.fast{'q' if self.fastq_export else 'a'}")

        result_queue = multiprocessing.Queue()
        writer_process = multiprocessing.Process(target=self.writer_worker, args=(result_queue, output_filepath_sequences, self.fastq_export, self.barcode_freq_dict, self.stats_dict))
        writer_process.start()

        # Run main logic
        with (
            concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor,
            open_read(self.input_file_path) as input_file
        ):
            # Track the futures
            futures: set[concurrent.futures.Future] = set()

            # Get the sequences of each read
            for i, fastq_seqs in enumerate(get_fastq_seqs(input_file)):
                # Keep no more than max_pending futures in memory
                if len(futures) >= max_pending:
                    # Wait for one future to complete before adding a new one
                    done, futures = concurrent.futures.wait(futures, return_when=concurrent.futures.FIRST_COMPLETED)

                    # As completed, immediately handle the result
                    for future in done:
                        result_queue.put(future.result())

                # Process the record
                futures.add(executor.submit(self.process_seqs, fastq_seqs))

            # Wait for any remaining records to be processed and handle them
            for future in concurrent.futures.as_completed(futures):
                result_queue.put(future.result())

            self.total_records: int = i + 1
        
        # Signal writer process to finish
        result_queue.put(None)
        writer_process.join()


    def export_barcode_freq(self) -> None:
        '''Export the barcode occurrence/frequency dictionary'''
        output_filepath_barcode_freq = path.join(self.output_dir_path, f"{self.input_file_basename}.tsv")
        with open(output_filepath_barcode_freq, 'w') as output_file_barcode_freq:
            for barcode in self.barcode_freq_dict:
                output_file_barcode_freq.write(f"{barcode}\t{self.barcode_freq_dict[barcode]}\n")


    def summarize(self) -> str:
        '''Returns a summary with different statistics regarding the processed records'''
        output_string: str = ""
        output_string += f"total_records\t{self.total_records}\n"
        output_string += f"unique_barcode_count\t{len(self.barcode_freq_dict)}\n"
        for stat in self.stats_dict:
            output_string += f"{stat}\t{self.stats_dict[stat]}\n"
        return output_string


if __name__ == '__main__':
    args = get_args()
    validate_args(args)
    bcsx = BarcodeSeqExtractor(args.input_file, args.output_dir, args.quality)
    bcsx.extract_seqs(args.threads, args.max_pending * args.threads)
    bcsx.export_barcode_freq()
    print(bcsx.summarize(), end='')
