#!/usr/bin/env python

# CLI
import argparse

# Decompress files
import gzip

# Concurrency
import concurrent.futures
import multiprocessing
from itertools import islice

# Type annotation
from typing import NamedTuple, cast, Literal, IO, Any, TextIO, Iterator, Callable
from io import TextIOWrapper
from multiprocessing.managers import DictProxy

# Misc
from os import path

class SeqRecord(NamedTuple):
    '''Tuple to describe data encompassed in FASTA/Q records'''
    id: str
    seq: str
    qual: str | None


def get_args() -> argparse.Namespace:
    '''Defines the command line arguments using argparse'''
    parser = argparse.ArgumentParser(prog=path.basename(__file__), description="Call consensus on sequences in a FASTA/Q file on a per-header basis")
    parser.add_argument('-i', '--input-file', help="Input FASTA/Q containing sequences to call consensus on", required=True)
    parser.add_argument('-c', '--cluster-file', help="Allows for mapping multiple headers to a singular header. Input TSV file with two columns, where the first column contains a singular header (does not need to exist in the FASTA/Q file), and the second column contains a series of comma-delimited headers (should be in the FASTA/Q file) to be considered for consensus calling together")
    parser.add_argument('-o', '--output-dir', help="Directory for output files", required=True)
    parser.add_argument('-f', '--input-format', help="Force the program to register the input file as a specific format", choices=['auto', 'fasta', 'fastq'], default='auto')
    parser.add_argument('-t', '--num-threads', help="The number of cores to be used by the program", type=int, default=1)
    parser.add_argument('-u', '--chunk-size', help="The number of records per chunk to use for batch multiprocessing (lower values reduce memory usage, at the cost of computation time)", type=int, default=1000)
    parser.add_argument('-m', '--max-outstanding', help="Multiplier applied to the number of cores to determine the maximum amount of chunks stored in memory", type=int, default=2)
    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    '''Validate the user-provided arguments, raising an error for incompatible configurations'''
    pass


def open_read(file_path: str) -> IO[Any] | TextIOWrapper:
    '''Takes an input file path and returns the file handle using the appropriate open() function, based on whether or not the file is compressed (ends with '.gz') or not.'''
    if file_path.endswith(".gz"):
        return gzip.open(file_path, 'rt')
    return open(file_path, 'r')


def get_fasta_records(fasta_file_handle: IO[Any] | TextIO) -> Iterator[SeqRecord]:
    '''Iteratively retrieve the header and nucleotide sequence from a FASTA file'''
    fasta_record: list[str] = list()
    for i, line in enumerate(fasta_file_handle):
        if i > 0 and line[0] == '>':
            yield SeqRecord(*cast(tuple[str, str], tuple(fasta_record)), None)
            fasta_record.clear()
        if line[0] == '>':
            fasta_record.append(line.strip('\n')[1:])
            fasta_record.append("")
        else:
            fasta_record[1] += line.strip('\n')
    yield SeqRecord(*cast(tuple[str, str], tuple(fasta_record)), None)


def get_fastq_records(fastq_file_handle: IO[Any] | TextIO) -> Iterator[SeqRecord]:
    '''Iteratively retrieve the header, nucleotide sequence, and quality score sequence from a FASTQ file'''
    fastq_record: list[str] = list()
    for i, line in enumerate(fastq_file_handle):
        if i > 0 and i % 4 == 0:
            yield SeqRecord(*fastq_record)
            fastq_record.clear()
        if i % 4 == 0 or i % 4 == 1 or i % 4 == 3:
            fastq_record.append(line.strip('\n')[(i % 4 == 0):])
    yield SeqRecord(*fastq_record)


def auto_determine_input_file_format(input_file_path: str) -> Literal['fasta', 'fastq'] | None:
    '''Automatically determine the input file format from its extension. Returns None if it cannot be determined'''
    # Ignore .gz extensions
    if input_file_path.endswith('.gz'):
        input_file_path = input_file_path.rstrip('.gz')

    # Check for the extension in these predefined sets
    _, extension = path.splitext(input_file_path)
    valid_input_format_dict: dict[str, set[str]] = {
        'fasta': {'.fa', '.fna', '.fasta'},
        'fastq': {'.fq', '.fastq'}
    }

    # Return appropriate file type if found
    for file_format in valid_input_format_dict:
        if extension in valid_input_format_dict[file_format]:
            return cast(Literal['fasta', 'fastq'], file_format)
    return None


def build_cluster_map(cluster_file_path: str) -> dict[str, str]:
    '''
    Iterate through and construct a dictionary from a header clustering file (TSV) with the following format:
        Column 1: Collapsed header (does not need to exist in the FASTA/Q file), will be used as the header in the final output file
        Column 2: Comma-delimited series of headers (that should be in the FASTA/Q file) to be considered for consensus calling together
    The resulting dictionary maps each header from the second column to the corresponding header in the first column.
    '''
    cluster_map: dict[str, str] = dict()

    with open_read(cluster_file_path) as cluster_file_handle:
        for cluster_line in cluster_file_handle:
            cluster_rep, cluster_members = cluster_line.strip('\n').split('\t')
            for cluster_member in cluster_members:
                if cluster_member in cluster_map:
                    raise KeyError("Cluster members must only appear in one cluster across the file")
                cluster_map[cluster_member] = cluster_rep

    return cluster_map


def merge_freq_dicts(primary_dict: dict[str, int], secondary_dict: dict[str, int]) -> dict[str, int]:
    '''Takes two frequency dictionaries as input (i.e., where dict values are integers) and full joins them, summing values that appear in both'''
    for key in secondary_dict:
        if key in primary_dict:
            primary_dict[key] += secondary_dict[key]
        else:
            primary_dict[key] = secondary_dict[key]
    return primary_dict


class SeqFreqDictBuilder:
    def __init__(self, input_file_path: str, input_file_format: Literal['fasta', 'fastq']):
        '''TODO'''
        self.input_file_path = input_file_path

        self.manager = multiprocessing.Manager()
        self.nt_seq_freq_dict: DictProxy[str, int] = self.manager.dict()

        # Determine appropriate iterator
        self.record_generator: Callable[[IO[Any] | TextIO], Iterator[SeqRecord]]
        match input_file_format:
            case 'fasta':
                self.record_generator = get_fasta_records
            case 'fastq':
                self.record_generator = get_fastq_records
        if self.record_generator is None:
            raise TypeError("Input file format must be 'fasta' or 'fastq'")


    @staticmethod
    def process_chunk(chunked_records: list[SeqRecord]) -> dict[str, int]:
        '''Count occurrences of each sequence in a chunk'''
        seq_freq_dict: dict[str, int] = dict()
        for record in chunked_records:
            if record.seq not in seq_freq_dict:
                seq_freq_dict[record.seq] = 1
            else:
                seq_freq_dict[record.seq] += 1
        return seq_freq_dict


    @staticmethod
    def merger_worker(result_queue: multiprocessing.Queue, primary_freq_dict: dict[str, int]):
        '''Writer process that consumes results from the queue and merges'''
        while True:
            freq_dict_chunk: list[dict[str, int]] = result_queue.get()
            if freq_dict_chunk is None:
                break
            for freq_dict in freq_dict_chunk:
                merge_freq_dicts(primary_freq_dict, freq_dict)


    def get_file_chunked_records(self, input_file_path: str, chunk_size: int) -> Iterator[list[SeqRecord]]:
        '''Iteratively yield chunks of records from a FASTA/Q file'''
        if self.record_generator is None:
            raise TypeError("Input file format must be 'fasta' or 'fastq'")

        with open_read(input_file_path) as input_file_handle:
            record_generator = self.record_generator(input_file_handle)
            while True:
                chunk: list[SeqRecord] = list(islice(record_generator, chunk_size))
                if not chunk:
                    break
                yield chunk


    def build_seq_freq_dict(self, num_threads: int, chunk_size: int, max_outstanding: int) -> DictProxy[str, int]:
        '''Iterate through the file, constructing a dictionary where the keys are sequences'''
        result_queue = multiprocessing.Queue()
        merger_process = multiprocessing.Process(target=self.merger_worker, args=(result_queue, self.nt_seq_freq_dict))
        merger_process.start()

        with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
            futures: set[concurrent.futures.Future] = set()
            
            # Submit chunked records for processing
            for chunk in self.get_file_chunked_records(self.input_file_path, chunk_size):
                futures.add(executor.submit(self.process_chunk, chunk))

                # Wait for outstanding futures
                if len(futures) >= max_outstanding:
                    done, _ = concurrent.futures.wait(futures, return_when=concurrent.futures.FIRST_COMPLETED)

                    # Merge processed chunk data
                    for future in done:
                        result_queue.put(future.result())
                        futures.remove(future)

            # Merge remaining data
            for future in concurrent.futures.as_completed(futures):
                result_queue.put(future.result())
        
        result_queue.put(None)
        merger_process.join()

        return self.nt_seq_freq_dict


class ConsensusCaller:
    def __init__(self, input_file_path: str, output_dir_path: str):
        self.manager = multiprocessing.Manager()


if __name__ == '__main__':
    args = get_args()
    validate_args(args)

    # Determine input file format
    if args.input_format == 'auto':
        input_file_format = auto_determine_input_file_format(args.input_file)
        if input_file_format is None:
            raise argparse.ArgumentError(cast(Any, args.input_format), "Input file format could not be automatically detected")

    if args.cluster_file is not None:
        cluster_map: dict[str, str] = build_cluster_map(args.cluster_file)
    
    sfdb = SeqFreqDictBuilder(args.input_file, input_file_format)
    nt_seq_freq_dict: dict[str, int] = sfdb.build_seq_freq_dict(args.num_threads, args.chunk_size, (args.num_threads * args.max_outstanding))

    for seq in nt_seq_freq_dict:
        print(f"{nt_seq_freq_dict[seq]} {len(seq)}")
