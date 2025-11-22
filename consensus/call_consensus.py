#!/usr/bin/env python

# BGMP 2025-26

# CLI
import argparse

# Decompress files
import gzip

# Sorting
import tempfile
import subprocess

# Concurrency
import concurrent.futures
import multiprocessing

# Translation
from Bio.Seq import Seq

# Type annotation
from typing import NamedTuple, cast, Literal, IO, Any, TextIO, Iterator, Callable
from io import TextIOWrapper
from multiprocessing.managers import DictProxy

# Random sampling
import random

# Consensus calling
from dataclasses import dataclass, field
from collections import defaultdict
from operator import attrgetter
from statistics import mean

# STDERR access
import sys

# Misc
from os import path

class SeqRecord(NamedTuple):
    '''Tuple to describe data encompassed in FASTA/Q records'''
    id: str
    seq: str
    qual: str | None


@dataclass(slots=True)
class AggregateSeqData:
    '''Dataclass to describe various summarized attributes of records with identical sequences in a cluster for consensus calling'''
    # The absolute frequency of this sequence in the cluster
    abs_freq: int = 0
    # The relative frequency of this sequence in the cluster
    rel_freq: int | float = -1
    # The mean quality score across all records with this sequence in the cluster
    avg_quals: list[int | float] = field(default_factory=list)
    # The mean quality score of the list of means
    avg_qual: int | float = -1


def unit_interval(value: str):
    '''TODO'''
    float_value: float = float(value)
    if float_value < 0 or float_value > 1:
        raise argparse.ArgumentTypeError(f"Value must be within the interval [0, 1]")
    return float_value


def get_args() -> argparse.Namespace:
    '''Defines the command line arguments using argparse'''
    parser = argparse.ArgumentParser(prog=path.basename(__file__), description="Call consensus on sequences in a FASTA/Q file on a per-header basis")
    subparsers = parser.add_subparsers(dest="command", required=True)
    
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument('-i', '--input-file', help="Input FASTA/Q containing sequences to call consensus on", required=True)
    common.add_argument('-o', '--output-file', help="Output sequence file", required=True)
    common.add_argument('-f', '--input-format', help="Force the program to register the input file as a specific format", choices=['auto', 'fasta', 'fastq'], default='auto')
    common.add_argument('-t', '--num-threads', help="The number of cores to be used by the program", type=int, default=1)

    sort_parser = subparsers.add_parser("sort", help="Sort a FASTA/Q file by header, optionally renaming headers based on a provided clustering file; headers that are members of clusters will be renamed to the cluster representative, and sorted accordingly", parents=[common])
    sort_parser.add_argument('-c', '--cluster-file', help="Allows for mapping multiple headers to a singular header. Input TSV file with two columns, where the first column contains a singular header (does not need to exist in the FASTA/Q file), and the second column contains a series of comma-delimited headers (should be in the FASTA/Q file) to be considered for consensus calling together")

    consensus_parser = subparsers.add_parser("consensus", help=f"Call consensus for sequences in a sorted FASTA/Q file (see the '{path.basename(__file__)} sort' command) on a per-header basis", parents=[common])
    consensus_parser.add_argument('-s', '--sample-size', help="For clusters with more sequences than this specified size, use a random sample of sequences from the cluster instead", type=int, default=500)
    consensus_parser.add_argument('-r', '--random-seed', help="Set a random seed for random sampling (reproducibility)")
    consensus_parser.add_argument('-p', '--plurality-min', help="A value p in the interval [0, 1]. Sequences must occur at a relative frequency greater than p across the cluster to be considered for consensus. 0: No minimum relative frequency; 1: Sequence must comprise entire cluster", type=unit_interval, default=0.5)
    consensus_parser.add_argument('-m', '--collision-margin', help="A value m in the interval [0, 1]. The difference in relative frequency between the sequences with the highest and second-highest frequencies must be greater than m in order to not be considered a collision", type=unit_interval, default=0.1)

    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    '''Validate the user-provided arguments, raising an error or warnings for incompatible configurations'''
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


def auto_determine_input_file_format(input_file_path: str) -> Literal['fasta', 'fastq']:
    '''Automatically determine the input file format from its extension. Raises an error if it cannot be determined'''
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
    raise Exception("Input format was unable to be automatically determined. Use --input-format to specify the correct format.")


def get_record_generator(input_file_format: Literal['fasta', 'fastq']) -> Callable[[IO[Any] | TextIO], Iterator[SeqRecord]]:
    '''Get the appropriate record Generator function for the provided file format'''
    match input_file_format:
        case 'fasta':
            return get_fasta_records
        case 'fastq':
            return get_fastq_records


class ClusterSorter:
    def __init__(self, input_file_path: str, cluster_file_path: str | None, output_sorted_file: str):
        '''TODO'''
        self.input_file_path = input_file_path
        self.input_file_format = auto_determine_input_file_format(self.input_file_path)
        self.output_sorted_file = output_sorted_file
        
        self.record_generator = get_record_generator(self.input_file_format)

        self.cluster_map: dict[str, str] | None = None
        if cluster_file_path is not None:
            self.cluster_map = self.build_cluster_map(cluster_file_path)


    @staticmethod
    def build_cluster_map(cluster_file_path: str) -> dict[str, str]:
        '''
        Iterate through and construct a dictionary from a header clustering file (TSV) with the following format:
            Column 1: Collapsed header (does not need to exist in the FASTA/Q file), will be used as the header in the final output file
            Column 2: Comma-delimited series of headers (that should be in the FASTA/Q file) to be considered for consensus calling together
        The resulting dictionary maps each header from the second column to the corresponding header in the first column.
        '''
        cluster_map: dict[str, str] = dict()

        with open_read(cluster_file_path) as cluster_file_handle:
            for i, cluster_line in enumerate(cluster_file_handle):
                cluster_rep, cluster_members = cluster_line.strip('\n').split('\t')
                cluster_members = cluster_members.split(',')
                for cluster_member in cluster_members:
                    if cluster_member in cluster_map:
                        raise KeyError(f"Cluster members must only appear in one cluster across the file: Duplicate '{cluster_member}' in cluster {i}")
                    cluster_map[cluster_member] = cluster_rep
        return cluster_map


    def get_cluster_rep(self, member: str) -> str:
        '''
        Get the corresponding cluster representative from the cluster map
        Returns the same string if there is no cluster map or if it is not a member of a cluster
        '''
        if self.cluster_map is None or member not in self.cluster_map:
            return member
        return self.cluster_map[member]


    def sort_file(self, num_threads: int):
        '''
        Create a temp copy of the input file with the headers renamed to the cluster representative as defined in the cluster map, and sorts by header
        Converts multi-line FASTA formats to one-line FASTA format
        '''
        # Use the unit separator character as a delimiter, permitting whitespace in FASTA/Q headers
        delim: str = '\x1F'
        with (
            open_read(self.input_file_path) as input_file_handle,
            tempfile.NamedTemporaryFile('w+') as temp_file_handle,
            open(self.output_sorted_file, 'w') as output_sorted_file
        ):
            for seq_record in self.record_generator(input_file_handle):
                header = self.get_cluster_rep(seq_record.id)
                match self.input_file_format:
                    case 'fasta':
                        temp_file_handle.write(f">{header}{delim}{seq_record.seq}\n")
                    case 'fastq':
                        temp_file_handle.write(f"@{header}{delim}{seq_record.seq}{delim}+{delim}{seq_record.qual}\n")
            # Ensure the file is complete (clear write buffer)
            temp_file_handle.flush()
            # Point to the beginning of the header-renamed file
            temp_file_handle.seek(0)

            # Run Unix sort on the file, sorting by the first field,
            # then replace the unit separators with newlines using sed
            p1 = subprocess.Popen(
                ["sort", "-t", delim, "-k1,1", f"--parallel={int(num_threads)}"],
                stdin=temp_file_handle,
                stdout=subprocess.PIPE,
                text=True,
                env={"LC_ALL": "C"}
            )

            p2 = subprocess.Popen(
                ["sed", f"s/{delim}/\\n/g"],
                stdin=p1.stdout,
                stdout=output_sorted_file,
                text=True
            )

            p2.communicate()


class ConsensusCaller:
    def __init__(self, input_file_path: str, output_file_path: str):
        '''TODO'''
        self.input_file_path = input_file_path
        self.output_file_path = output_file_path

        # Determine file format
        self.input_file_format = auto_determine_input_file_format(input_file_path)

        # Determine the appropriate generator function to extract records
        self.record_generator = get_record_generator(self.input_file_format)
        
        self.nt_seq_freq_dict: dict[str, int] = dict() #self.build_seq_freq_dict_from_file(input_file_path, self.record_generator, False)
        self.aa_seq_freq_dict: dict[str, int] = dict() #self.build_seq_freq_dict_from_file(input_file_path, self.record_generator, True)
    

    @staticmethod
    def build_seq_freq_dict_from_file(input_file_path: str, record_generator: Callable[[IO[Any] | TextIO], Iterator[SeqRecord]], protein: bool) -> dict[str, int]:
        '''TODO'''
        seq_freq_dict: dict[str, int] = dict()
        
        with open_read(input_file_path) as input_file_handle:
            for seq_record in record_generator(input_file_handle):
                seq: str = seq_record.seq
                if protein:
                    seq = str(Seq(seq).translate())
                if seq not in seq_freq_dict:
                    seq_freq_dict[seq] = 1
                else:
                    seq_freq_dict[seq] += 1
        
        return seq_freq_dict
    

    def get_cluster_sample(self, input_file_handle: IO[Any] | TextIO, max_sample_size: int) -> Iterator[list[SeqRecord]]:
        '''
        Yield each group of records that share headers, containing up to max_sample_size records (reservoir random sampling)
        Values of max_sample_size that are <1 remove the limit (i.e., yield entire clusters)
        ''' 
        record_generator: Iterator[SeqRecord] = self.record_generator(input_file_handle)
        sample: list[SeqRecord] = list()
        header: str | None = None
        i: int = 0
        for seq_record in record_generator:
            if header is not None and seq_record.id != header:
                header = seq_record.id
                yield sample
                sample.clear()
                i = 0
            if i < max_sample_size or max_sample_size < 1:
                sample.append(seq_record)
            else:
                j: int = random.randint(0, i)
                if j < max_sample_size:
                    sample[j] = seq_record
            header = seq_record.id
            i += 1
        yield sample


    @staticmethod
    def phred_offset(char: str, phred_offset: int=33) -> int:
        '''TODO'''
        return ord(char) - phred_offset


    @staticmethod
    def qual_mean(qual_seq: str) -> int | float:
        '''Returns the mean quality score value of a quality sequence'''
        return sum(map(ConsensusCaller.phred_offset, qual_seq)) / len(qual_seq)
    

    @staticmethod
    def call_consensus_for_cluster(cluster: list[SeqRecord], plurality_min: float, collision_margin: float, nt_seq_freq_dict: dict[str, int] | None, aa_seq_freq_dict: dict[str, int] | None) -> tuple[SeqRecord, float, int] | None:
        '''TODO'''
        # Collect data on a per-sequence basis within the cluster, populate dictionary
        cluster_data_by_seq: dict[str, AggregateSeqData] = defaultdict(AggregateSeqData)
        for seq_record in cluster:
            cluster_data_by_seq[seq_record.seq].abs_freq += 1
            # if seq_record.qual is not None:
            #     cluster_data_by_seq[seq_record.seq].avg_quals.append(ConsensusCaller.qual_mean(seq_record.qual))
        
        # Calculate the relative frequency for each sequence, storing the top two records with the highest relative frequency
        top1: tuple[str | None, AggregateSeqData] = (None, AggregateSeqData())
        top2: tuple[str | None, AggregateSeqData] = (None, AggregateSeqData())
        for k, v in cluster_data_by_seq.items():
            rel_freq = cluster_data_by_seq[k].abs_freq / len(cluster)
            print(cluster_data_by_seq[k].abs_freq, len(cluster), rel_freq)
            cluster_data_by_seq[k].rel_freq = rel_freq
            if rel_freq > top1[1].rel_freq:
                top2 = top1
                top1 = (k, v)
            elif rel_freq > top2[1].rel_freq:
                top2 = (k, v)
        
        # For plurality_min of 1, we want to consider a rel_freq of 1 to be 'higher', for comparison purposes
        if top1[1].rel_freq == 1:
            top1[1].rel_freq = float('inf')

        # Return the highest-plurality record
        if top1[1].rel_freq > plurality_min and top1[1].rel_freq - top2[1].rel_freq > collision_margin:
            return SeqRecord(seq_record.id, cast(str, top1[0]), None), min(top1[1].rel_freq, 1), len(cluster)


    def call_consensuses(self, max_sample_size: int, plurality_min: float, collision_margin: float):
        '''TODO'''
        with (
            open_read(self.input_file_path) as input_file_handle,
            open(self.output_file_path, 'w') as output_file_handle
        ):
            output_file_handle.write(f"barcode\tconsensus_sequence_nt\tconsensus_sequence_aa\trelative_freq\tcluster_size\n")
            for cluster in self.get_cluster_sample(input_file_handle, max_sample_size):
                # Get the result, continue if there is no consensus
                cluster_result: tuple[SeqRecord, float, int] | None = self.call_consensus_for_cluster(cluster, plurality_min, collision_margin, self.nt_seq_freq_dict, self.aa_seq_freq_dict)
                if cluster_result is None:
                    continue
                consensus_seq_record, plurality, cluster_size = cluster_result
                output_file_handle.write(f"{consensus_seq_record.id}\t{consensus_seq_record.seq}\t{str(Seq(consensus_seq_record.seq).translate())}\t{plurality:.3f}\t{cluster_size}\n")

                
if __name__ == '__main__':
    args = get_args()
    validate_args(args)

    match args.command:
        case 'sort':
            cs = ClusterSorter(args.input_file, args.cluster_file, args.output_file)
            cs.sort_file(args.num_threads)
        case 'consensus':
            random.seed(args.random_seed)
            cc = ConsensusCaller(args.input_file, args.output_file)
            cc.call_consensuses(args.sample_size, args.plurality_min, args.collision_margin)

    # for seq in cc.nt_seq_freq_dict:
    #     print(f"{cc.nt_seq_freq_dict[seq]} {len(seq)}")
