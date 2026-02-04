#!/usr/bin/env python

# BGMP 2025-26

# CLI
import argparse

# Decompress files
import gzip

# Sorting
import tempfile
import subprocess

# Translation
from Bio.Seq import Seq

# Type annotation
from typing import NamedTuple, cast, Literal, IO, Any, TextIO, Iterator, Iterable, Callable
from io import TextIOWrapper, StringIO
from collections.abc import Sequence

# Random sampling
import random

# Multiprocessing
from multiprocessing import Pool

# Consensus calling
from dataclasses import dataclass, field
from collections import defaultdict, Counter
from statistics import mean
from math import log, exp
import editdistance

# STDERR access
import sys

# Timing
import time

# Misc
from os import path

class SeqRecord(NamedTuple):
    '''Tuple to describe data encompassed in FASTA/Q records'''
    id: str
    seq: str
    qual: str | None

@dataclass()
class SeqConsensusData:
    '''TODO'''
    # The cumulative per-base quality score error probability exponents
    # (must be exponentiated with e to derive effective error probability)
    qual_seq: list[int | float]
    # The absolute frequency of this sequence in the cluster
    abs_freq: int = 1
    # Score
    score: int | float = 0


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
    common.add_argument('-o', '--output-file', help="Output sequence file. If omitted, print to STDOUT")
    common.add_argument('-l', '--output-log', help="File for logging diagnostic information. If omitted, print to STDERR")
    common.add_argument('-f', '--input-format', help="Force the program to register the input file as a specific format", choices=['auto', 'fasta', 'fastq'], default='auto')
    common.add_argument('-t', '--num-threads', help="The number of cores to be used by the program", type=int, default=1)

    sort_parser = subparsers.add_parser("sort", help="Sort a FASTA/Q file by header, optionally renaming headers based on a provided clustering file; headers that are members of clusters will be renamed to the cluster representative, and sorted accordingly", parents=[common])
    sort_parser.add_argument('-c', '--cluster-file', help="Allows for mapping multiple headers to a singular header. Input TSV file with two columns, where the first column contains a singular header (does not need to exist in the FASTA/Q file), and the second column contains a series of comma-delimited headers (should be in the FASTA/Q file) to be considered for consensus calling together")

    consensus_parser = subparsers.add_parser("consensus", help=f"Call consensus for sequences in a sorted FASTA/Q file (see the '{path.basename(__file__)} sort' command) on a per-header basis (id-cluster)", parents=[common])
    consensus_parser.add_argument('-s', '--sample-size', help="For id-clusters with more sequences than this specified size, use a random sample of sequences from the cluster instead", type=int, default=500)
    consensus_parser.add_argument('-r', '--random-seed', help="Set a random seed for random sampling (reproducibility)")
    consensus_parser.add_argument('-g', '--gap-penalty', help="A value g such that when normalizing scores by length of sequences, add g to the divisor for every gap in the aligned sequence", type=unit_interval, default=0.1)
    consensus_parser.add_argument('-y', '--lclust-similarity', help="A value y in the interval [0, 1]. Represents the sequence similarity used when clustering by Levenshtein distance (l-cluster), as a fraction of the length of the longer of the two sequences; 1 = do not calculate distances", type=unit_interval, default=0.95)
    consensus_parser.add_argument('-p', '--lclust-plurality-min', help="A value p in the interval [0, 1]. l-cluster sequences must comprise a relative frequency greater than p across the id-cluster for the id-cluster to be considered for consensus. 0: No minimum relative frequency; 1 = l-cluster must comprise entire id-cluster", type=unit_interval, default=0.5)
    consensus_parser.add_argument('-m', '--lclust-collision-margin', help="A value m in the interval [0, 1]. The difference in relative frequency between the l-clusters with the highest and second-highest frequencies must be greater than m in order to not be considered a collision", type=unit_interval, default=0.1)

    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    '''Validate the user-provided arguments, raising an error or warnings for incompatible configurations'''
    pass


def render_args(mode, args: argparse.Namespace) -> str:
    '''Render the run mode (e.g. 'sort', 'consensus') and all set arguments as a string, separated by newlines'''
    arg_summary: str = ""
    arg_summary += f"MODE: {mode}\n"
    for parameter, argument in vars(args).items():
        arg_summary += f"\t{parameter}: {argument}\n"
    return arg_summary


def open_read(file_path: str) -> IO[Any] | TextIOWrapper:
    '''Takes an input file path and returns the file handle using the appropriate open() function, based on whether or not the file is compressed (ends with '.gz') or not.'''
    if file_path.endswith(".gz"):
        return gzip.open(file_path, 'rt')
    return open(file_path, 'r')


def open_write(file_path: str, mode='w', default_file_handle=sys.stdout):
    '''TODO'''
    if file_path == '-' or file_path is None:
        return default_file_handle
    return open(file_path, mode)


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
            open_write(self.output_sorted_file) as output_sorted_file
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
        if top1[1].rel_freq > plurality_min:
            if top1[1].rel_freq - top2[1].rel_freq > collision_margin:
                return SeqRecord(seq_record.id, cast(str, top1[0]), None), min(top1[1].rel_freq, 1), len(cluster)
            else:
                # Attempt to resolve collision
                pass # Not implemented yet


    @staticmethod
    def msa(sequences: Iterable[str]) -> Iterator[SeqRecord]:
        '''Perform multiple sequence alignment using MAFFT'''
        # Cast the list of sequences into FASTA format
        # (redundant for FASTA input files but seems like the most valid approach, need to do this anyway for FASTQ inputs)
        fasta_str: str = ""
        for i, seq in enumerate(sequences):
            fasta_str += f">S{i}\n{seq}\n"

        msa_result = subprocess.Popen(
            ["mafft", "--preservecase", "--auto", "-"],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            encoding='utf-8',
            text=True
        )
        out, _ = msa_result.communicate(input=fasta_str)

        return get_fasta_records(StringIO(out))


    @staticmethod
    def cluster_seqs_by_ldist(seqs: Sequence[str], similarity: float) -> list[list[int]]:
        """
        Cluster unique sequences by Levenshtein distance, where each cluster contains sequences within the specified distance
        The returned list contains a lists of indices, where each index corresponds to the index of the sequence in the input data
        """
        # Get the number of sequences
        # Define an empty matrix with the length of the number of sequences,
        # where each index is a list of indices that correspond to sequences that are within the max_dist threshold
        # This represents a graph where nodes are sequences and edges connect nearby sequences
        num_seqs: int = len(seqs)
        adj: list[list[int]] = [[] for _ in range(num_seqs)]

        # Populate the matrix by checking the distances of each sequence from each other;
        # For both sequences, add the corresponding sequence to the appropriate list
        for i in range(num_seqs):
            for j in range(i + 1, num_seqs):
                # Convert the similarity into the maximum distance between two sequences, based on the longer sequence
                max_dist = (1 - similarity) * max(len(seqs[i]), len(seqs[j]))
                if editdistance.eval(seqs[i], seqs[j]) <= max_dist:
                    adj[i].append(j)
                    adj[j].append(i)

        # Build the clusters from connected nodes
        # Initialize a list to track which sequences have already been visited
        visited: list[bool] = [False] * num_seqs
        clusters: list[list[int]] = []

        # Iterate through all nodes; for each unvisited node,
        # add its neighbors (recursively) to the active cluster
        for i in range(num_seqs):
            if not visited[i]:
                stack = [i]
                comp = []
                visited[i] = True
                while stack:
                    u = stack.pop()
                    comp.append(u)
                    for v in adj[u]:
                        if not visited[v]:
                            visited[v] = True
                            stack.append(v)
                clusters.append(comp)

        return clusters


    @staticmethod
    def call_consensus_for_cluster_weighted(cluster: list[SeqRecord], lcluster_similarity: float, lcluster_plurality_min: float, lcluster_collision_margin: float, gap_penalty: int | float) -> tuple[SeqRecord | None, float | None, int, float]:
        '''TODO'''
        # Time this cluster
        time_start = time.perf_counter()

        # Get the frequency of each sequence in the cluster
        seq_data_dict: dict[str, SeqConsensusData] = dict()
        for seq_record in cluster:
            if seq_record.seq not in seq_data_dict:
                seq_data_dict[seq_record.seq] = SeqConsensusData(qual_seq=[0] * len(seq_record.seq))
            else:
                seq_data_dict[seq_record.seq].abs_freq += 1
            if seq_record.qual is not None:
                for i, qual_char in enumerate(seq_record.qual):
                    seq_data_dict[seq_record.seq].qual_seq[i] += (-ConsensusCaller.phred_offset(qual_char) / 10) * log(10)
        
        # LEVENSHTEIN DISTANCE STEP --------------------------------------

        # Use Levenshtein dist to determine if enough of the sequences are drastically different (and count that as a collision)
        seq_data_seqs: list[str] = list(seq_data_dict.keys())

        # Group seqs by the given similarity
        # Ultimately produce a list relative frequencies of each lcluster (as implemented, unassociated with sequences)
        rel_lcluster_freqs: list[int | float]
        # If similarity argument is 1, no actual lclustering is required
        if lcluster_similarity < 1:
            lclusters: list[list[int]] = ConsensusCaller.cluster_seqs_by_ldist(seq_data_seqs, lcluster_similarity)
        
            # Track the total number of reads in each cluster
            lcluster_read_counts: list[int] = list()
            for lcluster in lclusters:
                cluster_sum = 0
                for seq_data_index in lcluster:
                    cluster_sum += seq_data_dict[seq_data_seqs[seq_data_index]].abs_freq
                lcluster_read_counts.append(cluster_sum)
                
            assert sum(lcluster_read_counts) == len(cluster), "Cluster read count mismatch"

            # Convert to relative frequencies
            rel_lcluster_freqs = [f / len(cluster) for f in lcluster_read_counts]
        else:
            rel_lcluster_freqs = list()
            for seq in seq_data_dict:
                rel_lcluster_freqs.append(seq_data_dict[seq].abs_freq / len(cluster))

        # Get the highest two relative frequency clusters
        top_1: tuple[int | None, float] = None, float('-inf')
        top_2: tuple[int | None, float] = None, float('-inf')
        for i, rel_lcluster_freq in enumerate(rel_lcluster_freqs):
            if rel_lcluster_freq > top_1[1]:
                top_2 = top_1
                top_1 = i, rel_lcluster_freq
            elif rel_lcluster_freq > top_2[1]:
                top_2 = i, rel_lcluster_freq

        # For plurality_min of 1, we want to consider a relative frequency of 1 to be 'higher', for comparison purposes
        if top_1[1] == 1:
            top_1 = top_1[0], float('inf')

        # Return None if the relative frequency does not meet the threshold
        if top_1[1] < lcluster_plurality_min:
            return None, None, len(cluster), time.perf_counter() - time_start
        
        # Return None if the difference between the firstmost and secondmost highest relative frequencies is too small
        # i.e., we are not confident which groups of reads are biologically correct
        if top_1[1] - top_2[1] < lcluster_collision_margin:
            return None, None, len(cluster), time.perf_counter() - time_start

        # MSA STEP -------------------------------------------------------

        # Perform multiple sequence alignment
        alignments = list(ConsensusCaller.msa(seq_data_dict))

        # Majority call on a per-base basis
        # Initialize a list (of the length of the alignment) of empty Counters
        column_counts: list[Counter] = [Counter() for _ in range(len(alignments[0].seq))]
        # Determine the frequency of each base at each alignment position
        for alignment in alignments:
            for i, char in enumerate(alignment.seq):
                column_counts[i][char] += 1

        assert len(seq_data_dict) == len(alignments), "Alignment count mismatch"

        for i, seq in enumerate(seq_data_dict):
            score: int | float = 0
            # Since dictionaries are ordered, and MAFFT yields ordered output,
            # we can assume that the alignments remain in the corresponding order
            seq_alignment = alignments[i].seq
            seq_qual = seq_data_dict[seq].qual_seq
            qual_index_offset: int = 0
            for j, aligned_base in enumerate(seq_alignment):
                if aligned_base == '-':
                    qual_index_offset += 1
                else:
                    score += (column_counts[j][aligned_base] / column_counts[j].total()) * (1 - exp(seq_qual[j - qual_index_offset]))
            # Normalize score by length (e.g., mean score per base)
            # Normalization includes a penalty for gaps
            seq_data_dict[seq].score = score / (len(seq_qual) + (qual_index_offset * gap_penalty))
        
        # Get the highest scoring sequence (a.k.a. the consensus sequence)
        consensus_seq = max(seq_data_dict, key=lambda k: seq_data_dict[k].score)

        return SeqRecord(seq_record.id, consensus_seq, None), seq_data_dict[consensus_seq].score, len(cluster), time.perf_counter() - time_start
    

    def generate_cluster_worker_input(self, input_file_handle: IO[Any] | TextIOWrapper, max_sample_size: int, lcluster_similarity: float, lcluster_plurality_min: float, lcluster_collision_margin: float, msa_gap_penalty: int | float) -> Iterator[tuple[list[SeqRecord], float, float, float, int | float]]:
        '''
        Generates the input for the consensus calling, packing it into a tuple (to pass to multiple processes)
        '''
        for cluster in self.get_cluster_sample(input_file_handle, max_sample_size):
            yield cluster.copy(), lcluster_similarity, lcluster_plurality_min, lcluster_collision_margin, msa_gap_penalty


    @staticmethod
    def wrapper_call_consensus_for_cluster_weighted(call_consensus_input: tuple[list[SeqRecord], float, float, float, int | float]) -> tuple[SeqRecord | None, float | None, int, float]:
        '''
        Wrapper function for consensus calling function that unpacks the input tuple.
        '''
        return ConsensusCaller.call_consensus_for_cluster_weighted(*call_consensus_input)


    def call_consensuses(self, num_threads: int, max_sample_size: int, lcluster_similarity: float, lcluster_plurality_min: float, lcluster_collision_margin: float, msa_gap_penalty: int | float):
        '''TODO'''
        with (
            Pool(processes=num_threads) as pool,
            open_read(self.input_file_path) as input_file_handle,
            open_write(self.output_file_path) as output_file_handle
        ):
            output_file_handle.write(f"barcode\tconsensus_sequence_nt\tconsensus_sequence_aa\tscore\n")
            for i, cluster_result in enumerate(pool.imap_unordered(ConsensusCaller.wrapper_call_consensus_for_cluster_weighted, self.generate_cluster_worker_input(input_file_handle, max_sample_size, lcluster_similarity, lcluster_plurality_min, lcluster_collision_margin, msa_gap_penalty), chunksize=50)):
                print(f"id-cluster\t{i + 1}", end='', file = sys.stderr)
                consensus_seq_record, score, cluster_size, elapsed_time = cluster_result
                if consensus_seq_record is None:
                    print(f"\tCollision", end='', file=sys.stderr)
                else:
                    print(f"\tConsensus", end='', file=sys.stderr)
                    output_file_handle.write(f"{consensus_seq_record.id}\t{consensus_seq_record.seq}\t{str(Seq(consensus_seq_record.seq[:len(consensus_seq_record.seq) - (len(consensus_seq_record.seq) % 3)]).translate())}\t{score:.5f}\n")
                    output_file_handle.flush()
                print(f"\t{cluster_size}\t{elapsed_time:.6f}", file=sys.stderr, flush=True)

                
if __name__ == '__main__':
    args = get_args()
    validate_args(args)

    match args.command:
        case 'sort' as mode:
            print(render_args(mode, args), file=sys.stderr, flush=True)
            cs = ClusterSorter(args.input_file, args.cluster_file, args.output_file)
            cs.sort_file(args.num_threads)
        case 'consensus' as mode:
            print(render_args(mode, args), file=sys.stderr, flush=True)
            if args.random_seed is None:
                random_seed = random.randint(10000, 99999)
            random.seed(random_seed)
            cc = ConsensusCaller(args.input_file, args.output_file)            
            cc.call_consensuses(args.num_threads, args.sample_size, args.lclust_similarity, args.lclust_plurality_min, args.lclust_collision_margin, args.gap_penalty)
    