import Bio
import csv
import hashlib
from Bio import SeqIO
from Bio.Align import PairwiseAligner
import multiprocessing
import sys
import time
import os

# -----------------------------------------------------------------------------
# FUNCTIONS
# -----------------------------------------------------------------------------

def getFastaSeqs(filename):
    """
    Opens a FASTA file and returns a dictionary with sequence as value and ID as key.
    Removes the leading 'M' from the sequence.
    """
    fastaseqs = {}
    try:
        with open(filename) as handle:
            for seqrec in SeqIO.parse(handle, "fasta"):
                fastaseqs[str(seqrec.id).split(".")[0]] = str(seqrec.seq[1:])  # remove leading M
    except FileNotFoundError:
        print(f"Error: FASTA file not found: {filename}")
        sys.exit(1)
    return fastaseqs

def genMutantName(data_in):
    """
    Generates a mutant name based on the reference sequence, read sequence, and ID.
    Follows the naming convention from parse_SAM_alignments_MP.pdf:
      - For low numbers (<6) of mutations (and non-negative), list each mutation.
      - For higher numbers, use a hash of the read sequence.
    Returns: (BC, fail_flag, pct_ident, mutations, mutant_name, seq)
    """
    BC, IDalign, refseq, seq = data_in  # Unpack the Barcode, IDalign, refseq, and sequence
    fail_flag = 1
    pct_ident = 0
    mutations = -10000
    mutant_name = ''

    # Initialize the pairwise aligner
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 1.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -2.0
    aligner.extend_gap_score = -0.5

    try:
        # Perform alignment
        alignments = aligner.align(refseq, seq)
        if alignments:
            alignment = alignments[0]
            fail_flag = 0

            # Calculate percent identity using the aligned strings
            matches = sum(1 for i in range(len(alignment[0])) if alignment[0][i] == alignment[1][i])
            pct_ident = float(matches) / len(refseq) if len(refseq) > 0 else 0.0

            sizeread = len(seq.strip())
            sizeref = len(refseq.strip())
            if sizeread <= sizeref:
                mutations = int(sizeref - matches)
            else:
                mutations = -1 * int(sizeref - matches + (sizeread - sizeref)) - 10000

            # Use detailed mutation naming only if few mutations (<6) and non-negative
            if mutations < 6 and mutations >= 0:
                read_seq_align = str(alignment[1])
                ref_seq_align = str(alignment[0])
                ref_insertion_count = 0
                mutation_parts = [IDalign]
                for i in range(len(read_seq_align)):
                    if read_seq_align[i] != ref_seq_align[i]:
                        if read_seq_align[i] == '-':
                            mutation_parts.append(f"_{ref_seq_align[i]}{i + 1 - ref_insertion_count}X")
                        elif ref_seq_align[i] == '-':
                            mutation_parts.append(f"_X{i + 1 - ref_insertion_count}{read_seq_align[i]}")
                            ref_insertion_count += 1
                        else:
                            mutation_parts.append(f"_{ref_seq_align[i]}{i + 1 - ref_insertion_count}{read_seq_align[i]}")
                mutant_name = "".join(mutation_parts)
            else:
                # For many mutations, use a hash of the read sequence
                hash_object = hashlib.sha256(seq.encode('utf-8'))
                hex_dig = hash_object.hexdigest()
                mutant_name = IDalign + '_' + str(hex_dig)
        else:
            raise Exception("No alignment found")
    except Exception as e:
        print(f"Alignment error for BC {BC}: {e}")
        return BC, fail_flag, pct_ident, mutations, mutant_name, seq

    return BC, fail_flag, pct_ident, mutations, mutant_name, seq

# -----------------------------------------------------------------------------
# MAIN EXECUTION
# -----------------------------------------------------------------------------

if __name__ == '__main__':
    # -------------------------------------------------------------------------
    # Input/Output files (passed as command-line arguments)
    # -------------------------------------------------------------------------
    if len(sys.argv) != 5:
        print("Usage: python parse_sam.py <protein_seq_file> <BC_info_file> <sam_file> <analysis_dir>")
        sys.exit(1)

    protein_seq_file = sys.argv[1]    # FASTA file of protein sequences
    BC_info_file = sys.argv[2]          # CSV or FASTA file with good barcodes and translation info
    sam_file = sys.argv[3]            # SAM file from BBMap
    analysis_dir = sys.argv[4]        # Output directory

    mutID_BC_out = os.path.join(analysis_dir, 'C5seqs_mutID_all.csv')
    mutID_info_out = os.path.join(analysis_dir, 'C5seqs_mutID_info_all.csv')

    # -------------------------------------------------------------------------
    # Variables and Data Structures
    # -------------------------------------------------------------------------
    BClist = []                    # List of good barcodes
    maptransdict = {}              # BC -> translated sequence (cut at stop codon)
    seq_to_info_dict = {}          # BC -> "IDalign cigar" string (from SAM)
    seq_to_mutID_dict = {}         # seq -> mutant name (mutID)
    # Unique mutant info (aggregated by mutID)
    mutID_info_dict = {}           # key: mutID; value: dict with keys: IDalign, mutations, seq, pct_ident, BCs (set), BCcode
    # For parallel processing results
    seq_to_info2_dict = [{} for _ in range(2)]  # index 0: mutations, index 1: pct_ident

    NUM_PROCS = multiprocessing.cpu_count()   # Use all available CPUs
    BC_read_count = {}           # BC -> number of reads (from SAM)
    BC_cigar_dict = {}           # BC -> CIGAR string (from SAM)

    # Additional dictionaries from MP code:
    seq_to_BC_dict = {}          # seq -> set of BCs that produced this sequence
    seq_to_BCcode_dict = {}      # seq -> BC code (e.g., "NoGT_noStop" or "NoGT_Stop")

    # -------------------------------------------------------------------------
    # Load reference protein sequences
    # -------------------------------------------------------------------------
    print("Loading reference protein sequences...")
    proseqs = getFastaSeqs(protein_seq_file)

    # -------------------------------------------------------------------------
    # Load good barcodes from BC_info_file (if CSV, assume BC is in first column)
    # -------------------------------------------------------------------------
    print("Loading good barcodes...")
    try:
        # Try CSV approach first; if file is FASTA, this list will be constructed from headers.
        ext = os.path.splitext(BC_info_file)[1].lower()
        if ext in ['.fa', '.fasta']:
            BClist = [rec.id for rec in SeqIO.parse(BC_info_file, "fasta")]
        else:
            with open(BC_info_file) as infile:
                BClist = [line.split(',')[0] for line in infile if not line.startswith('"')]
        print(f"{len(BClist)} good barcodes loaded")
    except FileNotFoundError:
        print(f"Error: Barcode info file not found: {BC_info_file}")
        sys.exit(1)

    # -------------------------------------------------------------------------
    # Load translated sequences and build BC-related dictionaries
    # -------------------------------------------------------------------------
    print("Loading translated sequences and BC info...")
    prematurestop = []
    count_noGT_noStop = 0
    count_noGT_Stop = 0
    zero_length_trans_seq = 0

    ext = os.path.splitext(BC_info_file)[1].lower()
    try:
        if ext in ['.fa', '.fasta']:
            # Parse as FASTA: header is assumed to be the barcode; sequence is the translation.
            for seqrec in SeqIO.parse(BC_info_file, "fasta"):
                BC = seqrec.id
                seqt = str(seqrec.seq).replace(' ', '').replace('\t', '').replace('\n', '')
                stop_index = seqt.find('*')
                if stop_index == -1:
                    seqtcut = seqt
                    count_noGT_noStop += 1
                    BC_code = "NoGT_noStop"
                    prematurestop.append(BC)
                else:
                    seqtcut = seqt[:stop_index]
                    count_noGT_Stop += 1
                    BC_code = "NoGT_Stop"
                maptransdict[BC] = seqtcut
                if seqtcut in seq_to_BC_dict:
                    seq_to_BC_dict[seqtcut].add(BC)
                else:
                    seq_to_BC_dict[seqtcut] = {BC}
                if seqtcut not in seq_to_BCcode_dict:
                    seq_to_BCcode_dict[seqtcut] = BC_code
                if len(seqtcut) == 0:
                    zero_length_trans_seq += 1
        else:
            # Parse as CSV (expects at least three columns: BC, ..., sequence)
            with open(BC_info_file) as transfile:
                reader = csv.reader(transfile)
                for row in reader:
                    if len(row) < 3:
                        continue
                    BC = row[0]
                    seqt = row[2].replace(' ', '').replace('\t', '').replace('\n', '')
                    stop_index = seqt.find('*')
                    if stop_index == -1:
                        seqtcut = seqt
                        count_noGT_noStop += 1
                        BC_code = "NoGT_noStop"
                        prematurestop.append(BC)
                    else:
                        seqtcut = seqt[:stop_index]
                        count_noGT_Stop += 1
                        BC_code = "NoGT_Stop"
                    maptransdict[BC] = seqtcut
                    if seqtcut in seq_to_BC_dict:
                        seq_to_BC_dict[seqtcut].add(BC)
                    else:
                        seq_to_BC_dict[seqtcut] = {BC}
                    if seqtcut not in seq_to_BCcode_dict:
                        seq_to_BCcode_dict[seqtcut] = BC_code
                    if len(seqtcut) == 0:
                        zero_length_trans_seq += 1

        total_BCs = len(maptransdict)
        if total_BCs == 0:
            print("Warning: No translated sequences loaded from BC info file.")
        else:
            print(f"{total_BCs} translated sequences loaded, of which {len(prematurestop)} have a premature stop codon")
            print(f"No GT, no stop: {count_noGT_noStop} ({round(count_noGT_noStop / total_BCs * 100, 1)}%)")
            print(f"No GT, has stop: {count_noGT_Stop} ({round(count_noGT_Stop / total_BCs * 100, 1)}%)")
            print(f"Zero-length sequences: {zero_length_trans_seq} ({round(zero_length_trans_seq / total_BCs * 100, 1)}%)")
            print(f"{len(maptransdict)} BCs loaded")
    except FileNotFoundError:
        print(f"Error: Translation file not found: {BC_info_file}")
        sys.exit(1)

    # -------------------------------------------------------------------------
    # Load SAM alignment data: populate seq_to_info_dict, BC_cigar_dict, and BC_read_count
    # -------------------------------------------------------------------------
    print("Loading SAM alignment data...")
    try:
        with open(sam_file) as samfile:
            for line in samfile:
                if line.startswith('@'):
                    continue
                listWords = line.split("\t")
                BC = listWords[0].strip()
                if '_part_2' in BC:
                    BC = BC.split('_part_2')[0]
                IDalign = listWords[2].split(".")[0]
                cigar = listWords[5]
                # Store CIGAR string and read counts
                BC_cigar_dict[BC] = cigar
                BC_read_count[BC] = BC_read_count.get(BC, 0) + 1
                if BC in maptransdict:
                    seq = maptransdict[BC]
                    # For a given translated sequence, store the aligned ID and cigar info
                    seq_to_info_dict[seq] = IDalign + " " + cigar
        print(f"{len(seq_to_info_dict)} aligned IDs matched to sequences")
    except FileNotFoundError:
        print(f"Error: SAM file not found: {sam_file}")
        sys.exit(1)

    # -------------------------------------------------------------------------
    # Parallel mutant name generation with progress tracking
    # -------------------------------------------------------------------------
    print("Starting parallel mutant name generation...")
    # Prepare data for parallel processing.
    tasks = []
    for BC in BClist:
        if BC in maptransdict:
            seq = maptransdict[BC]
            if seq in seq_to_info_dict:
                IDalign = seq_to_info_dict[seq].split(" ")[0]
                if IDalign in proseqs:
                    tasks.append((BC, IDalign, proseqs[IDalign], seq))
    total_tasks = len(tasks)
    print(f"{total_tasks} sequences will have mutant names generated.")

    mutID_BC_data = []   # List of rows for C5seqs_mutID_all.csv (one row per BC)
    # Dictionary for unique mutant info (aggregating across BCs)
    # Key: mutant_name, Value: dict with keys: IDalign, mutations, seq, pct_ident, BCs (set), BCcode.
    mutID_info_dict = {}

    pool = multiprocessing.Pool(NUM_PROCS)
    results = []
    for task in tasks:
        results.append(pool.apply_async(genMutantName, (task,)))
    pool.close()

    start_time = time.time()
    completed = 0
    total_sequences = len(results)
    for res in results:
        BC, fail_flag, pct_ident, mutations, mutant_name, seq = res.get()
        completed += 1
        if fail_flag == 0:
            seq_to_mutID_dict[seq] = mutant_name
            seq_to_info2_dict[0][seq] = mutations
            seq_to_info2_dict[1][seq] = pct_ident

            cigar = BC_cigar_dict.get(BC, "<unknown>")
            read_count = BC_read_count.get(BC, 0)
            mutID_BC_data.append([BC, mutant_name, mutations, cigar, read_count])

            # Update unique mutant info
            BCs_for_seq = seq_to_BC_dict.get(seq, set())
            BCcode_for_seq = seq_to_BCcode_dict.get(seq, "<unknown>")
            IDalign = seq_to_info_dict[seq].split(" ")[0]
            if mutant_name not in mutID_info_dict:
                mutID_info_dict[mutant_name] = {
                    'IDalign': IDalign,
                    'mutations': mutations,
                    'seq': seq,
                    'pct_ident': pct_ident,
                    'BCs': set(),  # will collect all BCs producing this mutant
                    'BCcode': BCcode_for_seq
                }
            mutID_info_dict[mutant_name]['BCs'].update(BCs_for_seq)
        # Progress tracker
        progress = (completed / total_sequences) * 100
        elapsed_time = time.time() - start_time
        remaining_time = (elapsed_time / completed) * (total_sequences - completed) if completed > 0 else 0
        print(f"Progress: {progress:.2f}% | Elapsed Time: {elapsed_time:.2f}s | Remaining Time: {remaining_time:.2f}s", end='\r')
    print(f"\nProgress: 100.00% | Elapsed Time: {elapsed_time:.2f}s | Remaining Time: 0.00s")

    # -------------------------------------------------------------------------
    # Write output files
    # -------------------------------------------------------------------------
    print("Writing output files...")
    try:
        # Write the BC-to-mutant file (one row per BC)
        with open(mutID_BC_out, 'w', newline='') as outfi:
            writer = csv.writer(outfi)
            writer.writerow(['BC', 'mutID', 'mutations', 'cigar', 'reads'])
            for row in mutID_BC_data:
                writer.writerow(row)

        # Write the unique mutant info file (one row per unique mutID)
        with open(mutID_info_out, 'w', newline='') as outfi2:
            writer2 = csv.writer(outfi2)
            # Header: mutID, IDalign, numBCs, mutations, seq, pct_ident, BCs, BCcode
            writer2.writerow(['mutID', 'IDalign', 'numBCs', 'mutations', 'seq', 'pct_ident', 'BCs', 'BCcode'])
            for mutID, info in mutID_info_dict.items():
                numBCs = len(info['BCs'])
                # Convert set of BCs into a space-separated string
                BCs_str = " ".join(sorted(info['BCs']))
                writer2.writerow([mutID, info['IDalign'], numBCs, info['mutations'], info['seq'],
                                  info['pct_ident'], BCs_str, info['BCcode']])
        print("Output files written successfully.")
    except IOError as e:
        print(f"Error writing to file: {e}")
        sys.exit(1)

    print("Done")
    print("--- %s seconds ---" % (time.time() - start_time))
