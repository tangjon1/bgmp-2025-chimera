import regex
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import collections
import csv
import gzip
from concurrent.futures import ThreadPoolExecutor

# NEW
import argparse
import os
parser = argparse.ArgumentParser(prog=os.path.basename(__file__), description="Extract and count barcode sequences")
parser.add_argument('-i', '--input', help="Input files (number of files must be equal between input and output)", nargs='+')
parser.add_argument('-o', '--output', help="Output files (number of files must be equal between input and output)", nargs='+')
args = parser.parse_args()
input_filenames = args.input
output_filenames = args.output

for i in range(len(input_filenames)):#

    #input_file = "test.monomer.skera.fastq"
    #output_file = "outMP.test.monomer.skera.fasta"

    input_file = input_filenames[i]
    output_file = output_filenames[i]

    max_threads = 80  # You can change this to your desired number of threads

    d = collections.defaultdict(list)
    dnum = dict()
    outputseqs_noN = []

    motif = r'((GTCGCTGCCGAACAGC)(........................)(AGGAGAAGAGCGCACG)){e<4}'  # e<4 = less than 4 errors #sequences flanking the barcode

    #with open(input_file) as handle:
    #    records = list(SeqIO.parse(handle, "fastq"))

    #with open(input_file) as handle:
    with gzip.open(input_file, "rt") as handle:
        records = list(SeqIO.parse(handle, "fastq"))

    count = 0
    match_count = 0
    count_no_startseq = 0
    count_no_endseq_found = 0
    rc_count = 0

    print(str(len(records)), " total records")  # tells you total records found.


    def process_record(record):
        global count, match_count, count_no_startseq, count_no_endseq_found, rc_count
        count += 1  # count read
        rc_flag = 0
        seq = str(record.seq).upper()
        match = regex.search(motif, seq, regex.BESTMATCH)

        if match is None:  # check other strand, maybe it's reversed?
            seq = str(record.seq.reverse_complement()).upper()
            match = regex.search(motif, seq, regex.BESTMATCH)
            rc_flag = 1
            rc_count += 1

        if match is not None:

            match_count += 1
            barcode = match.group(3)
            bc_index = seq.find(barcode)
            sequence = seq[0:bc_index]

            d[barcode].append(sequence)  # store as a dictionary key = barcode, value = list of sequences

            if not barcode in dnum:  # new barcode not seen before
                dnum[barcode] = 1
            else:  # possible collision
                dnum[barcode] = dnum[barcode] + 1

            startseq_index = sequence.find("CATATG")
            if startseq_index != -1:
                endseq_index = sequence.find("GACGACCGCACGCTGCTG")
                if endseq_index != -1:
                    seq_between_ndei_kpni = sequence[startseq_index + 6:endseq_index]
                    if not 'N' in seq_between_ndei_kpni:
                        if seq_between_ndei_kpni != "":
                            outputseqs_noN.append(SeqRecord(id=barcode, seq=Seq(seq_between_ndei_kpni), description=""))
                else:
                    count_no_endseq_found += 1
            else:
                count_no_startseq += 1


    with ThreadPoolExecutor(max_threads) as executor:
        executor.map(process_record, records)

    print(str(rc_count), " times tried RC ", str(len(records)))
    print(str(len(d)), " BCs found out of ", str(len(records)))
    print(str(count_no_startseq), " no start (NdeI) site out of ", str(match_count))
    print(str(count_no_endseq_found), " no End  site found (but has start) out of ", str(match_count - count_no_startseq))
    print(str(len(outputseqs_noN)), " no ambiguous sites out of ",
          str(match_count - count_no_startseq - count_no_endseq_found))

    count_out3 = 0
    with open(output_file.replace(".fasta", "bc_stats.csv"), "w") as outfile:
        writer = csv.writer(outfile)
        for k, v in dnum.items():
            writer.writerow([k, str(v)])
            count_out3 += 1

    count_out3b = 0
    with open(output_file.replace(".fasta", "bc_stats_for_starcode.tsv"), "w") as outfile:
        writer = csv.writer(outfile, delimiter='\t')
        for k, v in dnum.items():
            writer.writerow([k, str(v)])
            count_out3b += 1

    count_out4 = 0
    with open(output_file.replace(".fasta", "bc_list.csv"), "w") as outfile:
        writer = csv.writer(outfile)
        for k, v in d.items():
            for ss in v:
                writer.writerow([k, str(ss)])
                count_out4 += 1

    output_file_x = open(output_file.replace(".fasta", "_noN.fasta"), 'w')
    SeqIO.write(outputseqs_noN, output_file_x, "fasta")
    output_file_x.close()

    print(str(count_out3), " write out bc_stats.csv")
    print(str(count_out4), " write out bc_list.csv")
    print(str(len(outputseqs_noN)), " write out fasta no N")

