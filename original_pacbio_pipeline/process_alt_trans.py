from Bio import SeqIO
from Bio.Seq import Seq
import csv
import sys

def read_and_translate_fasta(file_path):
    results = []
    for record in SeqIO.parse(file_path, "fasta"):
        barcode = record.id
        dna_seq = str(record.seq)
        protein_seq = str(record.seq.translate())
        results.append((barcode, dna_seq, protein_seq))
    return results

def write_to_csv(csv_file_path, data):
    with open(csv_file_path, 'w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        for row in data:
            csv_writer.writerow(row)

def main(input_fasta, output_csv):
    data = read_and_translate_fasta(input_fasta)
    write_to_csv(output_csv, data)
    print(f"Conversion complete. Data written to {output_csv}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python process_alt_trans.py <input_fasta> <output_csv>")
        sys.exit(1)

    input_fasta = sys.argv[1]
    output_csv = sys.argv[2]
    main(input_fasta, output_csv)
