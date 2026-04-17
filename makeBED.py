import argparse
import csv
import os
from pathlib import Path

def write_bed(input_file, output_file):
    with open(input_file, 'r', newline='') as f_in, open(output_file, 'w') as f_out:
        reader = csv.DictReader(f_in, delimiter='\t')
        for row in reader:
            chrom = row['CHROM']
            start = int(row['POS']) - 1  # Convert to 0-based indexing
            end = int(row['POS'])
            name = f"{row['REF']}|{row['ALT']}|{row['ANNOTATION']}|{row['REGION']}|{row['GENE']}|{row['PROTEIN']}"
            f_out.write(f"{chrom}\t{start}\t{end}\t{name}\n")

    with open(output_file, 'r+') as f:
        content = f.read().rstrip()
        f.seek(0)
        f.write(content)
        f.truncate()

def main():
    parser = argparse.ArgumentParser(description='Convert compiled variant txt files to BED files for IGV viewing.')

    parser.add_argument('-i1', '--input-all', help='stringent_compiled.txt (all variants)', required=True)
    parser.add_argument('-i2', '--input-high', help='final_stringent_compiled.txt (2+ caller variants)', required=True)
    parser.add_argument('-o1', '--output-all', help='Output BED for all variants', required=True)
    parser.add_argument('-o2', '--output-high', help='Output BED for high-confidence variants', required=True)

    args = parser.parse_args()

    write_bed(args.input_all, args.output_all)
    write_bed(args.input_high, args.output_high)

if __name__ == '__main__':
    main()
