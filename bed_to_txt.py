
#this script is for combining all your results from a collection of evolved yeast
#at this point, you have verified all mutations in IGV and are ready to combine the results into one file :)
#USAGE: "python3 bed_to_txt.py -d /path/to/bed/files -o /path/to/output.txt"


import argparse
import csv
import os
from pathlib import Path

KNOWN_SUFFIXES = ['_all_variants', '_highConfidenceVars']

def sample_name_from_path(bed_path):
    stem = Path(bed_path).stem
    for suffix in KNOWN_SUFFIXES:
        if stem.endswith(suffix):
            return stem[:-len(suffix)]
    return stem

def parse_bed(bed_path, writer):
    sample = sample_name_from_path(bed_path)
    with open(bed_path, 'r') as f:
        for line in f:
            if not line.strip():
                continue
            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = fields[2]  # end column is the 1-based POS
            ref, alt, annotation, region, gene, protein = fields[3].split('|')
            writer.writerow({
                'sample': sample,
                'CHROM': chrom,
                'POS': pos,
                'REF': ref,
                'ALT': alt,
                'ANNOTATION': annotation,
                'REGION': region,
                'GENE': gene,
                'PROTEIN': protein,
            })

def main():
    parser = argparse.ArgumentParser(description='Convert BED files back to a combined variant txt file.')
    parser.add_argument('-d', '--directory', default='.', help='Directory containing BED files (default: current directory)')
    parser.add_argument('-o', '--output', required=True, help='Output txt file')
    parser.add_argument('--pattern', default='*.bed', help='Glob pattern for BED files (default: *.bed)')
    args = parser.parse_args()

    bed_files = sorted(Path(args.directory).glob(args.pattern))
    if not bed_files:
        print(f"No BED files found in {args.directory} matching '{args.pattern}'")
        return

    fieldnames = ['sample', 'CHROM', 'POS', 'REF', 'ALT', 'ANNOTATION', 'REGION', 'GENE', 'PROTEIN']

    with open(args.output, 'w', newline='') as f_out:
        writer = csv.DictWriter(f_out, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for bed_file in bed_files:
            parse_bed(bed_file, writer)

    print(f"Wrote {len(bed_files)} file(s) to {args.output}")

if __name__ == '__main__':
    main()
