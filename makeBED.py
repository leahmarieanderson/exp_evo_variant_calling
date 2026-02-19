import argparse
import os
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description='Convert a txt file to BED file that we will use to view in IGV.')

    parser.add_argument('-i', '--input', help='Input txt file', required=True)
    parser.add_argument('-o1', '--output-high', help='High confidence BED', required=True)
    parser.add_argument('-o2', '--output-medium', help='Medium confidence BED', required=True)
    parser.add_argument('-o3', '--output-low', help='Low confidence BED', required=True)

    args = parser.parse_args()

    with open(args.input, 'r') as f_in, open(args.output_high, 'w') as f_out_high, open(args.output_medium, 'w') as f_out_medium, open(args.output_low, 'w') as f_out_low:
        header = next(f_in)
        for line in f_in:
            if line.strip():
                fields = line.strip().split()
                if len(fields) >= 3:
                    chrom = fields[0]
                    start = int(fields[1]) - 1  # Convert to 0-based indexing
                    end = int(fields[1])
                    name = fields[0] + '_' + fields[1] + '_' + fields[4]
                    if int(fields[7]) == 3:
                        f_out_high.write(f"{chrom}\t{start}\t{end}\t{name}\n")
                    elif int(fields[7]) == 2:
                        f_out_medium.write(f"{chrom}\t{start}\t{end}\t{name}\n")
                    elif int(fields[7]) == 1:
                        f_out_low.write(f"{chrom}\t{start}\t{end}\t{name}\n")

   # Remove trailing whitespace/newlines from the output files
    with open(args.output_high, 'r+') as f:
        content = f.read().rstrip()  # remove trailing whitespace/newlines
        f.seek(0)
        f.write(content)
        f.truncate()
    with open(args.output_medium, 'r+') as f:
        content = f.read().rstrip()  # remove trailing whitespace/newlines
        f.seek(0)
        f.write(content)
        f.truncate()
    with open(args.output_low, 'r+') as f:
        content = f.read().rstrip()  # remove trailing whitespace/newlines
        f.seek(0)
        f.write(content)
        f.truncate()
if __name__ == '__main__':
    main()

