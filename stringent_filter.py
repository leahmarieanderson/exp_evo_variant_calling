import csv
import argparse
import os

# filter values for a samtools called file
ST_QUAL_THRES = 70
ST_DP_THRES = 35

# filter values for a freebayes called file
FB_QUAL_THRES = 20
FB_DP_THRES = 10

# filter values for a lofreq called file
LOFREQ_QUAL_THRES = 20
LOFREQ_DP_THRES = 20

# use unique filters on the txt file (which represents our simpified vcf file)
def caller_filter(caller_name, input_file, output_file):
    caller_QUAL_THRES = 0
    caller_DP_THRES = 0
    # set filters to corresponding values according to their caller_name
    if caller_name == "samtools":
        caller_QUAL_THRES = ST_QUAL_THRES
        caller_DP_THRES = ST_DP_THRES
    elif caller_name == "freebayes":
        caller_QUAL_THRES = FB_QUAL_THRES
        caller_DP_THRES = FB_DP_THRES
    elif caller_name == "lofreq":
        caller_QUAL_THRES = LOFREQ_QUAL_THRES
        caller_DP_THRES = LOFREQ_DP_THRES
    
    # start writing 
    with open(input_file, 'r') as cleaned_file, open(output_file, 'w', newline='') as outfile:
        # Initialize CSV reader and writer
        reader = csv.reader(cleaned_file, delimiter='\t')
        fieldnames = next(reader)  # Read header line

        # Specify the columns you want to keep
        selected_columns = ['CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'ANNOTATION', 'REGION', 'PROTEIN']

        # Filter the fieldnames to keep only the selected columns
        filtered_fieldnames = [name for name in fieldnames if name in selected_columns]

        # Create a CSV writer with the filtered fieldnames
        writer = csv.writer(outfile, delimiter='\t')
        writer.writerow(filtered_fieldnames)  # Write filtered header to the output file

        for row in reader:        
            # Create a dictionary for easier access
            row_dict = dict(zip(fieldnames, row))

            # Extract QUAL and DP values
            qual = float(row_dict['QUAL'])
            info = row_dict['INFO']
            anno = row_dict['ANNOTATION']
            dp = None
            
            # filter based on values that we decided worked best for samtools
            if caller_name == "samtools":
                mq = None
                ref_for = None # DP4[0] in the vcf
                ref_rev = None # DP4[1] in the vcf
                alt_for = None # DP4[2] in the vcf
                alt_rev = None # DP4[3] in the vcf

                # Parse DP from INFO field since INFO field contains many variables
                for entry in info.split(';'):
                    
                    if entry.startswith('DP='): # get DP
                        dp = float(entry.split('=')[1])
                    
                    if entry.startswith('DP4='):
                        values = entry[4:]
                        ref_for = float(values.split(',')[0]) # get DP4[0] 
                        ref_rev = float(values.split(',')[1]) # get DP4[1]
                        alt_for = float(values.split(',')[2]) # get DP4[2]
                        alt_rev = float(values.split(',')[3]) # get DP4[3]
                    
                    if entry.startswith('MQ='): # get MQ
                        mq = float(entry.split('=')[1])

                # if the annotation is non-coding, make the filter more stringent
                if anno == "non-coding":
                    # Apply EXTRA stringent filters since it is non-coding
                    if (all(val is not None for val in [dp, mq, ref_for, ref_rev, alt_for, alt_rev]) 
                        and qual >= caller_QUAL_THRES * 2 and dp >= caller_DP_THRES * 2 
                        and mq > 30 and (alt_for + alt_rev) > 4 
                        and ((ref_for + alt_for)/(ref_for + ref_rev + alt_for + alt_rev)) > 0.01 
                        and ((ref_rev + alt_rev)/(ref_for + ref_rev + alt_for + alt_rev)) > 0.01
                        ):
                        # Create a filtered row with only the selected columns
                        filtered_row = [row_dict[col] for col in filtered_fieldnames]
                        writer.writerow(filtered_row)
                else: # just apply regular stringent filter based on the type of caller was used
                    if (all(val is not None for val in [dp, mq, ref_for, ref_rev, alt_for, alt_rev])  # all variables aren't None
                        and qual >= caller_QUAL_THRES and dp >= caller_DP_THRES  # greater than or equal to our QUAL and DP thresholds
                        and mq > 30 and (alt_for + alt_rev) >= 4  # Mapping Quality is higher than 30 and alt read depth greater than 4 on both forward and reverse combined
                        and ((ref_for + alt_for)/(ref_for + ref_rev + alt_for + alt_rev)) > 0.01 # percentage of forward read depth greater than 1 percent
                        and ((ref_rev + alt_rev)/(ref_for + ref_rev + alt_for + alt_rev)) > 0.01 # percentage of reverse read depth greater than 1 percent
                        ):
                        # Create a filtered row with only the selected columns
                        filtered_row = [row_dict[col] for col in filtered_fieldnames]
                        writer.writerow(filtered_row)

            # else if caller name is "freebayes", then we will filter on values we decided is best for freebayes
            elif caller_name == "freebayes":
                mqm = None # Mapping Quality Mean of Alt Alleles
                mqmr = None # Mapping Quality Mean of Reference Alleles
                srf = None # Number of reference observations on the forward strand, same as samtools DP[0] just different name
                srr = None # Number of reference observations on the reverse strand
                saf = None # Number of alternate observations on the forward strand
                sar = None # Number of alternate observations on the reverse strand
                
                # Parse DP from INFO field since INFO field contains many variables
                for entry in info.split(';'):
                    
                    if entry.startswith('DP='): # get DP
                        dp = float(entry.split('=')[1])
                    
                    elif entry.startswith('MQM='): # get MQM
                        num_reads = entry.split('=')[1].split(',') # Get value(s) after 'MQM=' there will be more than one value if it is multi-allelic
                        mqm = sum(float(read_num) for read_num in num_reads) # convert all values into floats and sum them up
                    
                    elif entry.startswith('MQMR='): # get MQMR
                        num_reads = entry.split('=')[1].split(',') # Get value(s) after 'MQMR=' there will be more than one value if it is multi-allelic
                        mqmr = sum(float(read_num) for read_num in num_reads)

                    elif entry.startswith('SAF='): # get SAF
                        num_reads = entry.split('=')[1].split(',') # Get value(s) after 'SAF=' there will be more than one value if it is multi-allelic
                        saf = sum(float(read_num) for read_num in num_reads)
                    
                    elif entry.startswith('SAR='): # get SAR
                        num_reads = entry.split('=')[1].split(',') # Get value(s) after 'SAR=' there will be more than one value if it is multi-allelic
                        sar = sum(float(read_num) for read_num in num_reads)
                    
                    elif entry.startswith('SRF='): # get SRF
                        num_reads = entry.split('=')[1].split(',') # Get value(s) after 'SRF=' there will be more than one value if it is multi-allelic
                        srf = sum(float(read_num) for read_num in num_reads)

                    elif entry.startswith('SRR='): # get SRR
                        num_reads = entry.split('=')[1].split(',') # Get value(s) after 'SRR=' there will be more than one value if it is multi-allelic
                        srr = sum(float(read_num) for read_num in num_reads)
                
                # if the annotation is non-coding, make the filter more stringent
                if anno == "non-coding":
                     # Apply EXTRA stringent filters since it is non-coding
                    if (all(val is not None for val in [dp, mqm, mqmr, saf, sar, srf, srr]) 
                        and qual >= caller_QUAL_THRES * 2 and dp >= caller_DP_THRES * 2 
                        and mqm > 30 and (saf + sar) > 4 
                        and ((srf + saf)/ dp) > 0.01 
                        and ((srr + sar)/ dp) > 0.01
                        ):
                        # Create a filtered row with only the selected columns
                        filtered_row = [row_dict[col] for col in filtered_fieldnames]
                        writer.writerow(filtered_row)
                else: # just apply regular stringent filter based on the type of caller was used
                    if (all(val is not None for val in [dp, mqm, mqmr, saf, sar, srf, srr]) 
                        and qual >= caller_QUAL_THRES and dp >= caller_DP_THRES
                        and mqm > 30 and (saf + sar) > 4 
                        and ((srf + saf)/ dp) > 0.01 
                        and ((srr + sar)/ dp) > 0.01
                        ):
                        # Create a filtered row with only the selected columns
                        filtered_row = [row_dict[col] for col in filtered_fieldnames]
                        writer.writerow(filtered_row)
            
            # else if caller name is "lofreq", then we will filter on values we decided is best for lofreq    
            elif caller_name == "lofreq":
                ref_for = None # DP4[0] in the vcf
                ref_rev = None # DP4[1] in the vcf
                alt_for = None # DP4[2] in the vcf
                alt_rev = None # DP4[3] in the vcf

                # Parse DP from INFO field since INFO field contains many variables
                for entry in info.split(';'):
                    
                    if entry.startswith('DP='): # get DP
                        dp = float(entry.split('=')[1])
                    
                    if entry.startswith('DP4='):
                        values = entry[4:]
                        ref_for = float(values.split(',')[0]) # get DP4[0] 
                        ref_rev = float(values.split(',')[1]) # get DP4[1]
                        alt_for = float(values.split(',')[2]) # get DP4[2]
                        alt_rev = float(values.split(',')[3]) # get DP4[3]
                    
                    if entry.startswith('MQ='): # get MQ
                        mq = float(entry.split('=')[1])

                # if the annotation is non-coding, make the filter more stringent
                if anno == "non-coding":
                    # Apply EXTRA stringent filters since it is non-coding
                    if (all(val is not None for val in [dp, ref_for, ref_rev, alt_for, alt_rev]) 
                        and qual >= caller_QUAL_THRES * 2 and dp >= caller_DP_THRES * 2 
                        and (alt_for + alt_rev) > 4 
                        and ((ref_for + alt_for)/(ref_for + ref_rev + alt_for + alt_rev)) > 0.01 
                        and ((ref_rev + alt_rev)/(ref_for + ref_rev + alt_for + alt_rev)) > 0.01
                        ):
                        # Create a filtered row with only the selected columns
                        filtered_row = [row_dict[col] for col in filtered_fieldnames]
                        writer.writerow(filtered_row)
                else: # just apply regular stringent filter based on the type of caller was used
                    if (all(val is not None for val in [dp, ref_for, ref_rev, alt_for, alt_rev])  # all variables aren't None
                        and qual >= caller_QUAL_THRES and dp >= caller_DP_THRES  # greater than or equal to our QUAL and DP thresholds
                        and (alt_for + alt_rev) >= 4  # Mapping Quality is higher than 30 and alt read depth greater than 4 on both forward and reverse combined
                        and ((ref_for + alt_for)/(ref_for + ref_rev + alt_for + alt_rev)) > 0.01 # percentage of forward read depth greater than 1 percent
                        and ((ref_rev + alt_rev)/(ref_for + ref_rev + alt_for + alt_rev)) > 0.01 # percentage of reverse read depth greater than 1 percent
                        ):
                        # Create a filtered row with only the selected columns
                        filtered_row = [row_dict[col] for col in filtered_fieldnames]
                        writer.writerow(filtered_row)

# filters the annotated_vcf.txt files and converts them into csv format
def filter_vcf(input_file):
    # variables so we can set filter thresholds appropriately
    caller_name = ""

    # Getting rid of all the comments about INFO
    with open(input_file, 'r') as infile:
        # Read through each line to determine whether it is data or just a comment
        input_lines = infile.readlines()
        cleaned_lines = []
        for line in input_lines:
            # First, find the type of caller used for the input file
            if line.startswith("##"):
                # Set up filter according to the caller used
                if "samtools" in line: 
                    caller_name = "samtools"
                elif "freebayes" in line: 
                    caller_name = "freebayes"
                elif "lofreq" in line:
                    caller_name = "lofreq"
            else: # just append the line (which should be data and not a comment) into the list
                if(line.startswith("#CHROM")): # just remove the '#' character from the start of our header line
                    cleaned_lines.append(line[1:])
                else:
                    cleaned_lines.append(line)
            
        # Create a new txt file that holds just the headers and data
        with open("temp.txt", "w") as temp_file:
            temp_file.writelines(cleaned_lines)

    # Create a new name for output file
    csv_name = input_file.replace("_annotated_vcf.txt", "_condensed.csv")

    # Filter the temp file we created based on caller used to call the variants
    caller_filter(caller_name, "temp.txt", csv_name)
    
    # remove temp.txt
    if os.path.exists("temp.txt"):
        os.remove("temp.txt")
    # return name of converted file
    return csv_name

# helper function to isolate sample name from file name
def find_second_underscore(s):
    first_index = s.find('_')
    if first_index == -1:
        return -1  # No underscores found
    second_index = s.find('_', first_index + 1)
    return second_index

# function that helps convert CHROM column strings to ints for sorting
def chromosome_conversion(chrom_number):
	chrom_conv = {'I':1, 'II':2, 'III':3, 'IV':4, 'V':5, 'VI':6, 'VII':7, 'VIII':8, 'IX':9, 'X':10, 
				'XI':11, 'XII':12, 'XIII':13, 'XIV':14, 'XV':15, 'XVI':16, 'Mito':17, 'mitochondrion':17, 'M':17, '0M':17}
	if chrom_number.startswith('chr'):
		chrom_number = chrom_number[3:]
	try:
		if int(chrom_number) in chrom_conv.values():
			return int(chrom_number)
	except ValueError:			
		return chrom_conv[chrom_number]

# given a csv, sort based on the CHROM column and POS column
def sort_csv(csv_name):
    with open(csv_name, 'r', newline='') as infile:
        reader = csv.DictReader(infile, delimiter='\t')
        rows = list(reader)

        # Sort the rows based on 'CHROM' and 'POS'
        rows_sorted = sorted(rows, key=lambda x: (chromosome_conversion(x['CHROM']), int(x['POS'])))

        # Write the sorted rows to a new tab-delimited CSV file
        final_result_name = csv_name.replace('all_condensed.csv','final_stringent_compiled.csv')
        with open(final_result_name, 'w', newline='') as outfile:
            fieldnames = reader.fieldnames
            writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            writer.writerows(rows_sorted)


def main(all_file_names):
    # convert annotated txt files into csv files
    converted_files = []
    for txtfile in all_file_names:
        csv_name = filter_vcf(txtfile)
        converted_files.append(csv_name)

    # Dictionary to track the count of each row with key-value pair {row : count_of_duplicates}
    row_and_dupe_count = {}

    # open all csvs in the converted_files and count number of duplicates using dictionary
    for file_name in converted_files:
            with open(file_name, 'r', newline='') as file:
                csv_reader = csv.reader(file, delimiter = '\t')
                header = next(csv_reader)  # Read the header

                # Count the occurrences of each row
                for row in csv_reader:
                    row_tuple = tuple(row)
                    
                    if row_tuple in row_and_dupe_count:
                        row_and_dupe_count[row_tuple] += 1
                    else:
                        row_and_dupe_count[row_tuple] = 1

    sample_name_end_index = find_second_underscore(converted_files[0])
    sample_name = converted_files[0][:sample_name_end_index]

    temp = sample_name + '_all_condensed.csv' # make a temp csv file name

    # Open the output file in write mode
    with open(temp, 'w', newline='') as outfile:
        seen_rows = set()
        csv_writer = csv.writer(outfile, delimiter = '\t')
        header_written = False

        for file_name in converted_files:
            with open(file_name, 'r', newline='') as infile:
                csv_reader = csv.reader(infile, delimiter='\t')
                header = next(csv_reader)  # Read the header

                # Add a new header for the duplicate count
                if not header_written:
                    csv_writer.writerow(header + ['NUM_OCCURANCES'])
                    header_written = True

                # Go through each row and write 
                for row in csv_reader:
                    row_tuple = tuple(row)
                    
                    if row_tuple not in seen_rows and row_tuple in row_and_dupe_count:
                        seen_rows.add(row_tuple)
                        csv_writer.writerow(row + [row_and_dupe_count[row_tuple]])

    sort_csv(temp) # sort the combined csv file
    # get rid of the intermediate csv file
    if os.path.exists(temp):
        os.remove(temp)

    print("Combined CSV saved to " + sample_name + "_final_stringent_compiled.csv")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine multiple CSV files into one.")
    parser.add_argument("all_file_names", type=str, nargs='+', help="Names of files to filter and combine (assumed to be in the current directory).")
    
    args = parser.parse_args()
    main(args.all_file_names)
