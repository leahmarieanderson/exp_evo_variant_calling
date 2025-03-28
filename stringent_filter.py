import csv
import argparse
import os

# filter values for a samtools called file
ST_QUAL_THRES = 70
ST_DP_THRES = 35
ST_Non_Coding_QUAL_THRES = 175
ST_Telomere_QUAL_THRES = 200

# filter values for a freebayes called file
FB_QUAL_THRES = 20
FB_DP_THRES = 10
FB_Non_Coding_QUAL_THRES = 600
FB_Telomere_QUAL_THRES = 650

# filter values for a lofreq called file
LOFREQ_QUAL_THRES = 20
LOFREQ_DP_THRES = 20
LOFREQ_Non_Coding_QUAL_THRES = 40
LOFREQ_Telomere_QUAL_THRES = 60

non_gff_annonations = ["missense", "intergenic", "synonymous", "5'-upstream", "nonsense", "indel-frameshift", "indel-inframe", "intron"]

# use unique filters on the txt file (which represents our simpified vcf file)
def caller_filter(caller_name, input_file, output_file):
    caller_QUAL_THRES = 0
    caller_DP_THRES = 0
    caller_NC_QUAL_THRES = 0
    caller_TELOMERE_QUAL_THRES = 0
    # set filters to corresponding values according to their caller_name
    if caller_name == "samtools":
        caller_QUAL_THRES = ST_QUAL_THRES
        caller_DP_THRES = ST_DP_THRES
        caller_NC_QUAL_THRES = ST_Non_Coding_QUAL_THRES
        caller_TELOMERE_QUAL_THRES = ST_Telomere_QUAL_THRES
    elif caller_name == "freebayes":
        caller_QUAL_THRES = FB_QUAL_THRES
        caller_DP_THRES = FB_DP_THRES
        caller_NC_QUAL_THRES = FB_Non_Coding_QUAL_THRES
        caller_TELOMERE_QUAL_THRES = FB_Telomere_QUAL_THRES
    elif caller_name == "lofreq":
        caller_QUAL_THRES = LOFREQ_QUAL_THRES
        caller_DP_THRES = LOFREQ_DP_THRES
        caller_NC_QUAL_THRES = LOFREQ_Non_Coding_QUAL_THRES
        caller_TELOMERE_QUAL_THRES = LOFREQ_Telomere_QUAL_THRES
    
    # start writing 
    with open(input_file, 'r') as cleaned_file, open(output_file, 'w', newline='') as outfile:
        # Initialize CSV reader and writer
        reader = csv.reader(cleaned_file, delimiter='\t')
        fieldnames = next(reader)  # Read header line

        # Specify the columns you want to keep
        selected_columns = ['CHROM', 'POS', 'REF', 'ALT', 'ANNOTATION', 'REGION', 'PROTEIN', "QUAL"]

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

                # if the annotation is not in non_gff_annotations, make the filter more stringent
                if anno not in non_gff_annonations:
                    # If telomere, then use the telomere qual value threshold
                    if anno == "telomere":
                        if (all(val is not None for val in [dp, mq, ref_for, ref_rev, alt_for, alt_rev]) 
                            and qual >= caller_TELOMERE_QUAL_THRES and dp >= caller_DP_THRES * 2 
                            and mq > 30 and (alt_for + alt_rev) > 4 
                            and ((ref_for + alt_for)/(ref_for + ref_rev + alt_for + alt_rev)) > 0.01 
                            and ((ref_rev + alt_rev)/(ref_for + ref_rev + alt_for + alt_rev)) > 0.01
                            ):
                            # Create a filtered row with only the selected columns
                            filtered_row = [row_dict[col] for col in filtered_fieldnames]
                            writer.writerow(filtered_row)
                    else: 
                        # Apply stringent filters since it is non-coding
                        if (all(val is not None for val in [dp, mq, ref_for, ref_rev, alt_for, alt_rev]) 
                            and qual >= caller_NC_QUAL_THRES and dp >= caller_DP_THRES * 2 
                            and mq > 30 and (alt_for + alt_rev) > 4 
                            and ((ref_for + alt_for)/(ref_for + ref_rev + alt_for + alt_rev)) > 0.01 
                            and ((ref_rev + alt_rev)/(ref_for + ref_rev + alt_for + alt_rev)) > 0.01
                            ):
                            # Create a filtered row with only the selected columns
                            filtered_row = [row_dict[col] for col in filtered_fieldnames]
                            writer.writerow(filtered_row)
                else: # it's coding, just apply regular stringent filter based on the type of caller was used
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
                
                # if the annotation is not in non_gff_annotations, make the filter more stringent
                if anno not in non_gff_annonations:
                    # If telomere, then use the telomere qual value threshold
                    if anno == "telomere":
                        if (all(val is not None for val in [dp, mqm, saf, sar, srf, srr]) 
                            and qual >= caller_TELOMERE_QUAL_THRES and dp >= caller_DP_THRES * 2 
                            and mqm > 30 and (saf + sar) > 4 
                            and ((srf + saf)/ dp) > 0.01 
                            and ((srr + sar)/ dp) > 0.01
                            ):
                            # Create a filtered row with only the selected columns
                            filtered_row = [row_dict[col] for col in filtered_fieldnames]
                            writer.writerow(filtered_row)
                    else: 
                        # Apply stringent filters since it is non-coding
                        if (all(val is not None for val in [dp, mqm, saf, sar, srf, srr]) 
                            and qual >= caller_NC_QUAL_THRES and dp >= caller_DP_THRES * 2 
                            and mqm > 30 and (saf + sar) > 4 
                            and ((srf + saf)/ dp) > 0.01 
                            and ((srr + sar)/ dp) > 0.01
                            ):
                            # Create a filtered row with only the selected columns
                            filtered_row = [row_dict[col] for col in filtered_fieldnames]
                            writer.writerow(filtered_row)
                else: # it's coding, just apply regular stringent filter based on the type of caller was used
                    if (all(val is not None for val in [dp, mqm, saf, sar, srf, srr]) 
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

                # if the annotation is not in non_gff_annotations, make the filter more stringent
                if anno not in non_gff_annonations:
                    # If telomere, then use the telomere qual value threshold
                    if anno == "telomere":
                        if (all(val is not None for val in [dp, ref_for, ref_rev, alt_for, alt_rev]) 
                            and qual >= caller_TELOMERE_QUAL_THRES and dp >= caller_DP_THRES * 2 
                            and (alt_for + alt_rev) > 4 
                            and ((ref_for + alt_for)/(ref_for + ref_rev + alt_for + alt_rev)) > 0.01 
                            and ((ref_rev + alt_rev)/(ref_for + ref_rev + alt_for + alt_rev)) > 0.01
                            ):
                            # Create a filtered row with only the selected columns
                            filtered_row = [row_dict[col] for col in filtered_fieldnames]
                            writer.writerow(filtered_row)
                    else: 
                        # Apply stringent filters since it is non-coding
                        if (all(val is not None for val in [dp, ref_for, ref_rev, alt_for, alt_rev]) 
                            and qual >= caller_NC_QUAL_THRES and dp >= caller_DP_THRES * 2 
                            and (alt_for + alt_rev) > 4 
                            and ((ref_for + alt_for)/(ref_for + ref_rev + alt_for + alt_rev)) > 0.01 
                            and ((ref_rev + alt_rev)/(ref_for + ref_rev + alt_for + alt_rev)) > 0.01
                            ):
                            # Create a filtered row with only the selected columns
                            filtered_row = [row_dict[col] for col in filtered_fieldnames]
                            writer.writerow(filtered_row)
                else: # it's coding, just apply regular stringent filter based on the type of caller was used
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
    csv_name = input_file.replace("_annotated_vcf.txt", "_condensed.txt")

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
        final_result_name = csv_name.replace('all_condensed.txt','final_stringent_compiled.txt')
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
    
    # Input CSV files
    csv_files = {
        "samtools": "",
        "freebayes": "",
        "lofreq": ""
    }

    for file in converted_files:
        if "samtools" in file:
            csv_files["samtools"] = file
        elif "freebayes" in file:
            csv_files["freebayes"] = file
        elif "lofreq" in file:
            csv_files["lofreq"] = file
        else:
            print(f"Warning: No match found for file {file}")

    # Dictionary to store merged data
    variant_dict = {}

    # Read each CSV file and store values
    for source, file in csv_files.items():
        with open(file, "r") as f:
            reader = csv.DictReader(f,  delimiter="\t")
            reader.fieldnames = [name.strip() for name in reader.fieldnames]
            for row in reader:
                key = (row["CHROM"], row["POS"], row["REF"], row["ALT"], row["ANNOTATION"], row["REGION"], row["PROTEIN"])

                if key not in variant_dict:
                    variant_dict[key] = {
                        "NUM_OCCURRENCES": 0,
                        "QUAL_samtools": None,
                        "QUAL_freebayes": None,
                        "QUAL_lofreq": None
                    }
                
                # Update count and QUAL value for the specific tool
                variant_dict[key]["NUM_OCCURRENCES"] += 1
                variant_dict[key][f"QUAL_{source}"] = row["QUAL"]

    sample_name_end_index = find_second_underscore(converted_files[0])
    sample_name = converted_files[0][:sample_name_end_index]

    temp = sample_name + '_all_condensed.txt' # make a temp csv file name
    # Write the combined CSV output
    with open(temp, "w", newline="") as f:
        fieldnames = ["CHROM", "POS", "REF", "ALT", "ANNOTATION", "REGION", "PROTEIN", 
                    "NUM_OCCURRENCES", "QUAL_samtools", "QUAL_freebayes", "QUAL_lofreq"]
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        
        for (CHROM, POS, REF, ALT, ANNOTATION, REGION, PROTEIN), values in variant_dict.items():
            writer.writerow({
                "CHROM": CHROM, "POS": POS, "REF": REF, "ALT": ALT, 
                "ANNOTATION": ANNOTATION, "REGION": REGION, "PROTEIN": PROTEIN,
                "NUM_OCCURRENCES": values["NUM_OCCURRENCES"],
                "QUAL_samtools": values["QUAL_samtools"],
                "QUAL_freebayes": values["QUAL_freebayes"],
                "QUAL_lofreq": values["QUAL_lofreq"]
            })

    sort_csv(temp) # sort the combined csv file
    # get rid of the intermediate csv file
    if os.path.exists(temp):
        os.remove(temp)

    print("Combined CSV saved to " + sample_name + "_final_stringent_compiled.txt")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine multiple CSV files into one.")
    parser.add_argument("all_file_names", type=str, nargs='+', help="Names of files to filter and combine (assumed to be in the current directory).")
    
    args = parser.parse_args()
    main(args.all_file_names)
