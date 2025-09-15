# Batch Submitter for CNV pipeline
# To run this script, do 
# python3 CNV_batch_submit.py path/to/directory ancestor_name
# ex. python3 CNV_batch_submit.py /net/dunham/vol2/Zilong/updating_pipeline_2024 AB0reseq_S34
# make sure that the path/to/directory argument contains the directory which has the results from the SNP pipeline

import os
import argparse
import shutil

def createOutputDirectories(work_dir_path):
	# Make a directory for the transposon outputs if it doesn't exist
	CNV_dir = os.path.join(work_dir_path, "CNV_graphs")
	os.makedirs(CNV_dir, exist_ok=True)

	# Make the error and output directories if they don't exist already
	outputs_dir = os.path.join(work_dir_path, "outputs")
	errors_dir = os.path.join(work_dir_path, "errors")
	os.makedirs(outputs_dir, exist_ok=True)
	os.makedirs(errors_dir, exist_ok=True)

def batchSubmit(work_dir_path, ancestor_name):
    # open sample list
    with open(f"{work_dir_path}/CNV_graphs/sample_list.txt", "r") as WD_samples:
        # iterate through the sample list; run CNV pipeline
        for sample in WD_samples:
            sample = sample.strip()

            # submit each individual sample job to cluster with qsub
            cmd = f"qsub -N 'CNV_{sample}' -wd '{work_dir_path}' -o '{work_dir_path}/outputs' -e '{work_dir_path}/errors' CNV_runner.sh {sample} 1000 1"

            if ancestor_name != "" and sample != ancestor_name:
                cmd += f" {ancestor_name}"

            os.system(cmd)

def main():
	parser = argparse.ArgumentParser(description="Batch Submit jobs to run CNV scripts on samples from WorkDirectory")
	parser.add_argument("work_dir_path", help="Path to working directory, should contain WorkDirectory with all your .bam files.")
	parser.add_argument("ancestor_name", nargs="?", default="", help="name of fastq folder to submit (default: "")")
	args = parser.parse_args()
	work_dir_path = args.work_dir_path
	ancestor_name = args.ancestor_name
	createOutputDirectories(work_dir_path)
	# get the sample list from directories inside of WorkDirectory
	os.system(f"cd {work_dir_path}/WorkDirectory && ls -d */ | sed 's:/*$::' > {work_dir_path}/CNV_graphs/sample_list.txt")
	batchSubmit(work_dir_path, ancestor_name)

if __name__ == "__main__":
    main()
