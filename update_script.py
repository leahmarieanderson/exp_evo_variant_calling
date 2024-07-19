def get_user_input(prompt, default=None):
    user_input = input(f"\n{prompt} \nDefault is :[{default}]\n Just press enter for default value : \n")
    return user_input if user_input else default

def get_script_path():
    while True:
        script_path = input("What is the path to your align.sh script?\nExample: /net/dunham/vol2/Zilong/updating_pipeline_2024/exp_evo_variant_calling/align.sh\n")
        if script_path.endswith("align.sh"):
            return script_path
        else:
            print("Error: The path must end with 'align.sh'. Please try again.")

def main():
    # Get script path from user
    script_path = get_script_path()

    # Get new values from the user
    wd_option = get_user_input("Path of directory where your fastq directory is located", "/net/dunham/vol2/Zilong/updating_pipeline_2024")
    mfree_option = get_user_input("How much memory to allocate to job?", "8G")
    jobname_option = get_user_input("What would you like the job name to be?", "cluster_job")
    h_rt_option = get_user_input("How much runtime to allocate? (HH:MM:SS)", "36:0:0")
    FOLDER_option = get_user_input("What is the name of the directory with your fastq files? \n ex. if fastq files are found in path /net/dunham/vol2/Zilong/updating_pipeline_2024/fastq , then we would just put 'fastq'", "fastq")
    SEQID_option = get_user_input("Give a project name for the bam headers : ", "mutation")
    REF_option = get_user_input("Path to your reference genome. \n ex. /net/dunham/vol2/Zilong/updating_pipeline_2024/genomes/sacCer3.fasta", "/net/dunham/vol2/Zilong/updating_pipeline_2024/genomes/sacCer3.fasta")
    # ANNOTATE_option = get_user_input("Path to get certain annotation scripts: ", "AB0reseq_S34")
    SCRIPTS_option = get_user_input("Path to your directory containing annotation_final.py", "/net/dunham/vol2/Zilong/updating_pipeline_2024/exp_evo_variant_calling/")

    # Read the existing script
    with open(script_path, "r") as file:
        lines = file.readlines()

    # Update specific lines
    updated_lines = []
    for line in lines:
        if line.startswith("#$ -wd"):
            updated_lines.append(f"#$ -wd {wd_option}\n")
        elif line.startswith("#$ -o"):
            updated_lines.append(f"#$ -o {wd_option}/outputs/\n")
        elif line.startswith("#$ -e"):
            updated_lines.append(f"#$ -e {wd_option}/errors/\n")
        elif line.startswith("#$ -l mfree"):
            updated_lines.append(f"#$ -l mfree={mfree_option}\n")
        elif line.startswith("#$ -l h_rt"):
            updated_lines.append(f"#$ -l h_rt={h_rt_option}\n")
        elif line.startswith("#$ -N"):
            updated_lines.append(f"#$ -N {jobname_option}\n")
        elif line.startswith("FOLDER="):
            updated_lines.append(f"FOLDER={FOLDER_option}\n")
        elif line.startswith("DIR="):
            updated_lines.append(f"DIR={wd_option}\n")
        elif line.startswith("SEQID="):
            updated_lines.append(f"SEQID={SEQID_option} # Project name and date for bam header\n")
        elif line.startswith("REF="):
            updated_lines.append(f"REF={REF_option} # Reference genome\n")
        # commenting out ANNOTATE for now since I don't k
        #elif line.startswith("ANNOTATE="):
        #    updated_lines.append(f"ANNOTATE={}\n")
        elif line.startswith("SCRIPTS="):
            updated_lines.append(f"SCRIPTS={SCRIPTS_option} # Path of annotation_final.py directory\n")
        else:
            updated_lines.append(line)

    # Write the updated script back to the file
    with open(script_path, "w") as file:
        file.writelines(updated_lines)

    print(f"Script '{script_path}' updated successfully.")

if __name__ == "__main__":
    main()
