import os
import glob
import subprocess
import sys

exp_evo_path = os.path.dirname(os.path.realpath(__file__)) # this is our current directory

work_dir_path = os.path.dirname(exp_evo_path) # this is the work directory we will put for the bash -wd directive 

def find_second_last_underscore(s):
    last_index = s.find('_R1_001')
    return last_index

def submit_job(script_name, sample_name, ancestor_name):
    # update script to set job name to sample_name
    with open(script_name, "r") as file:
        lines = file.readlines()

    # Print job submission details
    print(f"Submitting job for sample: {sample_name} with ancestor: {ancestor_name}")
    
    # Run the command to submit the job
    subprocess.run(['qsub', '-N', "MD" + sample_name, script_name, sample_name, ancestor_name])

def get_sample_names(directory):
    # Get all R1 sample files
    R1files = glob.glob(os.path.join(directory, '*R1*'))
    # Get all R2 sample files
    R2files = glob.glob(os.path.join(directory, '*R2*'))
    
    # Ensure both lists have the same number of files
    if len(R1files) != len(R2files):
        raise ValueError("The number of R1 files does not match the number of R2 files.")
    
    # Sort both lists
    R1files.sort()
    R2files.sort()
    
    # Ensure that each R1 file has a corresponding R2 file with the same prefix
    for R1, R2 in zip(R1files, R2files):
        R1_prefix = os.path.basename(R1).split('R1')[0]
        R2_prefix = os.path.basename(R2).split('R2')[0]
        
        if R1_prefix != R2_prefix:
            raise ValueError(f"Prefix mismatch between files: {R1} and {R2}")

    sample_name_list = []
    for file_path in R1files:
        # Extract the sample name from the file name
        file_name = os.path.basename(file_path)
        sample_name_end_index = find_second_last_underscore(file_name)
        sample_name = file_name[:sample_name_end_index]
        sample_name_list.append(sample_name)

    return sample_name_list

def multi_qsub(script_name, directory, ancestor_name):
    # directory should be the path to your fastq directory
    samples = get_sample_names(directory)
    
    # check if the ancestor name is even in the directory
    if ancestor_name not in samples:
        print(ancestor_name, " not found in ", directory, ". Exiting Program...")
        sys.exit(1)

    # Submit a job for each file
    for sample_name in samples:
        if sample_name != ancestor_name: # only qsub for with the samples, not ancestor
            submit_job(script_name, sample_name, ancestor_name)

def print_script_var(script_vars):
    # Print all script_variables
    print("\nCurrent Script Variables:")
    for var in script_vars:
        print(var)

def remove_last_word(s):
    # Strip the string of leading and trailing whitespace
    s = s.strip()
    
    # Find the position of the last whitespace
    last_whitespace_index = s.rfind(' ')
    
    # If there's no whitespace, return an empty string (i.e., the string itself is a single word)
    if last_whitespace_index == -1:
        return ''
    
    # Return the string up to the last whitespace
    return s[:last_whitespace_index + 1]

def remove_last_var(s):
    # Strip the string of leading and trailing whitespace
    s = s.strip()
    
    # Find the position of the first equals
    index = s.find('=')
    
    # If there's no equals, return an empty string (i.e., the string itself is a single word)
    if index == -1:
        return ''
    
    # Return the string up to the last whitespace
    return s[:index + 1]
    
def find_index_of_substring(strings, substring):
    for index, s in enumerate(strings):
        if substring in s:
            return index
    return -1  # Return -1 if the substring is not found

def main():

    # CHANGE THIS TO YOUR ALIGN SCRIPT
    script_name = 'align.sh'
    
    # Read the existing script
    with open(script_name, "r") as script:
        script_lines = script.readlines()

    # Initialize lists and counters
    auto_gen_bash_settings = []
    auto_gen_script_vars = []
    bash_settings = []
    script_variables = []

    # Iterate through script_lines and categorize lines
    for line in script_lines:
        # form an array of the currenty job directives in the script
        if line.startswith("#$ -"):
            bash_settings.append(line)
        # form an array of the current script variables 
        elif line.startswith("FOLDER="):
            script_variables.append(line.strip())
        elif line.startswith("DIR="):
            script_variables.append(line.strip())
        elif line.startswith("SEQID="):
            script_variables.append(line.strip())
        elif line.startswith("REF="):
            script_variables.append(line.strip())
        elif line.startswith("SCRIPTS="):
            script_variables.append(line.strip())
    
    # auto generate the bash settings (SGE Directives) using what we know about our own file structure.
    auto_gen_bash_settings.append("#$ -S /bin/bash\n")
    auto_gen_bash_settings.append(f"#$ -wd {work_dir_path}\n")
    auto_gen_bash_settings.append(f"#$ -o {work_dir_path}/outputs/\n")
    auto_gen_bash_settings.append(f"#$ -e {work_dir_path}/errors/\n")
    auto_gen_bash_settings.append("#$ -l mfree=8G\n")
    auto_gen_bash_settings.append("#$ -l h_rt=36:0:0\n")

    # check current bash settings to our auto generated bash settings
    # if different, then prompt users to keep old settings or auto generated one. 
    if set(bash_settings) != set(auto_gen_bash_settings):
        # print current bash settings in script
        print("\nCurrent Bash Settings for scripts:")
        for setting in bash_settings:
            print(setting)
        # print new bash settings that use your directories
        print("\n**NEW** Bash Settings for scripts:")
        for setting in auto_gen_bash_settings:
            print(setting)
        # prompt user to either use the new settings that are auto generated
        while True:
            user_input = input("Change your SGE Directives to fit your directories? (y/n) : ")
            if user_input.lower() == 'y':
                bash_settings = auto_gen_bash_settings
                break
            elif user_input.lower() == 'n':
                # just keep our current SGE directives and continue
                break
            else: 
                print("This is not a valid answer")
    
    # else the SGE Directives are the same so we just ignore and continue forward     

    print_script_var(script_variables)

    auto_gen_script_vars.append("FOLDER=fastq")
    auto_gen_script_vars.append(f"DIR={work_dir_path}")
    # Prompt user for what they want the SEQID to be
    user_input_SEQID = input("What would you like your SEQID to be? (this could be your project name): ")
    auto_gen_script_vars.append(f"SEQID={user_input_SEQID} # Project name and date for bam header")
    auto_gen_script_vars.append("REF=${DIR}/exp_evo_variant_calling/genomes/sacCer3.fasta # Reference genome")
    auto_gen_script_vars.append("SCRIPTS=${DIR}/exp_evo_variant_calling # Path of annotation_final.py directory")

    if set(script_variables) != set(auto_gen_script_vars):
        # print current bash settings in script
        print("\nCurrent Variables Settings for scripts:")
        for variable in script_variables:
            print(variable)
        # print new bash settings that use your directories
        print("\n**NEW** Variables Settings for scripts:")
        for variable in auto_gen_script_vars:
            print(variable)
        # prompt user to either use the new settings that are auto generated
        while True:
            user_input = input("Change your variable paths to fit your directories? (y/n) : ")
            if user_input.lower() == 'y':
                script_variables = auto_gen_script_vars
                break
            elif user_input.lower() == 'n':
                # just keep our current script variables and continue
                break
            else: 
                print("This is not a valid answer")



    # get index of each options from the arrays bash_settings and script variables
    shell_index = find_index_of_substring(bash_settings, ('#$ -S'))
    wd_index = find_index_of_substring(bash_settings, ('#$ -wd'))
    out_index = find_index_of_substring(bash_settings, ('#$ -o'))
    err_index = find_index_of_substring(bash_settings, ('#$ -e'))
    mfree_index = find_index_of_substring(bash_settings, ('#$ -l mfree'))
    h_rt_index = find_index_of_substring(bash_settings, ('#$ -l h_rt'))
    FOLDER_index = find_index_of_substring(script_variables, ('FOLDER='))
    DIR_index = find_index_of_substring(script_variables, ('DIR='))
    SEQID_index = find_index_of_substring(script_variables, ('SEQID='))
    REF_index = find_index_of_substring(script_variables, ('REF='))
    SCRIPTS_index = find_index_of_substring(script_variables, ('SCRIPTS='))

    # Get options to respective strings and use .strip() to remove leading and trailing whitespace
    shell_option = bash_settings[shell_index].strip()
    wd_option = bash_settings[wd_index].strip()
    out_option = bash_settings[out_index].strip()
    err_option = bash_settings[err_index].strip()
    mfree_option = bash_settings[mfree_index].strip()
    h_rt_option = bash_settings[h_rt_index].strip()
    FOLDER_option = script_variables[FOLDER_index].strip()
    DIR_option = script_variables[DIR_index].strip()
    SEQID_option = script_variables[SEQID_index].strip()
    REF_option = script_variables[REF_index].strip()
    SCRIPTS_option = script_variables[SCRIPTS_index].strip()

    # open align script to read through all lines
    with open(script_name, "r") as file:
        lines = file.readlines()
    
    # Update specific lines
    updated_lines = []
    for line in lines:
        if line.startswith("#$ -S"):
            updated_lines.append(f"{shell_option}\n")
        elif line.startswith("#$ -wd"):
            updated_lines.append(f"{wd_option}\n")
        elif line.startswith("#$ -o"):
            updated_lines.append(f"{out_option}\n")
        elif line.startswith("#$ -e"):
            updated_lines.append(f"{err_option}\n")
        elif line.startswith("#$ -l mfree"):
            updated_lines.append(f"{mfree_option}\n")
        elif line.startswith("#$ -l h_rt"):
            updated_lines.append(f"{h_rt_option}\n")
        elif line.startswith("FOLDER="):
            updated_lines.append(f"{FOLDER_option}\n")
        elif line.startswith("DIR="):
            updated_lines.append(f"{DIR_option}\n")
        elif line.startswith("SEQID="):
            updated_lines.append(f"{SEQID_option}\n")
        elif line.startswith("REF="):
            updated_lines.append(f"{REF_option}\n")
        elif line.startswith("SCRIPTS="):
            updated_lines.append(f"{SCRIPTS_option}\n")
        else:
            updated_lines.append(line)

    # Write the updated script back to the file
    with open(script_name, "w") as file:
        file.writelines(updated_lines)

    print(f"Script '{script_name}' updated.")

    fastq_dir = DIR_option[4:] + "/" + FOLDER_option[7:] + "/"
    # prompts user if they want to qsub all samples(except ancestor) in their ${DIR}/${FOLDER} i.e. your fastq folder
    user_submit = input(f"Would you like to qsub all samples in {fastq_dir} ? (y/n) : ")
    if user_submit.lower() == 'y':
        ancestor_name = input("What is the name of the ancestor? : ")
        # check if ancestor directory even exists in WorkDirectory

        if os.path.isdir(DIR_option[4:] + "/WorkDirectory"):
            multi_qsub(script_name, fastq_dir, ancestor_name)
        else:
            # Ancestor directory does not exist in WorkDirectory, prompt user to "qsub align.sh ancestor" first before batch submitting
            print(f"{ancestor_name} not found in {DIR_option[4:]}/WorkDirectory, cannot batch submit samples without relevant ancestor files \n")
            print(f"Align the ancestor files first with the command 'qsub align.sh {ancestor_name}' then use batch submit python script after that job is completed")

    else: 
        print("No qsub has been selected. Script completed successfully.")

if __name__ == '__main__':
    main()
