def print_bash_settings(bash_settings):
    # Print all bash_settings
    print("\nCurrent Bash Settings:")
    for setting in bash_settings:
        print(setting)

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

def remove_line_marker(arr):
    new_arr = []
    for s in arr:

        # Strip the string of leading and trailing whitespace
        s = s.strip()
        
        # Find the position of the first right bracket
        index = s.find(']')

        comment_index = s.find('#')

        if index == -1:
            raise RuntimeError("unable to find line marker to remove. Crashing Program.")
        
        if comment_index > 5: #removes real comment, reminder that the string looks like '[1] #$ -S /bin/bash' so we have to check after the 5th character
            new_arr.append(s[index + 2:comment_index])
        else: # there are no comments or there are '#' at the beginning like with bash code
            new_arr.append(s[index + 2:])
    return new_arr
    
def find_index_of_substring(strings, substring):
    for index, s in enumerate(strings):
        if substring in s:
            return index
    return -1  # Return -1 if the substring is not found

def main():
    script_name = 'align.sh'
    
    # Read the existing script
    with open(script_name, "r") as script:
        script_lines = script.readlines()

    # Initialize lists and counters
    bash_settings = []
    script_variables = []
    bash_line_num = 0
    script_var_line_num = 0

    # Iterate through script_lines and categorize lines
    for line in script_lines:
        if line.startswith("#$ -"):
            bash_settings.append(f"[{bash_line_num}] {line}")
            bash_line_num += 1
        elif line.startswith("FOLDER="):
            script_variables.append(f"[{script_var_line_num}] {line}")
            script_var_line_num += 1
        elif line.startswith("DIR="):
            script_variables.append(f"[{script_var_line_num}] {line}")
            script_var_line_num += 1
        elif line.startswith("SEQID="):
            script_variables.append(f"[{script_var_line_num}] {line}")
            script_var_line_num += 1
        elif line.startswith("REF="):
            script_variables.append(f"[{script_var_line_num}] {line}")
            script_var_line_num += 1
        elif line.startswith("SCRIPTS="):
            script_variables.append(f"[{script_var_line_num}] {line}")
            script_var_line_num += 1
    
    # allow user to change bash settings
    print_bash_settings(bash_settings)
    while True:
        user_input = input("Enter the number of the line to change (type 'c' to finish changes): ")
        if user_input.lower() == 'c':
            print_bash_settings(bash_settings)
            break
        elif user_input.isdigit() and int(user_input) < bash_line_num:
            user_changes = input("What would you like to change this to? : ").replace(" ","") # get input and remove whitespaces from input
            if(bash_settings[int(user_input)].find('=') == -1): # no equals character so we replace word after last whitespace
                bash_settings[int(user_input)] = remove_last_word(bash_settings[int(user_input)]) + user_changes 
            else: # there is an equals so we have to replace word after equals 
                bash_settings[int(user_input)] = remove_last_var(bash_settings[int(user_input)]) + user_changes 
            print_bash_settings(bash_settings)
        else:
            print("This is not a valid number")

    print_script_var(script_variables)
    while True:
        user_input = input("Enter the number of the line to change (type 'c' to finish changes): ")
        if user_input.lower() == 'c':
            print_script_var(script_variables)
            break
        elif user_input.isdigit() and int(user_input) < script_var_line_num: # script variables always have equals character
            user_changes = input("What would you like to change this to? : ").replace(" ","") # get input and remove whitespaces from input
            
            script_variables[int(user_input)] = remove_last_var(script_variables[int(user_input)]) + user_changes 
            print_script_var(script_variables)
        else:
            print("This is not a valid number")

    bash_settings = remove_line_marker(bash_settings)
    script_variables = remove_line_marker(script_variables)

    # get index of each options from the arrays bash_settings and script variables
    shell_index = find_index_of_substring(bash_settings, ('#$ -S'))
    wd_index = find_index_of_substring(bash_settings, ('#$ -wd'))
    out_index = find_index_of_substring(bash_settings, ('#$ -o'))
    err_index = find_index_of_substring(bash_settings, ('#$ -e'))
    mfree_index = find_index_of_substring(bash_settings, ('#$ -l mfree'))
    h_rt_index = find_index_of_substring(bash_settings, ('#$ -l h_rt'))
    jobname_index = find_index_of_substring(bash_settings, ('#$ -N'))
    FOLDER_index = find_index_of_substring(script_variables, ('FOLDER='))
    DIR_index = find_index_of_substring(script_variables, ('DIR='))
    SEQID_index = find_index_of_substring(script_variables, ('SEQID='))
    REF_index = find_index_of_substring(script_variables, ('REF='))
    SCRIPTS_index = find_index_of_substring(script_variables, ('SCRIPTS='))

    # Get options to respective strings
    shell_option = bash_settings[shell_index]
    wd_option = bash_settings[wd_index]
    out_option = bash_settings[out_index]
    err_option = bash_settings[err_index]
    mfree_option = bash_settings[mfree_index]
    h_rt_option = bash_settings[h_rt_index]
    jobname_option = bash_settings[jobname_index]
    FOLDER_option = script_variables[FOLDER_index]
    DIR_option = script_variables[DIR_index]
    SEQID_option = script_variables[SEQID_index]
    REF_option = script_variables[REF_index]
    SCRIPTS_option = script_variables[SCRIPTS_index]

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
        elif line.startswith("#$ -N"):
            updated_lines.append(f"{jobname_option}\n")
        elif line.startswith("FOLDER="):
            updated_lines.append(f"{FOLDER_option}\n")
        elif line.startswith("DIR="):
            updated_lines.append(f"{DIR_option}\n")
        elif line.startswith("SEQID="):
            updated_lines.append(f"{SEQID_option} # Project name and date for bam header\n")
        elif line.startswith("REF="):
            updated_lines.append(f"{REF_option} # Reference genome\n")
        elif line.startswith("SCRIPTS="):
            updated_lines.append(f"{SCRIPTS_option} # Path of annotation_final.py directory\n")
        else:
            updated_lines.append(line)

    # Write the updated script back to the file
    with open('align.sh', "w") as file:
        file.writelines(updated_lines)

    print("Script 'align.sh' updated successfully.")

if __name__ == "__main__":
    main()
