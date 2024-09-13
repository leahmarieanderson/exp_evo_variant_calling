# Exp_evo_variant_calling
The experimental evolution variant calling pipeline is a set of bash and python scripts that process sample data in the form of fastq files into annotated vcf files containing potential variants from one's experiments. 
## Installation
First, ssh onto the GS cluster. 
```php
ssh username@nexus.gs.washington.edu
```
Next, go to your directory which holds the your `fastq` directory.
In this example, this would be the `experiment1` directory
```bash
└── experiment1
    └── fastq
        ├── ancestor_R1_001.fastq.gz
        ├── ancestor_R2_001.fastq.gz
        ├── sample1_R1_001.fastq.gz
        ├── sample1_R2_001.fastq.gz
        ├── sample2_R1_001.fastq.gz
        └── sample2_R2_001.fastq.gz
```
So we would want to use the `cd` command to go to our experiment1 directory.
```php
cd path/to/experiment1
```
***IMPORTANT***  
If you haven't already, make sure to create an `errors` directory, `outputs` directory, and the `genomes` directory. The `errors` and `output` directory can be empty but your `genomes` directory must contain your reference genome and any other additional files that it comes with. 
***DOUBLE IMPORTANT***  
These directory names: `errors`, `outputs`, `genomes`, and `fastq` are hard coded into the scripts so if the directory names are different, the scripts will likely fail. (`experiment1` is not hard coded so this would just be the name of your directory that contains your fastq directory)
```bash
└── experiment1
    ├── errors (ADD THIS)
    ├── fastq
    │   ├── ancestor_R1_001.fastq.gz
    │   ├── ancestor_R2_001.fastq.gz
    │   ├── sample1_R1_001.fastq.gz
    │   ├── sample1_R2_001.fastq.gz
    │   ├── sample2_R1_001.fastq.gz
    │   └── sample2_R2_001.fastq.gz
    ├── genomes (ADD THIS)
    │   ├── sacCer3.fasta (REFERENCE GENOME)
    │   ├── sacCer3.dict
    │   └── ...
    └── outputs (ADD THIS)
```
You can make these directories using the `mkdir` command. Just make sure you are in the the correct directory first (in our example, this would be the `experiment1` directory).
```php
mkdir errors outputs
```
You'll likely want to use the `cp` command to copy the reference genome folder from some other source if you don't have it already.
```php
cp -r path/to/reference/genome/folder . 
```
Your directory tree should now look like what was previously shown above. 

If the Github is set to public, then you can simply go click on the green `<> Code` button on the Github and copy the HTTPS option's link. Then you could go back to your terminal and use `git clone`.
```php
git clone https://github.com/leahmarieanderson/exp_evo_variant_calling.git 
```
If the Github repository is private, You'll likely be your github username and password or personal access token. If that's the case, you can just `cp` from one of the directories on the GS cluster which has the code instead. 
```php
cp -r /net/dunham/vol2/Zilong/updating_pipeline_2024/exp_evo_variant_calling .
```
Now your directory tree should something like this:
```bash
└── experiment1
    ├── errors
    ├── exp_evo_variant_calling
    │   ├── align_ANC.sh
    │   ├── align_samtools.sh
    │   ├── align.sh
    │   ├── annotation_final.py
    │   ├── batch_submit.py
    │   ├── README.md
    │   └── stringent_filter.py
    ├── fastq
    │   ├── ancestor_R1_001.fastq.gz
    │   ├── ancestor_R2_001.fastq.gz
    │   ├── sample1_R1_001.fastq.gz
    │   ├── sample1_R2_001.fastq.gz
    │   ├── sample2_R1_001.fastq.gz
    │   └── sample2_R2_001.fastq.gz
    ├── genomes
    │   ├── sacCer3.fasta
    │   ├── sacCer3.dict
    │   └── ...
    └── outputs
```
## Usage
We have two main steps in our pipeline:
1. Align the ancestor strain and process its bam and vcf files.
2. Align and annotate the evolved strain's variants and put the results into comprehensive output files 

### Align Ancestor strain
First, we need to create the bam files as well as the samtools and freebayes vcfs for our ancestor fastq files.   
So, we need to get on the cluster. We can do so by running this command
```php
ssh grid-head3.gs.washington.edu
```
Next, go into the exp_evo_variant_calling directory
```php
cd path/to/exp_evo_variant_calling
```
Here, you will need to change some of the hard coded variables in the align scripts.  
Use `vim` or `nano` to open `align_ANC.sh` 
```php
nano align_ANC.sh
OR
vim align_ANC.sh
```
After you run this command, your terminal should open up a text editor and show you the code of the bash script.  
you only need to change a few lines of code. 
```php
#!/bin/bash
#$ -S /bin/bash
#$ -wd path/to/working/directory (CHANGE THIS)
#$ -o path/to/outputs (CHANGE THIS)
#$ -e path/to/errors (CHANGE THIS)
#$ -l mfree=8G
#$ -l h_rt=36:0:0
```
First, you'll likely see these lines of code at the very top of your script. These lines set up the job scheduler when you do a submit a job on the cluster later on. You'll want to change the `-wd`, `-o`, and `-e` lines to match your own directory paths. A reminder on how to get your path is to simply `cd` into the directory of your choice, then just type `pwd` into the terminal to get your current path.     
For example, given this directory tree, here's what the code should look like:
```bash
└── experiment1
    ├── errors
    ├── exp_evo_variant_calling
    ├── fastq
    ├── genomes
    └── outputs
-----------------------script below--------------------------------------
#!/bin/bash
#$ -S /bin/bash
#$ -wd path/to/experiment1
#$ -o path/to/experiment1/outputs
#$ -e path/to/experiment1/errors
#$ -l mfree=8G
#$ -l h_rt=36:0:0
```
Next, we want to change some of the file paths that are being used in the script to match your own. So given this directory tree and reference genome:
```bash
└── experiment1
    ├── errors
    ├── exp_evo_variant_calling
    ├── fastq
    ├── genomes
    │   ├── sacCer3.fasta (REFERENCE GENOME)
    │   ├── sacCer3.dict
    │   └── ...
    └── outputs
-----------------------script below--------------------------------------
FOLDER=fastq
ANC=$1
DIR=path/to/experiment1 (CHANGE THIS to be the same as the -wd line from before)
WORKDIR=${DIR}/WorkDirectory # Where files will be created
SEQDIR=${DIR}/${FOLDER} # Location of Fastqs
SEQID=experiment1 # Project name and date for bam header (CHANGE THIS to any name for your project)
REF=${DIR}/genomes/sacCer3.fasta # Reference genome (CHANGE THIS sacCer3.fasta to be your reference genome or keep it if that is your reference genome)
```
Now you can exit from your text editor and we can move on to submitting your job. Remember to save your changes before you exit!  


Next, you would want to submit a job for aligning your ancestor. You can do this by using `qsub` with our `align_ANC.sh` script.  
This is an example fastq directory. Note that `anc_AB_R1_001.fastq.gz` and `anc_AB_R2_001.fastq.gz` are our ancestor fastq files.
```bash
├── fastq
        ├── anc_AB_R1_001.fastq.gz
        ├── anc_AB_R2_001.fastq.gz
        ├── sample1_R1_001.fastq.gz
        ├── sample1_R2_001.fastq.gz
        ├── sample2_R1_001.fastq.gz
        └── sample2_R2_001.fastq.gz
```
We just want to use the prefix before `_R1_001.fastq.gz` as our argument.
So we would want to use this command to submit a job to align our ancestor.
```php
qsub align_ANC.sh anc_AB
```
The script will automatically check the `fastq` directory for a name that matches what we give it and use those fastqs for our bams and vcfs.  
After you finished qsub-ing, you can check on the status of your job by using the `qstat -u username` command in the terminal.  

After the job has finished running, you should have a new directory called `WorkDirectory` and your directory tree should look like:
```bash
└── experiment1
    ├── errors
    ├── exp_evo_variant_calling
    ├── fastq
    ├── genomes
    ├── outputs
    └── WorkDirectory
        └── anc_AB
            ├── dup_metrics
            ├── anc_AB_comb_R1R2.RG.MD.realign.sort.bam
            ├── anc_AB_comb_R1R2.RG.MD.realign.sort.bam.bai
            ├── anc_AB_comb_freebayes_BCBio.vcf
            └── anc_AB_samtools_AB.vcf
```
### Align and Annotate Evolved Sample
Now we want to align and annotate the evolved sample. We have two methods of doing this. 
#### Individual job submissions
Similarly to what we did for our `align_ANC.sh` script, we will need to `nano` or `vim` into an align script of your choice and change the respective bash variables and script variables. 
In our example, we will choose the `align.sh` script. So we will type the command:
```php
vim align.sh
OR
nano align.sh
```
As an example, given this directory tree, we will change our script accordingly.  
```bash
└── experiment1
    ├── errors
    ├── exp_evo_variant_calling
    ├── fastq
        ├── anc_AB_R1_001.fastq.gz
        ├── anc_AB_R2_001.fastq.gz
        ├── sample1_R1_001.fastq.gz (THIS IS OUR EVOLVED STRAIN TO ALIGN)
        ├── sample1_R2_001.fastq.gz (THIS IS OUR EVOLVED STRAIN TO ALIGN)
        ├── sample2_R1_001.fastq.gz
        └── sample2_R2_001.fastq.gz
    ├── genomes
    ├── outputs
    └── WorkDirectory
-----------------------script below--------------------------------------
#!/bin/bash
#$ -S /bin/bash
#$ -wd path/to/experiment1 (CHANGE THIS)
#$ -o path/to/experiment1/outputs (CHANGE THIS)
#$ -e path/to/experiment1/errors (CHANGE THIS)
#$ -l mfree=8G
#$ -l h_rt=36:0:0
#$ -N sample1 (CHANGE THIS)
```
We will only need to change the `-wd`, `o`, `-e`, and `-N` to appropriately named paths and names. (`-N` is simply the name of the job when we submit)  
After this, we will need to change the script variables.
```php
FOLDER=fastq
SAMPLE=$1 # Passed sample prefix (ex: Sample-01)
ANC=$2
DIR=/path/to/experiment1 (CHANGE THIS)
WORKDIR=${DIR}/WorkDirectory # Where files will be created
SEQDIR=${DIR}/${FOLDER} # Location of Fastqs
SEQID=experiment1 # Project name and date for bam header (CHANGE THIS)
REF=${DIR}/genomes/sacCer3.fasta # Reference genome (CHANGE THIS IF NOT YOUR REFERENCE GENOME)
ANNOTATE=${DIR}/genomes # Location of custom annotation scripts
SCRIPTS=${DIR}/exp_evo_variant_calling # Path of annotation_final.py directory
ANCBAM=${WORKDIR}/${ANC}/${ANC}_comb_R1R2.RG.MD.realign.sort.bam
VCFDIR=${WORKDIR}/${ANC}/
```
You'll only need to change the `DIR`, `SEQID`, and possibly `REF` if you have a different reference genome.  
Lastly, just like how we did a job submission with `align_ANC.sh`, we will do the same thing with our `align.sh` script.  
This command is a little different than our align_ANC.sh qsub but with one additional argument. It's the same as before with the format  
`qsub align_script evolved_name ancestor_name`
So with our example from before: 
```bash
├── fastq
        ├── anc_AB_R1_001.fastq.gz (Ancestor name: anc_AB)
        ├── anc_AB_R2_001.fastq.gz
        ├── sample1_R1_001.fastq.gz (Evolved strain name: sample1)
        ├── sample1_R2_001.fastq.gz 
        ├── sample2_R1_001.fastq.gz
        └── sample2_R2_001.fastq.gz
```
Our command would be:
```php
qsub align.sh sample1 anc_AB
```
After that, you can check on your job status by using `qstat -u username`.
### Batch submit
You can submit multiple samples for alignment and annotating as long as they all come from the same ancestor and  
they have both R1 and R2 fastqs in one directory.  

First, cd into your `exp_evo_variant_calling` directory
```php
cd path/to/exp_evo_variant_calling
```
Next, you'll want to run the `batch_submit` python script. (on default, it will call the `align.sh` script)
```php
python3 batch_submit.py
```
After running the program, you will be prompted with various lines that you can edit in the `align.sh` script.   
For example: 
```
Current Bash Settings:
[0] #$ -S /bin/bash

[1] #$ -wd /net/dunham/vol2/Zilong/updating_pipeline_2024

[2] #$ -o /net/dunham/vol2/Zilong/updating_pipeline_2024/outputs/

[3] #$ -e /net/dunham/vol2/Zilong/updating_pipeline_2024/errors/

[4] #$ -l mfree=8G

[5] #$ -l h_rt=36:0:0

Enter the number of the line to change (type 'c' to finish changes):
```
You can then enter the line number of the bash setting you want to change or just type "c" to continue if you don't want to make changes.  
After continuing, you can also change your script variables.  
```
Current Script Variables:
[0] FOLDER=fastq

[1] DIR=/net/dunham/vol2/Zilong/updating_pipeline_2024

[2] SEQID=experiment1 # Project name and date for bam header

[3] REF=${DIR}/genomes/sacCer3.fasta # Reference genome

[4] SCRIPTS=${DIR}/exp_evo_variant_calling # Path of annotation_final.py directory

Enter the number of the line to change (type 'c' to finish changes): 
```
Enter the line number of the variable you want to change or just type "c" to continue.   
NOTE: If you continue past this point, modifications will be finalized so if you accidentally  
overwrote a line you didn't intend to overwrite, just use `ctr + c` to interrupt your session and exit without saving changes.

After finalizing script modifications, you will be shown these prompts:
```
Script 'align.sh' updated.
Would you like to qsub all samples in path/to/FOLDER ? (y/n) : 
```
Double check that all samples in your folder have 2 fastq files and that they all use the same ancestor   
and that the `path/to/FOLDER` is correct. Enter "y" to continue.  
Example:
```
What is the name of the ancestor? : anc_AB
```
Enter the name of the ancestor (this should be the ancestor of all samples in your `FOLDER`) and this should start submitting jobs for all samples in the `FOLDER`.

## Align.sh Outputs
After your `align.sh` jobs have been completed, you will have a few new directories and files in your `WorkDirectory`.  
In this example, the ancestor is `anc_AB` and the sample that was submitted for alignment and annotation is `sample1`.
```bash
└── WorkDirectory
        ├── anc_AB
        └── sample1  *NEW*
            ├── bams
            ├── dup_metrics
            ├── intermediate_vcfs
            └── results
```
The directory that we are most interested in is the `results` directory. 
```
└── results
        ├── sample1_final_stringent_compiled.csv       
        ├── sample1_freebayes_BCBio_AncFiltered_annotated_vcf.txt
        ├── sample1_freebayes_BCBio_AncFiltered_condensed.csv
        ├── sample1_lofreq_AncFiltered_annotated_vcf.txt
        ├── sample1_lofreq_AncFiltered_condensed.csv
        ├── sample1_samtools_AB_AncFiltered_annotated_vcf.txt
        └── sample1_samtools_AB_AncFiltered_condensed.csv
```

- The `freebayes_BCBio_AncFiltered_annotated_vcf.txt`, `lofreq_AncFiltered_annotated_vcf.txt`, and `samtools_AB_AncFiltered_annotated_vcf.txt` are annotated files of each variant caller which has all the ancestor mutations already filtered out. This means that these files contains only the variants that were found over the course of your experiment.

- The `freebayes_BCBio_AncFiltered_condensed.csv`, `lofreq_AncFiltered_condensed.csv`, and `samtools_AB_AncFiltered_condensed.csv` are the same files as the ones mentioned above but we applyed a unique set of filter conditions for each file based on their specific variant caller and we condense the columns. The `_condensed.csv` files only show the columns that are relevant for our analysis. Such columns like `CHROM`, `POS`, `REF`, `ALT`, `ANNOTATION`, `REGION`, and `PROTEIN`. For the filter, a file's particular variant caller changes the filter conditions for `QUAL`,`DP`, and number of reads on the ref and alt alleles. We would keep any variants that pass the specified threshold for `QUAL`, `DP`, etc. 

- For example, we set our Samtools filter to have a default `QUAL` threshold of 75, so any variants under 75 for the `QUAL` would not make it into the `samtools_AB_AncFiltered_condensed.csv`. Our Freebayes filter on the other hand has a `QUAL` threshold of 20.

- The `final_stringent_compiled.csv` file is a combination of the `freebayes_BCBio_AncFiltered_condensed.csv`, `lofreq_AncFiltered_condense.csv`, and `samtools_AB_AncFiltered_condensed.csv` that has been sorted, removed duplicates, and added an additional column `NUM_OCCURANCES` that counts the number of times this variant has shown between the different variant callers. The higher this value, the more reliable this variant is since it means it was called by more variant callers. 

Below is a pipeline of how these file output files are generated. 
```
                                                (filter based on thresholds, keep only important columns)             (Append all vcfs together)
freebayes_BCBio_AncFiltered_annotated_vcf.txt ─────────[filter]───freebayes_BCBio_AncFiltered_condensed.csv────┐
                                                                                                               │ 
lofreq_AncFiltered_annotated_vcf.txt ──────────────────[filter]─────lofreq_AncFiltered_condensed.csv───────────┼───>  final_stringent_compiled.csv
                                                                                                               │
samtools_AB_AncFiltered_annotated_vcf.txt ─────────────[filter]─────samtools_AB_AncFiltered_condensed.csv──────┘
```

The other directory `intermediate_vcfs` contains the vcfs that were used to create the files in the `results` directory, some of these files are vcfs that have not been run through our annotation script and thus are without the `ANNOTATIONS`, `REGION`, and `PROTEIN` column, but can provide value if you wanted to look at the raw outputs of samtools, freebayes, and lofreq. 