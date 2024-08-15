# Exp_evo_variant_calling
The experimental evolution variant calling pipeline is a set of bash and python scripts that process sample data in the form of fastq files into annotated vcf files containing potential variants from one's experiments. 
## Installation
First, ssh onto the GS cluster. 
```php
ssh username@nexus.gs.washington.edu
```
Next, go to your directory which holds the your fastq directory.
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
    │   └── README.md
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
-----------------------------------------------------------------
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
-----------------------------------------------------------------
FOLDER=fastq
ANC=$1
DIR=path/to/experiment1 (CHANGE THIS to be the same as the -wd line from before)
WORKDIR=${DIR}/WorkDirectory # Where files will be created
SEQDIR=${DIR}/${FOLDER} # Location of Fastqs
SEQID=experiment1 # Project name and date for bam header (CHANGE THIS to any name for your project)
REF=${DIR}/genomes/sacCer3.fasta # Reference genome (CHANGE THIS sacCer3.fasta to be your reference genome or keep it if that is your reference genome)
```
Now you can exit from your text editor and we can move on to submitting your job. Remember to save your changes before you exit!  


Next, you would want to submit a job for aligning your ancestor. You can do this by using `qsub` with our align_ANC.sh script.  
This is an example fastq directory. Note that `anc_AB_R1_001.fastq.gz` and `anc_AB_R2_001.fastq.gz` are our ancestor fastq files.
```bash
├── fastq
    │   ├── anc_AB_R1_001.fastq.gz
    │   ├── anc_AB_R2_001.fastq.gz
    │   ├── sample1_R1_001.fastq.gz
    │   ├── sample1_R2_001.fastq.gz
    │   ├── sample2_R1_001.fastq.gz
    │   └── sample2_R2_001.fastq.gz
```
We just want to use the prefix before `_R1_001.fastq.gz` as our argument.
So we would want to use this command to submit a job to align our ancestor.
```php
qsub align_ANC.sh anc_AB
```
The script will automatically check the `fastq` directory for a name that matches what we give it and use those fastqs for our bams and vcfs.



