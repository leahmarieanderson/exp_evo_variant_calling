# exp_evo_variant_calling
The experimental evolution variant calling pipeline is a set of bash and python scripts that process sample data in the form of fastq files into annotated vcf files containing potential variants from one's experiments. 

*Authors: Zilong Zeng and Leah Anderson, 2026*

## Installation
First, log onto the GS cluster: 
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
If the Github is set to public, then you can simply go click on the green `<> Code` button on the Github and copy the HTTPS option's link. Then you can go back to your terminal and use `git clone`.
```php
git clone https://github.com/leahmarieanderson/exp_evo_variant_calling.git 
```
If you don't have a personal access token to GitHub, you can just 'cp' from one of the directories on the GS cluster which has the code instead. 
```php
cp -r /net/dunham/vol2/Zilong/updating_pipeline_2024/exp_evo_variant_calling .
```
Now your directory tree should something like this:
```bash
└── experiment1
    ├── exp_evo_variant_calling
        ├── align.sh
        ├── annotation_final.py
        ├── batch_submit.py
        ├── README.md
        ├── stringent_filter.py
        ├──  ...
        ├──  yEvo_allcolor_ancestor
            └── master_ancestor_freebayes_BCBio.vcf
            └── master_ancestor_gatk_haplo.vcf
            └── master_ancestor_lofreq.vcf
        └── genomes
           └── sacCer3.fasta
           └── sacCer3.dict
           └── ...
    ├── fastq
        ├── ancestor_R1_001.fastq.gz
        ├── ancestor_R2_001.fastq.gz
        ├── sample1_R1_001.fastq.gz
        ├── sample1_R2_001.fastq.gz
        ├── sample2_R1_001.fastq.gz
        └── sample2_R2_001.fastq.gz

```

## Alignment Script Setup
Now we want to change the our alignment script to match our unique personal directory. The bash settings of the script upon initial installation will likely not be applicable for your directory. 

First, go into the `exp_evo_variant_calling` directory so we can access our scripts.
```
cd exp_evo_variant_calling
```

Next, run the `batch_submit.py` and follow the steps below. 

```
python3 batch_submit.py
```

It will then prompt you to modify the current bash settings to bash settings that would fit your local path.

```
Current Bash Settings for scripts:
#$ -S /bin/bash

#$ -wd /net/dunham/vol2/Leah/yEvo_echinocandins/yEvo_sequencing_260416

#$ -o /net/dunham/vol2/Leah/yEvo_echinocandins/yEvo_sequencing_260416/outputs/

#$ -e /net/dunham/vol2/Leah/yEvo_echinocandins/yEvo_sequencing_260416/errors/

#$ -l mfree=8G

#$ -l h_rt=36:0:0


**NEW** Bash Settings for scripts:
#$ -S /bin/bash

#$ -wd /net/dunham/vol2/user/experiment1

#$ -o /net/dunham/vol2/user/experiment1/outputs/

#$ -e /net/dunham/vol2/user/experiment1/errors/

#$ -l mfree=8G

#$ -l h_rt=36:0:0

Change your SGE Directives to fit your directories? (y/n) : 
```

Enter `y` to continue

```
Current Script Variables:
FOLDER=fastq
DIR=/net/dunham/vol2/Leah/yEvo_echinocandins/yEvo_sequencing_260416
SCRIPTS=${DIR}/exp_evo_variant_calling # Path of annotation_final.py directory
SEQID=april2026seq # Project name and date for bam header
REF=${DIR}/exp_evo_variant_calling/genomes/sacCer3.fasta # Reference genome
What would you like your SEQID to be? (this could be your project name): 
```

Here, you just enter the project name of your experiment, this is purely for the bam header. Every variable will be updated at the end of running the script. 

```
Current Variables Settings for scripts:
FOLDER=fastq
DIR=/net/dunham/vol2/Leah/yEvo_echinocandins/yEvo_sequencing_260416
SCRIPTS=${DIR}/exp_evo_variant_calling # Path of annotation_final.py directory
SEQID=april2026seq # Project name and date for bam header
REF=${DIR}/exp_evo_variant_calling/genomes/sacCer3.fasta # Reference genome

**NEW** Variables Settings for scripts:
FOLDER=fastq
DIR=/net/dunham/vol2/user/experiment1
SEQID=test # Project name and date for bam header
REF=${DIR}/exp_evo_variant_calling/genomes/sacCer3.fasta # Reference genome
SCRIPTS=${DIR}/exp_evo_variant_calling # Path of annotation_final.py directory
Change your variable paths to fit your directories? (y/n) 
```

Double check if the `FOLDER` and `DIR` variable are valid file paths for your experiment and that the `FOLDER` variable is the name of the directory which contains the fastqs that you want to call variants on. If they aren't, you can just manually change them in the script later. Enter `y` to continue.

```
Script 'align.sh' updated.
Would you like to qsub all samples in /net/dunham/vol2/user/experiment1/fastq/ ? (y/n) : 
```
Enter `n` here as we are just trying to set up your bash settings and script variables.

## Usage
We have two ways to run our pipeline depending on the type of your experiment:

### Yeast Evolution experiment (yEvo)
yEvo experiments will typically use one of the following ancestors for their experiment:

* YMD4612_Pink_Ancestor_S46
* YMD4613_Purple_Ancestor_S47
* YMD4614_Blue_Ancestor_S48
* YMD4615_Yellow_Ancestor_S49
* YMD4616_Black_Ancestor_S50
* YMD4617_Gray_Ancestor_S51
* YMD4618_Orange_Ancestor_S52
* YMD4619_White_Ancestor_S53

Rather than having to run the pipeline on each of these ancestors for your corresponding sample, the `yEvo_allcolor_ancestor` directory contains a combination of variants from all 8 color ancestors for each variant caller. 

For example, `master_ancestor_freebayes_BCBio.vcf` is the combination of all unique variants found in all 8 of the color ancestors using the freebayes variant caller. The `master_ancestor_gatk_haplo.vcf` is the combination of all unique variants found in all 8 of the color ancestors using the gatk haplo variant caller. Lastly, the `master_ancestor_lofreq.vcf` is the combination of all unique variants found in all 8 of the color ancestors using the lofreq variant caller. 

You can either run the pipeline on one sample or do a batch submit on all the samples in a given directory (this directory needs to contain the fastq files of your samples and nothing else)

### Single Sample Job Submission
For yEvo jobs, you would use this command:

*Note: do not include the whole fastq file name, just the prefix. For example: if your forward read file is "MD4406_S2_R1_001.fastq.gz", then your sample name is "MD4406_S2".*

```
qsub align_script.sh <SAMPLE> --yevo
```

### Batch Job Submission
To submit multiple yEvo jobs, we will utilize the `batch_submit.py` script.
Run the following command:
```
python3 batch_submit.py
```
Go through the steps like how you did in the Alignment Script Setup section and on the last question, enter `y`.
```
Would you like to qsub all samples in /net/dunham/vol2/user/experiment1/fastq/ ? (y/n) : y
```
Then you will be prompted if this is a yevo experiment, enter `y`.
```
Is this a yevo experiment? (y/n) : y
```
And you should see an output of multiple job submissions being qsubbed. 

## Non-yEvo experiments
For non-yevo experiments, the workflow looks like this. 

1. Align the ancestor strain and process its bam and vcf files.
2. Align the evolved strain's variants, apply filters based on the results from aligning ancestor strain and put the results into comprehensive output files 

### Align Ancestor strain
First, we need to create the bam files as well as the samtools and freebayes vcfs for our ancestor fastq files.   

Start by going into the exp_evo_variant_calling directory that you have cloned or copied:
```php
cd path/to/exp_evo_variant_calling
```

#### *IMPORTANT*: make sure that you have your ancestor fastqs in the same directory as what you put for the FOLDER variable in the align.sh script

Next, align your ancestor by submitting a job to the cluster with `qsub`:

```
$ qsub align.sh my_ancestor_strain
```
*Note: do not include the whole fastq file name, just the prefix. For example: if your forward read file is "MD4406_S2_R1_001.fastq.gz", then your sample name is "MD4406_S2".*

This may take up to 2 hours to complete, depending on the size of your fastq. You need to wait for the entire job to complete before moving onto the next step.
After you finished qsub-ing, you can check on the status of your job by using the `qstat -u username` command in the terminal.  

*Note: if you are running this pipeline on yEvo samples, the ancestor file we typically use is YMD4612_pink_S1*

After the job has finished running, you should have a new directory called `WorkDirectory` and your directory tree should look like:
```bash
└── experiment1
    ├── errors
    ├── exp_evo_variant_calling
    ├── fastq
    ├── genomes
    ├── outputs
    ├──  yEvo_allcolor_ancestor
    └── WorkDirectory
        └── my_ancestor_strain
            ├── dup_metrics
            ├── anc_AB_freebayes_BCBio.vcf
            ├── anc_AB_gatk_haplo.vcf
            ├── anc_AB_gatk_haplo.vcf.idx
            ├── anc_AB_R1R2_MD.sort.bam
            ├── anc_AB_R1R2_MD.sort.bam.bai
            ├── anc_AB_R1R2_MD.sort.bam.sbi
            ├── anc_AB_S1_R1_001_fastqc.html
            └── anc_AB_S1_R2_001_fastqc.html
```

### Ancestor Output Explanation

1. `dup_metrics` is a directory which contains the duplicate metrics of the ancestor strain. It may be interesting to look at the `READ_PAIR_DUPLICATES` column to see how many read pairs were flagged as duplicates, or look at the `PERCENTAGE_DUPLICATION` column to see what percent of the reads are duplicates.
2. `anc_AB_freebayes_BCBio.vcf` is the output of one of the 2 variant callers that we call on the ancestor (FreeBayes).
3. `anc_AB_gatk_haplo.vcf` is the output of one of the 2 variant callers that we call on the ancestor (FreeBayes).
4. `anc_AB_gatk_haplo.vcf.idx` is an index file for the GATK HaplotypeCaller output.
5. `anc_AB_R1R2_MD.sort.bam`, `anc_AB_R1R2_MD.sort.bam.bai`, and `anc_AB_R1R2_MD.sort.bam.sbi` are the typical `.bam` file outputs that we get.
6. `anc_AB_S1_R1_001_fastqc.html` and `anc_AB_S1_R2_001_fastqc.html` are Fastq Quality Control Check outputs which can tell us analytical statistics on each of the short reads. 

### Align and Annotate Evolved Sample
Now we want to align and annotate the evolved samples. We have two methods of doing this. 
#### Method 1: Individual job submissions

If you only have 1 or 2 samples to run, you can just submit the qsub jobs individually this way: 

```php
qsub -N jobname align.sh sample1 my_ancestor_strain
```

#### Remember that you only need to use the prefix, not the full file name.

#### Example:

For evolved sample reads `evolved_AB_R1_001.fastq.gz` and ancestor reads `anc_AB_R1_001.fastq.gz`, the command would look like:

```
qsub -N jobname align.sh evolved_AB anc_AB
```

### Method 2: Batch submit
You can submit multiple samples for alignment and annotating as long as:
1. The samples all come from the same ancestor
2. The samples have both R1 and R2 fastqs in the directory that is listed in the `SEQDIR` variable in the `align.sh` script.

You can do this by submitting the `batch_submit` python script again.

```php
python3 batch_submit.py
```
Follow the prompts as described above. If you've already run the script on the ancestor, you shouldn't have to rename your working directory paths again. (So in that case, you can enter "y" or "n" for changing path names, it won't matter.)

Finally, when the program asks if you'd like to qsub all fastq's in the working directory, enter "y".

Then you will see this:
```
What is the name of your ancestor?
```
Enter your ancestor name. If the ancestor and all necessary files are not present in the `WorkingDirectory` are not present, the program will exit.

If everything works successfully, each sample in your fastq folder will be submitted to the cluster as a qsub for the `align.sh` script.
You can check on your job status by using `qstat -u username`.


## Align.sh Outputs
After your `align.sh` jobs have been completed, you will have a few new directories and files in your `WorkDirectory`.  
In this example, the ancestor is `anc_AB` and the sample that was submitted for alignment and annotation is `sample1`.
```bash
└── WorkDirectory
        ├── my_ancestor_sample
        └── sample1  *NEW*
            ├── dup_metrics
            ├── sample1_all_variants.bed
            ├── sample1_final_stringent_compiled.txt <-- compiled results from all 3 callers
            ├── sample1_freebayes_BCBio_AncFiltered_annotated_vcf.txt
            ├── sample1_freebayes_BCBio_AncFiltered_condensed.txt
            ├── sample1_freebayes_BCBio_AncFiltered.vcf
            ├── sample1_freebayes_BCBio.vcf
            ├── sample1_gatk_haplo_AncFiltered_annotated_vcf.txt
            ├── sample1_gatk_haplo_AncFiltered_condensed.txt
            ├── sample1_gatk_haplo_AncFiltered.vcf
            ├── sample1_gatk_haplo.vcf
            ├── sample1_gatk_haplo.vcf.idx
            ├── sample1_highConfidenceVars.bed
            ├── sample1_lofreq_AncFiltered_annotated_vcf.txt
            ├── sample1_lofreq_AncFiltered_condensed.txt
            ├── sample1_lofreq_AncFiltered_reheadered.vcf
            ├── sample1_lofreq_AncFiltered.vcf
            ├── sample1_R1R2_MD.sort.bam
            ├── sample1_R1R2_MD.sort.bam.bai
            ├── sample1_R1R2_MD.sort.bam.sbi
            ├── sample1_S1_R1_001_fastqc.html
            ├── sample1_S1_R2_001_fastqc.html
            └── sample1_stringent_compiled.txt 
```
### Output File explanation

- The `gatk_haplo.vcf` and `freebayes_BCBio.vcf` are the raw output files from calling the specified variant callers. 

- The `freebayes_BCBio_AncFiltered.vcf`, `lofreq_AncFiltered.vcf`, and `gatk_haplo_AncFiltered.vcf` are files of each variant caller which has all the ancestor mutations already filtered out. This means that these files contains only the variants that were found over the course of your experiment.

- The `freebayes_BCBio_AncFiltered_annotated_vcf.txt`, `lofreq_AncFiltered_annotated_vcf.txt`, and `gatk_haplo_AncFiltered_annotated_vcf.txt` are annotated files of each variant caller which has all the ancestor mutations already filtered out. 

- The `freebayes_BCBio_AncFiltered_condensed.txt`, `lofreq_AncFiltered_condensed.txt`, and `gatk_haplo_AncFiltered_condensed.txt` are the same files as the ones mentioned above but we applyed a unique set of filter conditions for each file based on their specific variant caller and we condense the columns. The `_condensed.csv` files only show the columns that are relevant for our analysis. Such columns like `CHROM`, `POS`, `REF`, `ALT`, `ANNOTATION`, `REGION`, and `PROTEIN`. For the filter, a file's particular variant caller changes the filter conditions for `QUAL`,`DP`, and number of reads on the ref and alt alleles. We would keep any variants that pass the specified threshold for `QUAL`, `DP`, etc. 

- The `lofreq_AncFiltered_reheadered.vcf` is just the `lofreq_AncFiltered.vcf` file but reformatted to be compatible for data processing.

- The `final_stringent_compiled.txt` file contains the variants that were present in at least 2 out of 3 variant callers. This means that these are most likely real variants. 
  
- The `all_variants.bed` show all the variants that were called by our pipeline, while the `final_stringent_compiled_3caller.bed` file show the respective variants which are called by 2+ callers. It is recommended to use these bed files to check each variant in a genome alignment viewing software such as IGV: https://igv.org/

By opening the `final_stringent_compiled.txt` file in a program like Microsoft Excel, you can sort the called variants by quality score and/or number of occurrences across the different variant callers.

### IGV App BED File Tutorial
Now we want to view our results in IGV to verify our variants. We will be making use of the output bed files from our pipeline. 

1. Open IGV and on the top left corner tab, select "S. cerevisiae (sacCer3)" as the genome
2. Upload your `.bam` and `.bam.bai` files that you want to verify.
3. After this, go to the regions tab on the top and select "import regions". You'll want to import the `all_variants.bed` file.
4. After this, click back onto the regions tab at the top and select "region navigator". This will show you all your variants from your imported bed file.
5. Finally, just go through each variant by clicking on each line, and clicking on "view" to see the region in IGV. You can also choose to remove a particular variant if doesn't seem to be high quality by clicking on "remove". Afterwards, you can click back on the regions tab on the top, and export regions to get a bed file of all the variants that you have verified. 

