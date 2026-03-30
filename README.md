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
## Usage
We have two main steps in our pipeline:
1. Align the ancestor strain and process its bam and vcf files.
2. Align the evolved strain's variants, apply filters based on the results from aligning ancestor strain and put the results into comprehensive output files 

### Align Ancestor strain
First, we need to create the bam files as well as the samtools and freebayes vcfs for our ancestor fastq files.   

Start by going into the exp_evo_variant_calling directory that you have cloned or copied:
```php
cd path/to/exp_evo_variant_calling
```
Here, we will need to change some of the working directories in the align script. We can do this by using the `batch_submit.py` script:

```php
$ python3 batch_submit.py
```
Next you will see this popup:

```
Current Bash Settings for scripts:
#$ -S /bin/bash

#$ -wd /net/dunham/vol2/Zilong/updating_pipeline_2024

#$ -o /net/dunham/vol2/Zilong/updating_pipeline_2024/outputs/

#$ -e /net/dunham/vol2/Zilong/updating_pipeline_2024/errors/

#$ -l mfree=8G

#$ -l h_rt=36:0:0


**NEW** Bash Settings for scripts:
#$ -S /bin/bash

#$ -wd /net/dunham/vol2/Leah/labmeeting_250613

#$ -o /net/dunham/vol2/Leah/labmeeting_250613/outputs/

#$ -e /net/dunham/vol2/Leah/labmeeting_250613/errors/

#$ -l mfree=8G

#$ -l h_rt=36:0:0

Change your SGE Directives to fit your directories? (y/n) :
```
Then this:
```
Current Script Variables:
FOLDER=fastq
DIR=/net/dunham/vol2/Leah/yEvo_sequencing250520
SEQID=delmont # Project name and date for bam header
REF=${DIR}/exp_evo_variant_calling/genomes/sacCer3.fasta # Reference genome
SCRIPTS=${DIR}/exp_evo_variant_calling # Path of annotation_final.py directory
What would you like your SEQID to be? (this could be your project name):
```
Here, you will enter the name of your project, then press enter.
Next, you will see a prompt that looks something like this:

```
Current Variables Settings for scripts:
FOLDER=fastq
DIR=/net/dunham/vol2/Zilong/updating_pipeline_2024
SEQID=test # Project name and date for bam header
REF=${DIR}/exp_evo_variant_calling/genomes/sacCer3.fasta # Reference genome
SCRIPTS=${DIR}/exp_evo_variant_calling # Path of annotation_final.py directory

**NEW** Variables Settings for scripts:
FOLDER=fastq
DIR=/net/dunham/vol2/Leah/yEvo_sequencing250520
SEQID=test # Project name and date for bam header
REF=${DIR}/exp_evo_variant_calling/genomes/sacCer3.fasta # Reference genome
SCRIPTS=${DIR}/exp_evo_variant_calling # Path of annotation_final.py directory
Change your variable paths to fit your directories? (y/n) :
```
If you enter "y", then this code will rewrite the file paths in `align.sh` to match the directory you are currently working in. If you don't want the file paths to change, then enter "n".

Next, the program will ask if you want to run `align.sh` on all samples in the provided directory:

```
Script 'align.sh' updated.
Would you like to qsub all samples in /net/dunham/vol2/Leah/yEvo_sequencing250520/fastq/ ? (y/n) :
```
#### *IMPORTANT*: If you have not yet run `align.sh` on your ancestor, you should enter "n" before proceeding. 

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
    └── WorkDirectory
        └── my_ancestor_strain
            ├── dup_metrics
            ├── anc_AB_freebayes_BCBio_quality_filter.vcf
            ├── anc_AB_freebayes_BCBio.vcf
            ├── anc_AB_gatk_haplo_quality_filter.vcf
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
3. `anc_AB_freebayes_BCBio_quality_filter.vcf` is a post-processed vcf file of the freebayes output where we filter out some reads based on certain thresholds on mapping quality, read depth, and other metrics.
4. `anc_AB_gatk_haplo.vcf` is the output of one of the 2 variant callers that we call on the ancestor (FreeBayes).
5. `anc_AB_gatk_haplo_quality_filter.vcf` is a post-processed vcf file of the GATK HaplotypeCaller output where we filter out reads based on mapping quality, read depth, and other metrics.
6. `anc_AB_gatk_haplo.vcf.idx` is an index file for the GATK HaplotypeCaller output.
7. `anc_AB_R1R2_MD.sort.bam`, `anc_AB_R1R2_MD.sort.bam.bai`, and `anc_AB_R1R2_MD.sort.bam.sbi` are the typical `.bam` file outputs that we get.
8. `anc_AB_S1_R1_001_fastqc.html` and `anc_AB_S1_R2_001_fastqc.html` are Fastq Quality Control Check outputs which can tell us analytical statistics on each of the short reads. 

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
            ├── sample1_final_stringent_compiled.txt <-- compiled results from all 3 callers
            ├── sample1_freebayes_BCBio_AncFiltered.vcf
            ├── sample1_freebayes_BCBio_AncFiltered_annotated_vcf.txt
            ├── sample1_freebayes_BCBio_AncFiltered_condensed.txt
            ├── sample1_freebayes_BCBio.vcf
            ├── sample1_gatk_haplo_AncFiltered.vcf
            ├── sample1_gatk_haplo_AncFiltered_annotated_vcf.txt
            ├── sample1_gatk_haplo_AncFiltered_condensed.txt
            ├── sample1_gatk_haplo.vcf
            ├── sample1_gatk_haplo.vcf.idx
            ├── sample1_lofreq_AncFiltered.vcf
            ├── sample1_lofreq_AncFiltered_annotated_vcf.txt
            ├── sample1_lofreq_AncFiltered_condensed.txt
            ├── sample1_R1R2_MD.sort.bam
            ├── sample1_R1R2_MD.sort.bam.bai
            ├── sample1_R1R2_MD.sort.bam.sbi
            ├── sample1_S1_R1_001_fastqc.html
            └── sample1_S1_R2_001_fastqc.html
```
### Output File explanation

- The `gatk_haplo.vcf` and `freebayes_BCBio.vcf` are the raw output files from calling the specified variant callers. 

- The `freebayes_BCBio_AncFiltered.vcf`, `lofreq_AncFiltered.vcf`, and `gatk_haplo_AncFiltered.vcf` are files of each variant caller which has all the ancestor mutations already filtered out. This means that these files contains only the variants that were found over the course of your experiment.

- The `freebayes_BCBio_AncFiltered_annotated_vcf.txt`, `lofreq_AncFiltered_annotated_vcf.txt`, and `gatk_haplo_AncFiltered_annotated_vcf.txt` are annotated files of each variant caller which has all the ancestor mutations already filtered out. 

- The `freebayes_BCBio_AncFiltered_condensed.txt`, `lofreq_AncFiltered_condensed.txt`, and `gatk_haplo_AncFiltered_condensed.txt` are the same files as the ones mentioned above but we applyed a unique set of filter conditions for each file based on their specific variant caller and we condense the columns. The `_condensed.csv` files only show the columns that are relevant for our analysis. Such columns like `CHROM`, `POS`, `REF`, `ALT`, `ANNOTATION`, `REGION`, and `PROTEIN`. For the filter, a file's particular variant caller changes the filter conditions for `QUAL`,`DP`, and number of reads on the ref and alt alleles. We would keep any variants that pass the specified threshold for `QUAL`, `DP`, etc. 

- The `final_stringent_compiled.txt` file is a combination of the `freebayes_BCBio_AncFiltered_condensed.txt`, `lofreq_AncFiltered_condensed.txt`, and `gatk_haplo_AncFiltered_condensed.txt` that has been sorted, removed duplicates, and added an additional column `NUM_OCCURANCES` that counts the number of times this variant has shown between the different variant callers. The higher this value, the more reliable this variant is since it means it was called by more variant callers. 

It is recommended that each variant is then checked in a genome alignment viewing software such as IGV: https://igv.org/
By opening the `final_stringent_compiled.txt` file in a program like Microsoft Excel, you can sort the called variants by quality score and/or number of occurrences across the different variant callers.
