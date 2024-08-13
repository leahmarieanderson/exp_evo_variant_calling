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
These directory names: `errors`, `outputs`, `genomes`, and `fastq` are hard coded into the scripts so the directory names are different, the scripts will likely fail.
```bash
└── experiment1
    ├── errors (ADD THIS)
    ├── fastq
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
