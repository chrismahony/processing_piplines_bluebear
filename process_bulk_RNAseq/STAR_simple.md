# Aligning FASTQ data to reference genome using STAR - simple method


1. Run FASTQC on all fastq files and multiqc

```bash


#!/bin/bash
#SBATCH -n 40
#SBATCH -N 1
#SBATCH --mem 200G
#SBATCH --time 24:0:0
#SBATCH --mail-type ALL
#SBATCH --account=croftap-XXXX # Change to your RDS


module purge; module load bluebear
module load FastQC/0.11.9-Java-11

cd /path_to_fastqs/

for f in *.fastq.gz; do fastqc -o ./ $f ; done

module purge; module load bluebear
module load MultiQC/1.9-foss-2019b-Python-3.7.4

multiqc .

```


2. Remove adapter is (if require)

