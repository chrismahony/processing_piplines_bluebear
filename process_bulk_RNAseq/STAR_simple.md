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

```bash

#!/bin/bash
#SBATCH -n 40
#SBATCH -N 1
#SBATCH --mem 200G
#SBATCH --time 24:0:0
#SBATCH --mail-type ALL
#SBATCH --account=croftap-XXXX # Change to your RDS


module purge; module load bluebear
module load Trim_Galore/0.6.5-GCCcore-8.3.0-Java-11-Python-3.7.4

cd /path_to_fastqs/

#!/bin/bash

# Input directory containing the raw fastq.gz files
input_dir="/path/to/your/raw_fastqs"

# Output directory for the trimmed fastqs
output_dir="/rds/projects/c/croftap-fbmacro01/trimmed_fastqs"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Run trim_galore on all .fastq.gz files in the input directory
for f in "$input_dir"/*.fastq.gz; do
    echo "Processing $f..."
    trim_galore \
        -q 20 \
        --fastqc_args "--outdir $output_dir" \
        --illumina \
        --gzip \
        -o "$output_dir" \
        --length 20 "$f"
done

echo "Trimming and FastQC completed for all files."



```


3. Align to reference genome using STAR, for this you will need to create a ref genome in STAR format

```bash

#!/bin/bash
#SBATCH -n 40
#SBATCH -N 1
#SBATCH --mem 180000
#SBATCH --time 10:0:0
#SBATCH --mail-type ALL
#SBATCH --account=croftap-XXXX


set -e
module purge; module load bluebear
module load bear-apps/2022b
module load STAR/2.7.11a-GCC-12.2.0


STAR --runThreadN 60 \
--runMode genomeGenerate \
--genomeDir /rds/path/genome_dir \
--genomeFastaFiles ./GRCh38.primary_assembly.genome.fa \
--sjdbGTFfile ./gencode.v29.annotation.gtf


```



4. Now run STAR on fastq files

```bash

#!/bin/bash
#SBATCH -n 72
#SBATCH -N 1
#SBATCH --mem 180000
#SBATCH --time 68:0:0
#SBATCH --mail-type ALL
#SBATCH --account=croftap-XXXX


set -e
module purge; module load bluebear
module load bear-apps/2022b
module load STAR/2.7.11a-GCC-12.2.0


STAR --genomeDir /rds/path/genome_dir \
--runThreadN 60 \
--readFilesIn /rds/path/s1_1_read1.fastq.gz /rds/path/s1_1_read2.fastq.gz,/rds/path/s1_2_read1.fastq.gz /rds/path/s1_2_read2.fastq.gz \
--outFileNamePrefix /rds/path/sample1_ \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard 

```
