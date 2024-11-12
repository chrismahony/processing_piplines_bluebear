# Processing scATACseq data


1. Complete this step first: https://github.com/chrismahony/processing_piplines_bluebear/tree/main/dowloand_data_novogene


2. Then make a directory call something like 'all_fastqs'. Put all the fastq files from all samples that you want to process. Do not make sub directories, just put all files in the one dir 'all_fastqs'

Make a directory called something like 'count', so your project folder would look like:
/croftap-XXXX/
|-- my_project
| |--fastqs #downloaded fatsqs from novogene
| |--all_fastqs   #where you put all fastqs
| |--count   #where you will run the next steps and perform alingment


You may have a load of sub dirs, each with fastq files in them, if so, open a termonal, cd to where all the sub dirs are and run this to move all the fastq files out:


```bash

for dir in */; do
    mv "$dir"*.fastq* .
done

```



4. The copy the script below and save as something like: count.txt
!Important!
Be sure to amend the:

```bash

#SBATCH --account=croftap-XXXX

```

To a RDS folder you have access to.

You will also need to amend:


```bash
FASTQ_DIR #the folder you put all the fastq files
OUTPUT_DIR #where you wan to store the data
REF_PATH #where you ref genome is

```


Main script:


```bash


#!/bin/bash
#SBATCH -n 70                            # Number of cores
#SBATCH --mem=399G                       # Total memory
#SBATCH --time=99:0:0                    # Max runtime
#SBATCH --mail-type=ALL
#SBATCH --account=croftap-XXX

set -e
module purge; module load bluebear
module load CellRanger-ATAC/2.0.0

# Set paths
FASTQ_DIR="/rds/my_path/all_fastqs"  
OUTPUT_DIR="/rds/my_path/count/"      
REF_PATH="/rds/bear-apps/apps-data/CellRanger/refdata-cellranger-mm10-2.1.0"     

# Extract unique sample names (removing path and everything after the first underscore)
for sample_name in $(ls $FASTQ_DIR/*.fastq.gz | sed 's|.*/||; s|_\(.*\)||' | sort | uniq); do
    
    # Gather unique sample names related to the current sample
    sample_names=$(ls $FASTQ_DIR | grep "^${sample_name%_*}_" | sed 's|_S.*||' | sort | uniq | tr '\n' ',' | sed 's|,$||')

    sample_output_dir="$OUTPUT_DIR/$sample_name"

    # Run CellRanger count for the sample
    echo "Running cellranger-atac count for sample: $sample_name"
    cellranger-atac count --id=$sample_name \
                     --transcriptome=$REF_PATH \
                     --fastqs=$FASTQ_DIR \
                     --sample=$sample_names \
                     --localcores=8 \
                     --localmem=64

   done



```



