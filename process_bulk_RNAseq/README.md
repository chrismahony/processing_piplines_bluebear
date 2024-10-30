# Processing Bulk RNAseq data

1. Complete this first: https://github.com/chrismahony/processing_piplines_bluebear/tree/main/dowloand_data_novogene




2. Organise your data a logical way, e.g.



```bash


/croftap-XXXX/
|-- my_project
| |--fastqs   #where you put fastqs
| |--count   #where you will run the next steps and perform alingment



``` 



3.  Save this script in to your 'count' folder (or what ever you name it) and call it something link 'count1.txt'. Be sure to amend paths to where your data is as well as #SBATCH --account=croftap-labdata2 to the name of an RDS folder you have access to.


```bash

#!/bin/bash
#SBATCH -n 60
#SBATCH -N 1
#SBATCH --mem 299G
#SBATCH --time 48:0:0
#SBATCH --mail-type ALL
#SBATCH --account=croftap-labdata2  #chnage this to the name of an RDS folder you have permission to access
#SBATCH --array=0-8  # if you have 3 samples then this whould be 0-2

export MRO_DISK_SPACE_CHECK=disable
set -e

module purge; module load bluebear
module load FastQC/0.11.9-Java-11


#root dir is a director with folders for each sample, in each folder there are fastq files
ROOT_DIR="/rds/projects/c/croftap-XXXX/bulkRNAseq_X204SC14204168-Z01-F001/fastqs/X204SC14204168-Z01-F001/01.RawData"
FASTQ_DIRS=("$ROOT_DIR"/*/)  

#where your index is stored
BOWTIE_INDEX="/add_your_path/mm10/mm10"

cd /rds/projects/c/croftap-XXXX/bulkRNAseq_X204SC14204168-Z01-F001/

#create output dirs if they dont exist

if [ ! -d "./int_outs" ]; then
  mkdir ./int_outs
  echo "Created directory: int_outs"
fi

if [ ! -d "./final_outs" ]; then
  mkdir ./final_outs
  echo "Created directory: final_outs"
fi

if [ ! -d "./final_outs" ]; then
  mkdir ./fastqc_results
  echo "Created directory: fastqc_results"
fi


cd ./int_outs/


####start processing

SAMPLE_INDEX=$SLURM_ARRAY_TASK_ID

SAMPLE_DIR=${FASTQ_DIRS[$SAMPLE_INDEX]}

FASTQ_FILES=("$SAMPLE_DIR"*.fq.gz)

 #Check if there are any FASTQ files
if [ ${#FASTQ_FILES[@]} -eq 0 ]; then
    echo "No FASTQ files found in $SAMPLE_DIR"
    exit 1
fi

echo "Running FastQC on $SAMPLE_ID"
for fastq_file in "${FASTQ_FILES[@]}"; do
    fastqc -o "../fastqc_results" "$fastq_file"
done


module purge; module load bluebear
module load bear-apps/2022b
module load Bowtie2/2.5.1-GCC-12.2.0
module load SAMtools/1.17-GCC-12.2.0


SAMPLE_ID=$(basename "$SAMPLE_DIR")
SAM_FILE="${SAMPLE_ID}_unsorted.sam"

bowtie2 -x "$BOWTIE_INDEX" -U "${FASTQ_FILES[@]}" -S "$SAM_FILE"

samtools view -bS "$SAM_FILE" | samtools sort - > "${SAMPLE_ID}_sorted.bam"
samtools view -bq 1 "${SAMPLE_ID}_sorted.bam" > "../final_outs/unique_${SAMPLE_ID}_sorted.bam"

module purge; module load bluebear
module load MultiQC/1.9-foss-2019b-Python-3.7.4


if [ "$SAMPLE_INDEX" -eq "${#FASTQ_DIRS[@]}" ]; then
    echo "All FastQC runs complete. Running MultiQC."
    multiqc ../fastqc_results -o ../final_outs/multiqc_report
fi


if [ "$SAMPLE_INDEX" -eq 0 ]; then
    sbatch --dependency=afterok:$SLURM_JOB_ID ../count_all.txt
fi

``` 

4. Next save the following script in the same folder as the above script. Save it as count_all.txt. Be sure to amend paths to where your data is as well as #SBATCH --account=croftap-labdata2 to the name of an RDS folder you have access to

```bash 


#!/bin/bash
#SBATCH -n 50
#SBATCH -N 1
#SBATCH --mem 180000
#SBATCH --time 68:0:0
#SBATCH --mail-type ALL
#SBATCH --account=croftap-labdata2


set -e

module purge; module load bluebear
module load Subread/2.0.1-GCC-8.3.0

BAM_DIR="/rds/projects/c/croftap-XXXX/bulkRNAseq_X204SC14204168-Z01-F001/count/final_outs"
BAM_FILES=$(ls ${BAM_DIR}/*.bam)

featureCounts -t exon -g gene_name -a /add_your_path/output/GENOME/MM10.gtf \
    -o ./all_counts.txt \
    $BAM_FILES

```


5. Now open a terminal, navigate to your folder and submit the first script:


```bash


sbatch count1.txt


```


6. Once all jobs in this array have finihsed then count_all.txt will run to process the final counts. Make sure you save the second script as count_all.txt or it will not run!!


