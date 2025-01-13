# Processing bulk ATACseq

1. First, download and unpack your data: https://github.com/chrismahony/processing_piplines_bluebear/tree/main/dowloand_data_novogene


2. Organise your data a logical way, e.g.

/croftap-XXXX/
|-- my_project
| |--fastqs   #where you put fastqs
| |--count   #where you will run the next steps and perform alingment


Save this next script in to your 'count' folder (or what ever you name it) and call it something link 'count1.txt'.
Be sure to amend paths to where your data is as well as #SBATCH --account=croftap-labdata2 to the name of an RDS folder you have access to.

Also amend the number of arrays, if you have 9 samples then you need:

```bash
#SBATCH --array=0-8

```

if you have 5 samples then you need:

```bash
#SBATCH --array=0-4

```


3. Next, you need to run fastqc, align your fastq files to reference genome and filter for unique reads

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

mkdir -p int_outs
mkdir -p final_outs
mkdir -p fastqc_results

cd ./int_outs/


####start processing

SAMPLE_INDEX=$SLURM_ARRAY_TASK_ID

SAMPLE_DIR=${FASTQ_DIRS[$SAMPLE_INDEX]}

FASTQ_FILES=("$SAMPLE_DIR"*.{fq,fastq}.gz)

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




```

4. Next ned to remove mt reads and index bams



```bash

#!/bin/bash
#SBATCH -n 60
#SBATCH -N 1
#SBATCH --mem 299G
#SBATCH --time 48:0:0
#SBATCH --mail-type ALL
#SBATCH --account=croftap-labdata2  #chnage this to the name of an RDS folder you have permission to access

export MRO_DISK_SPACE_CHECK=disable
set -e

module purge; module load bluebear
module load bear-apps/2022b
module load Bowtie2/2.5.1-GCC-12.2.0
module load SAMtools/1.17-GCC-12.2.0


# Directory containing the BAM files
INPUT_DIR="/rds/projects/c/croftap-croftapcarcia/RNA_ATAC_JAN2025/count/final_outs/"

# Iterate through all BAM files starting with "unique_" and ending with ".bam"
for BAM_FILE in "$INPUT_DIR"/unique_*.bam; do
    # Extract the basename without extension
    BASENAME=$(basename "$BAM_FILE" .bam)
    
    # Index the BAM file if not already indexed
    if [ ! -f "$BAM_FILE.bai" ]; then
        echo "Indexing $BAM_FILE..."
        samtools index "$BAM_FILE"
    fi

    # Define the output file path
    OUTPUT_FILE="$INPUT_DIR/${BASENAME}_no_mt.bam"
    
    # Use samtools to remove mt reads and save to new file
    samtools view -h "$BAM_FILE" | grep -v 'chrM' | samtools view -b -o "$OUTPUT_FILE"

 done


```


NB. in one case I have many very small reads (~8bp in lenght) along with many other good quality reads (~150bp). I filtered out these contaminating reads using the following script:



```bash

#!/bin/bash
#SBATCH -n 60
#SBATCH -N 1
#SBATCH --mem 299G
#SBATCH --time 48:0:0
#SBATCH --mail-type ALL
#SBATCH --account=croftap-labdata2  #chnage this to the name of an RDS folder you have permission to access

export MRO_DISK_SPACE_CHECK=disable
set -e

module purge; module load bluebear
module load bear-apps/2022b
module load SAMtools/1.17-GCC-12.2.0


# Define the directory containing the BAM files
input_dir="/rds/projects/c/croftap-croftapcarcia/RNA_ATAC_JAN2025/count/final_outs/"

# Loop through all BAM files matching the pattern in the specified directory
for bam_file in "$input_dir"/*_no_mt.bam; do
    if [ -f "$bam_file" ]; then
                
        # Generate the output filename
        filtered_bam="${bam_file%.bam}_filtered_nomt_final.bam"
        
        # Filter reads shorter than 10 bp and save to the filtered BAM file
        samtools view -h "$bam_file" | \
        awk 'length($10) >= 10 || $1 ~ /^@/' | \
        samtools view -b -o "$filtered_bam"
        
        # Index the filtered BAM file
        samtools index "$filtered_bam"
        
        echo "Finished processing $bam_file. Output: $filtered_bam"
    else
        echo "No files matching *_no_mt.bam in $input_dir."
    fi
done


```






5. Then Call peaks from all BAM files


```bash


#!/bin/bash
#SBATCH -n 72
#SBATCH -N 1
#SBATCH --mem 499G
#SBATCH --time 99:0:0
#SBATCH --qos castles
#SBATCH --mail-type ALL
#SBATCH --account=croftap-stia-atac
#SBATCH --constraint=icelake

set -e

module purge; module load bluebear


cd /rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/Bulk_human_data/GSE128644_ATACseq
source /rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/kolias_paper_data/my-virtual-env-haswell/bin/activate


macs3 callpeak -t nomt_ATAC_merged.bam /rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/Bulk_human_data/GSE128644_ATACseq/ATAC_tnfa72h/nomt_ATAC_TNF_merged.bam /rds/projects/m/mahonyc-kitwong-runx1/ATAC-processing_CM/nomt_unique_SRR8758528_sorted.bam /rds/projects/m/mahonyc-kitwong-runx1/ATAC-processing_CM/nomt_unique_SRR8758527_sorted.bam \
-n ALL_PEAKS_alltimes_merged_rest --broad -q 0.1 --nomodel -g 2.8E9 \


deactivate


```



6. Next remove those lying in blacklist regtions and generate a count matrix for analysis in R

```bash

!/bin/bash
#SBATCH -n 72
#SBATCH -N 1
#SBATCH --mem 499G
#SBATCH --time 99:0:0
#SBATCH --mail-type ALL
#SBATCH --account=croftap-stia-atac


cd /rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/Bulk_human_data/GSE128644_ATACseq

set -e

module purge; module load bluebear

#remove blacklist regions from hg38 nal merged peaks
module load bear-apps/2021b
module load BEDTools/2.30.0-GCC-11.2.0

bedtools intersect -a ALL_PEAKS_alltimes_merged_rest_peaks.broadPeak -b /rds/projects/m/mahonyc-kitwong-runx1/ATAC-processing_CM/hg38.blacklist.bed -v > ALL_PEAKS_alltimes_merged_rest_noblacklist.broadPeak

module purge; module load bluebear
module load Subread/2.0.1-GCC-8.3.0

#convert broad peak to saf
awk -F $'\t' 'BEGIN {OFS = FS}{ $2=$2+1; peakid=" ALL_PEAKS_alltimes_merged_rest_peak_"++nr;  print peakid,$1,$2,$3,"."}' ALL_PEAKS_alltimes_merged_rest_noblacklist.broadPeak > ALL_PEAKS_alltimes_merged_rest_noblacklist.broadPeak.saf

#count reads/peak/sample
featureCounts -p -F SAF -a ALL_PEAKS_alltimes_merged_rest_noblacklist.broadPeak.saf --fracOverlap 0.2 -o all_timepoints_merged_peaks_macs_broad.counts nomt_ATAC_merged.bam /rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/Bulk_human_data/GSE128644_ATACseq/ATAC_tnfa72h/nomt_ATAC_TNF_merged.bam /rds/projects/m/mahonyc-kitwong-runx1/ATAC-processing_CM/nomt_unique_SRR8758528_sorted.bam /rds/projects/m/mahonyc-kitwong-runx1/ATAC-processing_CM/nomt_unique_SRR8758527_sorted.bam



```
