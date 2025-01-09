# Processing bulk ATACseq

1. First, download and unpack your data: https://github.com/chrismahony/processing_piplines_bluebear/tree/main/dowloand_data_novogene

2. Next, you need to aling your fastq files to reference genome, filter for unique reads and remove mt reads

```bash


#!/bin/bash
#SBATCH -n 20
#SBATCH -N 1
#SBATCH --mem 180000
#SBATCH --time 68:0:0
#SBATCH --mail-type ALL
#SBATCH --account=croftap-stia-atac


set -e
module purge; module load bluebear


module load bear-apps/2021b
module load SAMtools/1.15.1-GCC-11.2.0

samtools idxstats unique_SRR8758526_sorted.bam | cut -f1 | grep -v Mt | xargs samtools view --threads 7 -b unique_SRR8758526_sorted.bam > nomt_unique_SRR8758526_sorted.bam
samtools idxstats unique_SRR8758524_sorted.bam | cut -f1 | grep -v Mt | xargs samtools view --threads 7 -b unique_SRR8758524_sorted.bam > nomt_unique_SRR8758524_sorted.bam

samtools index nomt_unique_SRR8758526_sorted.bam
samtools index nomt_unique_SRR8758524_sorted.bam



```

3. Then Call peaks from all BAM files


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



4. Next remove those lying in blacklist regtions and generate a count matrix for analysis in R

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
