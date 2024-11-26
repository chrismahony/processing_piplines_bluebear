# How to process bulk ATACseq or ChIPseq

Both protocols will use the same approch


1. Complete this first: https://github.com/chrismahony/processing_piplines_bluebear/tree/main/dowloand_data_novogene

2. Organise your data a logical way, e.g.


```bash

/croftap-XXXX/
|-- my_project
| |--fastqs   #where you put fastqs
| |--count   #where you will run the next steps and perform alingment

```

3.  Save the below script into your 'count' folder (or what ever you name it) and call it something link 'count1.txt'.

 Be sure to amend paths to where your data is as well as #SBATCH --account=croftap-labdata2 to the name of an RDS folder you have access to.

```bash

#!/bin/bash
#SBATCH -n 60
#SBATCH -N 1
#SBATCH --mem 299G
#SBATCH --time 48:0:0
#SBATCH --mail-type ALL
#SBATCH --account=croftap-labdata2  #chnage this to the name of an RDS folder you have permission to access


set -e
module purge; module load bluebear
module load Bowtie2/2.4.4-GCC-10.3.0


# Ref genome can be downloaded for bowtie2 website, download and amend path to this

bowtie2 -x /path_to_ref_genome/GRCh38_noalt_as -U \
/path_to_fatsq1/sample1_EKRN230060415-1A_HVLNCDSX7_L2_1.fq.gz \
/path_to_fatsq2/sample1_EKRN230060415-1A_HVLNCDSX7_L2_2.fq.gz -S \
/path/count/sample1_unsorted.sam


# cd to where the output of bowtie2 will go
cd /path/count

module purge; module load bluebear
module load bear-apps/2021b
module load SAMtools/1.15.1-GCC-11.2.0


# Now sort and remove PCR duplicates
samtools view -bS sample1_unsorted.sam | samtools sort - >sample1_sorted.bam
samtools view -bq 1 sample1_sorted.bam > unique_sample1_sorted.bam

```

For ATACseq, modify the last part to remove reads in mt regions

```bash

5. For ATACseq, now remove peaks in mt chromosomes


```bash


module load bear-apps/2021b
module load SAMtools/1.15.1-GCC-11.2.0

samtools idxstats unique_sample1_sorted.bam | cut -f1 | grep -v Mt | xargs samtools view --threads 7 -b unique_sample1_sorted.bam > nomt_unique_sample1_sorted.bam

```




6. For both ATACseq and ChIPseq the next step is to call peaks:

```bash

#!/bin/bash
#SBATCH -n 20                       
#SBATCH -N 1                       
#SBATCH --mem 180000                 
#SBATCH --time 48:0:0
#SBATCH --mail-type ALL
#SBATCH --account=croftap-stia-atac


set -e

module purge; module load bluebear
module load bear-apps/2022a
module load MACS2/2.2.9.1-foss-2022a

mkdir output_macs2
mkdir tmp

macs2 callpeak -t /path1/unique_sample1_sorted.bam /path2/unique_sample2_sorted.bam /path3unique_sample3_sorted.bam \
-n all_peaks -g 2.8E9 \
--outdir ./output_macs2 \
--tempdir ./tmp

```


7. Next (only for ATACseq), remove blacklist regions and quanitfy reads in peaks

```bash

set -e
module purge; module load bluebear
module load bear-apps/2021b
module load BEDTools/2.30.0-GCC-11.2.0

#this uses peak set defined from all samples
bedtools intersect -a ALL_PEAKS_alltimes_merged_rest_peaks.broadPeak -b /rds/projects/m/mahonyc-kitwong-runx1/ATAC-processing_CM/hg38.blacklist.bed -v > ALL_PEAKS_alltimes_merged_rest_noblacklist.broadPeak

module purge; module load bluebear
module load Subread/2.0.1-GCC-8.3.0

#convert broad peak to saf
awk -F $'\t' 'BEGIN {OFS = FS}{ $2=$2+1; peakid=" ALL_PEAKS_alltimes_merged_rest_peak_"++nr;  print peakid,$1,$2,$3,"."}' ALL_PEAKS_alltimes_merged_rest_noblacklist.broadPeak > ALL_PEAKS_alltimes_merged_rest_noblacklist.broadPeak.saf

#count reads/peak/sample
featureCounts -p -F SAF -a ALL_PEAKS_alltimes_merged_rest_noblacklist.broadPeak.saf \
--fracOverlap 0.2 \
-o all_timepoints_merged_peaks_macs_broad.counts nomt_ATAC_merged.bam \
/path/nomt_unique_sample1_sorted.bam \
/path/nomt_unique_sample2_sorted.bam \
/path/nomt_unique_sample3_sorted.bam

```

8. The resulting count matrix can be read into R and analysed using DESeq2 as nroally with bulk RNAseq
