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



Main script:


```bash


#!/bin/bash
#SBATCH -n 40
#SBATCH -N 1
#SBATCH --mem 199G
#SBATCH --time 72:0:0
#SBATCH --mail-type ALL
#SBATCH --account=croftap-stia-atac
#SBATCH --output=./log/test_%A_%a.out
#SBATCH --error=./log/test_%A_%a.err
#SBATCH --array=0-6

export MRO_DISK_SPACE_CHECK=disable

set -e
mkdir log
module purge; module load bluebear
module load CellRanger-ATAC/2.0.0


sample_ID=("ATACT2518_044a-AK14653_HHYV7DRX5,ATACT2518_044b-AK14654_HHYV7DRX5,ATACT2518_044c-AK14655_HHYV7DRX5,ATACT2518_044d-AK14656_HHYV7DRX5" \
"ATACT2519_019LKa-AK8422_H7V5NDRX5,ATACT2519_019LKb-AK8423_H7V5NDRX5,ATACT2519_019LKc-AK8424_H7V5NDRX5,ATACT2519_019LKd-AK8421_H7V5NDRX5" \
"ATACT2519_019RKa-AK13919_H7VMWDRX5,ATACT2519_019RKb-AK13921_H7VMWDRX5,ATACT2519_019RKc-AK13920_H7VMWDRX5,ATACT2519_019RKd-AK13922_H7VMWDRX5" \
"ATACT2519_021LKa-AK2543_H7VMWDRX5,ATACT2519_021LKb-GA01_H7VMWDRX5,ATACT2519_021LKc-X127_H7VMWDRX5,ATACT2519_021LKd-X172_H7VMWDRX5" \
"ATACT2519_021RKa-AK12521_H7VMWDRX5,ATACT2519_021RKb-AK12522_H7VMWDRX5,ATACT2519_021RKc-AK12523_H7VMWDRX5,ATACT2519_021RKd-AK11699_H7VMWDRX5" \
"ATACT2519_022a-X002_HHYV7DRX5,ATACT2519_022b-AK7142_HHYV7DRX5,ATACT2519_022c-X087_HHYV7DRX5,ATACT2519_022d-X133_HHYV7DRX5" \
"ATACT2519_023a-A6_HHYV7DRX5,ATACT2519_023b-AK8610_HHYV7DRX5,ATACT2519_023c-AK8609_HHYV7DRX5,ATACT2519_023d-AK8611_HHYV7DRX5")

output_ID=("ATACT2518_044" "ATACT2519_019LK" "ATACT2519_019RK" "ATACT2519_021LK" "ATACT2519_021RK" "ATACT2519_022" "ATACT2519_023")

cellranger-atac count --id=${output_ID[$SLURM_ARRAY_TASK_ID]} \
                   --reference=/rds/projects/c/croftap-mapjag-batch5/scATAC/count/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                   --fastqs=/rds/projects/c/croftap-mapjagb11/scATACseq/all_fastqs_all_batches \
                   --sample=${sample_ID[$SLURM_ARRAY_TASK_ID]} \
                   --localcores=8 \
                   --localmem=64



```



5. Should take ~6 h per sample, could be longer if more cells
