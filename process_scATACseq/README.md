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



6. Once all sample have fiihsed the easiest way to combine samples it to use cell ranger atac aggr in an sbatc script as below:



```bash

#!/bin/bash
#SBATCH -n 110
#SBATCH -N 1
#SBATCH --mem 399G
#SBATCH --time 99:0:0
#SBATCH --mail-type ALL
#SBATCH --account=croftap-stia-atac
#SBATCH --constraint=sapphire


export MRO_DISK_SPACE_CHECK=disable

set -e

module purge; module load bluebear

module load CellRanger-ATAC/2.0.0

cellranger-atac aggr --id=all_aggr_mapjag \
                   --reference=/rds/projects/c/croftap-mapjag-batch5/scATAC/count/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                   --normalize=none \
		   --csv=libraries_mapjag.csv

```


With several samples, you are likly going to been nearly the max memory allowance to process them all in a resobale time, not the #SBATCH constrain added to specify a specific node that allows more cores to be selected:


```bash
#SBATCH --constraint=sapphire

```


the libaires file should look like this:


```bash

library_id,fragments,cells
ATAC2518_19,/rds/projects/c/croftap-mapjagb6/atac/count/ATAC2518_19/outs/fragments.tsv.gz,/rds/projects/c/croftap-mapjagb6/atac/count/ATAC2518_19/outs/singlecell.csv
ATAC2518_31,/rds/projects/c/croftap-mapjagb6/atac/count/ATAC2518_31/outs/fragments.tsv.gz,/rds/projects/c/croftap-mapjagb6/atac/count/ATAC2518_31/outs/singlecell.csv
ATAC2518_32,/rds/projects/c/croftap-mapjagb6/atac/count/ATAC2518_32/outs/fragments.tsv.gz,/rds/projects/c/croftap-mapjagb6/atac/count/ATAC2518_32/outs/singlecell.csv
ATAC2518_36,/rds/projects/c/croftap-mapjagb6/atac/count/ATAC2518_36/outs/fragments.tsv.gz,/rds/projects/c/croftap-mapjagb6/atac/count/ATAC2518_36/outs/singlecell.csv
ATAC2518_39,/rds/projects/c/croftap-mapjagb6/atac/count/ATAC2518_39/outs/fragments.tsv.gz,/rds/projects/c/croftap-mapjagb6/atac/count/ATAC2518_39/outs/singlecell.csv
ATAC2519_14,/rds/projects/c/croftap-mapjag-batch5/scATAC/count/ATAC2519_14/outs/fragments.tsv.gz,/rds/projects/c/croftap-mapjag-batch5/scATAC/count/ATAC2519_14/outs/singlecell.csv
ATAC2519_15,/rds/projects/c/croftap-mapjag-batch5/scATAC/count/ATAC2519_15/outs/fragments.tsv.gz,/rds/projects/c/croftap-mapjag-batch5/scATAC/count/ATAC2519_15/outs/singlecell.csv
ATAC2518_28,/rds/projects/c/croftap-mapjagb8/scATACseq/count/ATAC2518_028/outs/fragments.tsv.gz,/rds/projects/c/croftap-mapjagb8/scATACseq/count/ATAC2518_028/outs/singlecell.csv
ATAC2518_29,/rds/projects/c/croftap-mapjagb8/scATACseq/count/ATAC2518_029/outs/fragments.tsv.gz,/rds/projects/c/croftap-mapjagb8/scATACseq/count/ATAC2518_029/outs/singlecell.csv


```



7. The way cellrnager call peaks is pretty average and generally peaks need to be called using MACS2, fitlered for low quality peaks and then these used for cellrnager reanalyze

   To call peaks on your data using MACS2 a sbatch script can be run as below:

   -g should be set at 1.87e9 for mouse


```bash

#!/bin/bash
#SBATCH -n 20                       
#SBATCH -N 1                       
#SBATCH --mem 180000                 
#SBATCH --time 48:0:0
#SBATCH --qos castles
#SBATCH --mail-type ALL
#SBATCH --account=croftap-stia-atac


set -e

module purge; module load bluebear
module load bear-apps/2022a
module load MACS2/2.2.9.1-foss-2022a

mkdir output_macs2
mkdir tmp

macs2 callpeak -t /path1/possorted_bam.bam /path2/possorted_bam.bam /path3/possorted_bam.bam \
-n all_peaks -g 2.8E9 \
--outdir ./output_macs2 \
--tempdir ./tmp

```


8. Filter peaks and remove low quality ones


First, open the .NARROWPEAK file and delete all colnames and other information above peaks (if there are anything)

Then start an IGV session using Bluebear ondemand.

Select Hg38 as a genome.

Then File > Load from file > your.NARROWPEAK file
 
Next in IGV, search for some key genes, e.g. COL1A1 etc

Look out for peaks with poor Signal Value and qValue to dertmine appropriate cuttoffs

This can be hard so another option is reading the file into R and looking at histograms


```R

ALL_PEAKS_annie_peaks <- read.delim("/rds/projects/c/croftap-mapjag-batch5/scATAC/aggr/ALL_PEAKS_annie_peaks.narrowPeak", header=FALSE)

ggplot(ALL_PEAKS_annie_peaks, aes(x = V9)) +
  geom_histogram(
    bins = 100, 
    aes(y = ..density.., fill = ..count..), 
    color = "black", 
    alpha = 0.7
  ) +
  scale_fill_gradient(low = "blue", high = "red", name = "Peak Count") +
  labs(
    title = "Distribution of qValues in Peaks",
    x = "nCount Peaks",
    y = "Density"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )+theme_ArchR()


ALL_PEAKS_annie_peaks$V9 %>% median()/4

> [1] 3.212112

```

A cutoff of 3 for the qValue is resonable



```R


ggplot(ALL_PEAKS_annie_peaks, aes(x = V7)) +
  geom_histogram(
    bins = 100, 
    aes(y = ..density.., fill = ..count..), 
    color = "black", 
    alpha = 0.7
  ) +
  scale_fill_gradient(low = "blue", high = "red", name = "Peak Count") +
  labs(
    title = "Distribution of Singal strength in Peaks",
    x = "nCount Peaks",
    y = "Density"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )+theme_ArchR()


ALL_PEAKS_annie_peaks$V7 %>% median()/4

[1] 0.93735




```bash



```



9. Now run cell range re-analyze

```bash


#!/bin/bash
#SBATCH -n 40
#SBATCH -N 1
#SBATCH --mem 399G
#SBATCH --time 48:0:0
#SBATCH --qos castles
#SBATCH --mail-type ALL
#SBATCH --account=croftap-stia-atac


export MRO_DISK_SPACE_CHECK=disable

set -e

module purge; module load bluebear

module load CellRanger-ATAC/2.0.0

cellranger-atac reanalyze --id=all_aggr_mapjag_reanalyze \
                   --reference=/rds/projects/c/croftap-mapjag-batch5/scATAC/count/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
--fragments=/rds/projects/c/croftap-mapjag-batch5/scATAC/aggr/all_aggr_mapjag/outs/fragments.tsv.gz \
--peaks=/rds/projects/c/croftap-mapjag-batch5/scATAC/aggr/ALL_PEAKS_MAPJAG_peaks_f2.txt


```

10. Now ready for analysis in R




