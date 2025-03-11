# Processing scATACseq data


1. Complete this step first: https://github.com/chrismahony/processing_piplines_bluebear/tree/main/dowloand_data_novogene


2. Then make a directory call something like 'all_fastqs'. Put all the fastq files from all samples that you want to process. Do not make sub directories, just put all files in the one dir 'all_fastqs'

Make a directory called something like 'count', so your project folder would look like:

```bash

/croftap-XXXX/
|-- my_project
| |--fastqs #downloaded fatsqs from novogene
| |--all_fastqs   #where you put all fastqs
| |--count   #where you will run the next steps and perform alingment

```


You may have a load of sub dirs, each with fastq files in them, if so, open a termonal, cd to where all the sub dirs are and run this to move all the fastq files out:


```bash

for dir in */; do
    mv "$dir"*.fastq* .
done

```

Alternative to move all files up one level into a dir called all_fastqs

```bash

cd /rds/projects  # navigate to the folder there all the sub dirs are found with the fastq files

mkdir all_fastqs

find ./ -mindepth 1 -type f -exec mv "{}" all_fastqs/ \;

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

export MRO_DISK_SPACE_CHECK=disable

set -e


module purge; module load bluebear
module load CellRanger-ATAC/2.0.0


# Set paths
FASTQ_DIR="/rds/projects/c/croftap-mapjagatac12/X204SC25022048-Z01-F001/fastqs/X204SC25022048-Z01-F001/01.RawData/all_fastqs"  
OUTPUT_DIR="/rds/projects/c/croftap-mapjagatac12/X204SC25022048-Z01-F001/fastqs/X204SC25022048-Z01-F001/count/"      
REF_PATH="/rds/projects/c/croftap-mapjag-batch5/scATAC/count/refdata-cellranger-arc-GRCh38-2020-A-2.0.0"     

# Extract unique sample names (removing path and everything after the first underscore)
for sample_name in $(ls $FASTQ_DIR/*.fastq.gz | sed 's|.*/||; s|_\(.*\)||' | sort | uniq); do
    
    # Gather unique sample names related to the current sample
    sample_names=$(ls $FASTQ_DIR | grep "^${sample_name%_*}_" | sed 's|_S.*||' | sort | uniq | tr '\n' ',' | sed 's|,$||')

    sample_output_dir="$OUTPUT_DIR/$sample_name"

    # Run CellRanger atac count for the sample
    cellranger-atac count --id=$sample_name \
                   --reference=$REF_PATH \
                   --fastqs=$FASTQ_DIR \
                   --sample=$sample_names \
                   --localcores=8 \
                   --localmem=64


done




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
module load bear-apps/2022b
module load MACS3/3.0.1-foss-2022b

mkdir output_macs3
mkdir tmp

macs3 callpeak -t /path1/possorted_bam.bam /path2/possorted_bam.bam /path3/possorted_bam.bam \
-n all_peaks -g 2.8E9 \
--outdir ./output_macs3 \
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

```

Singlal strength cutoff of 1 is a little low. Typically I would apply a threshold of 1.5 - 2

Next apply these cutt off in bash.

```bash

cd /rds/projects/your_path/

# Remove all peaks in col 9 with value < 3 and pipe to test.txt
cat ALL_PEAKS_peaks.narrowPeak | awk '$9 > 3' > test.txt 

# Remove all peaks in col 7 with value < 1.5 and pipe to test2.txt
cat test.txt | awk '$7 > 1.5' > test2.txt

# Subset to columns 1,2,3
cat test2.txt | cut -f 1-3 > test3.txt 

# Pipe test3 to a narrowPeak file 
mv test3.txt ALL_peaks_filtered.narrowPeak 

# Remove all intermediate test files
rm test*.* 

```


For convienvce, this is written in a complete pipeline:


```bash

cat ALL_PEAKS_peaks.narrowPeak | awk '$9 > 3' > test.txt 
cat test.txt | awk '$7 > 1.5' > test2.txt
cat test2.txt | cut -f 1-3 > test3.txt 
mv test3.txt ALL_peaks_filtered.narrowPeak 
rm test*.*


```



Next for cell ranger reanalyze, the peaks need to be in the correct order. To do this take the origional peak file and subset based on the filterest file from above


```bash

grep -Fwf ALL_peaks_filtered.narrowPeak ALL_PEAKS_peaks.narrowPeak > ALL_PEAKS_peaks_filtered2.narrowPeak
cat ALL_PEAKS_peaks_filtered2.narrowPeak | cut -f 1-3 > ALL_PEAKS_peaks_filtered2.txt 

```

Now to remove all peaks that are in unlocalised scafolds or in mitochonrial chromosomes


e.g. the start of my fial looks like:

```bash


head ALL_PEAKS_peaks_filtered2.txt
GL000008.2      106     4279
GL000008.2      83069   84024
GL000008.2      103773  104088
GL000008.2      118521  119272
GL000008.2      125221  125479
GL000008.2      126683  127050
GL000008.2      127342  127926
GL000008.2      132424  133084
GL000008.2      147357  147723
GL000008.2      153414  153740


```

```bash

grep "^chr" ALL_PEAKS_peaks_filtered2.txt > filtered_peaks.txt
grep -v "^chrM" filtered_peaks.txt > filtered_peaks_final.txt
[mahonyc@bear-pg-login06 output_macs2]$ head filtered_peaks_final.txt
chr1    9954    10770
chr1    10824   11675
chr1    12159   12418
chr1    15615   15856
chr1    16051   16580
chr1    17450   17673
chr1    19049   19350
chr1    20723   21655
chr1    21756   22612
chr1    23411   23941


```

And to check length of file (number of peaks)

```bash

 wc -l filtered_peaks_final.txt

```


Finally, make sure the peak file is in dos format


```bash

dos2unix filtered_peaks_final.txt

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
--fragments=/rds/projects/my_path/aggr/all_aggr_mapjag/outs/fragments.tsv.gz \
--peaks=/rds/projects/my_path/aggr/filtered_peaks_final.narrow.txt


```

10. Now ready for analysis in R

https://github.com/chrismahony/processing_piplines_bluebear/tree/main/analyse_scATACseq


