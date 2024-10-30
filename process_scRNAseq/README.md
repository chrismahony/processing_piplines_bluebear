# Processing scRNAseq data

This is an approch to process scRNAseq data using cellrnager. 
<br>

The script will find all fastq files in a folder, extract their names an dprocess them with CellRanger (10x Genomics). 
<br>


1. Complete this step: https://github.com/chrismahony/processing_piplines_bluebear/tree/main/dowloand_data_novogene
 <br>

2. Then make a directory call something like 'all_fastqs'. Put all the fastq files from all samples that you want to process.
   Do not make sub directories, just put all files in the one dir 'all_fastqs'
<br>

3. Make a directory called something like 'count', so your project folder would look like:

```bash

/croftap-XXXX/
|-- my_project
| |--fastqs #downloaded fatsqs from novogene
| |--all_fastqs   #where you put all fastqs
| |--count   #where you will run the next steps and perform alingment

```

<br>
4. The copy the script below and save as something like: count.txt
<br>

!Important!
<br>
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

Note: this scripe will also run cellrnager aggr, so that you have a final loupe file for interactive analysis to send to people AFTER analysing the data in R/Python and then importaing meta.data into Loupe

Script:

```bash

#!/bin/bash
#SBATCH -n 70                            # Number of cores
#SBATCH --mem=399G                       # Total memory
#SBATCH --time=99:0:0                    # Max runtime
#SBATCH --mail-type=ALL
#SBATCH --account=croftap-XXX

set -e
module purge; module load bluebear
module load CellRanger/7.0.0

# Set paths
FASTQ_DIR="/rds/my_path/all_fastqs"  
OUTPUT_DIR="/rds/my_path/count/"      
REF_PATH="/rds/bear-apps/apps-data/CellRanger/refdata-cellranger-mm10-2.1.0"     

# Generate aggregation CSV
AGGR_CSV="$OUTPUT_DIR/aggregation.csv"
echo "sample_id,molecule_h5" > $AGGR_CSV

# Extract unique sample names (removing path and everything after the first underscore)
for sample_name in $(ls $FASTQ_DIR/*.fastq.gz | sed 's|.*/||; s|_\(.*\)||' | sort | uniq); do
    
    # Gather unique sample names related to the current sample
    sample_names=$(ls $FASTQ_DIR | grep "^${sample_name%_*}_" | sed 's|_S.*||' | sort | uniq | tr '\n' ',' | sed 's|,$||')

    sample_output_dir="$OUTPUT_DIR/$sample_name"

    # Run CellRanger count for the sample
    echo "Running cellranger count for sample: $sample_name"
    cellranger count --id=$sample_name \
                     --transcriptome=$REF_PATH \
                     --fastqs=$FASTQ_DIR \
                     --sample=$sample_names \
                     --localcores=8 \
                     --localmem=64

    # Add the sample and molecule_info.h5 path to the aggregation CSV
    molecule_h5="$sample_output_dir/outs/molecule_info.h5"
    echo "$sample_name,$molecule_h5" >> $AGGR_CSV
done

# Once all samples are processed, run CellRanger aggr
echo "Running cellranger aggr"
cellranger aggr --id=aggregated_samples \
                --csv=$AGGR_CSV

```

<br>
8. Open a terminal, naviate to the directory with count.txt
<br>

9. Submit:

```bash

sbatch count.txt

```

