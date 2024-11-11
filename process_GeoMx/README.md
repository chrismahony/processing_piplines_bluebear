# Processing GeoMx data

1. Dowonload data to an RDS folder

2. Are the files BCL files (i.e. they are not fastqs)?

If so then you will need to run a script like this to process them to fastqs

```bash

#!/bin/bash
#SBATCH -J BCL_convert			
#SBATCH -n 8                                          
#SBATCH -N 1                                           
#SBATCH --mem 64G                                    
#SBATCH -t 0-05:00                                     
#SBATCH --mail-type=FAIL,END                           
#SBATCH --account=mcmurraj-bta-geomx-data-storage       

set -e
module purge; module load bluebear 
module load bear-apps/2022a
module load bcl-convert/4.0.3-2el7.x86_64

#### ENSURE THAT NO-LANE-SPLITTING IS SET TO FALSE IN THE SAMPLE SHEET ####
bcl-convert --bcl-input-directory /rds/path/BCL_data \
--output-directory /rds/patyh/FASTQs \
--sample-sheet /rds/path/SampleSheet.csv #example uploaded in this folder

```




3. Once this is completed then you need to aling reads to ref genome

Script like this should work:


```bash

#!/bin/bash
#SBATCH -J BCL_convert			
#SBATCH -n 8                                          
#SBATCH -N 1                                           
#SBATCH --mem 64G                                    
#SBATCH -t 0-05:00                                     
#SBATCH --mail-type=FAIL,END                           
#SBATCH --account=mcmurraj-bta-geomx-data-storage       

set -e
module purge; module load bluebear 
module load bear-apps/2021b
module load GeoMxNGSPipeline/2.3.3.10 # Required package to run script

## Run the actual script ========================================================
### CHANGE THE NUMBER OF THREADS IN THE CONFIG.INI FILE TO 8
geomxngspipeline --in=/rds/path_to_fastqs/GB500923-AH_FASTQ \
--out=/rds/where_you_want_output/00_DCC-files \
--ini=/rds/where_this_ini_file_is_stored/Annie_WTA_20231207T1039_GNP_config.ini \  # example uploaded, need editing based on names of wells and 'scan widths' 
--save-interim-files=true \
--threads=8



```


