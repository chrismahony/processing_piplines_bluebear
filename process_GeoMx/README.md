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
#SBATCH --qos castles				# Request Castles Node
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
--sample-sheet /rds/path/SampleSheet.csv #exmaple uploaded in this folder

```
