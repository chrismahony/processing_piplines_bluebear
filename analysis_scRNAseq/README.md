# Generic analysis of scRNAseq data

This script will:
  -read in data
  -process using seurat
  -aggregate samples and proces again
  -cluster at various resolutiona and FindAllMarkers at these resolutions

<br>
Then you would need to read the data in and manually annotate etc.

<br>

<br>

1. Complete this step first: https://github.com/chrismahony/processing_piplines_bluebear/tree/main/process_scRNAseq
<br>

2. Next, copy the following script, save it as something like 'R1.txt'

```bash
library(Signac)
library(Seurat)

#devtools::install_github("chrismahony/sprokforlife") #install if required
library(sporkforlife)


#optain a list of dirs to samples and smaple names
data.10x = list()
dirs <- list.dirs("/rds/my_path/count", recursive = FALSE)  #path to where cellrnager fisished the count step
dirs <- dirs[-1] #remove the aggregated folder, would need to modify this step if you did not do it
dirs <- paste(dirs, "/outs/filtered_feature_bc_matrix/", sep="")

sample_names <- sub('/rds/my_path/count', '', dirs)
sample_names <- sub('/outs/filtered_feature_bc_matrix/', '', sample_names)


#this sporkforlife function will process, aggregate samples and cluster to the number of target_n_cluster using Seurat
#once the target number of clusters has been reached then Seurat::FindAllMarkers() will run for that resolution

process_scrna_data(dirs, sample_names, target_n_clusters = 5,
                               resolution_range = seq(0.05, 0.3, by = 0.5),
                               min_nFeature_RNA = 500, max_nFeature_RNA = 7000,
                               max_percent_mt = 10, n_dims=50)

```

<br>

Make sure to amend the relavent paths to files and to save etc.

4. Next to submit this and run through sbatch, you need a second script to activate R and run the script, copy this and save with something like Rsb.txt

```bash

#!/bin/bash
#SBATCH -n 70                            
#SBATCH --mem=399G                      
#SBATCH --time=99:0:0                    
#SBATCH --mail-type=ALL
#SBATCH --account=croftap-XXX   #amend this to your RDS

set -e

module purge; module load bluebear
module load bear-apps/2021a
module load R-bundle-Bioconductor/3.13-foss-2021a-R-4.1.0

Rscript R1.txt

```


5. Open a terminal, navigate to where you save these files and submit:

```bash

sbatch Rsb.txt

```
