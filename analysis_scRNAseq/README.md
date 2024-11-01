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
library(sporkforlife)


data.10x = list()
dirs <- list.dirs("/rds/my_path/count", recursive = FALSE)  #path to where cellrnager fisished the count step
dirs <- dirs[-1] #remove the aggregated folder, would need to modify this step if you did not do it
dirs <- paste(dirs, "/outs/filtered_feature_bc_matrix/", sep="")

sample_names <- sub('/rds/my_path/count', '', dirs)
sample_names <- sub('/outs/filtered_feature_bc_matrix/', '', sample_names)


for (i in 1:length(dirs)) {
    data.10x[[i]] <- Read10X(data.dir = dirs[[i]])
}

names(data.10x) <- sample_names 
  
scrna.list = list()
for (i in 1:length(data.10x)) {
    print(paste("Processing sample:", sample_names[i]))
    scrna.list[[i]] = CreateSeuratObject(counts = data.10x[[i]], min.cells=3, min.features=200, project=sample_names[i]);
    scrna.list[[i]][["percent.mt"]] = PercentageFeatureSet(object=scrna.list[[i]], pattern = "^MT-");
    scrna.list[[i]] =subset(scrna.list[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 10)
    scrna.list[[i]] =NormalizeData(object = scrna.list[[i]]);
    scrna.list[[i]] =ScaleData(object = scrna.list[[i]]);
    scrna.list[[i]] =FindVariableFeatures(object = scrna.list[[i]]);
    scrna.list[[i]] =RunPCA(object = scrna.list[[i]], verbose = FALSE);     
   }


rm(data.10x)
names(scrna.list) <- sample_names
anchors <- FindIntegrationAnchors(object.list = scrna.list, dims = 1:50)
aggr <- IntegrateData(anchorset = anchors, dims = 1:50)
aggr <- FindVariableFeatures(aggr)
aggr <- ScaleData(aggr, verbose = FALSE)
aggr <- RunPCA(aggr, verbose = FALSE)
aggr <- RunUMAP(aggr, dims = 1:50)
aggr <- FindNeighbors(aggr, dims = 1:30)
aggr <- FindClusters(aggr, resolution = c(0.01, 0.05, 0.1, 0.2, 0.3), graph.name = 'integrated_snn')

Idents(aggr)<-'integrated_snn_res.0.01'  
res0.01markers<-FindAllMarkers(aggr, only.pos = T)

Idents(aggr)<-'integrated_snn_res.0.05'
res0.05markers<-FindAllMarkers(aggr, only.pos = T)

Idents(aggr)<-'integrated_snn_res.0.1'
res0.1markers<-FindAllMarkers(aggr, only.pos = T)

Idents(aggr)<-'integrated_snn_res.0.2'
res0.2markers<-FindAllMarkers(aggr, only.pos = T)

Idents(aggr)<-'integrated_snn_res.0.3'
res0.3markers<-FindAllMarkers(aggr, only.pos = T)

save.image("/rds/my_path/count/analysis_init.RData")


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
