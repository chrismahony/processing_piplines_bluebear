# Analysis of scATACseq data in R



1. Ensure you have fully processed data including calling peaks and filtering low quality peaks

2. Next, install sporkforlife if you have not already done so as well as loading required packages

```R
library(Signac)
library(Seurat)
devtools::install_github("chrismahony/sprokforlife") #install if required
library(sporkforlife)

```



3. With sporkforlife you can process you object in a single wrapper function, e.g.:

```R

aggr <- process_scATAC_data(
  h5_file = "/rds/projects/c/croftap-mapjagb11/scATACseq/aggr/all_aggr_mapjag/outs/filtered_peak_bc_matrix.h5",
  metadata_file = "/rds/projects/c/croftap-mapjagb11/scATACseq/aggr/all_aggr_mapjag/outs/singlecell.csv",
  fragments_file = "/rds/projects/c/croftap-mapjagb11/scATACseq/aggr/all_aggr_mapjag/outs/fragments.tsv.gz",
  libraries_map_file = "/rds/projects/c/croftap-mapjagb11/scATACseq/aggr/libraries_mapjag.csv"
)


```


Note the other possible paramters for this function, that can be adjusted as required, inparticular pay attention to the genome used and adjust to mouse if required:


```R

process_scATAC_data <- function(h5_file, 
                                metadata_file, 
                                fragments_file, 
                                libraries_map_file, 
                                annotation_db = EnsDb.Hsapiens.v86, 
                                min_cells = 10, 
                                min_features = 200, 
                                tss_threshold = 3, 
                                nucleosome_threshold = 4, 
                                ncount_min = 3000, 
                                ncount_max = 50000, 
                                pct_reads_min = 15, 
                                tss_min = 2, 
                                resolutions = c(0.01, 0.05, 0.1, 0.2))


```



4. Processing mouse data would look like:


```R

library(EnsDb.Mmusculus.v79)


aggr <- process_scATAC_data(
  h5_file = "/rds/projects/c/croftap-mapjagb11/scATACseq/aggr/all_aggr_mapjag/outs/filtered_peak_bc_matrix.h5",
  metadata_file = "/rds/projects/c/croftap-mapjagb11/scATACseq/aggr/all_aggr_mapjag/outs/singlecell.csv",
  fragments_file = "/rds/projects/c/croftap-mapjagb11/scATACseq/aggr/all_aggr_mapjag/outs/fragments.tsv.gz",
  libraries_map_file = "/rds/projects/c/croftap-mapjagb11/scATACseq/aggr/libraries_mapjag.csv", annotation_db = EnsDb.Mmusculus.v79)
)


```



5. Next step would probally be to use a scRNAseq refeence to inform your clustering:


```R
#load in data
load("/my_path/mt10_20-nc-protein-10000.rdata")

# Subset (if required)
Idents(PBMC1) <- 'TYPE'
tissue <- subset(PBMC1, ident="Tissue")
tissue <- tissue %>% ScaleData()

transfer.anchors <- FindTransferAnchors(
  reference = tissue,
  query = aggr,
  reduction = 'cca',
  features = aggr@assays[["RNA"]]@var.features
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = tissue$global1,
  weight.reduction = aggr[['lsi']],
  dims = 2:30
)

aggr <- AddMetaData(object = aggr, metadata = predicted.labels)

# Assing to appropritate meta.data slot
aggr$global1 <- aggr$predicted.id


# Repeat (if required) for futher annotaitons
transfer.anchors <- FindTransferAnchors(
  reference = tissue,
  query = aggr,
  reduction = 'cca',
  features = aggr@assays[["RNA"]]@var.features
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = tissue$clusters2312,
  weight.reduction = aggr[['lsi']],
  dims = 2:30
)

aggr <- AddMetaData(object = aggr, metadata = predicted.labels)
aggr$clusters2312 <- aggr$predicted.id



```
