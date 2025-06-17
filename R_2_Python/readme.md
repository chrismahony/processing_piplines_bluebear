## Easily write out that data you need to convert Anndata obj to Seurat


1. Load what you need

```python

import scanpy as sc
import anndata
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import numpy as np
import os
import pandas as pd

amp2 = sc.read("/rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/AMP2_data/AMPII_annotated.h5ad")


```

2. Save the meta data from your anndata obj

```pytyhon

adata_fibs = amp2[amp2.obs['cell_type_AMP2021'].isin(['Fibroblast'])]
meta_fibs=adata_fibs.obs
meta_fibs_df = pd.DataFrame(meta_fibs) 
meta_fibs_df.to_csv('/rds/path/AMP2_data/meta_data_fibs.csv') 

```


3. Delete everything from your anndata obj and save the obj

```python

del adata_fibs.var
del adata_fibs.obs
del adata_fibs.uns
del adata_fibs.obsm
del adata_fibs.layers
results_file = '//rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/AMP2_data/adata_fibs_no_meta.h5ad'  # the file that will store the analysis results
adata_fibs.write(results_file)

```

4. Now in R, convert the object, read it in and generate the Seurat object

```R

library(Seurat)
library(SeuratDisk)
options(bitmapType='cairo')
Convert("/rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/AMP2_data/adata_Myeloid_no_meta.h5ad", dest = "h5seurat", overwrite = F)
amp2_Myeloid <- LoadH5Seurat("/rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/AMP2_data/adata_Myeloid_no_meta.h5seurat")

meta_data <- read_csv("/rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/AMP2_data/meta_Myeloid_df.csv")

meta_data <-as.data.frame(meta_data)
library(tidyverse)
meta_data<- meta_data %>% remove_rownames %>% column_to_rownames(var="X1")
amp2_Myeloid<-AddMetaData(amp2_Myeloid, meta_data)
amp2_Myeloid <- amp2_Myeloid %>%
    NormalizeData() %>%
    ScaleData() %>%
    FindVariableFeatures() %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30, verbose = FALSE)

DimPlot(amp2_Myeloid, group.by = "cell_lineage")

```



