# This is a wlakthrough of generating a GRN from scRNAseq data and scATACseq data

This is from data where the RNAseq and ATACseq has been performed in seperate experiments.

First the scRNAseq/scATACseq needs to be processed seperatly, QC'ed and annotated. So that you have two seperate objects.

1. Setup:

This is using R v4.4.1 and Seurat v5.

Now install scMEGA and dependicies: https://github.com/CostaLab/scMEGA?tab=MIT-2-ov-file. 

First, rownames (genes) in your scRNAseq object must all be in upper case for scMEGA to work if you are owkring with HUMAN data then you can ignore this step). You can update you object like this:

```R

# scRNAseq obj
counts <- GetAssayData(obj.rna, assay = "RNA",slot = "counts")
rownames(counts) <- toupper(rownames(counts))
obj.rna <- CreateSeuratObject(counts = counts, meta.data = obj.rna@meta.data)

# scATACseq obj
gene.activity <- obj.atac@assays$ACTIVITY@counts
rownames(gene.activity) <- toupper(rownames(gene.activity))
counts <- GetAssayData(obj.atac, assay = "ATAC",slot = "counts")
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c("-", "-"),
  genome = 'mm10',
  min.cells = 0,
  min.features = 0
)
obj.atac <- CreateSeuratObject(counts = counts, meta.data = obj.atac@meta.data, assay = "ATAC")


```

2. Next to begin the this workflow we will pair scRNAseq cells with scATACseq cells. You dont want treatment groups to be mixed up.

i.e. you dont want Control scRNAseq cells to be paired with scATACseq cells. So split the RNA and ATAC objectes into treatment groups.

```R

# scRNAseq objs
Idents(fibs) <- "condition"
fibroblasts_ctrl <- subset(fibs, idents = "control")
fibroblasts_LPS <- subset(fibs, idents = "LPS")


fibroblasts_ctrl <- fibroblasts_ctrl %>%
    ScaleData() %>%
    FindVariableFeatures()

fibroblasts_LPS <- fibroblasts_LPS %>%
    ScaleData() %>%
    FindVariableFeatures()

# scATACseq objs
Idents(fibroblasts_ATAC) <- "conditions"
fibroblasts_ATAC_ctrl <- subset(fibroblasts_ATAC, idents = c("control"))
fibroblasts_ATAC_ctrl <- FindTopFeatures(fibroblasts_ATAC_ctrl, min.cutoff = 'q5')
fibroblasts_ATAC_ctrl <- RunSVD(fibroblasts_ATAC_ctrl)

fibroblasts_ATAC_LPS <- subset(fibroblasts_ATAC, idents = "LPS")
fibroblasts_ATAC_LPS <- FindTopFeatures(fibroblasts_ATAC_LPS, min.cutoff = 'q5')
fibroblasts_ATAC_LPS <- RunSVD(fibroblasts_ATAC_LPS)


```

3. This workflow will take Seurat v5 objects to run.. If they are not v5, then update them:

```R

fibroblasts_ctrl <- UpdateSeuratObject(fibroblasts_ctrl)
fibroblasts_LPS <- UpdateSeuratObject(fibroblasts_LPS)

fibroblasts_ATAC_ctrl <- UpdateSeuratObject(fibroblasts_ATAC_ctrl)
fibroblasts_ATAC_LPS <- UpdateSeuratObject(fibroblasts_ATAC_LPS)
fibroblasts_ATAC <- UpdateSeuratObject(fibroblasts_ATAC)

```

4. Next create a list of the ATAC obj and the RNA objs and pull the gene activities from the main obj

```R

ATAC_objs <- list()
RNA_objs <- list()


ATAC_objs[[1]] <- fibroblasts_ATAC_DEXM_v5
ATAC_objs[[2]] <- fibroblasts_ATAC_LPSM_v5


RNA_objs[[1]] <- fibroblasts_RNA_DEXM_v5
RNA_objs[[2]] <- fibroblasts_RNA_LPSM_v5

gene.activities <- GeneActivity(coculture_ATAC_QC_v5)


```


5. Now begin work flow

```R

library(Seurat)
library(scMEGA)
library(ArchR)



```
