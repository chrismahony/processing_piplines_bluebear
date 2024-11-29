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

names(ATAC_objs) <- c("fibroblasts_ATAC_DEXM_v5", "fibroblasts_ATAC_LPSM_v5")
names(RNA_objs) <- c("fibroblasts_RNA_DEXM_v5", "fibroblasts_RNA_LPSM_v5")

# This wil take some time to run
gene.activities <- GeneActivity(coculture_ATAC_QC_v5)


```


5. Now begin work flow

```R

library(Seurat)
library(scMEGA)
library(ArchR)
library(sporkforlife)


#use wrapper functions from sporkforlife

co_embed_list <- coembed_data_function(ATAC_objs, RNA_objs, gene.activities)
all_coembed_merge <- process_co_embed_list(co_embed_list, RNA_objs, ATAC_objs)

```


6. Now set a trajectory and plot it

```R

obj.pair_all_coembed_merge <- AddTrajectory(object = all_coembed_merge, 
                          trajectory = c("Fb_DexM", "Fb_LPSM"),
                          group.by = "orig.ident", 
                          reduction = "pca",
                          dims = 1:3, 
                          use.all = FALSE)


DimPlot(obj.pair_all_coembed_merge)

obj.pair_all_coembed_merge <- obj.pair_all_coembed_merge[, !is.na(obj.pair_all_coembed_merge$Trajectory)]

TrajectoryPlot(object = obj.pair_all_coembed_merge, 
                    reduction = "umap",
                    continuousSet = "blueYellow",
                    size = 1,
                   addArrow = FALSE)

```


7. Next generate a genome for adding motifs and running chromVar

```R

library(JASPAR2020)
library(TFBSTools)
library(BiocParallel)
register(SerialParam())
library(BSgenome.Hsapiens.UCSC.hg38)

DefaultAssay(obj.pair_all_coembed_merge) <- "peaks"

pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)


keepBSgenomeSequences <- function(genome, seqnames)
{
    stopifnot(all(seqnames %in% seqnames(genome)))
    genome@user_seqnames <- setNames(seqnames, seqnames)
    genome@seqinfo <- genome@seqinfo[seqnames]
    genome
}

sequences_to_keep <- paste0("chr", c(1:22, "X", "Y"))
BSgenome.Hsapiens.UCSC.hg38_CM <- keepBSgenomeSequences(BSgenome.Hsapiens.UCSC.hg38, sequences_to_keep)

DefaultAssay(obj.pair_all_coembed_merge) <- "peaks"

gr <- granges(obj.pair_all_coembed_merge)
seq_keep <- seqnames(gr) %in% seqnames(BSgenome.Hsapiens.UCSC.hg38_CM) 
seq_keep <- as.vector(seq_keep)
feat.keep <- GRangesToString(grange = gr[seq_keep])
obj.pair_all_coembed_merge2 <- obj.pair_all_coembed_merge
obj.pair_all_coembed_merge2[['peaks']] <- subset(obj.pair_all_coembed_merge2[["peaks"]], features = feat.keep)

DefaultAssay(obj.pair_all_coembed_merge2) <- "peaks"



obj.pair_all_coembed_merge2 <- AddMotifs(
  object = obj.pair_all_coembed_merge2,
  genome = BSgenome.Hsapiens.UCSC.hg38_CM,
  pfm = pfm,
    assay = "peaks"
)

obj.pair_all_coembed_merge2 <- RunChromVAR(
  object = obj.pair_all_coembed_merge2,
  genome = BSgenome.Hsapiens.UCSC.hg38_CM,
    assay = "peaks"
)

register(MulticoreParam(40, progressbar = F))





```



