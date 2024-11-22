# Analysis of spatial niches  

Goal here is to identify which cells are in contact with other cells, then colapse counts from neibouring cells to cell in question. i.e. cell idenity is informed by the expression of genes of surounding cells


1.  Load required packages and function

```R
library(spatula)
library(data.table)
library(sfdct)
library(sf)
library(rapportools)
library(purrr)
library(dplyr)
library(spatstat.geom)
library(igraph)
library(Matrix)

source("/rds/projects/c/croftap-mapjagx1/analysis/niches_int/repeat_clean/functions.txt")  # Your path to functions.txt, this file is found in this repo with all functions


```




2.  Read in the tx files. Here I am using outputs from baysor, but default Xenium outputs will also work

```R

all_csv <- dir(path="/rds/projects/c/croftap-mapjagx1/analysis/analysis_initial_setup/baysor_outs", pattern=".csv")
cell_stat_csvs <- dir(path="/rds/projects/c/croftap-mapjagx1/analysis/analysis_initial_setup/baysor_outs", pattern="stats.csv")
tx_csv <- all_csv[!all_csv %in% cell_stat_csvs]
samples<-paste("/rds/projects/c/croftap-mapjagx1/analysis/analysis_initial_setup/baysor_outs/", tx_csv, sep="")
rm(cell_stat_csvs)
rm(all_csv)


samplefov_list <- gsub(".csv", "",tx_csv)

#length(tfiles_list) prev
txfile <- list()
txfile <- lapply(seq_along(samples), function(i) {
  fread(tfiles_list[[i]])[, 
    .(x, y, z_location, gene, cell, SampleID = donor, FOV = fov_name)
  ][, 
    SampleFOV := paste(SampleID, FOV, sep = "_")
  ][, 
    `:=`(x = x * 1, y = y * 1) # Multiplication by 1 is redundant unless coercing
  ] %>% 
    sf::st_as_sf(coords = c("x", "y"), remove = FALSE) %>% 
    setDT()
})


names <- gsub("/rds/projects/c/croftap-mapjagx1/analysis/analysis_initial_setup/baysor_outs", "", samples)
names <- gsub(".csv", "", names)
names <- gsub("/", "", names)
names(txfile) <- names


```

3. Next need a list of meta.data for cells in each FOV. This is extracted from my xenium seurat object in $orig.ident

```R

meta_list <- list()
for (i in 1:length(unique(xenium_chrissy@meta.data$orig.ident))){
   meta_list[[i]] <-   xenium_chrissy@meta.data %>% filter(orig.ident==unique(xenium_chrissy@meta.data$orig.ident)[[i]])
   geoms_f_f2 <- geoms_f_f[geoms_f_f$cell %in% rownames(meta_list[[i]]),]
  meta_list[[i]]$cell_centroid <- st_centroid(geoms_f_f2$geometry)
   meta_list[[i]]$cellID <- rownames(meta_list[[i]])
   }

names(meta_list) <- unique(xenium_chrissy@meta.data$orig.ident)


4. There wil be cells in the txfile that did not pass QC thresholds in the Seurat processing. These need to be reomoved from the txfile. Also from the geoms file (a list of each FOV with geometries for each cell boundry).

```R

#make sure cells that appear in all files are the same (in this case a subset of data has been taken)
txfile <- txfile[names(txfile) %in% names(meta_list)]

for (i in 1:length(txfile)){
txfile[[i]] <- txfile[[i]][txfile[[i]]$cell %in% meta_list[[i]]$cellID,]
geoms_f[[i]] <- geoms_f[[i]][geoms_f[[i]]$cell %in% txfile[[i]]$cell,]

}


```

5. Next need to make sure that the correct columns exist in all the data

```R

for (i in 1:length(txfile)){
  txfile[[i]]$z <- txfile[[i]]$z_location
        meta_list[[i]]$cell_center_x <- meta_list[[i]]$Pos_X
    meta_list[[i]]$cell_center_y <- meta_list[[i]]$Pos_Y
        meta_list[[i]]$bbox_tx <- meta_list[[i]]$cell_centroid
                meta_list[[i]]$SampleFOV <- meta_list[[i]]$orig.ident
                    meta_list[[i]]$SampleID <- meta_list[[i]]$donor
                              meta_list[[i]]$geometry    <-   geoms_f[[i]]$geometry
}

```

At the end of this you need to have the following:

```R

geoms_f[[1]] %>% head()
                         geometry           cell num_transcripts named_new clustered_neighbors
2  POLYGON ((2152.097 653.7633...  CR9745ee3ab-3             111   Myeloid             Myeloid
5  POLYGON ((2146.086 699.092,...  CR9745ee3ab-5             109   Myeloid             Myeloid
8  POLYGON ((2179.76 697.0851,... CR9745ee3ab-12              23   Myeloid             Myeloid
9  POLYGON ((2238.514 729.8842... CR9745ee3ab-39              60   Myeloid             Myeloid
10 POLYGON ((2159.99 677.6901,... CR9745ee3ab-72              27  B/Plasma             Myeloid
12 POLYGON ((2155.514 729.2162... CR9745ee3ab-15              22  B/Plasma             Myeloid



meta_list[[1]] %>% head()
                orig.ident nCount_RNA nFeature_RNA RNA_snn_res.0.05 RNA_snn_res.0.1 RNA_snn_res.0.2 RNA_snn_res.0.3
CR9745ee3ab-11 D2518_11_M8         61           38                1               1               1               1
CR9745ee3ab-12 D2518_11_M8         23           21                6               7               7               8
CR9745ee3ab-14 D2518_11_M8         38           27                0               1               1               7
CR9745ee3ab-15 D2518_11_M8         22           21                0               3               3               3
CR9745ee3ab-17 D2518_11_M8         62           32                1               1               1               7
CR9745ee3ab-20 D2518_11_M8         41           28                1               1               1               7
               RNA_snn_res.0.4 seurat_clusters RNA_snn_res.0.6 RNA_snn_res.0.7 RNA_snn_res.0.8    named RNA_snn_res.0.13
CR9745ee3ab-11               1               1               8               7               8     macs                1
CR9745ee3ab-12               9               7              10              11              14 CD1C_DCs                7
CR9745ee3ab-14               7               1               7               9              10     fibs                1
CR9745ee3ab-15               3               3               3               2               2   Bcells                3
CR9745ee3ab-17               7               1               7               9              10     macs                1
CR9745ee3ab-20               7               1               7               9              10     macs                1
               RNA_snn_res.0.15 RNA_snn_res.0.5      global13        named_sub_new         named_sub_new_named
CR9745ee3ab-11                1               1       Myeloid    SPP1+ Macrophages      SPP1+ Macrophages_macs
CR9745ee3ab-12                8               8 B cells/ pDCs                 pDCs               pDCs_CD1C_DCs
CR9745ee3ab-14                1               7       Myeloid                 macs                   macs_fibs
CR9745ee3ab-15                3               3  Plasma cells Plasma cells/Myeloid Plasma cells/Myeloid_Bcells
CR9745ee3ab-17                1               7       Myeloid    S1008A+ Monocytes      S1008A+ Monocytes_macs
CR9745ee3ab-20                1               7       Myeloid    S1008A+ Monocytes      S1008A+ Monocytes_macs
                     named_sub_new2 global_final    donor clustered_neighbors myeloid_tissue1    Pos_Y    Pos_X
CR9745ee3ab-11    SPP1+ Macrophages      Myeloid D2518_11             Myeloid    0.0251628949 591.1937 2273.561
CR9745ee3ab-12                 pDCs      Myeloid D2518_11             Myeloid   -0.0009651733 667.5464 2179.552
CR9745ee3ab-14                 macs         Fibs D2518_11             Myeloid   -0.0573312406 577.0861 2307.277
CR9745ee3ab-15 Plasma cells/Myeloid     B/Plasma D2518_11             Myeloid   -0.0410209066 723.5449 2155.716
CR9745ee3ab-17    S1008A+ Monocytes      Myeloid D2518_11             Myeloid   -0.0301220307 609.0027 2293.943
CR9745ee3ab-20    S1008A+ Monocytes      Myeloid D2518_11             Myeloid   -0.0944189068 580.1567 2336.530
                           cell_centroid         cellID cell_center_x cell_center_y                   bbox_tx   SampleFOV
CR9745ee3ab-11  POINT (2152.19 666.4005) CR9745ee3ab-11      2273.561      591.1937  POINT (2152.19 666.4005) D2518_11_M8
CR9745ee3ab-12 POINT (2142.048 689.0993) CR9745ee3ab-12      2179.552      667.5464 POINT (2142.048 689.0993) D2518_11_M8
CR9745ee3ab-14 POINT (2187.308 682.6872) CR9745ee3ab-14      2307.277      577.0861 POINT (2187.308 682.6872) D2518_11_M8
CR9745ee3ab-15 POINT (2244.506 735.8092) CR9745ee3ab-15      2155.716      723.5449 POINT (2244.506 735.8092) D2518_11_M8
CR9745ee3ab-17 POINT (2166.887 685.7883) CR9745ee3ab-17      2293.943      609.0027 POINT (2166.887 685.7883) D2518_11_M8
CR9745ee3ab-20 POINT (2161.721 728.9233) CR9745ee3ab-20      2336.530      580.1567 POINT (2161.721 728.9233) D2518_11_M8
               SampleID                       geometry
CR9745ee3ab-11 D2518_11 POLYGON ((2152.097 653.7633...
CR9745ee3ab-12 D2518_11 POLYGON ((2146.086 699.092,...
CR9745ee3ab-14 D2518_11 POLYGON ((2179.76 697.0851,...
CR9745ee3ab-15 D2518_11 POLYGON ((2238.514 729.8842...
CR9745ee3ab-17 D2518_11 POLYGON ((2159.99 677.6901,...
CR9745ee3ab-20 D2518_11 POLYGON ((2155.514 729.2162...



txfile[[1]] %>% head()
          x        y z_location    gene          cell SampleID FOV   SampleFOV                  geometry        z
1: 2153.843 664.6056   21.72517 CLEC10A CR9745ee3ab-3 D2518_11  M8 D2518_11_M8 POINT (2153.843 664.6056) 21.72517
2: 2145.275 686.4321   21.61662    VCAN CR9745ee3ab-5 D2518_11  M8 D2518_11_M8 POINT (2145.275 686.4321) 21.61662
3: 2146.390 687.2050   21.78013    VCAN CR9745ee3ab-5 D2518_11  M8 D2518_11_M8   POINT (2146.39 687.205) 21.78013
4: 2146.746 685.7428   21.42938  ADAM28 CR9745ee3ab-5 D2518_11  M8 D2518_11_M8 POINT (2146.746 685.7428) 21.42938
5: 2147.210 686.8466   21.87589    FGL2 CR9745ee3ab-5 D2518_11  M8 D2518_11_M8  POINT (2147.21 686.8466) 21.87589
6: 2148.957 685.1083   22.75763    FGL2 CR9745ee3ab-5 D2518_11  M8 D2518_11_M8 POINT (2148.957 685.1083) 22.75763

```

6. Next Run the follwing to process the data for spatial niche analysis

```R

ver <- process_voronoi_data(txfile, meta_list)
glass_gridded_ls <- process_voronoi_and_glass(ver, txfile)
res <- process_voronoi_and_gcmat(ver, txfile, glass_gridded_ls$cellgeoms_fov)
obj_ksweep <- process_voronoi_and_seurat(res, max_k = 10)  #adjust max_k as required

```


7. Now we can proceed to processing counts and harmonising for sample and FOV. Where K represent n neighbouring cells. This can be incireased in step 6.

```R

library(rlang)
library(harmony)
gc()
objH <- list()
for (i in 1:length(obj_ksweep)){

obj_ksweep[[i]]$metadata<-obj_ksweep[[i]]$metadata %>% 
    as.data.frame %>%
    dplyr::mutate(cellID = tileID)
objH[[i]] <- dimred_and_cluster(
  obj_ksweep[[i]], 
  do_harmony = TRUE, 
  wts = "wts",
  vars_use = c("SampleFOV", "SampleID"),
  resolution_clustering = c(0.1),
  theta = c(0,0),
  sigma = 0.2, 
  max.iter.harmony = 12,
  max.iter.cluster = 40,
  do_QC = TRUE,
  do_cluster_after = TRUE,
  do_cluster_before = TRUE,
  return_object = TRUE,
  do_umap_after = TRUE,
  do_umap_before = TRUE
)
}

```


8. UMAP plotting and checking clustering

```R

library(ggplot2)
library(ArchR)

for (i in 1:10){
print(cbind(objH[[i]]$Humap$embedding %>% as.data.frame(),objH[[i]]$metadata) %>% ggplot(aes(x=V1, y=V2, color=SampleFOV))+ geom_point(shape=".")+theme_ArchR()+NoLegend()+ggtitle(paste("itteration", i)))
}

```



9. Futher clustering and plotting

```R

bjH[[10]]$Humap$clusters <- RunModularityClustering(
      SNN = objH[[10]]$Humap$fgraph, 
      resolution = 0.2,
      print.output = FALSE
    )

# Further plotting
clusters <- objH[[10]]$Humap$clusters %>% as.data.frame()
colnames(clusters) <- "clust1"
clusters$clust1 %>% table()
objH[[10]]$metadata$clust1 <- as.character(clusters$clust1)


cbind(objH[[10]]$Humap$embedding %>% as.data.frame(),objH[[10]]$metadata) %>% ggplot(aes(x=V1, y=V2, color=clust1))+ geom_point(shape=".")+theme_ArchR()+NoLegend()

```


10. Check that the spatial niches represent something simialr to the tissue morphology. Plot a single sample and focus in on a tissue fragment.

```

objH[[10]]$metadata %>% as.data.frame() %>% filter(SampleID == "D2518_24") %>%  ggplot() + 
  geom_sf(aes(geometry = polygon_centroid, color=clust1), alpha = 0.7) + 
            ggtitle("Tissue regions")+facet_wrap("SampleID", ncol=4) +xlim(900,3600)+ylim(5500,8200)
```
![Example plot](https://github.com/chrismahony/processing_piplines_bluebear/edit/main/process_Xenium_CosMx_R/niches/Figure_gh copy.pdf)


11. Now add Cluster ID onto the geoms object and examine

```R

#filter geoms (some FOVs removed due to low numhber of cells so impossible to predict neighbours)
geoms_final <- geoms_f[names(geoms_f) %in% names(res)]

#extract meta.data and split into a list
meta_new <- objH[[10]]$metadata
meta_list <- list()
for (i in 1:length(unique(objH[[10]]$metadata$SampleFOV))){
  
  meta_list[[i]] <- meta_new %>% filter(SampleFOV == unique(objH[[10]]$metadata$SampleFOV)[[i]])
  
}

names(meta_list) <- names(res)


# Match cells on geoms obj
for (i in 1:length(geoms_final)){
   index <- match(geoms_final[[i]]$cell, meta_list[[i]]$cellID)
geoms_final[[i]]$clust_new <- meta_list[[i]]$clust1[index]
}


# Plotting
cols <- ArchR::paletteDiscrete(meta_new$clust1) %>% as.data.frame()
geoms_final[["D2518_36_Z3"]] %>% as.data.frame %>% 
    ggplot()+
     geom_sf(aes(geometry = geometry, fill = clust_new), alpha = 0.7,
              color = "black")+ scale_fill_manual(values = cols$.) +theme(panel.background = element_rect(fill='white'))+
  theme(axis.line = element_line(color='black'),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank())#+xlim(2100,2300)+ylim(4370,4550)

```



12. Its possible that your will need to subset out a cluster and subcluster to refine analysis

13. Add details to Seurat object and caclulate porportion of cells in each nich

```R






```


