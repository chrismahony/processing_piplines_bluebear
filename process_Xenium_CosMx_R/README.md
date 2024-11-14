# Processing Xenium/CosmX data using R


## Generate Count matrices and draw cell boundaries


1. Assume you have a dir with multiple .csvs to read in, each FOV is a tx count table for a FOV


```R

all_csv <- dir(path="/rds/projects/c/croftap-mapjagx1/analysis/analysis_initial_setup/baysor_outs", pattern=".csv")


#remove irrelevant files from list if needed
cell_stat_csvs <- dir(path="/rds/projects/c/croftap-mapjagx1/analysis/analysis_initial_setup/baysor_outs", pattern="stats.csv")
tx_csv <- all_csv[!all_csv %in% cell_stat_csvs]
samples<-paste("/rds/projects/c/croftap-mapjagx1/analysis/analysis_initial_setup/baysor_outs/", tx_csv, sep="")



#read in csvs and combine
data = list()
for (i in 1:10) {
    data[[i]] <- read.csv(samples[[i]]);
}


library(data.table)
full_tx_table<-rbindlist(data)


#filter out BLANKS, unassinge tx and Neg probes at this stage


# I created this column for downstream analysis
full_tx_table <- filter(all_tx_f, donor_FOV== unique(all_tx_f$donor_FOV)[[i]])


```

Take care if you have .csvs across mutiple donors



2. You might have a single csv tx table with all detected txs, if so read in and split based on FOV


```R

data = list()
for (i in 1:length(unique(full_tx_table$donor_FOV))) {
    data[[i]] <- full_tx_table %>% filter(donor_FOV == unique(full_tx_table$donor_FOV)[[i]])
}

names(data) <- unique(full_tx_table$donor_FOV)

```

Take care if you have a tx table across mutiple donors


3. Define these functions to draw cells:

```R

library(furrr)
library(sfdct)
library(sf)


cellgeoms_with_area_centroid_slow <- function(segfile) {
  # Step 1: Calculate transcripts per cell
  transcriptspercell <- furrr::future_map_dfr(
    .x = unique(segfile$cell),
    .f = ~ data.frame(
      cell = .x,
      num_transcripts = sum(segfile$cell == .x)
    ),
    .options = furrr_options(seed = TRUE)
  )

  # Step 2: Filter cells with more than 5 transcripts
  cellidx <- transcriptspercell$cell[transcriptspercell$num_transcripts > 5]

  # Step 3: Calculate geometry for each cell using triangulation
  segfile.new <- furrr::future_map_dfr(
    .x = cellidx,
    .f = function(.x) {
      res <- st_as_sf(segfile[segfile$cell == .x, c('x', 'y')], coords = c('x', 'y')) %>%
        st_union() %>%  # st_union is needed here
        ct_triangulate()
      resdf <- data.frame(cell = .x, geometry = res)
      return(resdf)
    },
    .options = furrr_options(seed = TRUE)
  )

  # Step 4: Finalize cell geometries by merging geometries for each cell
  cellgeoms_final <- segfile.new$geometry %>%
    furrr::future_map(purrr::reduce, st_union, .options = furrr_options(seed = TRUE)) %>%
    st_sfc() %>%
    as.data.frame()

  # Step 5: Add transcript count information to the geometry data
  cellgeoms_final <- cellgeoms_final %>%
    cbind(transcriptspercell[transcriptspercell$cell %in% cellidx, ])

  # Step 6: Calculate area and centroids for each geometry
  # Calculate the area of each geometry
  cellgeoms_final$area <- st_area(cellgeoms_final$geometry)

  # Calculate the centroid of each geometry
  centroids <- st_centroid(cellgeoms_final$geometry)

  # Extract X and Y coordinates of the centroid
  centroid_coords <- st_coordinates(centroids)
  cellgeoms_final$centroid_x <- centroid_coords[, 1]
  cellgeoms_final$centroid_y <- centroid_coords[, 2]

  # Return the final data frame with area and centroid columns
  return(cellgeoms_final)
}



cellgeoms_with_area_centroid_fast <- function(segfile) {
  # Step 1: Calculate transcripts per cell
  transcriptspercell <- furrr::future_map_dfr(
    .x = unique(segfile$cell),
    .f = ~ data.frame(
      cell = .x,
      num_transcripts = sum(segfile$cell == .x)
    ),
    .options = furrr_options(seed = TRUE)
  )

  # Step 2: Filter cells with more than 5 transcripts
  cellidx <- transcriptspercell$cell[transcriptspercell$num_transcripts > 5]

  # Step 3: Calculate and union geometry for each cell
  segfile.new <- furrr::future_map_dfr(
    .x = cellidx,
    .f = function(.x) {
      cell_data <- segfile[segfile$cell == .x, c('x', 'y')]
      # Create an sf object and triangulate
      triangles <- st_as_sf(cell_data, coords = c('x', 'y')) %>%
        st_combine() %>%
        ct_triangulate() %>%
        st_union()  # Combine triangles into one geometry per cell
      data.frame(cell = .x, geometry = st_geometry(triangles))
    },
    .options = furrr_options(seed = TRUE)
  )

  # Step 4: Add transcript count information to the geometry data
  cellgeoms_final <- dplyr::left_join(segfile.new, transcriptspercell, by = "cell")

  # Step 5: Calculate area and centroids for each geometry
  cellgeoms_final <- cellgeoms_final %>%
    mutate(
      area = st_area(geometry),
      centroid = st_centroid(geometry)
    )

  # Extract X and Y coordinates of the centroid
  centroid_coords <- st_coordinates(cellgeoms_final$centroid)
  cellgeoms_final$centroid_x <- centroid_coords[, 1]
  cellgeoms_final$centroid_y <- centroid_coords[, 2]

  # Remove the temporary centroid column
  cellgeoms_final$centroid <- NULL

  # Return the final data frame with area and centroid columns
  return(cellgeoms_final)
}






```

4. Darw cells, this will take several minutes for each FOV. For a full dataset, you will need to run through sbatch 

```R

geoms_all <- list()


#slow takes roughly 10x longer, but is better for plotting cells. The fast is much faster but will have extra lines going through the cells when you plot du to the tirangulaiton function (i have not found a wa yto get around this yet!)

for (i in 1:length(data)){
  
#Need a column called 'x' and 'y' and 'cell' to draw cells
data[[i]]$cell <- data[[i]]$cell_id
#data[[i]]$x <- data[[i]]$x_location
#data[[i]]$y <- data[[i]]$y_location

#draw cells
geoms_all[[i]] <- cellgeoms_with_area_centroid_slow(data[[i]])

}

```

Now you have a list of cells, boudry plots, area and #txs which can be plotted directly


5. Plotting to check on individual FOVs

```R

geoms_all[[1]] %>% as.data.frame %>% 
    ggplot()+
     geom_sf(aes(geometry = geometry, fill = area), alpha = 0.7,
              color = "black")+theme_minimal()


geoms_all[[1]] %>% as.data.frame %>% 
    ggplot()+
     geom_sf(aes(geometry = geometry, fill = num_transcripts), alpha = 0.7,
)+theme_minimal()



```


6. Plot whole tissue with all cells and look at area and num tx

```R

combined_df <- do.call(rbind, geoms_all)

combined_df %>% as.data.frame %>% 
    ggplot()+
     geom_sf(aes(geometry = geometry, fill = num_transcripts), alpha = 0.7,
)+theme_minimal()


combined_df %>% as.data.frame %>% 
    ggplot()+
     geom_sf(aes(geometry = geometry, fill = area), alpha = 0.7,
              color = "black")+theme_minimal()

```


OPTIONAL (experimental and needs further testing)

Darw neuclei booundaries based on tx found in neuclei and remove cells with no tx in neucli


```R

data_nuc <- list()
geoms_nuc <- list()

for (i in 1:3){
  data[[i]]$cell <- data[[i]]$cell_id
data_nuc[[i]] <- data[[i]] %>% filter(overlaps_nucleus == 1)
  geoms_nuc[[i]] <- cellgeoms_with_area_centroid_fast(data_nuc[[i]])

}


combined_df_nuc <- do.call(rbind, geoms_nuc)

combined_df_nuc %>% as.data.frame %>% 
    ggplot()+
     geom_sf(aes(geometry = geometry, fill = area), alpha = 0.7,
              color = "black")+theme_minimal()


colnames(combined_df_nuc)[colnames(combined_df_nuc) == "geometry"] <- "geometry_nuc"
colnames(combined_df_nuc)[colnames(combined_df_nuc) == "area"] <- "area_nuc"
colnames(combined_df)[colnames(combined_df) == "geometry"] <- "geometry_cell"


df_new <- merge(combined_df_nuc[c(1:2,4)], combined_df, by="cell")


combined_df %>% as.data.frame %>% 
    ggplot()+
     geom_sf(aes(geometry = geometry_cell), alpha = 0.7,
              color = "black")+theme_minimal()


df_new %>% as.data.frame %>% 
ggplot() +
  geom_sf(aes(geometry = geometry_cell), fill = "red", alpha = 0.5, color = "black") +
  geom_sf(aes(geometry = geometry_nuc), fill = "blue") +
  theme_minimal() +
  theme(legend.position = "none") +
  ggtitle("Plotting Two Geometry Columns")+xlim(2400,2800)+ylim(200,600)



```





7. Next make a Seurat object. Need to gof from tx tables to count tables to processed data merged together

First defefine useful functions



```R

make_count_mtx <- function(genes, cells, remove_bg = TRUE) {
    if (remove_bg) {
        idx <- cells != 0  # Logical indexing is faster than `which`
        genes <- genes[idx]
        cells <- cells[idx]
    }

    # Retrieve unique levels once for efficiency
    unique_genes <- unique(genes)
    unique_cells <- unique(cells)
    
    # Map genes and cells to integer indices
    gene_idx <- match(genes, unique_genes)
    cell_idx <- match(cells, unique_cells)
    
    # Create sparse matrix
    counts <- Matrix::sparseMatrix(
        i = gene_idx, 
        j = cell_idx, 
        x = rep(1, length(gene_idx)),  # Each entry is 1
        dims = c(length(unique_genes), length(unique_cells))
    )
    
    # Set row and column names
    rownames(counts) <- unique_genes
    colnames(counts) <- unique_cells
    
    return(counts)
}


process_seurat_data <- function(data, combined_df, num_datasets = 3) {
  # Create empty lists to store counts and Seurat objects
  counts <- list()
  seurats <- list()

  # Loop through each dataset to create count matrices and Seurat objects
  for (i in 1:num_datasets) {
    counts[[i]] <- make_count_mtx(genes = data[[i]]$gene, data[[i]]$cell_id_new)
    seurats[[i]] <- CreateSeuratObject(counts = counts[[i]], project = unique(data[[i]]$donor_FOV))
  }

  # Merge all Seurat objects
  all_merged <- merge(x = seurats[[1]], y = seurats[2:length(seurats)])

  # Add cell area metadata
  cell_area_vector <- setNames(combined_df$area, combined_df$cell)
  all_merged <- AddMetaData(all_merged, metadata = cell_area_vector, col.name = "cell_area")

  # Remove cells with NA in cell_area
  na_cells <- all_merged@meta.data[is.na(all_merged@meta.data$cell_area), ]
  all_merged <- all_merged[, !colnames(all_merged) %in% rownames(na_cells)]

  # Subset based on feature and count thresholds
  all_merged <- subset(all_merged, subset = nFeature_RNA > round(median(all_merged$nFeature_RNA)/2) &
                         nCount_RNA > round(median(all_merged$nCount_RNA)/2) &
                         cell_area > round(median(all_merged$cell_area)/4))

  # Normalize and scale data, find variable features, run PCA and UMAP
  all_merged <- all_merged %>% 
    NormalizeData(scale.factor = median(all_merged$nCount_RNA)) %>% 
    ScaleData() %>% 
    FindVariableFeatures() %>% 
    RunPCA() %>% 
    RunUMAP(dims = 1:40) %>% 
    FindNeighbors(dims=1:40) %>% 
    FindClusters(res=c(0.01, 0.05, 0.1, 0.2))

  # Plot UMAP
  DimPlot(all_merged, raster = FALSE) + NoLegend()

  return(all_merged)
}
```

8. Now process data

```R

# First prepare geoms for next step
for (i in 1:length(geoms_all)){
  geoms_all[[i]]$cell <- paste0(geoms_all[[i]]$cell, "_", i)
  }
combined_df <- do.call(rbind, geoms_all)


#Now ready to process
#this function uses 'donor_FOV' column in your data df
all_merged <- process_seurat_data(data = data, combined_df = combined_df, num_datasets = 3)


```

9. Recluster, examine clusters and assing cell ids as required


10. Now loot to see where clusters are in space

```R

index <- match(combined_df$cell, rownames(all_merged@meta.data))
combined_df$clusters <- all_merged@meta.data$RNA_snn_res.0.2[index]


combined_df %>% as.data.frame %>% 
    ggplot()+
     geom_sf(aes(geometry = geometry, fill = clusters), alpha = 0.7,
              color = "black")+theme_minimal()


```

11. Now plot the expression of a gene of interest

```R

gene_name <- 'SOX2'
gene_expression <- all_merged@assays$RNA@data[gene_name, ]
combined_df$gene_expression <- gene_expression[match(combined_df$cell, names(gene_expression))]

combined_df %>% 
  as.data.frame() %>% 
  ggplot() +
  geom_sf(aes(geometry = geometry, fill = gene_expression), alpha = 0.7, color = "black") +
  scale_fill_viridis_c(option = "C") +  # You can change this scale to adjust color mapping
  theme_minimal() +
  labs(title = paste("Expression of", gene_name))


```


