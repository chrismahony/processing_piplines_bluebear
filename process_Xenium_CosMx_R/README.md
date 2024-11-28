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


3. Load libraris and functions
```R

library(furrr)
library(sfdct)
library(sf)
source("/path/functions_init.txt")  # this can be found here: https://github.com/chrismahony/processing_piplines_bluebear/blob/main/process_Xenium_CosMx_R/functions_init.txt


```

4. Draw cells, this will take several minutes for each FOV. For a full dataset, you will need to run through sbatch 

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
  geoms_nuc[[i]] <- cellgeoms_with_area_centroid_slow(data_nuc[[i]])

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


df_new$nuc_ratio <- df_new$area_nuc/df_new$area
ggplot(df_new, aes(x = nuc_ratio)) +
  geom_histogram(binwidth = 0.01, fill = "skyblue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Column Name", x = "Column Name", y = "Frequency") +
  theme_minimal()+xlim(0,1)



```





7. Next make a Seurat object. Need to gof from tx tables to count tables to processed data merged together

Make sure you have the functions you need loaded



```R

source("/path/functions_init.txt")  # this can be found here: https://github.com/chrismahony/processing_piplines_bluebear/blob/main/process_Xenium_CosMx_R/functions_init.txt
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



12. Further plotting options for default segmentation using spatialdataplot


https://github.com/chrismahony/processing_piplines_bluebear/blob/main/process_Xenium_CosMx_R/Plotting_spatial_data.md



13. Proximity enrichment once cell labels finalised

 https://github.com/chrismahony/processing_piplines_bluebear/blob/main/process_Xenium_CosMx_R/proximity_analysis/readme.md


14. Labeling cells and subtyping

    TBC
