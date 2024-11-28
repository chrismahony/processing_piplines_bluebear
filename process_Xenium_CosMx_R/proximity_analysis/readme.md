# Calculate proximity enrichent of cells in x/y space


This uses a df with cell names, x/y corrds for centroids of the cells and a identity for each cell types

```R

library(devtools)
library(spatula)
library(furrr)
source("/my_path/functions_prox.txt")

all_x_y <- rbindlist(data_x_y)
 
#need to add a column to all_x_y with the ident for each cell. in this case mine is called 'named'
 
coloc_res_coarse<- coloc_all_types(
        index_type = unique(all_x_y$named),
        coords = all_x_y[, c("Centroid.X.µm", "Centroid.Y.µm")],
        y = all_x_y$named,
        compartments = NULL,
        max_dist = 40,
        nperm = 1000,
        parallel = TRUE
    )
 
library(data.table)
library(splitstackshape)
library(ggpubr)
library(rstatix)
library(DescTools)
 
plt_df<-coloc_res_coarse %>%
    subset(pval < 0.05) %>%
    dplyr::select(index_type, type, zscore) %>%
    spread(type, zscore, fill = 0) %>%
    column_to_rownames('index_type') %>%
    as.matrix
 
plt_df %>%
    Heatmap(col = circlize::colorRamp2(c(-20, 0, 20), c('blue', 'white', 'red')), border = T, cluster_columns = T, cluster_rows = T)
 
```
