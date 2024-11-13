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


cellgeoms_baysor<-function(segfile){
            transcriptspercell<-furrr::future_map_dfr(.x = unique(segfile$cell), 
                                                  .f = ~ data.frame(
                                                      cell = .x, 
                                                      num_transcripts = sum(segfile$cell == .x)
                                                  ), 
                                                  .options = furrr_options(seed = TRUE)
        )
        cellidx <- transcriptspercell$cell[transcriptspercell$num_transcripts > 5]
        segfile.new <- furrr::future_map_dfr(.x = cellidx, function(.x) { 
            res <- st_as_sf(segfile[segfile$cell == .x, c('x', 'y')], coords = c('x', 'y')) %>%
                st_union() %>% #dont remove the union. It is needed here.
                ct_triangulate()
            resdf <- data.frame(cell = .x, geometry = res)
            return(resdf)
        }, .options = furrr_options(seed = TRUE))
        
        cellgeoms_final<-segfile.new$geometry %>% 
            furrr::future_map(purrr::reduce, st_union, .options = furrr_options(seed = TRUE)) %>%
            st_sfc() %>%
            as.data.frame()
        
        cellgeoms_final<-cellgeoms_final %>%
            cbind(transcriptspercell[transcriptspercell$cell %in% cellidx, ])
        
        return(cellgeoms_final)
        
    
    
}




calculate_area_and_centroid <- function(df) {
  # Calculate the area of each geometry
  df$area <- st_area(df$geometry)
  
  # Calculate the centroid of each geometry
  centroids <- st_centroid(df$geometry)
  
  # Extract X and Y coordinates of the centroid
  centroid_coords <- st_coordinates(centroids)
  df$centroid_x <- centroid_coords[, 1]
  df$centroid_y <- centroid_coords[, 2]
  
  # Return the modified data frame with area and centroid columns
  return(df)
}




```

4. Darw cells, this will take several minutes for each FOV. For a full dataset, you will need to run through sbatch 

```R

geoms_all <- list()


for (i in 1:length(data)){
  
#Need a column called 'x' and 'y' and 'cell' to draw cells
data[[i]]$cell <- data[[i]]$cell_id_new
#data[[i]]$x <- data[[i]]$x_location
#data[[i]]$y <- data[[i]]$y_location

#draw cells
geoms_all[[i]] <- cellgeoms_baysor(data[[i]])

#add area and num txs
geoms_all[[i]] <- calculate_area_and_centroid(geoms_all[[i]])
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
