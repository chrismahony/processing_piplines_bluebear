# This is an overview of analysing bulk RNAseq


1. Complete this step: https://github.com/chrismahony/processing_piplines_bluebear/tree/main/process_bulk_RNAseq

<br>

2. Start an R session

3. Load the required packages

```R

library(DESeq2)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

 ````


4. Install sprokforlife and load the package. Here you will have several handy functions that act as wrappers to process bulk and scRNAseq data

```R

devtools::install_github("chrismahony/sprokforlife")
library(sporkforlife)

```


5. Read in processed counts file and make a meta data df

```R

all_counts <- read.delim("/rds/my_path/all_counts.txt", row.names=1, comment.char="#")   #read in counts

all_counts <- all_counts %>% dplyr::select(-c(Chr, Start,End, Strand, Length))    #remove non-essential columns

colnames(all_counts) <- c("C3a", "C3b", "C3c", "c5a", "c5b", "c5c", "fba", "fbb", "fbc")  #rename columns

#make meta data from count matrix
meta_data=colnames(all_counts)
meta_data<-as.data.frame(meta_data)

#add a condition column
meta_data$condition <- c(rep("C3", 3), rep("C5", 3), rep("fb", 3))
colnames(meta_data)<-c("sample", "condition")

#make factors
meta_data$sample<-as.factor(meta_data$sample)
meta_data$condition<-as.factor(meta_data$condition)

```




6. Now process the count matrix and examine PCA plots

```R

order <- c("fb", "C3", "C5")  #the order you want your conditions to appear

result <- process_deseq_PCA(counts = all_counts, 
                             meta_data = meta_data, 
                             relevel_condition = order, 
                             min_gene_count = 200)

dds <- result$dds
pca_plot <- result$pca_plot

``` 


7. If you want to run go terms, then you need the correct reference genome

```R

library(gsfisher)
annotation_gs <- fetchAnnotation(species="mm", ensembl_version=NULL, ensembl_host=NULL)

```



8. Run the DE analysis, plot the heatmap and chnage orders and colors as required


```R

custom_colors <- c("fb" = "red", "C3" = "blue", "C5" = "green")

#define order of how you want to plot the heatmap
custom_order <- c("fb", "C3", "C5")


run_deseq_contrast_heatmap(dds = dds, 
                           meta_data = meta_data, 
                           targetvar = "condition", 
                           p_value_threshold = 0.05, 
                           log2fc_threshold = 1, 
                           condition_colors = custom_colors,
                           condition_order = custom_order,
                           go_analysis = TRUE, 
                           organism = "mm", 
                           annotation_gs = annotation_gs)

```


9. You might need to plot the heatmap/goterms seperatly and customise a bit

```R
#Heatmap
Heatmap(scale_sub_vsd,
              top_annotation = top_annotation,
              col = colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red")),
              row_names_gp = gpar(fontsize = 4),
              cluster_columns = FALSE,
              cluster_rows = TRUE,
              show_row_names = FALSE,
              show_column_names = FALSE,
              border = TRUE)



#go terms
sampleEnrichmentDotplot(go.results.top, selection_col = "description", selected_genesets = unique(go.results.top$description), sample_id_col = "comparison", fill_var = "odds.ratio", maxl=50, title="Go term",rotate_sample_labels = T)


```

10. Bulk RNAseq done

