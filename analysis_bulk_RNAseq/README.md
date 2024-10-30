# This is an overview of analysing bulk RNAseq


1. Complete this step: https://github.com/chrismahony/processing_piplines_bluebear/tree/main/process_bulk_RNAseq

<br>

2. Start an R session

3. Load the required packages

```R

library(DESeq2)
library(gsfisher)
library(ComplexHeatmap)
 ````


4. Install sprokforlife and load the package. Here you will have several handy functions that act as wrappers to process bulk and scRNAseq data

```R

devtools::install_github("chrismahony/sprokforlife")
library(sporkforlife)

```
