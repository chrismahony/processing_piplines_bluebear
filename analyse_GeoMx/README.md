
1. Load packages

```R
#remotes::install_github("DavisLaboratory/GeoMXAnalysisWorkflow", build_vignettes = FALSE)  #install as required
library(GeoMXAnalysisWorkflow)
library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)

```


2. Read in in the .DCC data

```R

DCCFiles <- dir(file.path("/rds/projects/c/croftap-actacfbmac/geomix_nov2024/DCC_files"), pattern = ".dcc$",
                full.names = TRUE, recursive = TRUE)


# PKC for Human is uploaded to this repo
PKCFiles <- dir(file.path("/rds/projects/c/croftap-actacfbmac/geomix_nov2024/"), pattern = ".pkc$",
                                full.names = TRUE, recursive = TRUE)


SampleAnnotationFile <-
    dir("/rds/projects/c/croftap-actacfbmac/geomix_nov2024/Annie_run_3_20240903T1516/", pattern = "annotaitons.xlsx",
        full.names = TRUE, recursive = TRUE)



demoData <-
    readNanoStringGeoMxSet(dccFiles = DCCFiles,
                           pkcFiles = PKCFiles,
                           phenoDataFile = SampleAnnotationFile,
                               phenoDataDccColName = "Sample_ID",
                           phenoDataSheet = "annotaitons",#must match name of sheet in excel
                           protocolDataColNames = c("Tags", "Scan Name"),
                           experimentDataColNames = c("panel"))


```
   
3. Calculate QC metrics and plot thresholds




```R

library(knitr)
pkcs <- annotation(demoData)
modules <- gsub(".pkc", "", pkcs)
kable(data.frame(PKCs = pkcs, modules = modules))

#demoData <- shiftCountsOne(demoData, useDALogic = TRUE)

QC_params <-
    list(minSegmentReads = 1000, # Minimum number of reads (1000)
         percentTrimmed = 80,    # Minimum % of reads trimmed (80%)
         percentStitched = 80,   # Minimum % of reads stitched (80%)
         percentAligned = 75,    # Minimum % of reads aligned (80%)
         percentSaturation = 50, # Minimum sequencing saturation (50%)
         minNegativeCount = 10,   # Minimum negative control counts (10)
         maxNTCCount = 1000,     # Maximum counts observed in NTC well (1000)
         minNuclei = 10,         # Minimum # of nuclei estimated (100)
         minArea = 600)         # Minimum segment area (5000)
demoData <-
    setSegmentQCFlags(demoData, 
                      qcCutoffs = QC_params)        

# Collate QC Results
QCResults <- protocolData(demoData)[["QCFlags"]]
flag_columns <- colnames(QCResults)
QC_Summary <- data.frame(Pass = colSums(!QCResults[, flag_columns]),
                         Warning = colSums(QCResults[, flag_columns]))
QCResults$QCStatus <- apply(QCResults, 1L, function(x) {
    ifelse(sum(x) == 0L, "PASS", "WARNING")
})
QC_Summary["TOTAL FLAGS", ] <-
    c(sum(QCResults[, "QCStatus"] == "PASS"),
      sum(QCResults[, "QCStatus"] == "WARNING"))

library(ggplot2)

col_by <- "Tags"

# Graphical summaries of QC statistics plot function
QC_histogram <- function(assay_data = NULL,
                         annotation = NULL,
                         fill_by = NULL,
                         thr = NULL,
                         scale_trans = NULL) {
    plt <- ggplot(assay_data,
                  aes_string(x = paste0("unlist(`", annotation, "`)"),
                             fill = fill_by)) +
        geom_histogram(bins = 50) +
        geom_vline(xintercept = thr, lty = "dashed", color = "black") +
        theme_bw() + guides(fill = "none") +
        facet_wrap(as.formula(paste("~", fill_by)), nrow = 4) +
        labs(x = annotation, y = "Segments, #", title = annotation)
    if(!is.null(scale_trans)) {
        plt <- plt +
            scale_x_continuous(trans = scale_trans)
    }
    plt
}


QC_histogram(sData(demoData), "Trimmed (%)", col_by, 80)
QC_histogram(sData(demoData), "Stitched (%)", col_by, 80)
QC_histogram(sData(demoData), "Aligned (%)", col_by, 75)
QC_histogram(sData(demoData), "Saturated (%)", col_by, 50) +
    labs(title = "Sequencing Saturation (%)",
         x = "Sequencing Saturation (%)")

QC_histogram(sData(demoData), "Area", col_by, 600, scale_trans = "log10")
QC_histogram(sData(demoData), "Nuclei", col_by, 10)


```



