
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


4. Calculate Negative gemotric means (means of log or exponential data) of negative probes

```R

negativeGeoMeans <- 
    esBy(negativeControlSubset(demoData), 
         GROUP = "Module", 
         FUN = function(x) { 
             assayDataApply(x, MARGIN = 2, FUN = ngeoMean, elt = "exprs") 
         }) 
protocolData(demoData)[["NegGeoMean"]] <- negativeGeoMeans

# explicitly copy the Negative geoMeans from sData to pData
negCols <- paste0("NegGeoMean_", modules)
pData(demoData)[, negCols] <- sData(demoData)[["NegGeoMean"]]
for(ann in negCols) {
    plt <- QC_histogram(pData(demoData), ann, "Segment", 2, scale_trans = "log10")
    print(plt)
}

plt

```

5. Gene detection analysis

```R

col_by <- "Tags2"

pData(demoData) <- pData(demoData)[, !colnames(pData(demoData)) %in% negCols]

## Remove Flagged segments
passedQC <- demoData[, QCResults$QCStatus == "PASS"]

### Subsetting our dataset has removed samples which did not pass QC
dim(demoData)
dim(passedQC)
demoData <- passedQC

# Probe QC
## Set QC Flags
### Generally keep the qcCutoffs parameters unchanged. Set removeLocalOutliers to 
### FALSE if you do not want to remove local outliers
demoData <- setBioProbeQCFlags(demoData, 
                               qcCutoffs = list(minProbeRatio = 0.1, percentFailGrubbs = 20), 
                               removeLocalOutliers = F)

ProbeQCResults <- fData(demoData)[["QCFlags"]]

### Define QC table for Probe QC
qc_df <- data.frame(Passed = sum(rowSums(ProbeQCResults[, -1]) == 0),
                    Global = sum(ProbeQCResults$GlobalGrubbsOutlier),
                    Local = sum(rowSums(ProbeQCResults[, -2:-1]) > 0
                                & !ProbeQCResults$GlobalGrubbsOutlier))
## Exclude Outliers
### Subset object to exclude all that did not pass Ratio & Global testing
ProbeQCPassed <- subset(demoData, 
                        fData(demoData)[["QCFlags"]][,c("LowProbeRatio")] == F &
                          fData(demoData)[["QCFlags"]][,c("GlobalGrubbsOutlier")] == F)
dim(ProbeQCPassed)

demoData <- ProbeQCPassed 


# Create Gene-level Count Data
## Check how many unique targets the object has
length(unique(featureData(demoData)[["TargetName"]]))

### collapse to targets
target_demoData <- aggregateCounts(demoData)
dim(target_demoData)

# exprs(target_demoData)[1:5, 1:2]

# Limit of Quantification (LOQ)
### Define LOQ SD threshold and minimum value
cutoff <- 2
minLOQ <- 2

## Calculate LOQ per module tested
LOQ <- data.frame(row.names = colnames(target_demoData))
for(module in modules) {
  vars <- paste0(c("NegGeoMean_", "NegGeoSD_"),
                 module)
  if(all(vars[1:2] %in% colnames(pData(target_demoData)))) {
    LOQ[, module] <-
      pmax(minLOQ,
           pData(target_demoData)[, vars[1]] * 
             pData(target_demoData)[, vars[2]] ^ cutoff)
  }
}
pData(target_demoData)$LOQ <- LOQ

# Filtering
LOQ_Mat <- c()
for(module in modules) {
  ind <- fData(target_demoData)$Module == module
  Mat_i <- t(esApply(target_demoData[ind, ], MARGIN = 1,
                     FUN = function(x) {
                       x > LOQ[, module]
                     }))
  LOQ_Mat <- rbind(LOQ_Mat, Mat_i)
}
### ensure ordering since this is stored outside of the geomxSet
LOQ_Mat <- LOQ_Mat[fData(target_demoData)$TargetName, ]

## Segment Gene Detection
### Save detection rate information to pheno data
pData(target_demoData)$GenesDetected <- colSums(LOQ_Mat, na.rm = T)
pData(target_demoData)$GeneDetectionRate <- pData(target_demoData)$GenesDetected / nrow(target_demoData)

### Determine detection thresholds: 1%, 5%, 10%, 15%, >15%
pData(target_demoData)$DetectionThreshold <- 
  cut(pData(target_demoData)$GeneDetectionRate,
      breaks = c(0, 0.01, 0.05, 0.1, 0.15, 1),
      labels = c("<1%", "1-5%", "5-10%", "10-15%", ">15%"))

### stacked bar plot of different cut points (1%, 5%, 10%, 15%)
ggplot(pData(target_demoData),
       aes(x = DetectionThreshold)) +
  geom_bar(aes(fill = `Segment`)) +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Gene Detection Rate",
       y = "Segments, #",
       fill = "Segment Type")

```


6. Retain only GOI

```R

### cut percent genes detected at 1, 5, 10, 15
kable(table(pData(target_demoData)$DetectionThreshold,
            pData(target_demoData)$Segment))

## Filter (generally 5-10% is a good filer, but might need adjusting)
target_demoData <- target_demoData[, pData(target_demoData)$GeneDetectionRate >= .1] # In the vignette, they do 10% of Genes

dim(target_demoData)


## Gene Detection Rate
library(scales) # for percent

### Calculate detection rate:
LOQ_Mat <- LOQ_Mat[, colnames(target_demoData)]
fData(target_demoData)$DetectedSegments <- rowSums(LOQ_Mat, na.rm = T)
fData(target_demoData)$DetectionRate <-
  fData(target_demoData)$DetectedSegments / nrow(pData(target_demoData))

### Gene of interest detection table
goi <- c("PDCD1", "CD274", "IFNG", "CD8A", "CD68", "RORC",
         "KRT18", "NPHS1", "NPHS2", "CALB1", "CLDN8")
goi_df <- data.frame(
  Gene = goi,
  Number = fData(target_demoData)[goi, "DetectedSegments"],
  DetectionRate = percent(fData(target_demoData)[goi, "DetectionRate"]))

head(goi_df)



## Gene Filtering
### Plot detection rate:
plot_detect <- data.frame(Freq = c(1, 5, 10, 20, 30, 50))
plot_detect$Number <- unlist(lapply(c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5),
                function(x) {sum(fData(target_demoData)$DetectionRate >= x)}))
plot_detect$Rate <- plot_detect$Number / nrow(fData(target_demoData))
rownames(plot_detect) <- plot_detect$Freq

ggplot(plot_detect, aes(x = as.factor(Freq), y = Rate, fill = Rate)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = formatC(Number, format = "d", big.mark = ",")),
            vjust = 1.6, color = "black", size = 4) +
  scale_fill_gradient2(low = "orange2", mid = "lightblue",
                       high = "dodgerblue3", midpoint = 0.65,
                       limits = c(0,1),
                       labels = scales::percent) +
  theme_bw() +
  scale_y_continuous(labels = scales::percent, limits = c(0,1),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "% of Segments",
       y = "Genes Detected, % of Panel > LOQ")


### Subset to target genes detected in at least 10% of the samples.
###   Also manually include the negative control probe, for downstream use
negativeProbefData <- subset(fData(target_demoData), CodeClass == "Negative")
neg_probes <- unique(negativeProbefData$TargetName)
target_demoData <- target_demoData[fData(target_demoData)$DetectionRate >= 0.1 |
                    fData(target_demoData)$TargetName %in% neg_probes, ]
dim(target_demoData)

### retain only detected genes of interest
goi <- goi[goi %in% rownames(target_demoData)]

```


7. Normalization

```R

# Normalisation --------------------------------------------------------------
library(reshape2)  # for melt
library(cowplot)   # for plot_grid

head(pData(target_demoData))
# Graph Q3 value vs negGeoMean of Negatives
ann_of_interest <- "Segment"
Stat_data <- data.frame(row.names = colnames(exprs(target_demoData)),
             Segment = colnames(exprs(target_demoData)),
             Annotation = pData(target_demoData)[, ann_of_interest],
             Q3 = unlist(apply(exprs(target_demoData), 2,
                               quantile, 0.75, na.rm = T)),
             NegProbe = exprs(target_demoData)[neg_probes, ])
Stat_data_m <- melt(Stat_data, measure.vars = c("Q3", "NegProbe"),
                    variable.name = "Statistic", value.name = "Value")

plt1 <- ggplot(Stat_data_m,
               aes(x = Value, fill = Statistic)) +
  geom_histogram(bins = 40) + theme_bw() +
  scale_x_continuous(trans = "log2") +
  facet_wrap(~Annotation, nrow = 1) + 
  scale_fill_brewer(palette = 3, type = "qual") +
  labs(x = "Counts", y = "Segments, #")

plt2 <- ggplot(Stat_data,
               aes(x = NegProbe, y = Q3, color = Annotation)) +
  geom_abline(intercept = 0, slope = 1, lty = "dashed", color = "darkgray") +
  geom_point() + guides(color = "none") + theme_bw() +
  scale_x_continuous(trans = "log2") + 
  scale_y_continuous(trans = "log2") +
  theme(aspect.ratio = 1) +
  labs(x = "Negative Probe GeoMean, Counts", y = "Q3 Value, Counts")

plt3 <- ggplot(Stat_data,
               aes(x = NegProbe, y = Q3 / NegProbe, color = Annotation)) +
  geom_hline(yintercept = 1, lty = "dashed", color = "darkgray") +
  geom_point() + theme_bw() +
  scale_x_continuous(trans = "log2") + 
  scale_y_continuous(trans = "log2") +
  theme(aspect.ratio = 1) +
  labs(x = "Negative Probe GeoMean, Counts", y = "Q3/NegProbe Value, Counts")

btm_row <- plot_grid(plt2, plt3, nrow = 1, labels = c("B", ""),
                     rel_widths = c(0.43,0.57))
plot_grid(plt1, btm_row, ncol = 1, labels = c("A", "")) # Separation of Q3 & NegProbe (there's a "shift" of the Q3 away from the NegProbe)


## Perform Normalisation (data structure can hold both, but need to take note of the "toElt" argument)
# Q3 norm (75th percentile) for WTA/CTA  with or without custom spike-ins
target_demoData <- normalize(target_demoData,
                             norm_method = "quant", 
                             desiredQuantile = .75,
                             toElt = "q_norm")

# Background normalization for WTA/CTA without custom spike-in
target_demoData <- normalize(target_demoData ,
                             norm_method = "neg",
                             fromElt = "exprs",
                             toElt = "neg_norm")

## Visualise
### visualize the first 10 segments with each normalization method
boxplot(exprs(target_demoData)[,1:10],
        col = "#9EDAE5", main = "Raw Counts",
        log = "y", names = 1:10, xlab = "Segment",
        ylab = "Counts, Raw")

### Q3
boxplot(assayDataElement(target_demoData[,1:10], elt = "q_norm"),
        col = "#2CA02C", main = "Q3 Norm Counts",
        log = "y", names = 1:10, xlab = "Segment",
        ylab = "Counts, Q3 Normalized")

### Background
boxplot(assayDataElement(target_demoData[,1:10], elt = "neg_norm"),
        col = "#FF7F0E", main = "Neg Norm Counts",
        log = "y", names = 1:10, xlab = "Segment",
        ylab = "Counts, Neg. Normalized")

```

8. UMAP and tSNE plotting of samples

```R

# Unsupervised Analysis ----------
##UMAPS
# install.packages("umap")
library(umap)
library(Rtsne)

# update defaults for umap to contain a stable random_state (seed)
custom_umap <- umap::umap.defaults
custom_umap$random_state <- 42

# run UMAP
umap_out <- umap(t(log2(assayDataElement(target_demoData , elt = "q_norm"))), config = custom_umap)
pData(target_demoData)[, c("UMAP1", "UMAP2")] <- umap_out$layout[, c(1,2)]
ggplot(pData(target_demoData),
       aes(x = UMAP1, y = UMAP2, color = `Segment`)) +
  geom_point(size = 3) +
  theme_bw()

## tSNE
# run tSNE
set.seed(42) # set the seed for tSNE as well
tsne_out <- Rtsne(t(log2(assayDataElement(target_demoData , elt = "q_norm"))),
        perplexity = ncol(target_demoData)*.15)
pData(target_demoData)[, c("tSNE1", "tSNE2")] <- tsne_out$Y[, c(1,2)]
ggplot(pData(target_demoData),
       aes(x = tSNE1, y = tSNE2, color = `Segment`)) +
  geom_point(size = 3) +
  theme_bw()

```


9. Plotting gene expression

```R


library(pheatmap)  # for pheatmap
# create a log2 transform of the data for analysis
assayDataElement(object = target_demoData, elt = "log_q") <-
    assayDataApply(target_demoData, 2, FUN = log, base = 2, elt = "q_norm")

# create CV function
calc_CV <- function(x) {sd(x) / mean(x)}
CV_dat <- assayDataApply(target_demoData,
                         elt = "log_q", MARGIN = 1, calc_CV)
# show the highest CD genes and their CV values
sort(CV_dat, decreasing = TRUE)[1:5]
#>   CAMK2N1    AKR1C1      AQP2     GDF15       REN 
#> 0.5886006 0.5114973 0.4607206 0.4196469 0.4193216

# Identify genes in the top 3rd of the CV values
GOI <- names(CV_dat)[CV_dat > quantile(CV_dat, 0.8)]
pheatmap(assayDataElement(target_demoData[GOI, ], elt = "log_q"),
         scale = "row", 
         show_rownames = FALSE, show_colnames = FALSE,
         border_color = NA,
       breaks = seq(-3, 3, 0.05),
         color = colorRampPalette(c("purple3", "black", "yellow2"))(120),
         annotation_col = 
             pData(target_demoData)[, c("Segment", "Slide Name")])

pheatmap(assayDataElement(target_demoData[c("PTPRC", "COL1A1"), ], elt = "log_q"),
         scale = "row", 
         show_rownames = T, show_colnames = FALSE,
         border_color = NA,
       breaks = seq(-3, 3, 0.05),
         color = colorRampPalette(c("purple3", "black", "yellow2"))(120),
         annotation_col = 
             pData(target_demoData)[, c("Segment", "Slide Name")])


 pheatmap(assayDataElement(target_demoData[c("PTPRC", "COL1A1"), ], elt = "log_q"),
         scale = "row", 
         show_rownames = T, show_colnames = T,
         border_color = NA,
         clustering_method = "average",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         breaks = seq(-3, 3, 0.05),
         color = colorRampPalette(c("purple3", "black", "yellow2"))(120),
         annotation_col = 
             pData(target_demoData)[, c("Segment", "Slide Name")])


```

10. DE analysis and plotting volcanoplot

```R

results <- c()
for(status in c("LM", "PI")) {
    ind <- pData(target_demoData)$Conditon2 == status
    mixedOutmc <-
        mixedModelDE(target_demoData[, ind],
                     elt = "log_q",
                     modelFormula = ~ testRegion + (1 + testRegion | slide),
                     groupVar = "testRegion",
                     nCores = 1,
                     multiCore = F)
    
    # format results as data.frame
    r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
    tests <- rownames(r_test)
    r_test <- as.data.frame(r_test)
    r_test$Contrast <- tests
    
    # use lapply in case you have multiple levels of your test factor to
    # correctly associate gene name with it's row in the results table
    r_test$Gene <- 
        unlist(lapply(colnames(mixedOutmc),
                      rep, nrow(mixedOutmc["lsmeans", ][[1]])))
    r_test$Subset <- status
    r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
    r_test <- r_test[, c("Gene", "Subset", "Contrast", "Estimate", 
                         "Pr(>|t|)", "FDR")]
    results <- rbind(results, r_test)
}


library(ggrepel) 
# Categorize Results based on P-value & FDR for plotting
results$Color <- "NS or FC < 0.5"
results$Color[results$`Pr(>|t|)` < 0.05] <- "P < 0.05"
results$Color[results$FDR < 0.05] <- "FDR < 0.05"
results$Color[results$FDR < 0.001] <- "FDR < 0.001"
results$Color[abs(results$Estimate) < 0.5] <- "NS or FC < 0.5"
results$Color <- factor(results$Color,
                        levels = c("NS or FC < 0.5", "P < 0.05",
                                   "FDR < 0.05", "FDR < 0.001"))

# pick top genes for either side of volcano to label
# order genes for convenience:
results$invert_P <- (-log10(results$`Pr(>|t|)`)) * sign(results$Estimate)
top_g <- c()
for(cond in c("LM")) {
    ind <- results$Subset == cond
    top_g <- c(top_g,
               results[ind, 'Gene'][
                   order(results[ind, 'invert_P'], decreasing = TRUE)[1:15]],
               results[ind, 'Gene'][
                   order(results[ind, 'invert_P'], decreasing = FALSE)[1:15]])
}
top_g <- unique(top_g)
results <- results[, -1*ncol(results)] # remove invert_P from matrix

# Graph results
ggplot(results,
       aes(x = Estimate, y = -log10(`Pr(>|t|)`),
           color = Color, label = Gene)) +
    geom_vline(xintercept = c(0.5, -0.5), lty = "dashed") +
    geom_hline(yintercept = -log10(0.05), lty = "dashed") +
    geom_point() +
    labs(x = "Enriched in Tubules <- log2(FC) -> Enriched in Glomeruli",
         y = "Significance, -log10(P)",
         color = "Significance") +
    scale_color_manual(values = c(`FDR < 0.001` = "dodgerblue",
                                  `FDR < 0.05` = "lightblue",
                                  `P < 0.05` = "orange2",
                                  `NS or FC < 0.5` = "gray"),
                       guide = guide_legend(override.aes = list(size = 4))) +
    scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
    geom_text_repel(data = subset(results, Gene %in% top_g & `Pr(>|t|)` < 0.05),
                    size = 4, point.padding = 0.15, color = "black",
                    min.segment.length = .1, box.padding = .2, lwd = 2,
                    max.overlaps = 50) +
    theme_bw(base_size = 16) +
    theme(legend.position = "bottom") +
    facet_wrap(~Subset, scales = "free_y")


```
