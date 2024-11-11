
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




```
   


