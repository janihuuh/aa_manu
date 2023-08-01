library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(data.table)
library(Seurat)
library(cowplot)
library(ineq)
library(SingleCellExperiment)
library(DESeq2)
library(circlize)
library(ComplexHeatmap)
library(dbscan)
library(ggrastr)

## Set global ggplot2-themes
theme_set(theme_classic(base_size = 17))

## Run all fun_* codes
for(code in list.files("src/R/general/", "^fun", full.names = T, recursive = T)){
  message(code)
  source(code)
}

for(code in list.files("src/functional/", "^fun_helper", full.names = T, recursive = T)){
  message(code)
  source(code)
}
