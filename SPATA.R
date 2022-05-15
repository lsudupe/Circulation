### SCRIPT: SPATA method in Circulation data

### 10.04.22 Laura Sudupe , git @lsudupe

#https://themilolab.github.io/SPATA/articles/spata-compatibility.html

###### libraries
library("base")
base::install.packages("devtools")

if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils'))

devtools::install_github('cole-trapnell-lab/leidenbase')
devtools::install_github('cole-trapnell-lab/monocle3')

#devtools::install_github(repo = 'satijalab/seurat', ref = 'develop')

devtools::install_github(repo = "kueckelj/confuns")
devtools::install_github(repo = "theMILOlab/SPATA")

library(SPATA)
library(magrittr)

library(tidyverse)
library(patchwork)

library(monocle3)
library(Seurat)

#install.packages("remotes")
remotes::install_github("kueckelj/SPATA2")
library(SPATA2)

fibro <- readRDS("./objects/integrated/integrated.fb.rds")
control <- fibro
dpi3 <- fibro
dpi5_female <- fibro
dpi5_male <- fibro

###subset the data
control_s <- subset(fibro, subset = sample == "control")


###create spata object
control_spata <- transformSeuratToSpata(
  control_s,
  sample_name = "control",
  method = "spatial",
  coords_from = "pca",
  assay_name = "SCT",
  assay_slot = NULL,
  image_name = "control",
  gene_set_path = NULL,
  verbose = TRUE
)


#https://rdrr.io/github/kueckelj/SPATA2/man/transformSeuratToSpata.html


createSegmentation(fibro)

