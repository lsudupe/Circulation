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

fibro <- readRDS("./objects/integrated/integrated.fb.rds")


spata_obj <- initiateSpataObject_10X(input_path = "./data/infarto_5dpi_hembra/",
                                    sample_names= c("integrated"),
                                    gene_set_path = fibro)

https://rdrr.io/github/kueckelj/SPATA2/man/transformSeuratToSpata.html


fibro <- loadSpataObject("./objects/integrated/integrated.fb.rds")

createSegmentation(fibro)

