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

fibro <- readRDS("./objects/integrated/integrated.sct.rds")


###subset the data
control_s <- subset(fibro, subset = sample == "control")
dpi3 <- subset(fibro, subset = sample == "dpi3")
dpi5_female <- subset(fibro, subset = sample == "dpi5_female")
dpi5_male <- subset(fibro, subset = sample == "dpi5_male")

###create spata object
dpi3_spata <- transformSeuratToSpata(
  dpi3,
  sample_name = "dpi3",
  method = "spatial",
  coords_from = "pca",
  assay_name = "SCT",
  assay_slot = NULL,
  image_name = "dpi3",
  gene_set_path = NULL,
  verbose = TRUE
)

dpi5_male_spata <- transformSeuratToSpata(
  dpi5_male,
  sample_name = "dpi5_male",
  method = "spatial",
  coords_from = "pca",
  assay_name = "SCT",
  assay_slot = NULL,
  image_name = "dpi5_male",
  gene_set_path = NULL,
  verbose = TRUE
)

#https://rdrr.io/github/kueckelj/SPATA2/man/transformSeuratToSpata.html

dpi5_female_spata <- transformSeuratToSpata(
  dpi5_female,
  sample_name = "dpi5_female",
  method = "spatial",
  coords_from = "pca",
  assay_name = "SCT",
  assay_slot = NULL,
  image_name = "dpi5_female",
  gene_set_path = NULL,
  verbose = TRUE
)

saveRDS(dpi3_spata, "./objects/spata/dpi3_spata.rds")
saveRDS(dpi5_female_spata, "./objects/spata/dpi5_female_spata.rds")
saveRDS(dpi5_male_spata, "./objects/spata/dpi5_male_spata.rds")

dpi3_spata <- readRDS("./objects/spata/dpi3_spata.rds")
dpi5_female_spata <- readRDS("./objects/spata/dpi5_female_spata.rds")
dpi5_male_spata <- readRDS("./objects/spata/dpi5_male_spata.rds")


dpi5_male_spata <- SPATA2::createSegmentation(dpi5_male_spata)
SPATA2::plotSegmentation(object = dpi5_male_spata)
plotSurface(object = dpi3_spata, color_by = "segmentation")

SPATA2::createSegmentation(dpi3_spata)
SPATA2::plotSegmentation(object = dpi3_spata)

