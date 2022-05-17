## SCRIPT: individual analysis spatial samples CIRCULATION project

## 05.04.22 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)

###Read data
control <- readRDS("./objects/initial/control_noGFP.rds")
dpi3 <- readRDS("./objects/initial/dpi3_noGFP.rds")
dpi5_female <- readRDS("./objects/initial/dpi5_female_noGFP.rds")
dpi5_male <- readRDS("./objects/initial/dpi5_male_noGFP.rds")

control <- SCTransform(control, assay = "Spatial", verbose = FALSE) %>%
  RunPCA(assay = "SCT",npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE)%>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0.4)

dpi3 <- SCTransform(dpi3, assay = "Spatial", verbose = FALSE) %>%
  RunPCA(assay = "SCT",npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE)%>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0.4)

dpi5_female <- SCTransform(dpi5_female, assay = "Spatial", verbose = FALSE) %>%
  RunPCA(assay = "SCT",npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE)%>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0.4)

dpi5_male <- SCTransform(dpi5_male, assay = "Spatial", verbose = FALSE) %>%
  RunPCA(assay = "SCT",npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE)%>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0.4)

saveRDS(control, "./objects/individual/control.rds")
saveRDS(dpi3, "./objects/individual/dpi3.rds")
saveRDS(dpi5_female, "./objects/individual/dpi5_female.rds")
saveRDS(dpi5_male, "./objects/individual/dpi5_male.rds")
