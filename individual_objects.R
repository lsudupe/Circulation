## SCRIPT: individual object creation and first analaysis spatial samples CIRCULATION project

## 05.04.22 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(dplyr)

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


################################normalize gfp

control_gfp <- readRDS("./objects/initial/control_sp.rds")
dpi3_gfp <- readRDS("./objects/initial/dpi3_sp.rds")
dpi5_female_gfp <- readRDS("./objects/initial/dpi5_female_sp.rds")
dpi5_male_gfp <- readRDS("./objects/initial/dpi5_male_sp.rds")


control_gfp <- SCTransform(control_gfp, assay = "Spatial", verbose = FALSE) %>%
  RunPCA(assay = "SCT",npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE)%>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0.4)

dpi3_gfp <- SCTransform(dpi3_gfp, assay = "Spatial", verbose = FALSE) %>%
  RunPCA(assay = "SCT",npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE)%>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0.4)

dpi5_female_gfp <- SCTransform(dpi5_female_gfp, assay = "Spatial", verbose = FALSE) %>%
  RunPCA(assay = "SCT",npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE)%>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0.4)

dpi5_male_gfp <- SCTransform(dpi5_male_gfp, assay = "Spatial", verbose = FALSE) %>%
  RunPCA(assay = "SCT",npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE)%>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0.4)

saveRDS(control_gfp, "./objects/individual/control_gfp.rds")
saveRDS(dpi3_gfp, "./objects/individual/dpi3_gfp.rds")
saveRDS(dpi5_female_gfp, "./objects/individual/dpi5_female_gfp.rds")
saveRDS(dpi5_male_gfp, "./objects/individual/dpi5_male_gfp.rds")

