## SCRIPT: Integration spatial samples CIRCULATION project

## 05.04.22 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)

#Data--------------------------------------
control <- readRDS("./objects/initial/control_sp.rds")
infarto_3dpi <- readRDS("./objects/initial/3dpi_sp.rds")
infarto_5dpi_h <- readRDS("./objects/initial/5dpi_h_sp.rds")
infarto_5dpi_m <- readRDS("./objects/initial/5dpi_m_sp.rds")

#Take out GFP for the integration#########################
samples <- c(control, infarto_3dpi, infarto_5dpi_h, infarto_5dpi_m)
names(samples) <- c("control", "3dpi", "5dpi_h", "5dpi_m")

for (i in 1:length(samples)){
  a <- samples[[i]]
  counts <- as.data.frame(a@assays[["Spatial"]]@counts)
  counts <- counts[row.names(counts) != "GFP", , drop = FALSE]
  counts <- as.matrix(counts)
  a@assays[["Spatial"]]@counts <- counts

  data <- as.data.frame(a@assays[["Spatial"]]@data)
  data <- data[row.names(data) != "GFP", , drop = FALSE]
  data <- as.matrix(data)
  a@assays[["Spatial"]]@data <- data

  ##save object
  saveRDS(a,file = paste0("./objects/initial/",names(samples[i]),"_noGFP.rds"))
}

control <- readRDS("./objects/initial/control_noGFP.rds")
infarto_3dpi <- readRDS("./objects/initial/3dpi_noGFP.rds")
infarto_5dpi_h <- readRDS("./objects/initial/5dpi_h_noGFP.rds")
infarto_5dpi_m <- readRDS("./objects/initial/5dpi_m_noGFP.rds")

control@meta.data$sample <- "control"
infarto_3dpi@meta.data$sample <- "3dpi"
infarto_5dpi_h@meta.data$sample <- "5dpi_h"
infarto_5dpi_m@meta.data$sample <- "5dpi_m"

##########################################################

samples <- c(control, infarto_3dpi, infarto_5dpi_h, infarto_5dpi_m)
names(samples) <- c("control", "3dpi", "5dpi_h", "5dpi_m")

###Integration
###https://satijalab.org/seurat/articles/integration_introduction.html#performing-integration-on-datasets-normalized-with-sctransform-1
#list <- SplitObject(combined, split.by = c("sample")
samples <- lapply(X = samples, FUN = SCTransform, assay="Spatial")
features <- SelectIntegrationFeatures(object.list = samples, nfeatures = 3000)
samples <- PrepSCTIntegration(object.list = samples, anchor.features = features)


anchors <- FindIntegrationAnchors(object.list = samples, normalization.method = "SCT",
                                  anchor.features = features)
integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

saveRDS(integrated, "./objects/integrated/integrated.gfp.rds")
integrated <- readRDS("./objects/integrated/integrated.gfp.rds")


###Transformation
integrated <- RunPCA(integrated, assay = "integrated",npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE)%>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0.4)

saveRDS(integrated, "./objects/integrated/integrated.gfp.rds")

###Markers

#markers <- Seurat::FindAllMarkers(object = integrated, 
 #                                 assay = "integrated",
  #                                slot = "data",
   #                               verbose = TRUE, 
    #                              only.pos = TRUE)

#saveRDS(markers, "./results/clusters/integrated_markers.rds")
