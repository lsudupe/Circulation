## SCRIPT: cluster analysis spatial samples CIRCULATION project

## 05.04.22 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)

###Read data
combined <- readRDS("./objects/initial/combined_noGFP.rds")

###Integration
###https://satijalab.org/seurat/articles/integration_introduction.html#performing-integration-on-datasets-normalized-with-sctransform-1
list <- SplitObject(combined, split.by = "sample")
saveRDS(list, "./objects/initial/list_sp.rds")
list <- lapply(X = list, FUN = SCTransform, assay="Spatial")
features <- SelectIntegrationFeatures(object.list = list, nfeatures = 3000)
list <- PrepSCTIntegration(object.list = list, anchor.features = features)


anchors <- FindIntegrationAnchors(object.list = list, normalization.method = "SCT",
                                  anchor.features = features)
combined <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

saveRDS(combined, "./objects/integrated/integrated.rds")
integrated <- readRDS("./objects/integrated/integrated.rds")


###Transformation

integrated <- RunPCA(integrated, assay = "integrated",npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE)%>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0.4)

saveRDS(integrated, "./objects/integrated/integrated.sct.rds")

###Markers

markers <- Seurat::FindAllMarkers(object = integrated, 
                                  assay = "integrated",
                                  slot = "data",
                                  verbose = TRUE, 
                                  only.pos = TRUE)

saveRDS(markers, "./results/integrated_markers.rds")

###Plots

pdf(file.path("./results/clusters/",filename = "umap_r0.4_integrated.pdf"))
DimPlot(integrated, group.by = c("sample"), label = T) + ggtitle("UMAP_r0.4")
dev.off()

pdf(file.path("./results/clusters/",filename = "umap_r0.4_integrated_clusters.pdf"))
DimPlot(integrated, group.by = c("seurat_clusters"), label = T) + ggtitle("UMAP_r0.4")
dev.off()

pdf(file.path("./results/clusters/",filename = "spatial_integrated.pdf"))
SpatialPlot(integrated, group.by = c("seurat_clusters"),label = TRUE, combine = FALSE)
dev.off()





