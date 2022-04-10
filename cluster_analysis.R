## SCRIPT: cluster analysis spatial samples CIRCULATION project

## 05.04.22 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)

###Data
integrated <- readRDS("./objects/integrated/integrated.gfp.rds")

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


###Individual GFP plots
pdf(file.path("./results/gfp/",filename = "sham_gfp.pdf"))
SpatialPlot(object = control, features = c("GFP")) 
dev.off()

SpatialPlot(object = infarto_3dpi, features = c("GFP"))
SpatialPlot(object = infarto_5dpi_h, features = c("GFP"))
SpatialPlot(object = infarto_5dpi_m, features = c("GFP"))


pdf(file.path("./results/gfp/",filename = "3dpi_gfp.pdf"))
SpatialFeaturePlot(infarto_3dpi, features = "GFP",
                   min.cutoff = 0,
                   max.cutoff = 60,)
dev.off()

pdf(file.path("./results/gfp/",filename = "5dpi_h_gfp.pdf"))
SpatialFeaturePlot(infarto_5dpi_h, features = "GFP",
                   min.cutoff = 0,
                   max.cutoff = 80,)
dev.off()


pdf(file.path("./results/gfp/",filename = "5dpi_m_gfp.pdf"))
SpatialFeaturePlot(infarto_5dpi_m, features = "GFP",
                   min.cutoff = 0,
                   max.cutoff = 100,)
dev.off()




