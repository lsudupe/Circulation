### SCRIPT: Deconvolution using SPOTlight package in IBEX

## 12.09.21 Laura Sudupe , git @lsudupe

## https://github.com/MarcElosua/SPOTlight

#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")


#install.packages("Seurat")
install.packages("devtools")
install.packages("base")

devtools::install_github("https://github.com/MarcElosua/SPOTlight")
BiocManager::install("SPOTlight")

library("Seurat")
library("SPOTlight")
library("base")

##Data-----------------------------------------------------
integrated <- readRDS("./objects/integrated/integrated.sct.rds")
fibro.ref <- sc.ref <- readRDS("./objects/sc/Tokio/sc.combined.sct.fibro.rds")

##marker genes 
Seurat::Idents(object = fibro.ref) <- fibro.ref@meta.data$ident
cluster_markers_fibro <- Seurat::FindAllMarkers(object = fibro.ref, 
 assay = "integrated",
  slot = "data",
 verbose = TRUE, 
  only.pos = TRUE)
saveRDS(cluster_markers_fibro,"./results/spotlight/marker_genes/combined.fibro_markers.RDS")

sc_markers <- readRDS("./results/spotlight/marker_genes/combined.fibro_markers.RDS")

#Decomposition fibro.ref
spot_SCT <- spotlight_deconvolution(
  se_sc = fibro.ref,
  counts_spatial = combined@assays[["Spatial"]]@counts,
  clust_vr = "seurat_clusters", # Variable in sc_seu containing the cell-type annotation
  cluster_markers = sc_markers, # Dataframe with the marker genes
  cl_n = 100, # number of cells per cell type to use
  hvg = 3000, # Number of HVG to use
  ntop = NULL, # How many of the marker genes to use (by default all)
  transf = "uv", # Perform unit-variance scaling per cell and spot prior to factorzation and NLS
  method = "nsNMF", # Factorization method
  min_cont = 0 # Remove those cells contributing to a spot below a certain threshold 
)






