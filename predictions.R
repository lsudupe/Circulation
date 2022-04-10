## SCRIPT: Spatial object label predictions with sc data

## 08.04.22 Laura Sudupe , git @lsudupe


#libraries
source("packages.R")

##Data-----------------------------------------------------
integrated <- readRDS("./objects/integrated/integrated.gfp.rds")
sc.ref <- readRDS("./objects/sc/Tokio/sc.combined.sct.fibro.rds")

#renormalize both 
sc.ref <- SCTransform(sc.ref, ncells = 3000, verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)
integrated <- SCTransform(integrated, assay = "Spatial", verbose = FALSE) %>%
  RunPCA(verbose = FALSE)

#anchors
anchors <- FindTransferAnchors(reference = sc.ref, 
                              query = integrated,
                               normalization.method = "SCT")

#predictions
integrated.predictions.assay <- TransferData(anchorset = anchors, refdata = sc.ref@meta.data[["seurat_clusters"]],
                                         prediction.assay = TRUE,
                                         weight.reduction = "pcaproject", 
                                         dims = 1:30)

integrated[["predictions"]] <- integrated.predictions.assay
DefaultAssay(integrated) <- "predictions"

#Based on these prediction scores, we can also predict cell types whose location is spatially restricted
#We use the same methods based on marked point processes to define spatially variable features
#but use the cell type prediction scores as the “marks” rather than gene expression.

integrated <- FindSpatiallyVariableFeatures(integrated, assay = "predictions", selection.method = "markvariogram",
                                        features = rownames(integrated), r.metric = 5, slot = "data")
top.clusters.integrated <- SpatiallyVariableFeatures(integrated)
integrated[["top.clusters"]] <- top.clusters.integrated

##save predictions objects
saveRDS(integrated, file = "./objects/predictions/integrated.pre.rds")

pdf(paste("./results/predictions/integrated_FB.pdf",sep=""))
SpatialPlot(object = integrated, features = top.clusters.integrated)#, ncol = 3)
dev.off()

