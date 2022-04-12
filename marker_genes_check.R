## SCRIPT: marker genes check

## 28.11.21 Laura Sudupe , git @lsudupe

#Libraries-------------------------------------------------
source("packages.R")
features <- source("marker_genes.R")

###Read data

a<- readRDS("./objects/integrated/integrated.gfp.rds")

###Select the path for the results
paths <- "./results/marker_genes/"

###################################################3

###Plot and save heatmap

###################################################3
png(file.path(paths,filename = "Macrophagues feature.png"))
DoHeatmap(subset(a, downsample = 100), features = Macrophagues, size = 3)
dev.off()
png(file.path(paths,filename = "Schwann.cells feature.png"))
DoHeatmap(subset(a, downsample = 100), features = Schwann.cells, size = 3)
dev.off()
png(file.path(paths,filename ="DC feature.png"))
DoHeatmap(subset(a, downsample = 100), features = DC.like.cells, size = 3)
dev.off()
png(file.path(paths,filename ="Fibro feature.png"))
DoHeatmap(subset(a, downsample = 100), features = fibroblast, size = 3)
dev.off()
png(file.path(paths,filename ="Circulation fibro feature.png"))
DoHeatmap(subset(a, downsample = 100), features = circulation.fibroblast, size = 3)
dev.off()
png(file.path(paths,filename ="Cardiomyocites feature.png"))
DoHeatmap(subset(a, downsample = 100), features = cardiomyocite, size = 3)
dev.off()
png(file.path(paths,filename ="Endohtelial feature.png"))
DoHeatmap(subset(a, downsample = 100), features = endothelial, size = 3)
dev.off()
png(file.path(paths,filename ="Smooth muscle feature.png"))
DoHeatmap(subset(a, downsample = 100), features = smooth.muscle, size = 3)
dev.off()
png(file.path(paths,filename ="B cell feature.png"))
DoHeatmap(subset(a, downsample = 100), features = B.cell, size = 3)
dev.off()
png(file.path(paths,filename ="Machofage/monocytes feature.png"))
DoHeatmap(subset(a, downsample = 100), features = machofage.monocytes, size = 3)
dev.off()
png(file.path(paths,filename ="T cells feature.png"))
DoHeatmap(subset(a, downsample = 100), features = T.cell, size = 3)
dev.off()
png(file.path(paths,filename ="Neutrophils feature.png"))
DoHeatmap(subset(a, downsample = 100), features = neutrophils, size = 3)
dev.off()
png(file.path(paths,filename ="Dentritic feature.png"))
DoHeatmap(subset(a, downsample = 100), features = dentritic, size = 3)
dev.off()
png(file.path(paths,filename ="Granulocytes feature.png"))
DoHeatmap(subset(a, downsample = 100), features = granulocytes, size = 3)
dev.off()

###################################################3

###Plot and save doheatmap

###################################################3

#Density plots in clusters
library(GSEABase)
##create geneSet
fibro.sig <- GeneSet(fibroblast, setName="geneSetFibro")
circulation.sig <- GeneSet(circulation.fibroblast, setName="geneSetCirc")
geneSets <- GeneSetCollection(fibro.sig, circulation.sig)
##enrichment
a@assays[["integrated"]]@counts <- a@assays[["integrated"]]@data
matrix <- a@assays[["integrated"]]@counts
#matrix <- as.matrix(matrix)
ES <- enrichIt(obj = matrix, 
                      gene.sets = geneSets, 
                      groups = 1000, cores = 2)
##save meta
a <- AddMetaData(a, ES)
##fibro
pdf(file.path(paths,filename ="Fibro dens.pdf"))
print(ridgeEnrichment(a@meta.data, gene.set = "geneSetFibro", group = "seurat_clusters", facet = "sample", add.rug = TRUE))
dev.off()

pdf(file.path(paths,filename ="circ dens.pdf"))
print(ridgeEnrichment(a@meta.data, gene.set = "geneSetCirc", group = "seurat_clusters", facet = "sample", add.rug = TRUE))
dev.off()

#####check same of them spatially
DefaultAssay(a) <- "SCT"
  
for (i in circulation.fibroblast){
  #p1 <- SpatialFeaturePlot(integrated, features = i, alpha = 0.6, combine = FALSE)
  #fix.p1 <- scale_fill_continuous(type = "viridis")
  #p2 <- lapply(p1, function (x) x + fix.p1)
  
  pdf(paste("D2",i,".pdf",sep=""))
  #print(CombinePlots(p2))
  SpatialFeaturePlot(a, features = i, alpha = 0.6, combine = FALSE)
  dev.off()
}
