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

pdf(file.path(paths,filename ="Dentritic feature umap.pdf"))
FeaturePlot(a,features = dentritic)
dev.off()
pdf(file.path(paths,filename = "Macrophagues feature umap.pdf"))
FeaturePlot(a, features = Macrophagues)
dev.off()
pdf(file.path(paths,filename = "Schwann.cells feature umap.pdf"))
FeaturePlot(a, features = Schwann.cells)
dev.off()
pdf(file.path(paths,filename ="DC feature umap.pdf"))
FeaturePlot(a,features = DC.like.cells)
dev.off()
pdf(file.path(paths,filename ="Fibro feature umap.pdf"))
FeaturePlot(a, features = fibroblast)
dev.off()
pdf(file.path(paths,filename ="Circulation fibro feature umap.pdf"))
FeaturePlot(a, features = circulation.fibroblast)
dev.off()
pdf(file.path(paths,filename ="Cardiomyocites feature umap.pdf"))
FeaturePlot(a, features = cardiomyocite)
dev.off()
pdf(file.path(paths,filename ="Endohtelial feature umap.pdf"))
FeaturePlot(a,features = endothelial)
dev.off()
pdf(file.path(paths,filename ="Smooth muscle feature umap.pdf"))
FeaturePlot(a,features = smooth.muscle)
dev.off()
pdf(file.path(paths,filename ="B cell feature umap.pdf"))
FeaturePlot(a,features = B.cell)
dev.off()
pdf(file.path(paths,filename ="Machofage/monocytes feature umap.pdf"))
FeaturePlot(a,features = machofage.monocytes)
dev.off()
pdf(file.path(paths,filename ="T cells feature umap.pdf"))
FeaturePlot(a,features = T.cell)
dev.off()
pdf(file.path(paths,filename ="Neutrophils feature umap.pdf"))
FeaturePlot(a, features = neutrophils)
dev.off()
pdf(file.path(paths,filename ="Dentritic feature umap.pdf"))
FeaturePlot(a,features = dentritic)
dev.off()
pdf(file.path(paths,filename ="Granulocytes feature umap.pdf"))
FeaturePlot(a, features = granulocytes)
dev.off()

