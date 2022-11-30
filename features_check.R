### SCRIPT:  Read genes and create geneSet CIRCULATION proyect

## 29.11.22 Laura Sudupe , git @lsudupe

###### Libraries
library(openxlsx)
library("Seurat")
library("GSEABase")
library("ggplot2")
library(UCell)
set.seed(123)

###### Read the gene sets
genes <- read.xlsx("./Circulation/data/silvia.xlsx", colNames = TRUE)
verde <- na.omit(genes$verde)
morado <- na.omit(genes$morado)


######spatial
heart <- readRDS("./Circulation/objects/integrated/integrated.rds")
DefaultAssay(heart) <- "SCT"

###### Check intersect genes B and object
v.genes <- intersect(heart@assays[["SCT"]]@data@Dimnames[[1]], verde)
v.list <- as.list(v.genes)
names(v.list) <- as.list(v.genes)
m.genes <- intersect(heart@assays[["SCT"]]@data@Dimnames[[1]], morado)
m.list <- as.list(m.genes)
names(m.list) <- as.list(m.genes)

##### Check individually


for (i in m.list){
  library(BuenColors)
  color <- jdb_palette("brewer_spectra")
  pdf(paste("./Circulation/results/Adrian.genes/silvia.december/individual/morado/",i,"_morado.pdf",sep=""))
  print(SpatialFeaturePlot(heart, features = i,  combine = TRUE, ncol = 2))
  dev.off()
}


##### Add moduleScore
DefaultAssay(heart) <- "Spatial"
v.genes <- list(intersect(heart@assays[["SCT"]]@data@Dimnames[[1]], verde))
heart_module <- AddModuleScore(heart, features = v.genes, name = c("verde_"))
m.genes <- list(intersect(heart@assays[["SCT"]]@data@Dimnames[[1]], morado))
heart_module <- AddModuleScore(heart_module, features = m.genes, name = c("morado_"))

##### Add UCellScore

heart_Uscore <- AddModuleScore_UCell(heart, features = v.genes, name = "_verde")
heart_Uscore <- AddModuleScore_UCell(heart_Uscore, features = m.genes, name = "_morado")

#####

pdf(paste("./Circulation/results/Adrian.genes/silvia.december/Ucell/geneset/morado/_geneset_morado.pdf",sep=""))
print(SpatialFeaturePlot(heart_Uscore, features = c("signature_1_morado"),  combine = TRUE, ncol = 2))
dev.off()

##### blanco y rojo

bl <- colorRampPalette(c("navy","royalblue","lightskyblue"))(200)                      
re <- colorRampPalette(c("mistyrose", "red2","darkred"))(200)

p1 <- SpatialFeaturePlot(heart_module, features = c("verde_1"), combine = FALSE, ncol = 2)
fix.p1 <- scale_fill_gradientn(colours=c(bl,"white", re), na.value = "grey98",limits = c(0,0.23))
p2 <- lapply(p1, function (x) x + fix.p1)
  
pdf(paste("./Circulation/results/Adrian.genes/silvia.december/geneset/verde/geneset_verde_bw.pdf",sep=""))
print(CombinePlots(p2))
dev.off()

##### normal
b <- c(min(heart_Uscore@meta.data[["signature_1_morado"]]), max(heart_Uscore@meta.data[["signature_1_morado"]]))
p1 <- SpatialFeaturePlot(heart_Uscore, features = c("signature_1_morado"),  combine = FALSE, ncol = 2)
fix.p1 <- scale_fill_gradientn(colours=color,
                               breaks=b,
                               labels=c("Min","Max"),
                               limits = b)
p2 <- lapply(p1, function (x) x + fix.p1)

pdf(paste("./Circulation/results/Adrian.genes/silvia.december/Ucell/geneset/morado/geneset_morado.pdf",sep=""))
print(CombinePlots(p2))
dev.off()






