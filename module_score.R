### SCRIPT: addModuleScore() function from Seurat for enrichement analysis CIRCULATION data in IBEX

## 08.04.22 Laura Sudupe , git @lsudupe

## https://www.waltermuskovic.com/2021/04/15/seurat-s-addmodulescore-function/


###### libraries
library("Seurat")
library("GSEABase")
library("ggplot2")
library("escape")
library("dittoSeq")
library("dichromat")


###### Read the gene sets
#genes.velocity <- read.csv("./data/genes_velocity.csv", header = TRUE, sep = ",")
#B <- read.table("./data/clustB_signature.txt", header = FALSE)
#B <- as.vector(B$V1)
#D1 <- as.vector(genes.velocity$Dynamics_1)
#D2 <- as.vector(genes.velocity$Dynamics_2)
#D2[D2 == ""] <- NA 
#D2 <- D2[!is.na(D2)]

##create list
#sig.list <- list(B, D1, D2)
#names(sig.list) <- c("B", "D1", "D2")
#saveRDS(sig.list, "./data/sig.list.rds")

##our gene list
sig.list <- readRDS("./data/sig.list.rds")
B <- sig.list[1]
D1 <- sig.list[2]
D2 <- sig.list[3]
               
###### Download  data
#a <- readRDS("../data/circulation_spatial/integrated.gfp.b.rds")
fibro <- readRDS("./objects/integrated/integrated.fb.rds")
no.fibro <- readRDS("./objects/integrated/integrated.nofb.rds")

###### Check intersect genes B and object
b.fibro <- intersect(fibro@assays[["SCT"]]@data@Dimnames[[1]], B)
b.no.fibro <- intersect(no.fibro@assays[["SCT"]]@data@Dimnames[[1]], B)

########subset assay with B genes and create new object for the enrichment analysis
fibro.new.counts <- GetAssayData(fibro@assays[["SCT"]])[b.fibro,]
#fibro[["B"]] <- CreateAssayObject(counts = fibro.new.counts)
fibro.new <- CreateSeuratObject(counts =fibro.new.counts, assay = "RNA")

no.fibro.new.counts <- GetAssayData(no.fibro@assays[["SCT"]])[b.no.fibro,]
#no.fibro[["B"]] <- CreateAssayObject(counts = no.fibro.new.counts)
no.fibro.new <- CreateSeuratObject(counts =no.fibro.new.counts, assay = "RNA")

#objects <- c(fibro, no.fibro)
objects <- c(fibro, no.fibro)
names(objects) <- c("fibro", "no.fibro")

###############################################333
for (i in 1:length(objects)){
  ###matrix
  a <- objects[[i]]
  #a@assays[["SCT"]]@counts <- a@assays[["SCT"]]@data
  ###### addModuleScore score 
  a <- AddModuleScore(a, features = B, assay = "SCT", name = "geneSetB_")
  a <- AddModuleScore(a, features = D1, assay = "SCT", name = "geneSetD1_")
  a <- AddModuleScore(a, features = D2, assay = "SCT", name = "geneSetD2_")
  #extract module Scores
  d1 <- as.vector(a@meta.data[["geneSetD1_1"]])
  d2 <- as.vector(a@meta.data[["geneSetD2_1"]])
  #calculate relation d1/d2
  d1.d2 <- log10(d1/d2)
  #is.na(d1.d2) <-sapply(d1.d2, is.infinite)
  d1.d2[is.na(d1.d2)] = 0
  a@meta.data$relation_log.d1.d2 <- d1.d2
  ##save object
  saveRDS(a,file = paste0("./results/module_score/all_genes/",names(objects[i]),".rds"))
}

###################################################

fibro.new <- readRDS("./results/module_score/all_genes/fibro.rds")
no.fibro.new <- readRDS("./results/module_score/all_genes/no.fibro.rds")

####add image info
fibro.new@images <- fibro@images
no.fibro.new@images <- no.fibro.new@images
fibro.new@meta.data[["ident"]] <- fibro@meta.data[["ident"]]
no.fibro.new@meta.data[["ident"]] <- no.fibro@meta.data[["ident"]]
fibro.new@meta.data[["sample"]] <- fibro@meta.data[["sample"]]
no.fibro.new@meta.data[["sample"]] <- no.fibro@meta.data[["sample"]]
fibro.new@reductions <- fibro@reductions
no.fibro.new@reductions <- no.fibro@reductions

#Density plots

dens <- c("geneSetD1_1","geneSetD2_1", "relation_log.d1.d2","geneSetB_1")

##fibro
for (i in dens){
  pdf(file.path("./results/module_score/all_genes/",filename = paste(i,"dens.fibro.pdf",sep="")))
  print(dittoRidgePlot(fibro.new, i, group.by = "sample"))
  dev.off()
}

##no fibro
for (i in dens){
  pdf(file.path("./results/module_score/all_genes/",filename = paste(i,"dens.no.fibro.pdf",sep="")))
  print(dittoRidgePlot(no.fibro.new, i, group.by = "sample"))
  dev.off()
}


####Check best values to plot spatial STANDARIZED#############################3
meta <- fibro.new@meta.data
dpi3 <- meta[grepl("dpi3", meta[,6]),]
dpi3_b <- as.vector(dpi3$geneSetB)
dpi5_female <- meta[grepl("dpi5_female", meta[,6]),]
dpi5_male <- meta[grepl("dpi5_male", meta[,6]),]

library(robustHD)
library(data.table)
dpi3_b_z <- standardize(dpi3_b, centerFun = mean, scaleFun = sd)
hist(dpi3_b_z)

meta <- fibro.new@meta.data
meta.d1.d2 <- meta$d1.d2
meta.z.d1.d2 <- standardize(meta.d1.d2, centerFun = mean, scaleFun = sd)
hist(meta.z.d1.d2)
fibro.new@meta.data$d1.d2 <- meta.z.d1.d2

#Spatial plots
bl <- colorRampPalette(c("navy","royalblue","lightskyblue"))(200)                      
re <- colorRampPalette(c("mistyrose", "red2","darkred"))(200)
library("GiNA")
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

feature.list <- c("geneSetD1_1","geneSetD2_1","geneSetB_1")
b <- c(0,0.7)
feature.list <- c("d1.d2")
b <- c(-2,2)
feature.list <- c("relation_log.d1.d2")
b <- c(-2,2)

##fibro
for (i in feature.list){
  pdf(file.path("./results/module_score/all_genes/",filename = paste(i,"dens.fibro.z.pdf",sep="")))
  print(dittoRidgePlot(fibro.new, i, group.by = "sample"))
  dev.off()
}

##################################

for (i in feature.list){
  p1 <- SpatialFeaturePlot(fibro.new, features = i, combine = FALSE)
  fix.p1 <- scale_fill_gradientn(colours=c(bl,"white", re), na.value = "grey98",limits = b,breaks=b, labels = c("min", "max"))
  #fix.p1 <- scale_fill_gradientn(colours=myPalette(100),breaks=b, labels = c("min", "max"),limits = b)
  p2 <- lapply(p1, function (x) x + fix.p1)
  
  pdf(file.path("./results/module_score/all_genes/",filename = paste(i,"spatial.pdf",sep="")))
  print(CombinePlots(p2))
  dev.off()
}



