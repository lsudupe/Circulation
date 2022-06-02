## SCRIPT: individual object creation and first analaysis spatial samples CIRCULATION project

## 05.04.22 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(dplyr)
library("GSEABase")
library("AUCell")

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

######plot to check
DimPlot(control, group.by = c("seurat_clusters"), label = T) + ggtitle("UMAP_r0.4")
SpatialPlot(dpi5_male, group.by = c("seurat_clusters"),label = TRUE, combine = FALSE)

###### Read the gene sets###################
genes.velocity <- read.csv("./data/genes_velocity.csv", header = TRUE, sep = ",")
B <- read.table("./data/clustB_signature.txt", header = FALSE)
B <- as.vector(B$V1)
D1 <- as.vector(genes.velocity$Dynamics_1)
D2 <- as.vector(genes.velocity$Dynamics_2)
D2[D2 == ""] <- NA 
D2 <- D2[!is.na(D2)]

###### Create genesets
b.sig <- GeneSet(B, setName="geneSetB")
d1.sig <- GeneSet(D1, setName="geneSetD1")
d2.sig <- GeneSet(D2, setName="geneSetD2")


#Create geneSet for FB
fb.genes <- read.csv("./results/DE/top200_tokio_FB.csv")
fb.genes <- as.vector(fb.genes$gene)
fb.sig <- GeneSet(fb.genes, setName="geneSetFB")
geneSets <- GeneSetCollection(d1.sig, d2.sig,b.sig, fb.sig)

objects <- c(control, dpi3, dpi5_female, dpi5_male)
names(objects) <- c("control", "dpi3","dpi5_female","dpi5_male")

###############################################333
for (i in 1:length(objects)){
  ###matrix
  a <- objects[[i]]
  a@assays[["SCT"]]@counts <- a@assays[["SCT"]]@data
  matrix <- a@assays[["SCT"]]@counts
  matrix <- as.matrix(matrix)
  ###### AUC score 
  cells_rankings <- AUCell_buildRankings(matrix, nCores=1)#, plotStats=TRUE)
  #rankings <- getRanking(cells_rankings)
  cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings,aucMaxRank=nrow(cells_rankings)*0.25)
  #extract AUC values
  auc_per_cell_all <- as.data.frame(t(getAUC(cells_AUC)))
  #calculate relation d1/d2
  d1 <- as.vector(auc_per_cell_all$geneSetD1)
  d2 <- as.vector(auc_per_cell_all$geneSetD2)
  d1.d2 <- log10(d1/d2)
  is.na(d1.d2) <-sapply(d1.d2, is.infinite)
  d1.d2[is.na(d1.d2)] = 0
  auc_per_cell_all$relation_log.d1.d2 <- d1.d2
  ##save meta
  a <- AddMetaData(a, auc_per_cell_all)
  saveRDS(a,file = paste0("./results/individual/",names(objects[i]),".rds"))
  
}
#Genes in the gene sets NOT available in the dataset: 
#  geneSetD1: 	8 (14% of 57)
#geneSetD2: 	5 (12% of 43)
#geneSetB: 	31 (14% of 216)
#Genes in the gene sets NOT available in the dataset: 
#  geneSetD1: 	2 (4% of 57)
#geneSetB: 	11 (5% of 216)
#Genes in the gene sets NOT available in the dataset: 
#  geneSetD1: 	1 (2% of 57)
#geneSetB: 	7 (3% of 216)
#Genes in the gene sets NOT available in the dataset: 
#  geneSetD1: 	1 (2% of 57)
#geneSetB: 	6 (3% of 216)

control <- readRDS("./results/individual/control.rds")
dpi3 <- readRDS("./results/individual/dpi3.rds")
dpi5_female <- readRDS("./results/individual/dpi5_female.rds")
dpi5_male <- readRDS("./results/individual/dpi5_male.rds")




