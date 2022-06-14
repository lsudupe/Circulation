## SCRIPT: individual object AUCell analysis B,Da,Db,Da/Db and FB genesets CIRCULATION project

## 13.06.22 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(dplyr)
library("GSEABase")
library("AUCell")

###Read data
control <- readRDS("./objects/individual/control.rds")
dpi3 <- readRDS("./objects/individual/dpi3.rds")
dpi5_female <- readRDS("./objects/individual/dpi5_female.rds")
dpi5_male <- readRDS("./objects/individual/dpi5_male.rds")

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
d1.sig <- GeneSet(D1, setName="geneSetDA")
d2.sig <- GeneSet(D2, setName="geneSetDB")

#Create geneSet for FB
fb.genes <- read.csv("./results/DE/top200_tokio_FB.csv")
fb.genes <- as.vector(fb.genes$gene)
fb.sig <- GeneSet(fb.genes, setName="geneSetFB")

###### Create geneSet nine
marina_nine <- read.delim("./data/cluster9_markers.txt")
#Top100 marina
top100 <- marina_nine %>%
  #group_by(cluster) %>%
  top_n(n = 100,
        wt = avg_log2FC)

top100 <- cbind(genes = rownames(top100), top100)
rownames(top100) <- 1:nrow(top100)
top100 <- as.vector(top100$genes)

###### Create genesets
nine.sig <- GeneSet(top100, setName="geneSet9")

geneSets <- GeneSetCollection(d1.sig, d2.sig,b.sig, fb.sig, nine.sig)


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
  d1 <- as.vector(auc_per_cell_all$geneSetDA)
  d2 <- as.vector(auc_per_cell_all$geneSetDB)
  d1.d2 <- log10(d1/d2)
  is.na(d1.d2) <-sapply(d1.d2, is.infinite)
  d1.d2[is.na(d1.d2)] = 0
  auc_per_cell_all$ratio_dinamics <- d1.d2
  ##save meta
  a <- AddMetaData(a, auc_per_cell_all)
  saveRDS(a,file = paste0("./results/individual/",names(objects[i]),".enrich.rds"))
  
}

control <- readRDS("./results/individual/control.enrich.rds")
dpi3 <- readRDS("./results/individual/dpi3.enrich.rds")
dpi5_female <- readRDS("./results/individual/dpi5_female.enrich.rds")
dpi5_male <- readRDS("./results/individual/dpi5_male.enrich.rds")

#####ratio estand

b <- standardize(control@meta.data[["ratio_dinamics"]], centerFun = mean, scaleFun = sd)
control@meta.data[["ratio_stand"]] <- b

b <- standardize(dpi3@meta.data[["ratio_dinamics"]], centerFun = mean, scaleFun = sd)
dpi3@meta.data[["ratio_stand"]] <- b

b <- standardize(dpi5_female@meta.data[["ratio_dinamics"]], centerFun = mean, scaleFun = sd)
dpi5_female@meta.data[["ratio_stand"]] <- b

b <- standardize(dpi5_male@meta.data[["ratio_dinamics"]], centerFun = mean, scaleFun = sd)
dpi5_male@meta.data[["ratio_stand"]] <- b


saveRDS(control, "./results/individual/control.enrich.rds")
saveRDS(dpi3, "./results/individual/dpi3l.enrich.rds")
saveRDS(dpi5_female, "./results/individual/dpi5_female.enrich.rds")
saveRDS(dpi5_male, "./results/individual/dpi5_male.enrich.rds")