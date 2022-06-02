## SCRIPT: regress out the unwanted contribution spatial samples CIRCULATION project

## 02.05.22 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(dplyr)
library("GSEABase")

#Read sc data Tokio
tokio <- readRDS("./objects/sc/Tokio/sc.combined.sct.fibro.rds")


########################select genes for the FB signature

##diferential expression
tokio_de <- Seurat::FindAllMarkers(object = tokio)
saveRDS(tokio_de, "./results/DE/tokio_FB_noFB.rds")

#Read objects and select top
de <- readRDS("./results/DE/tokio_FB_noFB.rds")

#take our no-fb data
de <- de[!grepl("no-FB", de[,6]),]
#Filter
de_subset<- subset(de, p_val_adj < 0.05 & 0.5 < avg_log2FC)
#Top100
top.200.de <- de_subset %>%
 group_by(cluster) %>%
top_n(n = 200,
     wt = avg_log2FC)
#Save
write.csv(top.200.de,"./results/DE/top200_tokio_FB.rds", row.names = FALSE)

#Create geneSet for FB
fb.genes <- as.vector(top.200.de$gene)
fb.sig <- GeneSet(fb.genes, setName="geneSetFB")

###### Read the gene sets
genes.velocity <- read.csv("./data/genes_velocity.csv", header = TRUE, sep = ",")
B <- read.table("./data/clustB_signature.txt", header = FALSE)
B <- as.vector(B$V1)
D1 <- as.vector(genes.velocity$Dynamics_1)
D2 <- as.vector(genes.velocity$Dynamics_2)
D2[D2 == ""] <- NA 
D2 <- D2[!is.na(D2)]

####Check how many genes are the same in each signature
FB.B <- b.fibro <- intersect(fb.genes, B)
FB.D1 <- b.fibro <- intersect(fb.genes, D1)
FB.D2 <- b.fibro <- intersect(fb.genes, D2)

##########################################3

control <- readRDS("./results/individual/control.rds")
dpi3 <- readRDS("./results/individual/dpi3.rds")
dpi5_female <- readRDS("./results/individual/dpi5_female.rds")
dpi5_male <- readRDS("./results/individual/dpi5_male.rds")

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
  cells_AUC <- AUCell_calcAUC(fb.sig, cells_rankings,aucMaxRank=genes.porcentage)
  #extract AUC values
  auc_per_cell_all <- as.data.frame(t(getAUC(cells_AUC)))
  ##save meta
  a <- AddMetaData(a, auc_per_cell_all)
  saveRDS(a,file = paste0("./results/individual/",names(objects[i]),"FB.rds"))
  
}

