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

####AUCell score for each object?
###matrix
a <- objects[[i]]
a@assays[["SCT"]]@counts <- a@assays[["SCT"]]@data
matrix <- a@assays[["SCT"]]@counts
matrix <- as.matrix(matrix)
###### AUC score 
cells_rankings <- AUCell_buildRankings(matrix, nCores=1)#, plotStats=TRUE)
#rankings <- getRanking(cells_rankings)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings,aucMaxRank=genes.porcentage)
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


sctransform::

