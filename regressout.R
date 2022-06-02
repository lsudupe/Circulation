## SCRIPT: regress out the unwanted contribution spatial samples CIRCULATION project

## 02.05.22 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)

#Read sc data Tokio
tokio <- readRDS("./objects/sc/Tokio/sc.combined.sct.fibro.rds")

tokio_de <- Seurat::FindAllMarkers(object = tokio)

saveRDS(tokio_de, "./results/DE/tokio_FB_noFB.rds")

#Read objects and top
de <- readRDS("./results/DE/tokio_FB_noFB.rds")
#Filter
de_subset<- subset(de, p_val_adj < 0.05 & 0.5 < avg_log2FC)
#Top5
top.100.de <- de_subset %>%
 group_by(cluster) %>%
top_n(n = 100,
     wt = avg_log2FC)
#Save
write.csv(top100.de,"./results/genes/top100de_FB_noFB.csv", row.names = FALSE)



sctransform::

