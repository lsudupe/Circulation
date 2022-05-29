## SCRIPT: differential expression analysis individual segmentation spatial samples CIRCULATION project

## 29.06.22 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(robustHD)
library(data.table)


#Data---------------------------------
#control <- readRDS("./results/individual/control.rds")
dpi3 <- readRDS("./objects/individual/segmentation/dpi3.seg.rds")
dpi5_female <- readRDS("./objects/individual/segmentation/dpi5_female.seg.rds")
dpi5_male <- readRDS("./objects/individual/segmentation/dpi5_male.seg.rds")

#####Change labels to get only IZ, BZ, RZ

########################dpi3

Idents(dpi3) <- (dpi3@meta.data[["segmen"]])
#Levels: BZ_1 BZ_2_3 BZ_4 IZ RZ
dpi3 <- RenameIdents(object = dpi3, `BZ_1` = "BZ",
                                    `BZ_2_3` = "BZ",`BZ_4` = "BZ",`IZ` = "IZ",
                                    `RZ` = "RZ")
###Markers dpi3

dpi3_IZ_RZ <- Seurat::FindMarkers(object = dpi3, 
                                  assay = "SCT",ident.1 = "IZ", ident.2 = "RZ")
dpi3_IZ_BZ <- Seurat::FindMarkers(object = dpi3, 
                                  assay = "SCT",ident.1 = "IZ", ident.2 = "BZ")

saveRDS(dpi3_IZ_RZ, "./results/DE/dpi3_IZ_RZ.rds")
saveRDS(dpi3_IZ_BZ, "./results/DE/dpi3_IZ_BZ.rds")

write.csv(dpi3_IZ_RZ,"./results/DE/dpi3_IZ_RZ.csv", row.names = TRUE)
write.csv(dpi3_IZ_BZ,"./results/DE/dpi3_IZ_BZ.csv", row.names = TRUE)

########################dpi5_female
unique(dpi5_female@meta.data[["segmen"]])
Idents(dpi5_female) <- (dpi5_female@meta.data[["segmen"]])
#Levels: BZ_1 BZ_2_3 BZ_4 IZ RZ
dpi5_female <- RenameIdents(object = dpi5_female, `BZ_1` = "BZ",
                     `BZ_2` = "BZ",`BZ_3` = "BZ",`BZ_4` = "BZ",`IZ` = "IZ",
                     `RZ` = "RZ")
###Markers dpi5_female

dpi5_female_IZ_RZ <- Seurat::FindMarkers(object = dpi5_female, 
                                  assay = "SCT",ident.1 = "IZ", ident.2 = "RZ")
dpi5_female_IZ_BZ <- Seurat::FindMarkers(object = dpi5_female, 
                                  assay = "SCT",ident.1 = "IZ", ident.2 = "BZ")

saveRDS(dpi5_female_IZ_RZ, "./results/DE/dpi5_female_IZ_RZ.rds")
saveRDS(dpi5_female_IZ_BZ, "./results/DE/dpi5_female_IZ_BZ.rds")

write.csv(dpi5_female_IZ_RZ,"./results/DE/dpi5_female_IZ_RZ.csv", row.names = TRUE)
write.csv(dpi5_female_IZ_BZ,"./results/DE/dpi5_female_IZ_BZ.csv", row.names = TRUE)


########################dpi5_male
unique(dpi5_male@meta.data[["segmen"]])
Idents(dpi5_male) <- (dpi5_male@meta.data[["segmen"]])
#Levels: BZ_1 BZ_2_3 BZ_4 IZ RZ
dpi5_male <- RenameIdents(object = dpi5_male, `BZ_1` = "BZ",
                            `BZ_2` = "BZ",`BZ_3` = "BZ",`BZ_4` = "BZ",`IZ` = "IZ",
                            `RZ` = "RZ")
###Markers dpi5_male

dpi5_male_IZ_RZ <- Seurat::FindMarkers(object = dpi5_male, 
                                         assay = "SCT",ident.1 = "IZ", ident.2 = "RZ")
dpi5_male_IZ_BZ <- Seurat::FindMarkers(object = dpi5_male, 
                                         assay = "SCT",ident.1 = "IZ", ident.2 = "BZ")

saveRDS(dpi5_male_IZ_RZ, "./results/DE/dpi5_male_IZ_RZ.rds")
saveRDS(dpi5_male_IZ_BZ, "./results/DE/dpi5_male_IZ_BZ.rds")

write.csv(dpi5_male_IZ_RZ,"./results/DE/dpi5_male_IZ_RZ.csv", row.names = TRUE)
write.csv(dpi5_male_IZ_BZ,"./results/DE/dpi5_male_IZ_BZ.csv", row.names = TRUE)



#Read objects and top
#de <- readRDS("./results/integrated_label_markers.rds")
#Filter
#de_subset<- subset(de, p_val_adj < 0.05 & 0.5 < avg_log2FC)

#Top5
#top3.de <- de_subset %>%
 # group_by(cluster) %>%
  #top_n(n = 3,
   #     wt = avg_log2FC)

#Save
#write.csv(top3.de,"./results/genes/top3de.csv", row.names = FALSE)