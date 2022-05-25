## SCRIPT: individual analysis AUCell and segmentation spatial samples CIRCULATION project

## 05.04.22 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(robustHD)
library(data.table)


#Data---------------------------------
#control <- readRDS("./results/individual/control.rds")
dpi3 <- readRDS("./results/individual/dpi3.rds")
dpi5_female <- readRDS("./results/individual/dpi5_female.rds")
dpi5_male <- readRDS("./results/individual/dpi5_male.rds")



######dpi3
enrichment.meta <- read.csv( "./data/infarto_3dpi/Enrichment.csv")
enrichment.meta <- as.vector(enrichment.meta$Enrichment)

dpi3e@meta.data["segmen"] <- as.factor(enrichment.meta)
d1d2<- as.vector(dpi3@meta.data[["relation_log.d1.d2"]])


d1d2<- standardize(d1d2, centerFun = mean, scaleFun = sd)
dpi3@meta.data["d1d2"] <- d1d2

SpatialPlot(dpi3, group.by = c("segmen"),label = TRUE, combine = FALSE)

pdf(file.path("./results/individual/dpi3",filename = paste("d1d2.boxplot.pdf",sep="")))
dpi3@meta.data%>% 
  ggplot(aes(x=d1d2, y= segmen, fill=segmen)) + 
  geom_boxplot() +  
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  xlim(-0.5, 0.5) +
  theme(plot.title = element_text(hjust=0.5, face="bold"))
dev.off()

######dpi5_male
enrichment.meta <- read.csv( "./data/infarto_5dpi_macho/Enrichment.csv")
enrichment.meta <- as.vector(enrichment.meta$Enrichment)

dpi5_male@meta.data["segmen"] <- as.factor(enrichment.meta)
d1d2 <- as.vector(dpi5_male@meta.data[["relation_log.d1.d2"]])


d1d2<- standardize(d1d2, centerFun = mean, scaleFun = sd)
dpi5_male@meta.data["d1d2"] <- d1d2

SpatialPlot(dpi5_male, group.by = c("segmen"),label = TRUE, combine = FALSE)

pdf(file.path("./results/individual/dpi5_male",filename = paste("d1d2.boxplot.pdf",sep="")))
dpi5_male@meta.data%>% 
  ggplot(aes(x=d1d2, y= segmen, fill=segmen)) + 
  geom_boxplot() +  
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  xlim(-0.5, 0.5) +
  theme(plot.title = element_text(hjust=0.5, face="bold"))
dev.off()
