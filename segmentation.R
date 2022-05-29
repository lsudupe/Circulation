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
enrichment.meta <- read.csv( "./data/dpi3/Enrichment.csv")
enrichment.meta <- as.vector(enrichment.meta$Enrichment)

dpi3@meta.data["segmen"] <- as.factor(enrichment.meta)
d1d2<- as.vector(dpi3@meta.data[["relation_log.d1.d2"]])
d1d2<- standardize(d1d2, centerFun = mean, scaleFun = sd)
dpi3@meta.data["d1d2"] <- d1d2

saveRDS(dpi3, "./objects/individual/segmentation/dpi3.seg.rds")

######dpi5_female
enrichment.meta <- read.csv( "./data/dpi5_female/enrichment_2.csv")
enrichment.meta <- as.vector(enrichment.meta$Enrichment)

dpi5_female@meta.data["segmen"] <- as.factor(enrichment.meta)
d1d2<- as.vector(dpi5_female@meta.data[["relation_log.d1.d2"]])
d1d2<- standardize(d1d2, centerFun = mean, scaleFun = sd)
dpi5_female@meta.data["d1d2"] <- d1d2

saveRDS(dpi5_female, "./objects/individual/segmentation/dpi5_female.seg.rds")

######dpi5_male
enrichment.meta <- read.csv( "./data/dpi5_male/Enrichment.csv")
enrichment.meta <- as.vector(enrichment.meta$Enrichment)

dpi5_male@meta.data["segmen"] <- as.factor(enrichment.meta)
d1d2<- as.vector(dpi5_male@meta.data[["relation_log.d1.d2"]])
d1d2<- standardize(d1d2, centerFun = mean, scaleFun = sd)
dpi5_male@meta.data["d1d2"] <- d1d2

saveRDS(dpi5_male, "./objects/individual/segmentation/dpi5_male.seg.rds")

#####################


#Create a custom color scale
library(RColorBrewer)
myColors <- brewer.pal(7,"Paired")
#myColors <- as.vector(names(myColors))
seg <- as.factor(c("IZ","BZ_1","BZ_2","BZ_3","BZ_2_3","BZ_4","RZ"))
names(myColors) <- levels(seg)
#colScale <- scale_colour_manual(name = "grp",values = myColors)

nombres <- c("IZ","BZ_1","BZ_2","BZ_3","BZ_2_3","BZ_4","RZ")
colors <- c("coral4", "skyblue", "red", "darkgreen", "yellow","khaki2", "darkolivegreen", "gray89")
names(colors)  <- nombres

####################dpi3
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

####################dpi5_female

pdf(file.path("./results/individual/dpi5_male",filename = paste("prueba..pdf",sep="")))
dpi3@meta.data%>% 
  ggplot(aes(x=geneSetB, y= segmen, fill=segmen)) + 
  geom_boxplot(aes(fill=segmen)) +  
  scale_fill_manual(values =colors ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  #xlim(-0.5, 0.5) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) 
dev.off()

pdf(file.path("./results/individual/dpi5_male",filename = paste("spatial.seg.pdf",sep="")))
SpatialPlot(dpi5_male, group.by = c("segmen"),label = TRUE, cols = myColors,combine = FALSE)
dev.off()

####################dpi5_male

pdf(file.path("./results/individual/dpi5_male",filename = paste("prueba..pdf",sep="")))
  dpi5_male@meta.data%>% 
  ggplot(aes(x=d1d2, y= segmen, fill=segmen)) + 
  geom_boxplot() +  
  scale_color_manual(colScale) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  xlim(-0.5, 0.5) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) 
dev.off()

