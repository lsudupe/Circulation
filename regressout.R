## SCRIPT: regress out the unwanted contribution spatial samples CIRCULATION project

## 02.05.22 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library("GSEABase")
library("AUCell")

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
write.csv(top.200.de,"./results/DE/top200_tokio_FB.csv", row.names = FALSE)
fb.genes <- read.csv("./results/DE/top200_tokio_FB.csv")

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


####Read the data
#control <- readRDS("./objects/individual/segmentation/dpi5_male.seg.rds")

dpi3 <- readRDS("./objects/individual/segmentation/dpi3.seg.rds")
dpi5_female <- readRDS("./objects/individual/segmentation/dpi5_female.seg.rds")
dpi5_male <- readRDS("./objects/individual/segmentation/dpi5_male.seg.rds")


objects <- c(dpi3, dpi5_female, dpi5_male )
names(objects) <- c("dpi3", "dpi5_female", "dpi5_male")

###############################################333
for (i in 1:length(objects)){
  ###matrix
  a <- objects[[i]]
  meta <- a@meta.data
  lm <- lm(meta$geneSetB ~ meta$geneSetFB, data =meta)
  residuals <- lm$residuals
  a@meta.data[["residualsB"]] <- residuals
  lm <- lm(meta$geneSetD1 ~ meta$geneSetFB, data =meta)
  residuals <- lm$residuals
  a@meta.data[["residualsD1"]] <- residuals
  lm <- lm(meta$geneSetD2 ~ meta$geneSetFB, data =meta)
  residuals <- lm$residuals
  a@meta.data[["residualsD2"]] <- residuals
  saveRDS(a,file = paste0("./results/results_regressout/",names(objects[i]),".B.rds"))
  #feature_plot
  #p1 <- SpatialFeaturePlot(a, features = "nCount_SCT", combine = FALSE)
  #fix.p1 <- scale_fill_continuous(limits = c(0,25000), 
  #                     	      breaks = c(0,25000),
  #                              type = "viridis")
  #p2 <- lapply(p1, function (x) x + fix.p1)
  
  #pdf(paste(names(objects[i]),"count.dens.pdf",sep=""))
  #pdf("count.dens.pdf")
  #print(CombinePlots(p2))
  #dev.off()
}         

###################################################PLOTS

dpi3 <- readRDS("./results/results_regressout/dpi3.B.rds")
dpi5_female <- readRDS("./results/results_regressout/dpi5_female.B.rds")
dpi5_male <- readRDS("./results/results_regressout/dpi5_male.B.rds")


###################Spatial plots
install.packages("colorspace")
#https://cran.r-project.org/web/packages/colorspace/vignettes/colorspace.html#Installation
library(colorspace)
q4 <- colorspace::sequential_hcl(7, palette = "Plasma")

feature.list <- c("geneSetB")
b <- c(0,0.8)

#########B

#for (i in feature.list){
  p1 <- SpatialFeaturePlot(dpi3, features = "geneSetB", combine = FALSE)
  #fix.p1 <- scale_fill_gradientn(colours=c(bl,"white", re), na.value = "grey98",limits = b,breaks=b, labels = c("min", "max"))
  fix.p1 <- scale_fill_gradientn(colours=q4,breaks=b, labels = c("min", "max"),limits = b)
  p2 <- lapply(p1, function (x) x + fix.p1)
  
#  pdf(file.path("./results_regressout/",filename = paste(i,"spatial.resB.dpi3.pdf",sep="")))
  print(CombinePlots(p2))
#  dev.off()
#}







