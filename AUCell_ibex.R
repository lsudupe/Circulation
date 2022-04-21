### SCRIPT: AUCell method enrichement analysis Tokio data in IBEX

## 03.03.22 Laura Sudupe , git @lsudupe

#https://bioconductor.org/packages/release/bioc/vignettes/AUCell/inst/doc/AUCell.html
#https://rdrr.io/bioc/AUCell/man/GeneSet-methods.html

###### libraries
library("Seurat")
library("GSEABase")
library("AUCell")
library("ggplot2")
library("escape")
library("dittoSeq")
library("dichromat")

##CHOOSE THE auxMaxRank
genes.porcentage <-  210#porcentage depending number of genes in spatial

###### Read the gene sets
genes.velocity <- read.csv("./data/genes_velocity.csv", header = TRUE, sep = ",")
B <- read.table("./data/clustB_signature.txt", header = FALSE)
B <- as.vector(B$V1)
D1 <- as.vector(genes.velocity$Dynamics_1)
D2 <- as.vector(genes.velocity$Dynamics_2)
D2[D2 == ""] <- NA 
D2 <- D2[!is.na(D2)]

###### Create genesets
#b.sig <- GeneSet(B, setName="geneSetB")
d1.sig <- GeneSet(D1, setName="geneSetD1")
d2.sig <- GeneSet(D2, setName="geneSetD2")
geneSets <- GeneSetCollection(d1.sig, d2.sig)

###### Download  data
#a <- readRDS("../data/circulation_spatial/integrated.gfp.b.rds")
fibro <- readRDS("./objects/integrated/integrated.gfp.fb.rds")
no.fibro <- readRDS("./objects/integrated/integrated.gfp.nofb.rds")

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

objects <- c(fibro.new, no.fibro.new)
names(objects) <- c("fibro", "no.fibro")

###############################################333
for (i in 1:length(objects)){
          ###matrix
          a <- objects[[i]]
          a@assays[["RNA"]]@counts <- a@assays[["RNA"]]@data
          matrix <- a@assays[["RNA"]]@counts
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
          ##################ssGSEA#############3
          ###### Enrichment
          #i.ES <- enrichIt(obj = matrix, 
          #                 gene.sets = geneSets, 
          #                 groups = 1000, cores = 4)
          ##save meta
          #a <- AddMetaData(a, i.ES)
          #d1 <- as.vector(a@meta.data[["geneSetD1"]])
          #d2 <- as.vector(a@meta.data[["geneSetD2"]])
          #d1.d2 <- log10(d1/d2)
          #d1.d2[is.na(d1.d2)] = 0
          #a@meta.data[["relation_log(d1.d2)"]] <- d1.d2
          ##save object
          saveRDS(a,file = paste0(names(objects[i]),"B.rds"))
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

fibro.new <- readRDS("./fibro.rds")
no.fibro.new <- readRDS("./no.fibroB.rds")

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

dens <- c("geneSetD1","geneSetD2", "relation_log.d1.d2")

##fibro
for (i in dens){
  pdf(paste(i,"dens.fibro.all.B.pdf",sep=""))
  print(dittoRidgePlot(fibro.new, i, group.by = "sample"))
  dev.off()
}

##no fibro
for (i in dens){
  pdf(paste(i,"dens.nofibro.B.pdf",sep=""))
  print(dittoRidgePlot(no.fibro.new, i, group.by = "sample"))
  dev.off()
}

#Spatial plots
bl <- colorRampPalette(c("navy","royalblue","lightskyblue"))(200)                      
re <- colorRampPalette(c("mistyrose", "red2","darkred"))(200)

feature.list <- c("geneSetD1","geneSetD2","relation_log.d1.d2")
b <- c(0.4,0.65)
feature.list <- c("relation_log.d1.d2")
l <- c(-0.3,0.5)


for (i in feature.list){
  p1 <- SpatialFeaturePlot(fibro.new, features = i, combine = FALSE)
  fix.p1 <- scale_fill_gradientn(colours=c(bl,"white", re), na.value = "grey98",limits = l,breaks=l, labels = c("min", "max"))
  #fix.p1 <- scale_fill_gradientn(colours=myPalette(100),breaks=b, labels = c("min", "max"),limits = b)
  p2 <- lapply(p1, function (x) x + fix.p1)
  
  pdf(paste(i,"spatial.all.pdf",sep=""))
  print(CombinePlots(p2))
  dev.off()
}

p1 <- SpatialFeaturePlot(fibro.new, features = "relation_log(d1.d2)")
fix.p1 <- scale_fill_gradientn(colours=c(bl,"white", re), na.value = "grey98",limits = c(-4,3))
#fix.p1 <- scale_fill_continuous(type = "viridis")
p2 <- lapply(p1, function (x) x + fix.p1)

pdf("d1.d2.spatial.ssGSEA.pdf")
print(CombinePlots(p2))
dev.off()

#Density plots in clusters

##fibro
for (i in dens){
  pdf(paste(i,"clusterdens.fibro.pdf",sep=""))
  print(ridgeEnrichment(fibro.new@meta.data, gene.set = i, group = "ident", facet = "sample", add.rug = TRUE))
  dev.off()
}

##no fibro
for (i in dens){
  pdf(paste(i,"clusterdens.nofibro.pdf",sep=""))
  print(ridgeEnrichment(no.fibro.new@meta.data, gene.set = i, group = "ident", facet = "sample", add.rug = TRUE))
  dev.off()
}

library(dplyr)
library(ggplot2)
library(hrbrthemes)

fibro.new.meta <- fibro.new@meta.data
fibro.new.control <- fibro.new.meta[grepl("1.control", fibro.new.meta[,4]),]
fibro.new.3dpi <- fibro.new.meta[grepl("2.3dpi", fibro.new.meta[,4]),]
fibro.new.5.dpi.female <- fibro.new.meta[grepl("3.5dpi_female", fibro.new.meta[,4]),]
fibro.new.5.dpi.male <- fibro.new.meta[grepl("4.5dpi_male", fibro.new.meta[,4]),]


# "geneSetB","geneSetD1","geneSetD2","relation_log.d1.d2")


####boxplot

pdf("boxplot.geneSetB.pdf")
fibro.new.meta%>% 
ggplot(aes(x=sample, y= geneSetB, fill=ident)) + 
geom_boxplot() +  
theme_classic() +
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
theme(plot.title = element_text(hjust=0.5, face="bold"))
dev.off()

####scaterplot

pdf("scaterplott.5dpi_male.d1.d22.pdf")
ggplot(fibro.new.5.dpi.male, aes(x=geneSetD2, y=relation_log.d1.d2, shape=ident, color=ident)) +
  geom_point() + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)
dev.off()

library("GiNA")
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

# "Postn","Aspn"
p1 <- SpatialFeaturePlot(object = fibro.new, 
                         features = c("Cthrc1"),
                         combine = FALSE) 
fix.p1 <- scale_fill_gradientn(colors=myPalette(100),
                               limits = c(0,4),
                               breaks=c(0,4),
                               labels=c("Min","Max"),)
p2 <- lapply(p1, function (x) x + fix.p1)

pdf("c.pdf")
print(CombinePlots(p2))
dev.off()

SpatialPlot(fibro.new,features = c("Postn"))

pdf("spatial_integrated.pdf")
SpatialDimPlot(object = fibro.new ,group.by = c("ident"), combine = FALSE)
dev.off()

