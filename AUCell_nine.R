### SCRIPT: AUCell method enrichement analysis Tokio data in IBEX

## 03.03.22 Laura Sudupe , git @lsudupe

#https://bioconductor.org/packages/release/bioc/vignettes/AUCell/inst/doc/AUCell.html
#https://rdrr.io/bioc/AUCell/man/GeneSet-methods.html

install.packages("GINA")

###### libraries
library("Seurat")
library("GSEABase")
library("AUCell")
library("ggplot2")
library("escape")
library("dittoSeq")
library("dichromat")
library("dplyr")

##CHOOSE THE auxMaxRank
genes.porcentage <-  1500#porcentage depending number of genes in spatial

###### Read the gene sets
marina_nine <- read.delim("./data/cluster9_markers.txt")
#Top100 marina
top100 <- marina_nine %>%
  #group_by(cluster) %>%
  top_n(n = 100,
        wt = avg_log2FC)

top100 <- top100$ as.vector(nine[,7])
top100 <- cbind(genes = rownames(top100), top100)
rownames(top100) <- 1:nrow(top100)
top100 <- as.vector(top100[,1])


circ <- readRDS( "./results/circulation_markers.rds")
nine <- circ[grepl("9", circ[,6]),]

#Top5
top200 <- nine %>%
  #group_by(cluster) %>%
  top_n(n = 200,
        wt = avg_log2FC)

nine <- as.vector(nine[,7])

#Save
write.csv(nine,"./results/genes/top100nine.csv", row.names = FALSE)
nine <- read.csv("./results/genes/top100nine.csv")
nine <- as.vector(nine[,1])

####Check which nine genes are same as B
ninexb <- intersect(nine, B)
ninexd1 <- intersect(nine, D1)
ninexd2 <- intersect(nine, D2)

###### Take out B genes from nine vector
nine.no.b <- nine [! nine %in% ninexb]
###### Create genesets
nine.sig <- GeneSet(nine.no.b, setName="geneSet9")
#nine.sig <- GeneSet(nine, setName="geneSet9")

###### Download  data
a <- readRDS("./objects/integrated/integrated.sct.rds")
fibro <- readRDS("./objects/integrated/integrated.fb.rds")
no.fibro <- readRDS("./objects/integrated/integrated.nofb.rds")

###### Check intersect genes B and object
b.fibro <- intersect(fibro@assays[["SCT"]]@data@Dimnames[[1]], nine)
b.no.fibro <- intersect(no.fibro@assays[["SCT"]]@data@Dimnames[[1]], nine)

########subset assay with B genes and create new object for the enrichment analysis
fibro.new.counts <- GetAssayData(fibro@assays[["SCT"]])[b.fibro,]
#fibro[["B"]] <- CreateAssayObject(counts = fibro.new.counts)
fibro.new <- CreateSeuratObject(counts =fibro.new.counts, assay = "RNA")

no.fibro.new.counts <- GetAssayData(no.fibro@assays[["SCT"]])[b.no.fibro,]
#no.fibro[["B"]] <- CreateAssayObject(counts = no.fibro.new.counts)
no.fibro.new <- CreateSeuratObject(counts =no.fibro.new.counts, assay = "RNA")

#objects <- c(fibro, no.fibro)
objects <- c(fibro, no.fibro)
names(objects) <- c("fibro", "no.fibro")

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
          cells_AUC <- AUCell_calcAUC(nine.sig, cells_rankings,aucMaxRank=genes.porcentage)
          #extract AUC values
          auc_per_cell_all <- as.data.frame(t(getAUC(cells_AUC)))
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
          saveRDS(a,file = paste0("./results/Aucell/nine_genes/",names(objects[i]),".nine.all.rds"))
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

fibro.new <- readRDS("./results/Aucell/nine_genes/fibro.nine.noB.rds")
no.fibro.new <- readRDS("./results/Aucell/nine_genes/no.fibro.nine.noB.rds")
a <- readRDS("./results/Aucell/nine_genes/NA.nine.all.rds")

#Density plots

dens <- c("geneSet9")

##fibro
for (i in dens){
  pdf(file.path("./results/Aucell/nine_genes/",filename = paste(i,"dens.fibro.nine.noB.pdf",sep="")))
  print(dittoRidgePlot(fibro.new, i, group.by = "sample"))
  dev.off()
}

##no fibro
for (i in dens){
  pdf(file.path("./results/Aucell/nine_genes/",filename = paste(i,"dens.no.fibro.nine.noB.pdf",sep="")))
  print(dittoRidgePlot(no.fibro.new, i, group.by = "sample"))
  dev.off()
}

#Spatial plots
bl <- colorRampPalette(c("navy","royalblue","lightskyblue"))(200)                      
re <- colorRampPalette(c("mistyrose", "red2","darkred"))(200)
library("GiNA")
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

feature.list <- c("geneSet9")
b <- c(0.1,0.36)

for (i in feature.list){
  p1 <- SpatialFeaturePlot(fibro.new, features = i, combine = FALSE)
  #fix.p1 <- scale_fill_gradientn(colours=c(bl,"white", re), na.value = "grey98",limits = b,breaks=b, labels = c("min", "max"))
  fix.p1 <- scale_fill_gradientn(colours=myPalette(100),breaks=b, labels = c("min", "max"),limits = b)
  p2 <- lapply(p1, function (x) x + fix.p1)
  
  pdf(file.path("./results/Aucell/nine_genes/",filename = paste(i,"spatial.nine.noB.fibro.pdf",sep="")))
  print(CombinePlots(p2))
  dev.off()
}

for (i in feature.list){
  p1 <- SpatialFeaturePlot(no.fibro.new, features = i, combine = FALSE)
  #fix.p1 <- scale_fill_gradientn(colours=c(bl,"white", re), na.value = "grey98",limits = b,breaks=b, labels = c("min", "max"))
  fix.p1 <- scale_fill_gradientn(colours=myPalette(100),breaks=b, labels = c("min", "max"),limits = b)
  p2 <- lapply(p1, function (x) x + fix.p1)
  
  pdf(file.path("./results/Aucell/nine_genes/",filename = paste(i,"spatial.nine.noB.no.fibro.pdf",sep="")))
  print(CombinePlots(p2))
  dev.off()
}

for (i in feature.list){
  p1 <- SpatialFeaturePlot(a, features = i, combine = FALSE)
  #fix.p1 <- scale_fill_gradientn(colours=c(bl,"white", re), na.value = "grey98",limits = b,breaks=b, labels = c("min", "max"))
  fix.p1 <- scale_fill_gradientn(colours=myPalette(100),breaks=b, labels = c("min", "max"),limits = b)
  p2 <- lapply(p1, function (x) x + fix.p1)
  
  pdf(file.path("./results/Aucell/nine_genes/",filename = paste(i,"spatial.nine.all.pdf",sep="")))
  print(CombinePlots(p2))
  dev.off()
}

#Density plots in clusters

##fibro
for (i in dens){
  pdf(file.path("./results/Aucell/nine_genes/",filename = paste(i,"dens.cluster.nine.fibro.pdf",sep="")))
  print(ridgeEnrichment(fibro.new@meta.data, gene.set = i, group = "ident", facet = "sample", add.rug = TRUE))
  dev.off()
}

##no fibro
for (i in dens){
  pdf(file.path("./results/Aucell/nine_genes/",filename = paste(i,"dens.cluster.nine.no.fibro.pdf",sep="")))
  print(ridgeEnrichment(no.fibro.new@meta.data, gene.set = i, group = "ident", facet = "sample", add.rug = TRUE))
  dev.off()
}

##all
for (i in dens){
  pdf(file.path("./results/Aucell/nine_genes/",filename = paste(i,"dens.cluster.nine.all.pdf",sep="")))
  print(ridgeEnrichment(a@meta.data, gene.set = i, group = "seurat_clusters", facet = "sample", add.rug = TRUE))
  dev.off()
}

library(dplyr)
library(ggplot2)
library(hrbrthemes)

fibro.new.meta <- fibro.new@meta.data
no.fibro.new.meta <- no.fibro.new@meta.data

# "geneSetB","geneSetD1","geneSetD2","relation_log.d1.d2")


####boxplot

pdf(file.path("./results/Aucell/nine_genes/",filename = paste("fibro.boxplot.nine.pdf")))
fibro.new.meta%>% 
ggplot(aes(x=sample, y= nine, fill=ident)) + 
geom_boxplot() +  
theme_classic() +
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
theme(plot.title = element_text(hjust=0.5, face="bold"))
dev.off()

####scaterplot

pdf(file.path("./results/Aucell/all_genes/",filename = paste("scaterplot.3dpi.d1.d2.pdf")))
ggplot(fibro.new.3dpi, aes(x=geneSetD1, y=geneSetD2, shape=ident, color=ident)) +
  geom_point() + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)
dev.off()

#####Check marker genes
genes <- c("Cthrc1", "Postn","Aspn")
#SpatialPlot(fibro.new,features = c("Postn"))
DefaultAssay(fibro.new) <- "SCT"
for (i in genes){
  p1 <- SpatialFeaturePlot(object = fibro.new, 
                         features = i,
                         combine = FALSE) 
  fix.p1 <- scale_fill_gradientn(colors=myPalette(100),
                               limits = c(0,4),
                               breaks=c(0,4),
                               labels=c("Min","Max"),)
  p2 <- lapply(p1, function (x) x + fix.p1)

  pdf(file.path("./results/Aucell/all_genes/",filename = paste(i,"spatial.all.pdf",sep="")))
  print(CombinePlots(p2))
  dev.off()
}





