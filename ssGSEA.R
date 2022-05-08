### SCRIPT: ssGSEA method enrichement analysis Circulation data in IBEX

## 08.04.22 Laura Sudupe , git @lsudupe

#https://ncborcherding.github.io/vignettes/escape_vignette.html


###### libraries
library("Seurat")
library("escape")
#install.packages("GSEABase")
library("GSEABase")
library("dittoSeq")
library("ggplot2")
library("dichromat")

##CHOOSE THE auxMaxRank
#genes.porcentage <-  210#porcentage depending number of genes in spatial

###### Read the gene sets
genes.velocity <- read.csv("./data/genes_velocity.csv", header = TRUE, sep = ",")
B <- read.table("./data/clustB_signature.txt", header = FALSE)
B <- as.vector(B$V1)
D1 <- as.vector(genes.velocity$Dynamics_1)
D2 <- as.vector(genes.velocity$Dynamics_2)
D2[D2 == ""] <- NA 
D2 <- D2[!is.na(D2)]

###### Create genesets
b.sig <- GeneSet(B, setName="geneSetB")
d1.sig <- GeneSet(D1, setName="geneSetD1")
d2.sig <- GeneSet(D2, setName="geneSetD2")
geneSets <- GeneSetCollection(d1.sig, d2.sig, b.sig)

###### Download  data
#a <- readRDS("../data/circulation_spatial/integrated.gfp.b.rds")
fibro <- readRDS("./objects/integrated/integrated.fb.rds")
no.fibro <- readRDS("./objects/integrated/integrated.nofb.rds")

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
  ##################ssGSEA#############3
  ###### Enrichment
  i.ES <- enrichIt(obj = matrix, 
                   gene.sets = geneSets, 
                   groups = 1000, cores = 4)
  ##save meta
  a <- AddMetaData(a, i.ES)
  d1 <- as.vector(a@meta.data[["geneSetD1"]])
  d2 <- as.vector(a@meta.data[["geneSetD2"]])
  d1.d2 <- log10(d1/d2)
  d1.d2[is.na(d1.d2)] = 0
  a@meta.data[["d1.d2"]] <- d1.d2
  ##save object
  saveRDS(a,file = paste0("./results/ssGSEA/all_genes/",names(objects[i]),".rds"))
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

fibro.new <- readRDS("./results/ssGSEA/all_genes/fibro.rds")
no.fibro.new <- readRDS("./results/ssGSEA/all_genes/no.fibro.rds")

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

dens <- c("geneSetD1","geneSetD2", "d1.d2","geneSetB")

##fibro
for (i in dens){
  pdf(file.path("./results/ssGSEA/all_genes/",filename = paste(i,"dens.fibro.pdf",sep="")))
  print(dittoRidgePlot(fibro.new, i, group.by = "sample"))
  dev.off()
}

##no fibro
for (i in dens){
  pdf(file.path("./results/ssGSEA/all_genes/",filename = paste(i,"dens.no.fibro.pdf",sep="")))
  print(dittoRidgePlot(no.fibro.new, i, group.by = "sample"))
  dev.off()
}


####Check best values to plot spatial STANDARIZED#############################3
meta <- fibro.new@meta.data
dpi3 <- meta[grepl("dpi3", meta[,6]),]
dpi3_b <- as.vector(dpi3$geneSetB)
dpi5_female <- meta[grepl("dpi5_female", meta[,6]),]
dpi5_male <- meta[grepl("dpi5_male", meta[,6]),]

library(robustHD)
library(data.table)
dpi3_b_z <- standardize(dpi3_b, centerFun = mean, scaleFun = sd)
hist(dpi3_b_z)

meta <- fibro.new@meta.data
meta.d1.d2 <- meta$d1.d2
meta.z.d1.d2 <- standardize(meta.d1.d2, centerFun = mean, scaleFun = sd)
hist(meta.z.d1.d2)
fibro.new@meta.data$d1.d2 <- meta.z.d1.d2

#Spatial plots
bl <- colorRampPalette(c("navy","royalblue","lightskyblue"))(200)                      
re <- colorRampPalette(c("mistyrose", "red2","darkred"))(200)
library("GiNA")
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

feature.list <- c("geneSetD1","geneSetD2","geneSetB")
b <- c(1000,4000)
feature.list <- c("d1.d2")
b <- c(-2,2)
feature.list <- c("relation_log.d1.d2")
b <- c(-0.1,0.2)

##fibro
for (i in feature.list){
  pdf(file.path("./results/ssGSEA/all_genes/",filename = paste(i,"dens.fibro.z.pdf",sep="")))
  print(dittoRidgePlot(fibro.new, i, group.by = "sample"))
  dev.off()
}

##################################

for (i in feature.list){
  p1 <- SpatialFeaturePlot(fibro.new, features = i, combine = FALSE)
  fix.p1 <- scale_fill_gradientn(colours=c(bl,"white", re), na.value = "grey98",limits = b,breaks=b, labels = c("min", "max"))
  #fix.p1 <- scale_fill_gradientn(colours=myPalette(100),breaks=b, labels = c("min", "max"),limits = b)
  p2 <- lapply(p1, function (x) x + fix.p1)
  
  pdf(file.path("./results/ssGSEA/all_genes/",filename = paste(i,"spatial.pdf",sep="")))
  print(CombinePlots(p2))
  dev.off()
}

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