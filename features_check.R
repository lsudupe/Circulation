### SCRIPT:  Read genes and create geneSet CIRCULATION proyect

## 29.11.22 Laura Sudupe , git @lsudupe

###### Libraries
library(openxlsx)
library("Seurat")
library("GSEABase")
library("ggplot2")
library(UCell)
set.seed(123)

###### Read the gene sets
genes <- read.xlsx("./Circulation/data/Silvia_clean.xlsx", colNames = TRUE)
verde <- na.omit(genes$verde)
morado <- na.omit(genes$morado)

###### Read data
#Variables---------------------------------
DIR_ROOT <- file.path(getwd())
DIR_DATA <- file.path(DIR_ROOT, "./Circulation/objects/individual/figures/")

#Data--------------------------------------
samples <- dir(path = DIR_DATA)
#samples <- samples[! samples %in% c("V11B18-363_A1")]
#samples <- samples[samples %in% c("control","dpi3","dpi5_female","dpi5_male")]

# Add area info in each sample
lista <- c(samples)
names(lista) <- samples
prueba <-c()

for (i in lista[-1]){
  a <- lista[i]
  sam <- readRDS(paste0(DIR_DATA, i,"/", i,".rds"))
  ## add area data
  area_a <- read.csv(paste0(DIR_ROOT, "/Circulation/data/", i, "/Enrichment.csv"))
  area_a <- as.vector(area_a$Enrichment)
  sam@meta.data["area"] <- as.factor(area_a)
  ## add object to list
  prueba[[length(prueba) + 1]] <- sam
}

control <- readRDS(paste0(DIR_DATA, "control/control.rds"))
prueba[[length(prueba) + 1]] <- control
names(prueba) <- c("dpi3","dpi5_female", "dpi5_male", "control")

control <-prueba[[4]]
dpi3 <-prueba[[1]]
dpi5_female <-prueba[[2]]
dpi5_male <-prueba[[3]]

control@meta.data[["orig.ident"]] <- "control"
dpi3@meta.data[["orig.ident"]] <- "dpi3"
dpi5_female@meta.data[["orig.ident"]] <- "dpi5_female"
dpi5_male@meta.data[["orig.ident"]] <- "dpi5_male"

##Merge them
combined <- merge(control, y = c(dpi3,dpi5_female,dpi5_male), 
                  add.cell.ids = c("control", "dpi3", "dpi5_female", "dpi5_male"), project = "infart")

saveRDS(combined, paste0(DIR_DATA, "combined.rds"))
combined <- readRDS(paste0(DIR_DATA, "combined.rds"))


#### Add enrichment score per list objetc, segmentation plots, boxpplots

objects <- c(control, dpi3, dpi5_female, dpi5_male)
names(objects) <- c("control", "dpi3","dpi5_female","dpi5_male")
enrich <- c()

###############################################333
for (i in 1:length(objects)){
  ##matrix
  a <- objects[[i]]
  DefaultAssay(a) <- "SCT"
  
  ## Check intersect genes B and object
  v.genes <- list(intersect(a@assays[["SCT"]]@data@Dimnames[[1]], verde))
  m.genes <- list(intersect(a@assays[["SCT"]]@data@Dimnames[[1]], morado))
  
  ## Add UCellScore
  a <- AddModuleScore_UCell(a, features = v.genes, name = "_verde")
  a <- AddModuleScore_UCell(a, features = m.genes, name = "_morado")
  
  ## Calculate ratio
  verde <- as.vector(a$signature_1_verde)
  morado <- as.vector(a$signature_1_morado)
  ratio <- log10(morado/verde)
  is.na(ratio) <-sapply(ratio, is.infinite)
  ratio[is.na(ratio)] = 0
  a$ratio <- ratio
  
  ## Save XXXXXXXX
  saveRDS(a,file = paste0(DIR_DATA, names(objects[i]),"/enrich/", names(objects[i]),".enrich.rds"))
  enrich[[length(enrich) + 1]] <- a
  
}
names(enrich) <- c("control", "dpi3","dpi5_female","dpi5_male")



##########################################################
######spatial
heart <- readRDS("./Circulation/objects/integrated/integrated.rds")
DefaultAssay(heart) <- "SCT"

###### Check intersect genes B and object
v.genes <- intersect(heart@assays[["SCT"]]@data@Dimnames[[1]], verde)
v.list <- as.list(v.genes)
names(v.list) <- as.list(v.genes)
m.genes <- intersect(heart@assays[["SCT"]]@data@Dimnames[[1]], morado)
m.list <- as.list(m.genes)
names(m.list) <- as.list(m.genes)

##### Check individually


for (i in m.list){
  library(BuenColors)
  color <- jdb_palette("brewer_spectra")
  pdf(paste("./Circulation/results/Adrian.genes/silvia.december/individual/morado/",i,"_morado.pdf",sep=""))
  print(SpatialFeaturePlot(heart, features = i,  combine = TRUE, ncol = 2))
  dev.off()
}


##### Add moduleScore
DefaultAssay(heart) <- "Spatial"
v.genes <- list(intersect(heart@assays[["SCT"]]@data@Dimnames[[1]], verde))
heart_module <- AddModuleScore(heart, features = v.genes, name = c("verde_"))
m.genes <- list(intersect(heart@assays[["SCT"]]@data@Dimnames[[1]], morado))
heart_module <- AddModuleScore(heart_module, features = m.genes, name = c("morado_"))

##### Add UCellScore

heart_Uscore <- AddModuleScore_UCell(heart, features = v.genes, name = "_verde")
heart_Uscore <- AddModuleScore_UCell(heart_Uscore, features = m.genes, name = "_morado")

##### Calculate relation d1/d2
verde <- as.vector(heart_Uscore$signature_1_verde)
morado <- as.vector(heart_Uscore$signature_1_morado)
ratio <- log10(morado/verde)
is.na(ratio) <-sapply(ratio, is.infinite)
ratio[is.na(ratio)] = 0
heart_Uscore$ratio <- ratio

#####

pdf(paste("./Circulation/results/Adrian.genes/silvia.december/Ucell/geneset/ratio/ratio.pdf",sep=""))
print(SpatialFeaturePlot(heart_Uscore, features = c("ratio"),  combine = TRUE, ncol = 2))
dev.off()

##### blanco y rojo

bl <- colorRampPalette(c("navy","royalblue","lightskyblue"))(200)                      
re <- colorRampPalette(c("mistyrose", "red2","darkred"))(200)

l <- c(min(heart_Uscore@meta.data[["ratio"]]), max(heart_Uscore@meta.data[["ratio"]]))
p1 <- SpatialFeaturePlot(heart_Uscore, features = c("ratio"), combine = FALSE, ncol = 2)
fix.p1 <- scale_fill_gradientn(colours=c(bl,"white", re),
                               breaks=b,
                               labels=c("Min","Max"),
                               na.value = "grey98",
                               limits = l)
p2 <- lapply(p1, function (x) x + fix.p1)
  
pdf(paste("./Circulation/results/Adrian.genes/silvia.december/Ucell/geneset/ratio/ratio_bw.pdf",sep=""))
print(CombinePlots(p2))
dev.off()

##### normal
b <- c(min(heart_Uscore@meta.data[["ratio"]]), max(heart_Uscore@meta.data[["ratio"]]))
p1 <- SpatialFeaturePlot(heart_Uscore, features = c("ratio"),  combine = FALSE, ncol = 2)
fix.p1 <- scale_fill_gradientn(colours=color,
                               breaks=b,
                               labels=c("Min","Max"),
                               limits = b)
p2 <- lapply(p1, function (x) x + fix.p1)

pdf(paste("./Circulation/results/Adrian.genes/silvia.december/Ucell/geneset/ratio/ratio_.pdf",sep=""))
print(CombinePlots(p2))
dev.off()








