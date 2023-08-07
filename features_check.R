### SCRIPT:  Read genes and create geneSet CIRCULATION proyect

## 29.11.22 Laura Sudupe , git @lsudupe

###### Libraries
library(openxlsx)
library("Seurat")
library("GSEABase")
library("ggplot2")
library(UCell)
library(RColorBrewer)
library(dplyr)
set.seed(123)

###### Read the gene sets
genes <- read.xlsx("./Circulation/data/Silvia_clean.xlsx", colNames = TRUE)
verde <- na.omit(genes$verde)
morado <- na.omit(genes$morado)

fb.genes <- read.csv("./Circulation/results/DE/top200_tokio_FB.csv")
fb.g <- na.omit(fb.genes$gene)

nine <- read.csv("./Circulation/results/genes/top100nine.csv")
nine.g <- na.omit(nine$x)

genes.velocity <- read.csv("./Circulation/data/genes_velocity.csv", header = TRUE, sep = ",")
B <- read.table("./Circulation/data/clustB_signature.txt", header = FALSE)
B <- as.vector(B$V1)
DA <- as.vector(genes.velocity$Dynamics_1)
DB <- as.vector(genes.velocity$Dynamics_2)
DB[DB == ""] <- NA 
DB <- DB[!is.na(DB)]

###### Read data
#Variables---------------------------------
DIR_ROOT <- file.path(getwd())
DIR_DATA <- file.path(DIR_ROOT, "./Circulation/objects/individual/figures/")
DIR_RES <-  file.path(DIR_ROOT, "./Circulation/results/Adrian.genes/silvia.december/segm/")

#Data--------------------------------------
samples <- dir(path = DIR_DATA)
#samples <- samples[! samples %in% c("V11B18-363_A1")]
#samples <- samples[samples %in% c("control","dpi3","dpi5_female","dpi5_male")]

# Add area info in each sample
lista <- c(samples)
names(lista) <- samples
lista <- lista[-1]
prueba <-c()

for (i in lista){
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
  fb.genes <- list(intersect(a@assays[["SCT"]]@data@Dimnames[[1]], fb.g))
  nine.genes <- list(intersect(a@assays[["SCT"]]@data@Dimnames[[1]], nine.g))
  b.genes <- list(intersect(a@assays[["SCT"]]@data@Dimnames[[1]], B))
  da.genes <- list(intersect(a@assays[["SCT"]]@data@Dimnames[[1]], DA))
  db.genes <- list(intersect(a@assays[["SCT"]]@data@Dimnames[[1]], DB))
  
  
  ## Add UCellScore
  a <- AddModuleScore_UCell(a, features = v.genes, name = "_verde")
  a <- AddModuleScore_UCell(a, features = m.genes, name = "_morado")
  a <- AddModuleScore_UCell(a, features = m.genes, name = "_fb")
  a <- AddModuleScore_UCell(a, features = nine.genes, name = "_nine")
  a <- AddModuleScore_UCell(a, features = b.genes, name = "_b")
  a <- AddModuleScore_UCell(a, features = da.genes, name = "_da")
  a <- AddModuleScore_UCell(a, features = db.genes, name = "_db")
  
  
  ## Calculate ratio verde/morado
  v <- as.vector(a$signature_1_verde)
  m <- as.vector(a$signature_1_morado)
  ratio <- log10(m/v)
  is.na(ratio) <-sapply(ratio, is.infinite)
  ratio[is.na(ratio)] = 0
  a$ratio <- ratio
  
  ## Calculate ratio da/db
  da <- as.vector(a$signature_1_da)
  db <- as.vector(a$signature_1_db)
  ratio.dadb <- log10(da/db)
  is.na(ratio.dadb) <-sapply(ratio.dadb, is.infinite)
  ratio.dadb[is.na(ratio.dadb)] = 0
  a$ratio.dadb <- ratio.dadb
  
  ## Regress out FB values
  meta <- a@meta.data
  lm <- lm(meta$signature_1_verde ~ meta$signature_1_fb, data =meta)
  residuals <- lm$residuals
  a@meta.data[["residuals_verde"]] <- residuals
  
  lm <- lm(meta$signature_1_morado ~ meta$signature_1_fb, data =meta)
  residuals <- lm$residuals
  a@meta.data[["residuals_morado"]] <- residuals
  
  ## Save
  saveRDS(a,file = paste0(DIR_DATA, names(objects[i]),"/enrich/", names(objects[i]),".enrich_2.rds"))
  enrich[[length(enrich) + 1]] <- a
  
}
names(enrich) <- c("control", "dpi3","dpi5_female","dpi5_male")

#####################seg.PLOTS#####################
## Create a custom color scale
nombres <- c("IZ","BZ_1","BZ_2","BZ_3","BZ_2_3","BZ_4","RZ")
colors <- c("coral4", "skyblue", "red", "darkgreen", "yellow","khaki2", "darkolivegreen", "gray89")
names(colors)  <- nombres

## Objects
new <- enrich[-1]

for (i in 1:length(new)){
  ## Object
  a <- new[[i]]
  
  ## Order area as factor
  a@meta.data[["area"]] <- factor(a@meta.data[["area"]],
                                          levels = c("IZ","BZ_1","BZ_2","BZ_3","BZ_2_3","BZ_4","RZ"),ordered = TRUE)
 
  ## boxplot
  # verde
  pdf(file.path(paste0(DIR_RES, names(new[i]),"/", names(new[i]), "boxplot_verde.pdf")))
  print(a@meta.data%>% 
    ggplot(aes(x=signature_1_verde, y= area, fill=area)) + 
    geom_boxplot() +  
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    #xlim(-0.5, 0.5) +
    theme(plot.title = element_text(hjust=0.5, face="bold")))
  dev.off() 
  pdf(file.path(paste0(DIR_RES, names(new[i]),"/", names(new[i]), "boxplot_verde_residuals.pdf")))
  print(a@meta.data%>% 
          ggplot(aes(x=residuals_verde, y= area, fill=area)) + 
          geom_boxplot() +  
          theme_classic() +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
          #xlim(-0.5, 0.5) +
          theme(plot.title = element_text(hjust=0.5, face="bold")))
  dev.off() 
  
  # morado
  pdf(file.path(paste0(DIR_RES, names(new[i]),"/", names(new[i]), "boxplot_morado.pdf")))
  print(a@meta.data%>% 
    ggplot(aes(x=signature_1_morado, y= area, fill=area)) + 
    geom_boxplot() +  
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    #xlim(-0.5, 0.5) +
    theme(plot.title = element_text(hjust=0.5, face="bold")))
  dev.off()
  pdf(file.path(paste0(DIR_RES, names(new[i]),"/", names(new[i]), "boxplot_morado_residuals.pdf")))
  print(a@meta.data%>% 
          ggplot(aes(x=residuals_morado, y= area, fill=area)) + 
          geom_boxplot() +  
          theme_classic() +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
          #xlim(-0.5, 0.5) +
          theme(plot.title = element_text(hjust=0.5, face="bold")))
  dev.off() 
  
  pdf(file.path(paste0(DIR_RES, names(new[i]),"/", names(new[i]), "boxplot_ratio.pdf")))
  print(a@meta.data%>% 
    ggplot(aes(x=ratio, y= area, fill=area)) + 
    geom_boxplot() +  
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    #xlim(-0.5, 0.5) +
    theme(plot.title = element_text(hjust=0.5, face="bold")))
  dev.off()
  
  ## spatial plots
  library(BuenColors)
  color <- jdb_palette("brewer_spectra")
  # verde
  b <- c(min(a@meta.data[["signature_1_verde"]]), max(a@meta.data[["signature_1_verde"]]))
  p1 <- SpatialFeaturePlot(a, features = c("signature_1_verde"),  combine = FALSE, ncol = 1)
  fix.p1 <- scale_fill_gradientn(colours=color,
                                 breaks=b,
                                 labels=c("Min","Max"),
                                 limits = b)
  p2 <- lapply(p1, function (x) x + fix.p1)
  
  pdf(paste0("./Circulation/results/Adrian.genes/silvia.december/Ucell/geneset/verde/",names(new[i]), "/", paste(names(new[i]), "_enrich_verde.pdf",sep="")))
  print(CombinePlots(p2))
  dev.off()
  
  b <- c(min(a@meta.data[["residuals_verde"]]), max(a@meta.data[["residuals_verde"]]))
  p1 <- SpatialFeaturePlot(a, features = c("residuals_verde"),  combine = FALSE, ncol = 1)
  fix.p1 <- scale_fill_gradientn(colours=color,
                                 breaks=b,
                                 labels=c("Min","Max"),
                                 limits = b)
  p2 <- lapply(p1, function (x) x + fix.p1)
  
  pdf(paste0("./Circulation/results/Adrian.genes/silvia.december/Ucell/geneset/verde/",names(new[i]), "/", paste(names(new[i]), "enrich_verde_residuals.pdf",sep="")))
  print(CombinePlots(p2))
  dev.off()

  # morado
  b <- c(min(a@meta.data[["signature_1_morado"]]), max(a@meta.data[["signature_1_morado"]]))
  p1 <- SpatialFeaturePlot(a, features = c("signature_1_morado"),  combine = FALSE, ncol = 1)
  fix.p1 <- scale_fill_gradientn(colours=color,
                                 breaks=b,
                                 labels=c("Min","Max"),
                                 limits = b)
  p2 <- lapply(p1, function (x) x + fix.p1)
  
  pdf(paste0("./Circulation/results/Adrian.genes/silvia.december/Ucell/geneset/morado/",names(new[i]), "/", paste(names(new[i]), "enrich_morado.pdf",sep="")))
  print(CombinePlots(p2))
  dev.off()
  
  b <- c(min(a@meta.data[["residuals_morado"]]), max(a@meta.data[["residuals_morado"]]))
  p1 <- SpatialFeaturePlot(a, features = c("residuals_morado"),  combine = FALSE, ncol = 1)
  fix.p1 <- scale_fill_gradientn(colours=color,
                                 breaks=b,
                                 labels=c("Min","Max"),
                                 limits = b)
  p2 <- lapply(p1, function (x) x + fix.p1)
  
  pdf(paste0("./Circulation/results/Adrian.genes/silvia.december/Ucell/geneset/morado/",names(new[i]), "/", paste(names(new[i]), "enrich_morado_residuals.pdf",sep="")))
  print(CombinePlots(p2))
  dev.off()
  
  # ratio
  b <- c(min(a@meta.data[["ratio"]]), max(a@meta.data[["ratio"]]))
  p1 <- SpatialFeaturePlot(a, features = c("ratio"),  combine = FALSE, ncol = 1)
  fix.p1 <- scale_fill_gradientn(colours=color,
                                 breaks=b,
                                 labels=c("Min","Max"),
                                 limits = b)
  p2 <- lapply(p1, function (x) x + fix.p1)
  
  pdf(paste0("./Circulation/results/Adrian.genes/silvia.december/Ucell/geneset/ratio/",names(new[i]), "/", paste(names(new[i]), "ratio.pdf",sep="")))
  print(CombinePlots(p2))
  dev.off()
}


###################################################
##### Check individually


for (i in m.list){
  library(BuenColors)
  color <- jdb_palette("brewer_spectra")
  pdf(paste("./Circulation/results/Adrian.genes/silvia.december/individual/morado/",i,"_morado.pdf",sep=""))
  print(SpatialFeaturePlot(heart, features = i,  combine = TRUE, ncol = 2))
  dev.off()
}



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








