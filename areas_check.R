### SCRIPT:  Calcagno paper areas CIRCULATION proyect

## 20.08.23 Laura Sudupe , git @lsudupe

###### Libraries

library(UCell)
library(RColorBrewer)
library(dplyr)
library(Seurat)
library(ggplot2)


# spatial data
DIR_ROOT <- file.path(getwd())
DIR_DATA <- file.path(DIR_ROOT, "./Circulation/objects/individual/figures/")
combined <- readRDS(paste0(DIR_DATA, "combined.rds"))

# areas genes

rz_genes <- c("Ckm", "Ech1", "Cox7a1", "Cox8b", "Acadm", "Tcap", 
              "Hadhb", "Ckmt2", "Eno3", "Acaa2", "mt-Nd4l", "Idh2", 
              "Ldhb", "Hrc", "Acsl1", "Chchd10", "Hopx", "mt-Atp8", 
              "Fhl2", "Acadvl", "Etfb", "Eci1", "Mpc2", "Hadha", "Gm15543",
              "Etfa", "Hadh", "Ttn", "Pfkm", "Ptgds", "Mgst3", "Myoz2", 
              "Nipsnap2", "Pygm", "Dgat2", "Mpc1", "Acadl", "Etfdh",
              "Cpt1b", "Rpl3l", "Acat1", "Got1", "Fndc5", "Acss1", "Ndufv2", 
              "Phyh", "Ndufs1", "Pgam2", "Coq8a", "Uqcrc2")

bz1_genes <- c("Casq2", "Hbb-bt", "Cilp", "Myh7", "Hba-a2", "Hba-a1", "Nmrk2")

bz2_genes <- c("Nppa", "Sprr1a", "Postn", "Timp1", "Ccn2", "Col1a2", "Col1a1", "Fn1", "Rtn4", 
               "Col5a2", "Sparc", "Xirp2", "Col3a1", "Tnc", "Serpine1", "Cilp", "Bgn", "Col8a1", 
               "Mfap5", "Col15a1", "Fstl1", "Tgm2", "Hspb1", "Vim", "Serpinb1a", "Serpinh1", "Gpx3", 
               "Lox", "Cthrc1", "Col4a1", "Lgals1", "Ltbp2", "Fbn1", "S100a6", "Nppb", "S100a11", 
               "Tuba1a", "Actg1", "Col5a1", "Emp1", "Col4a2", "Anxa2", "Myl6", "Actn1", "Thbs4", 
               "Fbln2", "Tmsb10", "Thbs1", "Hba-a2", "Myh7")

iz_genes <- c("Spp1", "Lyz2", "Apoe", "Ftl1", "Col1a1", "Fn1", "Ctsb", "Col1a2", "Ctss", 
              "Lgals3", "Col3a1", "Tyrobp", "C1qb", "Tmsb10", "Sparc", "Pf4", "Fabp5", "Bgn", 
              "Eef1a1", "Mpeg1", "Gpnmb", "Actg1", "Sfrp2", "Timp1", "Grn", "Fcer1g", "Vim", 
              "Cfl1", "Cd74", "Hmox1", "C1qa", "Lgmn", "Serpinh1", "Ctsz", "S100a6", "Sh3bgrl3", 
              "Arpc1b", "C1qc", "Ctsl", "Mfap4", "Cd68", "Col5a2", "Wfdc17", "Postn", "Col5a1", 
              "S100a4", "Lgals1", "Ccl6", "Fbln2", "Eln")


## Add UCellScore
vector<- ScoreSignatures_UCell(combined@assays[["Spatial"]]@counts, features = list(rz_genes))
combined@meta.data[["signature_1_rz"]] <- as.vector(vector)

vector<- ScoreSignatures_UCell(combined@assays[["Spatial"]]@counts, features = list(bz1_genes))
combined@meta.data[["signature_1_bz1"]] <- as.vector(vector)

vector<- ScoreSignatures_UCell(combined@assays[["Spatial"]]@counts, features = list(bz2_genes))
combined@meta.data[["signature_1_bz2"]] <- as.vector(vector)

vector<- ScoreSignatures_UCell(combined@assays[["Spatial"]]@counts, features = list(iz_genes))
combined@meta.data[["signature_1_iz"]] <- as.vector(vector)



## plots
color <- brewer.pal(11,"Spectral")
color <- rev(color)

b <- c(min(combined@meta.data[["signature_1_rz"]]), max(combined@meta.data[["signature_1_rz"]]))
p1 <- SpatialFeaturePlot(combined, features = c("signature_1_rz"),  combine = FALSE, ncol = 1)
fix.p1 <- scale_fill_gradientn(colours=color,
                               breaks=b,
                               labels=c("Min","Max"),
                               limits = b)
p2 <- lapply(p1, function (x) x + fix.p1)

pdf(paste0("./Circulation/results/areas_check/rz.pdf",sep=""))
print(CombinePlots(p2))
dev.off()

b <- c(min(combined@meta.data[["signature_1_bz1"]]), max(combined@meta.data[["signature_1_rz"]]))
p1 <- SpatialFeaturePlot(combined, features = c("signature_1_bz1"),  combine = FALSE, ncol = 1)
fix.p1 <- scale_fill_gradientn(colours=color,
                               breaks=b,
                               labels=c("Min","Max"),
                               limits = b)
p2 <- lapply(p1, function (x) x + fix.p1)

pdf(paste0("./Circulation/results/areas_check/bz1.pdf",sep=""))
print(CombinePlots(p2))
dev.off()

b <- c(min(combined@meta.data[["signature_1_bz2"]]), max(combined@meta.data[["signature_1_rz"]]))
p1 <- SpatialFeaturePlot(combined, features = c("signature_1_bz2"),  combine = FALSE, ncol = 1)
fix.p1 <- scale_fill_gradientn(colours=color,
                               breaks=b,
                               labels=c("Min","Max"),
                               limits = b)
p2 <- lapply(p1, function (x) x + fix.p1)

pdf(paste0("./Circulation/results/areas_check/bz2.pdf",sep=""))
print(CombinePlots(p2))
dev.off()


b <- c(min(combined@meta.data[["signature_1_iz"]]), max(combined@meta.data[["signature_1_rz"]]))
p1 <- SpatialFeaturePlot(combined, features = c("signature_1_iz"),  combine = FALSE, ncol = 1)
fix.p1 <- scale_fill_gradientn(colours=color,
                               breaks=b,
                               labels=c("Min","Max"),
                               limits = b)
p2 <- lapply(p1, function (x) x + fix.p1)

pdf(paste0("./Circulation/results/areas_check/iz.pdf",sep=""))
print(CombinePlots(p2))
dev.off()

