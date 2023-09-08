### SCRIPT:  Calcagno paper CM score CIRCULATION proyect

## 04.09.23 Laura Sudupe , git @lsudupe

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


# Define the gene lists for each category
RZ_genes <- c(
  "Mhrt", "Myh6", "Tnnt2", "Ttn", "Pde4d", "Fhl2", "Pln", "Rbm20", "Ctnna3",
  "Kcnd2", "Pde4dip", "Atp2a2", "Cacna1c", "Ryr2", "Grm1", "Fhod3", "Fgf13",
  "Corin", "Nexn", "Ppargc1a", "Coro6", "Magi2", "Vegfa", "Ivns1abp", "Lmo7",
  "Cpeb3", "Ccdc141", "Ldb3", "Erbb4", "Ppip5k2", "Atcayos", "Obscn", "Mlip",
  "Tnnc1", "Rbm24", "Kcnj3", "Txlnb", "Esrrg", "Slc4a3", "Mybpc3", "Tcap",
  "Myocd", "Tnni3", "Rcan2", "Cacnb2", "Sorbs1", "Pcdh7", "Chrm2", "Tnni3k",
  "Palld", "Trdn", "Actn2", "Neat1", "Lrrc2", "Dmd", "Fgf14", "Prune2",
  "Hs3st5", "Angpt1", "Lmod2", "Rnf207", "Synpo2", "Sox6", "Cmya5", "Tacc2",
  "Park2", "Kcnn2", "Thrb", "Acacb", "Akap6", "Gpcpd1", "Nebl", "Cdh2",
  "Mylk3", "Vldlr", "Ank2", "Dmpk", "Psme4", "Smtn", "Tecrl", "Ppp1r3a"
)

BZ1_genes <- c(
  "Shroom3", "Sorbs2", "Celf2", "Pdlim5", "Ppm1e", "Rbpms", "Slc8a1", "Ankrd1",
  "Zfp697", "Grk5", "Cdh2", "Efna5", "Col4a5", "Magi1", "Airn", "Ryr2", "Gm12295",
  "Kcnq1", "Lpp", "Dmd", "Scn5a", "Lrrfip1", "Chrm2", "Ppp1r14c", "Pam", "Tox3",
  "Dysf", "Tlr4", "B4galt1", "Mical2", "Ror1", "Lhfp", "Clic5", "Ptprk", "Ehbp1",
  "Pfkp", "Prkg1", "Pkia", "Gm13481", "Akap6", "Fam198b", "Uck2", "Actn4", "Ank3",
  "Cacna1c", "Msrb3", "Cap2", "Wls", "Tln2", "Cacnb2", "Adamtsl5", "Nppb", "Nrp1",
  "Pde4b", "Phldb2", "Gfod1", "Zfp608", "Hip1", "Etv6", "Ophn1", "Utrn", "Pcdh7",
  "Pdlim3", "Zfp568", "Gk", "Dtna", "Arhgap23", "Pip5k1b", "Man1a", "Hecw2", "Srl",
  "Phactr2", "Pdgfc", "Dpysl3", "Sgcd", "Msn", "Meis1", "Kcng2", "Xpr1", "Clybl",
  "Col4a1", "Ankrd10", "Col4a2", "Picalm", "Tom1", "Hivep3", "Plce1", "Lbh", "Arhgap6",
  "Parm1", "Pcgf5", "03-Mar", "Strip2", "Ghr", "Pkp2", "Ctps", "Zbtb38", "Adam19",
  "Bre", "Cblb", "Wwc2", "Zfpm2", "Shb", "Kif5b", "Snap91", "Chd9", "Sgms1", "Stat3",
  "Osmr", "Snta1", "Strn", "Atp8a2", "Cdh13", "Prune2", "Rab30", "Exoc6b", "Mnat1",
  "2010111I01Rik", "Pkig", "Cp", "Iqsec2", "Hspg2", "Cpeb1", "Heg1", "Luzp1", "Mid1",
  "Snd1", "Bmpr1a", "Afap1l1", "Maml3", "Pls3", "Ppip5k1", "Tnni3k", "Rbms3", "Rnf150",
  "Gata4", "Pgm5", "Rcan2", "Jph2", "Tnks", "Tgfb2", "Slc9a7", "Large", "Cstf3", "Syk",
  "Slc25a12", "Lnx2", "Man1c1", "Phf21a", "Gphn", "Gipc2", "Lnpep", "Bcl2", "Gm20732",
  "Masp1", "Hk1", "Celf1", "Qk", "Nuak1", "Rasal2", "Mpp7", "Mettl9", "Myocd", "Plcb1",
  "Osbpl9", "Mbd5", "St5", "Magi3", "Mical3", "Itga7", "Exoc4", "Atxn7l1", "Bcl2l1",
  "Gyg", "Trdn", "Myo1d", "Nedd4l", "Srgap2", "Rbfox2", "Ncam1", "Samd4", "Plin2", "Rlf",
  "Memo1", "Adamts9", "Ctnna1", "Sertad2", "Ttll11", "Cog5", "Mllt3", "Faf1", "Mast2",
  "Zfr", "Ssh2", "Actn2", "Trim55", "Nedd9", "Fam49b", "Rsu1", "Nlk", "Rarb", "Prkca",
  "Trabd2b", "Smpx", "Kat2b", "Rbm20", "Macrod2", "Ambra1", "Fam189a2", "Rbms1", "Lclat1",
  "Itga9", "Tbc1d5", "Mypn", "Lrba", "Pdss1", "Wdr7", "Immp2l", "Rps6ka3", "2210408I21Rik",
  "Sorbs1", "Prdm16", "Elp4", "Ammecr1", "Lnx1", "Spata5", "Lyrm4", "Bmpr2", "Mmd", "Grb10",
  "Tbx20", "Ptpn2", "Diaph2", "Ppp3cc", "Plpp1", "Parva", "Dzip3", "Ldb3", "Ube2g1", "Dym",
  "Cbfa2t2", "Rab3ip", "Usp32", "Nppa", "AW554918", "Pja2", "Sh3rf2", "Enox2", "Ptpn9",
  "Park2", "Aff2", "Fhit", "Prkdc", "Mlip", "Pcsk5", "Ankrd28", "Pdzrn3", "Lncpint", "B3galt1",
  "Lrch1", "Arhgap44", "Arid5b", "Vti1a", "Rusc2", "Limch1", "Reln", "Thsd4", "Myh7", "Elmo1",
  "Eda", "Kctd8"
)

BZ2_genes <- c(
  "Xirp2", "Enah", "Sorbs2", "Flnc", "Cacna2d1", "Rras2", "Dmd", "Palld", "Ank3",
  "Nrap", "Usp53", "Erc2", "Nrxn3", "Uck2", "Tacc2", "Large", "Ankrd1", "Tlr4",
  "Alpk2", "Zfp697", "Rtn4", "Airn", "Myo18b", "Gbe1", "1110002E22Rik", "Nav2",
  "Ppm1e", "Akap6", "Inpp4b", "Actn2", "Smpx", "Serpine1", "Pdlim5", "Pde4dip",
  "Lmcd1", "Ttn", "Sox6", "Sorbs1", "Zak", "Camk2d", "Adk", "Atxn10", "Pkp2", "Trdn",
  "Dtna", "Fgf1", "Samd4", "Fhl1", "Gm13481", "Prune2", "Ctnna3", "2810474O19Rik",
  "Ldb3", "Myocd", "Sacs", "Pcdh7", "Fgf13", "Limch1", "Tgfb2", "Fhod3", "Slc8a1",
  "Cdh2", "Efna5", "Rcan2", "Clic5", "Usp28", "Ube2e2", "Lrrfip1", "Obscn", "Mypn",
  "Ankrd10", "Vcl"
)

#####Mapping of CM gene-set scorers to visium clusters
# Iterate over each gene list and compute the score
gene_lists <- list(RZ_genes, BZ1_genes, BZ2_genes)
names(gene_lists) <- c("RZ_genes", "BZ1_genes", "BZ2_genes")

for (gene_set_name in names(gene_lists)) {
  genes_in_spatial <- gene_lists[[gene_set_name]][gene_lists[[gene_set_name]] %in% rownames(combined)]
  gene_counts <- GetAssayData(combined, assay = "Spatial", slot = "counts")[genes_in_spatial,]
  gene_score <- Matrix::colSums(gene_counts) / combined@meta.data$nCount_Spatial * 10000
  combined <- AddMetaData(combined, gene_score, col.name = gene_set_name)
}

# Verify the addition of the CM scores to the metadata
head(combined@meta.data)

#####classification
# Find all markers in the combined object
combined <-PrepSCTFindMarkers(combined)
zone_markers <- FindAllMarkers(combined, only.pos = TRUE, logfc.threshold = 0.25)
# Get the top 5 markers for each cluster based on average log2 fold change
top_zone_markers <- zone_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
# Plot the heatmap for the top markers with specific cluster colors
pdf("./Circulation/results/cm_score/one.pdf")
print(DoHeatmap(combined, features = c(top_zone_markers$gene), 
          group.colors = c("#7570B3", "#BA3859", "#FF0000", "darkgrey", "darkgrey")) + 
  scale_fill_gradientn(colors = c("grey", "white", brewer.pal(8, "Dark2")[3])))
dev.off()
# Plot the heatmap for the top markers with a broader color palette
pdf("./Circulation/results/cm_score/two.pdf")
print(DoHeatmap(combined, features = c(top_zone_markers$gene), 
          group.colors = c(brewer.pal(8, "Set2")[1:8], brewer.pal(8, "Accent")[1:8])) + 
  scale_fill_gradientn(colors = c("grey", "white", brewer.pal(8, "Dark2")[3])))
dev.off()

#####CM Score
# For the CM score, you want to take into account all three subclasses: RZ_genes, BZ1_genes, and BZ2_genes.
# Therefore, we'll combine the genes from all these subclasses to generate the CM score.
subclass_genes <- c(RZ_genes, BZ1_genes, BZ2_genes) 

# Extract the genes specific to the CM class that are present in your dataset
genes_in_spatial_and_cluster_class <- subclass_genes[subclass_genes %in% rownames(combined)]
gene_counts_class <- GetAssayData(combined, assay = "Spatial", slot = "counts")[genes_in_spatial_and_cluster_class, ]

# Compute the CM score
class_threshold <- 5  # Adjust this threshold based on your data, if necessary
class_score <- Matrix::colSums(gene_counts_class) / combined@meta.data$nCount_Spatial * 10000

# Add the computed CM score as metadata
combined <- AddMetaData(combined, class_score, col.name = "CM_Score")

# Create a new assay with the previously computed subclass scores 
# Adjust column indexes based on your metadata to include the RZ_genes, BZ1_genes, and BZ2_genes scores
# Create a new assay with the subclass scores
gene_scores <- CreateAssayObject(t(combined@meta.data[,c("RZ_genes", "BZ1_genes", "BZ2_genes", "CM_Score")]))
combined[["genescores"]] <- gene_scores


#####CM or IZ?
# Set the identity class for the object
Idents(combined) <- "sample"
# Subset by specific samples, for example, "control" and "dpi3"
combined_control <- subset(combined, idents = "control")
combined_dpi3 <- subset(combined, idents = "dpi3")
combined_dpi5f <- subset(combined, idents = "dpi5_female")
combined_dpi5m <- subset(combined, idents = "dpi5_male")

# Create a histogram for the CM score
hist(combined$CM_Score, breaks = 30)

# Create a dataframe for ggplot
data <- data.frame(type = c(rep("control", length(combined_control$CM_Score)), 
                            rep("dpi3", length(combined_dpi3$CM_Score)), 
                            rep("dpi5_female", length(combined_dpi5f$CM_Score)), 
                            rep("dpi5_male", length(combined_dpi5m$CM_Score))),
                   value = c(combined_control$CM_Score, 
                             combined_dpi3$CM_Score, 
                             combined_dpi5f$CM_Score, 
                             combined_dpi5m$CM_Score))

# Generate the histogram using ggplo
p <- data %>%
  ggplot(aes(x = value, fill = type)) +
  geom_histogram(color = "#e9ecef", alpha = 0.6, position = 'identity') +
  scale_fill_manual(values = c("control" = "#FF3333", 
                               "dpi3" = "#00CCCC", 
                               "dpi5_female" = "#FF9999",  # Choose an appropriate color
                               "dpi5_male" = "#99CCFF"),   # Choose an appropriate color
                    name = "Data Type",  # Legend title
                    breaks = c("control", "dpi3", "dpi5_female", "dpi5_male"),
                    labels = c("Control", "DPI 3", "DPI 5 Female", "DPI 5 Male")) +
  labs(fill = "") + xlim(0, 1500) + theme_light() +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank())


# Display the plot
pdf("./Circulation/results/cm_score/histo_cm_score.pdf")
print(p)
dev.off()

#vilinplot
# Generate a violin plot for the CM score
p <- VlnPlot(combined, 
        features = "CM_Score", 
        pt.size = -1, 
        group.by = "seurat_clusters", 
        cols = c(brewer.pal(8, "Set2")[1:8], brewer.pal(8, "Accent")[1:8])) +
  theme_classic() +
  theme(legend.position = "none", title = element_text(size = 0))
pdf("./Circulation/results/cm_score/violin_cm_score.pdf")
print(p)
dev.off()

# Plot the CM_Score on the UMAP
DefaultAssay(combined) <- "Spatial"
combined <- FindVariableFeatures(combined)
combined <- ScaleData(combined)
combined <- RunPCA(combined)
combined <- RunUMAP(combined, dims = 1:25)
p1 <- DimPlot(combined, group.by = "seurat_clusters", label = TRUE) + NoLegend() + theme_minimal()
p2 <- FeaturePlot(combined, "CM_Score") + ggtitle("CM_score")

pdf("./Circulation/results/cm_score/dimplot_cm_score.pdf")
print(p1 | p2)
dev.off()

# ROC analysis
roc_results <- FindMarkers(combined, ident.1 = c(5,6), ident.2 = c(4), test.use = "roc")
print(roc_results)

#in my case, score <350 no RZ
#Idents(combined) <- "high_resolution"
#object_integrated_CMs <- subset(combined, idents = c(1:9,13)) #subseting CMs in new object

#####RZ or BZ?
#dataframe
data_bz <- data.frame(type = c(rep("control", length(combined_control$BZ1_genes)), 
                            rep("dpi3", length(combined_dpi3$BZ1_genes)), 
                            rep("dpi5_female", length(combined_dpi5f$BZ1_genes)), 
                            rep("dpi5_male", length(combined_dpi5m$BZ1_genes))),
                   value = c(combined_control$BZ1_genes, 
                             combined_dpi3$BZ1_genes, 
                             combined_dpi5f$BZ1_genes, 
                             combined_dpi5m$BZ1_genes))

# Generate the histogram using ggplo
p <- data_bz %>%
  ggplot(aes(x = value, fill = type)) +
  geom_histogram(color = "#e9ecef", alpha = 0.6, position = 'identity') +
  scale_fill_manual(values = c("control" = "#FF3333", 
                               "dpi3" = "#00CCCC", 
                               "dpi5_female" = "#FF9999",  # Choose an appropriate color
                               "dpi5_male" = "#99CCFF"),   # Choose an appropriate color
                    name = "Data Type",  # Legend title
                    breaks = c("control", "dpi3", "dpi5_female", "dpi5_male"),
                    labels = c("Control", "DPI 3", "DPI 5 Female", "DPI 5 Male")) +
  labs(fill = "") + xlim(0, 1500) + theme_light() +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank())


# Display the plot
pdf("./Circulation/results/cm_score/histo_bz_score.pdf")
print(p)
dev.off()

#vilinplot
# Generate a violin plot for the BZ score
p <- VlnPlot(combined, 
             features = "BZ1_genes", 
             pt.size = -1, 
             group.by = "seurat_clusters", 
             cols = c(brewer.pal(8, "Set2")[1:8], brewer.pal(8, "Accent")[1:8])) +
  theme_classic() +
  theme(legend.position = "none", title = element_text(size = 0))
pdf("./Circulation/results/cm_score/violin_bz_score.pdf")
print(p)
dev.off()

p1 <- DimPlot(combined, group.by = "seurat_clusters", label = TRUE) + NoLegend() + theme_minimal()
p2 <- FeaturePlot(combined, "BZ1_genes", cols = c("lightgrey", "darkred")) + ggtitle("BZ_genes")

pdf("./Circulation/results/cm_score/dimplot_bz_score.pdf")
print(p1 | p2)
dev.off()



#####BZ1 or BZ2?
#dataframe
data_bz2 <- data.frame(type = c(rep("control", length(combined_control$BZ2_genes)), 
                               rep("dpi3", length(combined_dpi3$BZ2_genes)), 
                               rep("dpi5_female", length(combined_dpi5f$BZ2_genes)), 
                               rep("dpi5_male", length(combined_dpi5m$BZ2_genes))),
                      value = c(combined_control$BZ2_genes, 
                                combined_dpi3$BZ2_genes, 
                                combined_dpi5f$BZ2_genes, 
                                combined_dpi5m$BZ2_genes))

# Generate the histogram using ggplo
p <- data_bz2 %>%
  ggplot(aes(x = value, fill = type)) +
  geom_histogram(color = "#e9ecef", alpha = 0.6, position = 'identity') +
  scale_fill_manual(values = c("control" = "#FF3333", 
                               "dpi3" = "#00CCCC", 
                               "dpi5_female" = "#FF9999",  # Choose an appropriate color
                               "dpi5_male" = "#99CCFF"),   # Choose an appropriate color
                    name = "Data Type",  # Legend title
                    breaks = c("control", "dpi3", "dpi5_female", "dpi5_male"),
                    labels = c("Control", "DPI 3", "DPI 5 Female", "DPI 5 Male")) +
  labs(fill = "") + xlim(0, 1500) + theme_light() +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank())


# Display the plot
pdf("./Circulation/results/cm_score/histo_bz2_score.pdf")
print(p)
dev.off()

#vilinplot
# Generate a violin plot for the BZ2 score
p <- VlnPlot(combined, 
             features = "BZ2_genes", 
             pt.size = -1, 
             group.by = "seurat_clusters", 
             cols = c(brewer.pal(8, "Set2")[1:8], brewer.pal(8, "Accent")[1:8])) +
  theme_classic() +
  theme(legend.position = "none", title = element_text(size = 0))
pdf("./Circulation/results/cm_score/violin_bz2_score.pdf")
print(p)
dev.off()

p1 <- DimPlot(combined, group.by = "seurat_clusters", label = TRUE) + NoLegend() + theme_minimal()
p2 <- FeaturePlot(combined, "BZ2_genes", cols = c("lightgrey", "darkred")) + ggtitle("BZ2_genes")

pdf("./Circulation/results/cm_score/dimplot_bz2_score.pdf")
print(p1 | p2)
dev.off()













