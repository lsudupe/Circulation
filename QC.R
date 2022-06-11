## SCRIPT: QC of the seurat spatial objects CIRCULATION project

## 29.03.22 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)

#Data--------------------------------------
control <- readRDS("./objects/initial/control_sp.rds")
dpi3 <- readRDS("./objects/initial/3dpi_sp.rds")
dpi5_female <- readRDS("./objects/initial/5dpi_h_sp.rds")
dpi5_male <- readRDS("./objects/initial/5dpi_m_sp.rds")

samples <- c(control, dpi3, dpi5_female, dpi5_male)
names(samples) <- c("control", "dpi3", "dpi5_female", "dpi5_male")

for (i in 1:length(samples)){
  a <- samples[[i]]
  #Add number of genes per UMI for each cell to metadata
  a$log10GenesPerUMI <- log10(a$nFeature_Spatial) / log10(a$nCount_Spatial)
  # Compute percent mito ratio
  a <- PercentageFeatureSet(a, "^mt-", col.name = "percent_mito")
  ##add sample column
  a@meta.data$sample <- names(samples[i])
  ##save object
  saveRDS(a,file = paste0("./objects/initial/",names(samples[i]),"_sp.rds"))
}

control <- readRDS("./objects/initial/control_sp.rds")
dpi3 <- readRDS("./objects/initial/dpi3_sp.rds")
dpi5_female <- readRDS("./objects/initial/dpi5_female_sp.rds")
dpi5_male <- readRDS("./objects/initial/dpi5_male_sp.rds")

control@meta.data[["sample"]] <- "control"
dpi3@meta.data[["sample"]] <- "dpi3"
dpi5_female@meta.data[["sample"]] <- "dpi5_female"
dpi5_male@meta.data[["sample"]] <- "dpi5_male"

##Merge them
combined <- merge(control, y = c(dpi3, dpi5_female, dpi5_male ), 
                     add.cell.ids = c("control","dpi3","dpi5_female", "dpi5_male"), project = "Circulation")

##Visualization
# Visualize the number of spots counts per sample
pdf(file.path("./results/QC",filename = "Number of spot per sample.pdf"))
combined@meta.data%>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NSpots")
dev.off()

# Visualize the distribution of genes detected per spot via boxplot
pdf(file.path("./results/QC",filename = "genes detected per spot boxplot.pdf"))
combined@meta.data %>% 
  ggplot(aes(x=sample, y=log10(nFeature_Spatial), fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Ngenes vs Npots")
dev.off()

# Visualize the distribution of genes detected per spot via histogram
pdf(file.path("./results/QC",filename = "genes detected per spot histogram.pdf"))
combined@meta.data %>% 
  ggplot(aes(color=sample, x=nFeature_Spatial, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)
dev.off()

# Visualize the distribution of mitochondrial gene expression detected per spot
pdf(file.path("./results/QC",filename = "mito percentage per spot.pdf"))
combined@meta.data %>% 
  ggplot(aes(color=sample, x=percent_mito, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)
dev.off()

# Visualize the number UMIs/transcripts per cell
pdf(file.path("./results/QC",filename = "number UMIs sati transcripts per cell.pdf"))
combined@meta.data %>% 
  ggplot(aes(color=sample, x=nCount_Spatial, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)
dev.off()


# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI (novelty score)
pdf(file.path("./results/QC",filename = "genes detected per UMI novelty score.pdf"))
combined@meta.data %>% 
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)
dev.off()

##########SPATIAL plots

feature.list <- c("nCount_Spatial", "nFeature_Spatial","percent_mito", "log10GenesPerUMI")

for (i in feature.list){
  p1 <- SpatialFeaturePlot(combined, features = i, combine = FALSE)
  fix.p1 <- scale_fill_continuous(type = "viridis")
  p2 <- lapply(p1, function (x) x + fix.p1)
  
  pdf(paste("./results/QC/",i,"spatial_no_range.pdf",sep=""))
  print(CombinePlots(p2))
  dev.off()
}

##specific plots
p1 <- SpatialFeaturePlot(combined, features = "nCount_Spatial",alpha = 0.5,combine = FALSE)
fix.p1 <- scale_fill_continuous(limits = c(0,25000), 
                                breaks = c(0,25000),
                                type = "viridis")
p2 <- lapply(p1, function (x) x + fix.p1)

pdf(paste("./results/QC/nCount_Spatial.pdf",sep=""))
CombinePlots(p2)
dev.off()

p1 <- SpatialFeaturePlot(combined, features = "nFeature_Spatial",alpha = 0.5,combine = FALSE)
fix.p1 <- scale_fill_continuous(limits = c(0,6000), 
                                breaks = c(0,6000),
                                type = "viridis")
p2 <- lapply(p1, function (x) x + fix.p1)

pdf(paste("./results/QC/nFeature_Spatial.pdf",sep=""))
CombinePlots(p2)
dev.off()

p1 <- SpatialFeaturePlot(combined, features = "percent_mito",alpha = 0.5,combine = FALSE)
fix.p1 <- scale_fill_continuous(limits = c(0,100), 
                                breaks = c(0,50,100),
                                type = "viridis")
p2 <- lapply(p1, function (x) x + fix.p1)

pdf(paste("./results/QC/percent_mito.pdf",sep=""))
CombinePlots(p2)
dev.off()


#####save combined
saveRDS(combined,"./objects/initial/combined.rds")

combined.meta <- combined@meta.data

#####Violin plots
a <- ggplot(combined.meta, aes(x=sample, y=nCount_Spatial, color=sample)) +
  geom_violin(trim=FALSE) +
  labs(title="nCount_Spatial circulation",x="sample", y = "nCount_Spatial")

b <- ggplot(combined.meta, aes(x=sample, y=nFeature_Spatial, color=sample)) +
  geom_violin(trim=FALSE) +
  labs(title="nFeature_Spatial circulation",x="sample", y = "nFeature_Spatial")

c <- ggplot(combined.meta, aes(x=sample, y=percent_mito, color=sample)) +
  geom_violin(trim=FALSE) +
  labs(title="percent_mito circulation",x="sample", y = "percent_mito")

pdf(paste("./results/QC/features_violin.pdf",sep=""))
cowplot::plot_grid(a,NULL,b,NULL,c, nrow = 5, ncol = 1, rel_heights = c(4,0.5,3,0.5,3))
dev.off()


###########count/feature correlation

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
pdf(paste("./results/QC/count_features_cor.pdf",sep=""))
combined.meta %>% 
  ggplot(aes(x=nCount_Spatial, y=nFeature_Spatial, color=percent_mito)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)
dev.off()

FeatureScatter(control, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial")
FeatureScatter(dpi3, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial")
FeatureScatter(dpi5_female, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial")
FeatureScatter(dpi5_male, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial")


#####Check GFP#############################################################
library("GiNA")
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

p1 <- SpatialFeaturePlot(object = combined, 
                         features = c("GFP"),
                         combine = FALSE) 
fix.p1 <- scale_fill_gradientn(colors=myPalette(100),
                               limits = c(0,300),
                               breaks=c(0,300),
                               labels=c("Min","Max"),)
p2 <- lapply(p1, function (x) x + fix.p1)

pdf(file.path("./results/gfp/",filename = "combined_gfp.pdf"))
print(CombinePlots(p2))
dev.off()

#Take out GFP for the integration#########################
samples <- c(control, dpi3, dpi5_female, dpi5_male)
names(samples) <- c("control", "dpi3", "dpi5_female", "dpi5_male")

for (i in 1:length(samples)){
  a <- samples[[i]]
  counts <- as.data.frame(a@assays[["Spatial"]]@counts)
  counts <- counts[row.names(counts) != "GFP", , drop = FALSE]
  counts <- as.matrix(counts)
  a@assays[["Spatial"]]@counts <- counts
  
  data <- as.data.frame(a@assays[["Spatial"]]@data)
  data <- data[row.names(data) != "GFP", , drop = FALSE]
  data <- as.matrix(data)
  a@assays[["Spatial"]]@data <- data
  
  ##save object
  saveRDS(a,file = paste0("./objects/initial/",names(samples[i]),"_noGFP.rds"))
}

control <- readRDS("./objects/initial/control_noGFP.rds")
dpi3 <- readRDS("./objects/initial/dpi3_noGFP.rds")
dpi5_female <- readRDS("./objects/initial/dpi5_female_noGFP.rds")
dpi5_male <- readRDS("./objects/initial/dpi5_male_noGFP.rds")

##Merge them
combined <- merge(control, y = c(dpi3, dpi5_female, dpi5_male ), 
                  add.cell.ids = c("control","dpi3","dpi5_female", "dpi5_male"), project = "Circulation")

#####save combined
saveRDS(combined,"./objects/initial/combined_noGFP.rds")
