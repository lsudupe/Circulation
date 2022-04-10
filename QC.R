## SCRIPT: QC of the seurat spatial objects CIRCULATION project

## 29.03.22 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)

#Data--------------------------------------
control <- readRDS("./objects/initial/control_sp.rds")
infarto_3dpi <- readRDS("./objects/initial/3dpi_sp.rds")
infarto_5dpi_h <- readRDS("./objects/initial/5dpi_h_sp.rds")
infarto_5dpi_m <- readRDS("./objects/initial/5dpi_m_sp.rds")

samples <- c(control, infarto_3dpi, infarto_5dpi_h, infarto_5dpi_m)
names(samples) <- c("control", "3dpi", "5dpi_h", "5dpi_m")

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
infarto_3dpi <- readRDS("./objects/initial/3dpi_sp.rds")
infarto_5dpi_h <- readRDS("./objects/initial/5dpi_h_sp.rds")
infarto_5dpi_m <- readRDS("./objects/initial/5dpi_m_sp.rds")

##Merge them
combined <- merge(control, y = c(infarto_3dpi, infarto_5dpi_h, infarto_5dpi_m ), 
                     add.cell.ids = c("control","3dpi","5dpi_h", "5dpi_m"), project = "Circulation")

##Visualization
# Visualize the number of spots counts per sample
png(file.path("./results/QC",filename = "Number of spot per sample.png"))
combined@meta.data%>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NSpots")
dev.off()

# Visualize the distribution of genes detected per spot via boxplot
png(file.path("./results/QC",filename = "genes detected per spot boxplot.png"))
combined@meta.data %>% 
  ggplot(aes(x=sample, y=log10(nFeature_Spatial), fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Ngenes vs Npots")
dev.off()

# Visualize the distribution of mitochondrial gene expression detected per spot
png(file.path("./results/QC",filename = "mito percentage per spot.png"))
combined@meta.data %>% 
  ggplot(aes(color=sample, x=percent_mito, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)
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
FeatureScatter(infarto_3dpi, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial")
FeatureScatter(infarto_5dpi_h, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial")
FeatureScatter(infarto_5dpi_m, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial")

####filter and save combined data

combined <- subset(combined, subset = nFeature_Spatial > 300 & nFeature_Spatial < 5500 & percent_mito < 50)
saveRDS(combined, "./objects/initial/combined_sp.rds")


control <- subset(control, subset = nFeature_Spatial > 300 & nFeature_Spatial < 5500 & percent_mito < 50)
infarto_3dpi <- subset(infarto_3dpi, subset = nFeature_Spatial > 300 & nFeature_Spatial < 5500 & percent_mito < 50)
infarto_5dpi_h <- subset(infarto_5dpi_h, subset = nFeature_Spatial > 300 & nFeature_Spatial < 5500 & percent_mito < 50)
infarto_5dpi_m <- subset(infarto_5dpi_m, subset = nFeature_Spatial > 300 & nFeature_Spatial < 5500 & percent_mito < 50)

saveRDS(control,"./objects/initial/control_sp.rds")
saveRDS(infarto_3dpi,"./objects/initial/3dpi_sp.rds")
saveRDS(infarto_5dpi_h,"./objects/initial/5dpi_h_sp.rds")
saveRDS(infarto_5dpi_m,"./objects/initial/5dpi_m_sp.rds")
