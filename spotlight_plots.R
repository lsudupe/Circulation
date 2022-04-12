### SCRIPT: Analyze deconvolution results

## 07.11.21 Laura Sudupe , git @lsudupe

#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")


#install.packages("Seurat")
#install.packages("devtools")
#install.packages("base")

#devtools::install_github("https://github.com/MarcElosua/SPOTlight")

library("Seurat")
library("SPOTlight")
library("base")
library("dplyr")
library("gt")
library("DT")
library("ggplot2")
library("data.table")

##Read the data
spatial <- readRDS("./objects/integrated/integrated.gfp.rds")
spotlight_ls <- readRDS("./objects/spotligth/spot_integrated_fibro.rds")

###Select the path for the results
paths <- "./results/spotlight/"

## separate list
nmf_mod <- spotlight_ls[[1]]
decon_mtrx <- spotlight_ls[[2]]

#each topic profile specifity
h <- NMF::coef(nmf_mod[[1]])
rownames(h) <- paste("Topic", 1:nrow(h), sep = "_")
topic_profile_plts <- SPOTlight::dot_plot_profiles_fun(
  h = h,
  train_cell_clust = nmf_mod[[2]])

pdf("topic profile.pdf")
topic_profile_plts[[2]] + ggplot2::theme(
  axis.text.x = ggplot2::element_text(angle = 90), 
  axis.text = ggplot2::element_text(size = 12))
dev.off()

#Visualization, join decomposition with metadata
rownames(decon_mtrx) <- colnames(spatial)
decon_mtrx_sub <- decon_mtrx[, colnames(decon_mtrx) != "res_ss"]

##########min contribution == 0
spatial_min_0 <- copy(spatial)
decon_mtrx_sub_min_0 <- copy(decon_mtrx_sub)
decon_mtrx_min0 <- copy(decon_mtrx)

decon_mtrx_sub_min_0[decon_mtrx_sub_min_0 < 0.0] <- 0
decon_mtrx_min0 <- cbind(decon_mtrx_sub_min_0, "res_ss" = decon_mtrx_min0[, "res_ss"])

decon_df_min0 <- decon_mtrx_min0 %>%
  data.frame() %>%
  tibble::rownames_to_column("barcodes")

spatial_min_0@meta.data <- spatial_min_0@meta.data %>%
  tibble::rownames_to_column("barcodes") %>%
  dplyr::left_join(decon_df_min0, by = "barcodes") %>%
  tibble::column_to_rownames("barcodes")

cell_types_all_0.0 <- colnames(decon_mtrx_min0)[which(colnames(decon_mtrx_min0) != "res_ss")]

pdf(file.path(paths,filename ="Fibro deco.pdf"))
print(ridgeEnrichment(spatial_min_0@meta.data, gene.set = "FB", group = "seurat_clusters", facet = "sample", add.rug = TRUE))
dev.off()

pdf(file.path(paths,filename ="no.FB deco.pdf"))
print(ridgeEnrichment(spatial_min_0@meta.data, gene.set = "no.FB", group = "seurat_clusters", facet = "sample", add.rug = TRUE))
dev.off()

###########predictions 0.0
pdf("FB.predictions.0.pdf")
SpatialPlot(object = spatial_min_0, features = cell_types_all_0.0)
dev.off()

##########min contribution == 0.1
spatial_min_0.1 <- copy(spatial)
decon_mtrx_sub_min_0.1 <- copy(decon_mtrx_sub)
decon_mtrx_min0.1 <- copy(decon_mtrx)

decon_mtrx_sub_min_0.1[decon_mtrx_sub_min_0.1 < 0.4] <- 0
decon_mtrx_min0.1 <- cbind(decon_mtrx_sub_min_0.1, "res_ss" = decon_mtrx_min0.1[, "res_ss"])

decon_df_min0.1 <- decon_mtrx_min0.1 %>%
  data.frame() %>%
  tibble::rownames_to_column("barcodes")

spatial_min_0.1@meta.data <- spatial_min_0.1@meta.data %>%
  tibble::rownames_to_column("barcodes") %>%
  dplyr::left_join(decon_df_min0.1, by = "barcodes") %>%
  tibble::column_to_rownames("barcodes")

cell_types_all_0.1 <- colnames(decon_mtrx_min0.1)[which(colnames(decon_mtrx_min0.1) != "res_ss")]

pdf(file.path(paths,filename ="Fibro deco.pdf"))
print(ridgeEnrichment(spatial_min_0.1@meta.data, gene.set = "FB", group = "seurat_clusters", facet = "sample", add.rug = TRUE))
dev.off()

pdf(file.path(paths,filename ="no.FB deco.pdf"))
print(ridgeEnrichment(spatial_min_0.1@meta.data, gene.set = "no.FB", group = "seurat_clusters", facet = "sample", add.rug = TRUE))
dev.off()

###########predictions 0.1
pdf("FB.predictions.0.1.pdf")
SpatialPlot(object = spatial_min_0.1, features = cell_types_all_0.0)
dev.off()

##########min contribution == 0.5
spatial_min_0.5 <- copy(spatial)
decon_mtrx_sub_min_0.5 <- copy(decon_mtrx_sub)
decon_mtrx_min0.5 <- copy(decon_mtrx)

decon_mtrx_sub_min_0.5[decon_mtrx_sub_min_0.5 < 0.5] <- 0
decon_mtrx_min0.5 <- cbind(decon_mtrx_sub_min_0.5, "res_ss" = decon_mtrx_min0.5[, "res_ss"])

decon_df_min0.5 <- decon_mtrx_min0.5 %>%
  data.frame() %>%
  tibble::rownames_to_column("barcodes")

spatial_min_0.5@meta.data <- spatial_min_0.5@meta.data %>%
  tibble::rownames_to_column("barcodes") %>%
  dplyr::left_join(decon_df_min0.5, by = "barcodes") %>%
  tibble::column_to_rownames("barcodes")

cell_types_all_0.5 <- colnames(decon_mtrx_min0.5)[which(colnames(decon_mtrx_min0.5) != "res_ss")]

###########predictions 0.5
pdf("FB.predictions.0.5.pdf")
SpatialPlot(object = spatial_min_0.5, features = cell_types_all_0.0)
dev.off()

source("scatterpie_function.R")

###DECONVOLUTION PLOTS

#pdf("deco_0.pdf")
#scatterpie_plot(se_obj = spatial_min_0,
 #                         cell_types_all = cell_types_all_0.0,
  #                        slice = "control",
   #                       pie_scale = 0.4)
#dev.off()




pdf("control.04.pdf")
spatial_scatterpie_modify(se_obj = spatial_min_0.1,
                          cell_types_all = cell_types_all_0.0,
                          img_path = "./control/tissue_lowres_image.png",
                          slice = "control",
                          pie_scale = 0.7)
dev.off()

pdf("dpi3.04.pdf")
spatial_scatterpie_modify(se_obj = spatial_min_0.1,
                          cell_types_all = cell_types_all_0.0,
                          img_path = "./dpi3/tissue_lowres_image.png",
                          slice = "dpi_3",
                          pie_scale = 0.7)
dev.off()

pdf("dpi5h.04.pdf")
spatial_scatterpie_modify(se_obj = spatial_min_0.1,
                          cell_types_all = cell_types_all_0.0,
                          img_path = "./dpi5_h/tissue_lowres_image.png",
                          slice = "dpi_5_h",
                          pie_scale = 0.7)
dev.off()

pdf("dpi5m.04.pdf")
spatial_scatterpie_modify(se_obj = spatial_min_0.1,
                          cell_types_all = cell_types_all_0.0,
                          img_path = "./dpi5_m/tissue_lowres_image.png",
                          slice = "dpi_5_m",
                          pie_scale = 0.7)
dev.off()


#pdf("deco_day7_0.5_BvsD1vsD2vsOthers.pdf")
#spatial_scatterpie_modify(se_obj = spatial_min_0.5,
 #                         cell_types_all = cell_types_all_0.0,
  #                        img_path = "./visium/day7_spatial/sample1/spatial/tissue_lowres_image.png",
   #                       pie_scale = 0.4)
#dev.off()








