### SCRIPT: Select negative control threshold and use in the other distributions

## 03.04.22 Laura Sudupe , git @lsudupe

###Packages
library("Seurat")
library("data.table")
library("dplyr")
library("dittoSeq")


############################GFP###############################3
###Read the data
GFP <- readRDS("objects/integrated/integrated.sct.GFP.rds")

####Plot the distributions
DefaultAssay(GFP) <- "Spatial"
i <- "GFP"
pdf(file.path("./results/threshold/",filename = paste(i,"threshold.sample.line.pdf",sep="")))
print(dittoRidgePlot(GFP, i, group.by = "sample",add.line = 10,min=0, max=100))#,split.by = "sample"))
dev.off()

###Where is my 99 threshold
data <- dittoRidgePlot(GFP, i, group.by = "seurat_clusters",split.by = "sample")
data <- as.data.frame(data[["data"]])
control <- data[grepl("control", data[,4]),]
counts.control <- as.vector(control$var.data)
counts.all <- as.vector(data$var.data)

quantile(counts.control, prob)
quantile <- quantile(counts.control, probs = c(0.99)) #171.3 es 99%, 96 es 95%

ggplot(data, aes(x=var.data, colour=sample)) +
  geom_density() +
  xlim(0, 100) +
  geom_vline(data=data, aes(xintercept=quantile,  colour=sample),
           linetype="dashed", size=1)

SpatialFeaturePlot(object = gfp, 
                   features = c("GFP"),
                   combine = FALSE)
########################################################################

###########################marker genes##############################3
integrated <- readRDS("objects/integrated/integrated.sct.rds")

fibroblast <- c("Col1a1","Col3a1","Ddr2","Ly6a","Thy1","S100a4","Tcl21",
                "Crispld2","Ogn","Islr","Pdgfra","Dcn","Gsn")
circulation.fibroblast <- c("Cthrc1","Ddah1","Fmod","Lox","Comp","Ptn",
                            "Gsn","Pdgfra","Ddr2","Tcf21","P4hb","Acta2",
                            "Postn","Ckap4","Fstl1")

DefaultAssay(integrated) <- "integrated"
pdf(file.path("./results/threshold",filename ="integrated.marker.genes.pdf"))
DoHeatmap(subset(integrated, downsample = 100), features = fibroblast, size = 3)
dev.off()

DefaultAssay(integrated) <- "integrated"
pdf(file.path("./results/threshold",filename ="integrated.circulation.genes.pdf"))
DoHeatmap(subset(integrated, downsample = 100), features = circulation.fibroblast, size = 3)
dev.off()


########################################################################

###########################DECO##############################3
spatial <- readRDS("objects/integrated/integrated.sct.rds")
spotlight_ls <- readRDS("./results/deco.ibex/fb_deco/spot_integrated_fibro.rds")


## separate list
nmf_mod <- spotlight_ls[[1]]
decon_mtrx <- spotlight_ls[[2]]

##############0.0 filter
#Visualization, join decomposition with metadata
rownames(decon_mtrx) <- colnames(spatial)
decon_mtrx_sub <- decon_mtrx[, colnames(decon_mtrx) != "res_ss"]


spatial_min_0 <- spatial
decon_mtrx_sub_min_0 <- decon_mtrx_sub
decon_mtrx_min0 <- decon_mtrx

decon_mtrx_sub_min_0[decon_mtrx_sub_min_0 < 0.5] <- 0
decon_mtrx_min0 <- cbind(decon_mtrx_sub_min_0, "res_ss" = decon_mtrx_min0[, "res_ss"])

decon_df_min0 <- decon_mtrx_min0 %>%
  data.frame() %>%
  tibble::rownames_to_column("barcodes")

spatial_min_0@meta.data <- spatial_min_0@meta.data %>%
  tibble::rownames_to_column("barcodes") %>%
  dplyr::left_join(decon_df_min0, by = "barcodes") %>%
  tibble::column_to_rownames("barcodes")
#####################


matrix <- as.data.frame(decon_mtrx)
matrix["samples"] <- data[,4]
matrix$res_ss <- NULL

control <- matrix[grepl("control", matrix[,3]),]
counts.control <- as.vector(control$no.FB)
counts.all <- as.vector(data$var.data)


quantile <- quantile(counts.control, probs = c(0.01)) #0.519 es 99% fb, 0.658 no fb 99, 0.481 no.fib 0.01%

i <- "no.FB"
pdf(file.path("./results/threshold/",filename = paste(i,"threshold.sample.line.deco.0.0.no.fibro.pdf",sep="")))
print(dittoRidgePlot(spatial_min_0, i, group.by = "sample",add.line = 0.481))#,min=0, max=100))#,split.by = "sample"))
dev.off()







