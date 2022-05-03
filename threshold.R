### SCRIPT: Select negative control threshold and use in the other distributions

## 03.04.22 Laura Sudupe , git @lsudupe

###Packages
library("Seurat")


###Read the data
GFP <- readRDS("objects/integrated/integrated.sct.GFP.rds")


####Plot the distributions
DefaultAssay(GFP) <- "Spatial"
i <- "GFP"
pdf(file.path("./results/threshold/",filename = paste(i,"threshold.sample.pdf",sep="")))
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
