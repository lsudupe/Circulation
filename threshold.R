### SCRIPT: Select negative control threshold and use in the other distributions

## 03.04.22 Laura Sudupe , git @lsudupe

###Packages
library("Seurat")


###Read the data
GFP <- readRDS("objects/integrated/integrated.sct.GFP.rds")

DefaultAssay(GFP) <- "Spatial"
i <- "GFP"
pdf(file.path("./results/threshold/",filename = paste(i,"threshold.pdf",sep="")))
print(dittoRidgePlot(GFP, i, group.by = "sample"))
dev.off()

SpatialFeaturePlot(object = gfp, 
                   features = c("GFP"),
                   combine = FALSE)
