### SCRIPT: Check and manipulate spatial data

## 10.04.22 Laura Sudupe , git @lsudupe

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

####SUBSET THE SPATIAL DATA

spatial@meta.data[["ident"]] <- spatial@meta.data[["seurat_clusters"]]
Seurat::Idents(object = spatial) <- spatial@meta.data[["ident"]]

spatial.fb <- subset(x = spatial, idents = c("1","4","5","6"))
spatial.no.fb <- subset(x = spatial, idents = c("0","2","3","7","8"))

spatial.fb@meta.data[["ident"]] <- spatial.fb@active.ident
spatial.no.fb@meta.data[["ident"]] <- spatial.no.fb@active.ident

saveRDS(spatial.fb, file = "./objects/integrated/integrated.gfp.fb.rds")
saveRDS(spatial.no.fb, file = "./objects/integrated/integrated.gfp.nofb.rds")
