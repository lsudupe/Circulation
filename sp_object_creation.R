## SCRIPT: Spatial seurat object creation with visium data CIRCULATION project

## 29.03.22 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)

#Data--------------------------------------
control <- Load10X_Spatial(data.dir = "./data/control/",
                           filename = "filtered_feature_bc_matrix.h5",
                           assay = "Spatial",
                           slice = "control",
                           filter.matrix = TRUE)

infarto_3dpi <- Load10X_Spatial(data.dir = "./data/infarto_3dpi/",
                           filename = "filtered_feature_bc_matrix.h5",
                           assay = "Spatial",
                           slice = "dpi_3",
                           filter.matrix = TRUE)

infarto_5dpi_h <- Load10X_Spatial(data.dir = "./data/infarto_5dpi_hembra/",
                                filename = "filtered_feature_bc_matrix.h5",
                                assay = "Spatial",
                                slice = "dpi_5_female",
                                filter.matrix = TRUE)

infarto_5dpi_m <- Load10X_Spatial(data.dir = "./data/infarto_5dpi_macho/",
                                  filename = "filtered_feature_bc_matrix.h5",
                                  assay = "Spatial",
                                  slice = "dpi_5_male",
                                  filter.matrix = TRUE)
###########################################################
#save objects
saveRDS(control,"./objects/initial/control_sp.rds")
saveRDS(infarto_3dpi,"./objects/initial/3dpi_sp.rds")
saveRDS(infarto_5dpi_h,"./objects/initial/5dpi_h_sp.rds")
saveRDS(infarto_5dpi_m,"./objects/initial/5dpi_m_sp.rds")


