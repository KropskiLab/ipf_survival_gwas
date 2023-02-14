library(Seurat)
library(dplyr)
library(ggplot2)
library(harmony)
library(SeuratDisk)
options(future.globals.maxSize = 40000 * 1024^2)

VlnPlot(ild_all, features = 'PCSK6', split.by='dataset')
Idents(ild_all) <- 'dataset'
ild_all <- subset(ild_all, idents = c('Yale/BWH'))
gc()

ild_all <- NormalizeData(ild_all)
ild_all <- FindVariableFeatures(ild_all, nfeatures = 3000, selection.method = "vst")
ild_all <- ScaleData(ild_all)
ild_all <- RunPCA(ild_all)
ild_all <- RunUMAP(ild_all, dims = 1:35, verbose = F)
DimPlot(ild_all, group.by = 'orig.ident')

library(SeuratDisk)
SaveH5Seurat(ild_all, file='~/Desktop/ild_all/yale.h5seurat')
ild_all <- LoadH5Seurat(file='~/Desktop/ild_all/yale.h5seurat')
vumc <- readRDS('~/Desktop/ild_all/sa_2020_reannotated.rds')
Idents(vumc) <- 'celltype_new'
gc()

ild.anchors <- FindTransferAnchors(reference = vumc, query = ild_all,
                                        dims = 1:35, reference.reduction = "pca")
predictions <- TransferData(anchorset = ild.anchors, refdata = vumc$celltype_new,
                            dims = 1:35)
ild_all <- AddMetaData(ild_all, metadata = predictions)
DimPlot(ild_all, group.by = 'predicted.id', cols = 'polychrome')
Idents(ild_all) <- 'predicted.id'


DotPlot(ild_all, features = c('PCSK6'), split.by = 'Status')
DotPlot(vumc, features = c('PCSK6'), split.by = 'Status')


DotPlot(ild_all, features = c('PCSK6'), cols = c('grey', 'black'))
DotPlot(vumc, features = c('PCSK6'), cols = c('grey', 'black'))


SaveH5Seurat(ild_all, file='~/Desktop/ild_all/yale_annotated.h5seurat')


vumc <- readRDS("~/Downloads/ild_all.rds")
VlnPlot(vumc, features = 'PCSK6', group.by='Diagnosis')
Idents(vumc) <- 'dataset'
vumc <- subset(vumc, idents = c('VUMC/TGen'))
Idents(vumc) <- 'Diagnosis'
vumc <- subset(vumc, idents = c('Control', 'IPF'))
gc()

vumc <- FindVariableFeatures(vumc, nfeatures = 3000, selection.method = "vst")
vumc <- ScaleData(vumc)
vumc <- RunPCA(vumc)
vumc <- RunHarmony(vumc, max.iter.harmony = 50, group.by.vars = 'Sample_Name')
vumc <- RunUMAP(vumc, reduction = 'harmony', dims = 1:35)
DimPlot(vumc, group.by = 'orig.ident')

sa_2020 <- readRDS("Desktop/ild_all/sa_2020_reannotated.rds")

ild.anchors <- FindTransferAnchors(reference = sa_2020, query = vumc,
                                   dims = 1:35, reference.reduction = "pca")
predictions <- TransferData(anchorset = ild.anchors, refdata = sa_2020$celltype_new,
                            dims = 1:35)
vumc <- AddMetaData(vumc, metadata = predictions)
DimPlot(vumc, group.by = 'predicted.id', cols = 'polychrome')
Idents(vumc) <- 'predicted.id'
VlnPlot(vumc, features = 'PCSK6', sort = T)



SaveH5Seurat(vumc, file='~/Desktop/ild_all/vumc_ipf.h5seurat')



vumc$predicted.id <- factor(x = vumc$predicted.id, levels = c("AT1", "AT2",  "Basal", "Ciliated", "KRT5-/KRT17+", "PNEC/Ionocyte", "Proliferating epithelial", "Secretory - MUC5B+", "Secretory - SCGB1A1+/SCGB3A2+", "Secretory - SCGB3A2+","Transitional", 
"Arteriole", "aCap", 'gCap', 'Peribronchiolar', 'Venule', 'Lymphatic', "Adventitial FB", "FB - HAS1+", "FB - WNT2+", "MyoFB", "MyoFB - activated", "Mesothelial", "Pericyte", "SMC", "Proliferating FB",
"Macrophage", "Monocyte", "moDC", "cDC", "pDC", "Mast", "NK", "CD4", "CD8","Treg",  "B cells", "Plasma", "Proliferating immune"))

DotPlot(vumc, features = 'PCSK6', group.by = 'predicted.id', cols = c('grey', 'black'))
gc()

yale <- LoadH5Seurat('~/Desktop/ild_all/yale_annotated.h5seurat')

yale$predicted.id <- factor(x = yale$predicted.id, levels = c("AT1", "AT2",  "Basal", "Ciliated", "KRT5-/KRT17+", "PNEC/Ionocyte", "Proliferating epithelial", "Secretory - MUC5B+", "Secretory - SCGB1A1+/SCGB3A2+", "Secretory - SCGB3A2+","Transitional", 
                                                              "Arteriole", "aCap", 'gCap', 'Peribronchiolar', 'Venule', 'Lymphatic', "Adventitial FB", "FB - HAS1+", "FB - WNT2+", "MyoFB", "MyoFB - activated", "Mesothelial", "Pericyte", "SMC", "Proliferating FB",
                                                              "Macrophage", "Monocyte", "moDC", "cDC", "pDC", "Mast", "NK", "CD4", "CD8","Treg",  "B cells", "Plasma", "Proliferating immune"))


DotPlot(yale, features = 'PCSK6', group.by = 'predicted.id', cols = c('grey', 'black'))
DotPlot(yale, features = 'PCSK6', group.by = 'predicted.id',cols = c('light blue', 'red')) + coord_flip()+theme(axis.text.x = element_text(angle = 90, hjust=1))


SaveH5Seurat(yale, file='~/Desktop/ild_all/yale_annotated.h5seurat', overwrite = T)
Convert('~/Desktop/ild_all/yale_annotated.h5seurat', dest = 'h5ad')
Convert('~/Desktop/ild_all/vumc_ipf.h5seurat', dest = 'h5ad')

