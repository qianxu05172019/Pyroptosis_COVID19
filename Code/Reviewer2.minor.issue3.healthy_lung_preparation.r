library(anndata)
library(Seurat)
library(patchwork)

setwd("Z:/Cailab/Qian_writting/pyroptosis_covid/healthy_lung_upper_airway")
ad <- read_h5ad('madissoon19_lung.processed.h5ad')
meta <- ad$obs
table(meta$CellType)
ad <- ad$X

set.seed(111)
macrophage <- rownames(meta)[meta$CellType %in% c('Macrophage_Dividing','Macrophage_MARCOneg','Macrophage_MARCOpos')]
macrophage <- sample(x = macrophage, size = 1000, replace = F)
ad_macrophage <- ad[rownames(ad) %in% macrophage,]
ad_macrophage <- as.matrix(ad_macrophage)
ad_macrophage <- t(ad_macrophage)
ad_macrophage <- CreateSeuratObject(ad_macrophage,project = 'Healthy macrophage')
saveRDS(ad_macrophage,file = 'Healthymacrophage.rds')

rm(ad_macrophage)

set.seed(111)
NK <- rownames(meta)[meta$CellType %in% c('NK_Dividing','NK')]
NK <- sample(x = NK, size = 1000, replace = F)
ad_NK <- ad[rownames(ad) %in% NK,]
ad_NK <- as.matrix(ad_NK)
ad_NK <- t(ad_NK)
ad_NK <- CreateSeuratObject(ad_NK,project = 'Healthy NK')
saveRDS(ad_NK,file = 'HealthyNK.rds')
dim(ad_NK)
rm(ad_NK)

set.seed(111)
T_CD4 <- rownames(meta)[meta$CellType %in% c('T_CD4')]
T_CD4 <- sample(x = T_CD4, size = 1000, replace = F)
ad_T_CD4 <- ad[rownames(ad) %in% T_CD4,]
ad_T_CD4 <- as.matrix(ad_T_CD4)
ad_T_CD4 <- t(ad_T_CD4)
ad_T_CD4 <- CreateSeuratObject(ad_T_CD4,project = 'Healthy T_CD4')
saveRDS(ad_T_CD4,file = 'HealthyT_CD4.rds')
rm(ad_T_CD4)

set.seed(111)
T_CD8_CytT <- rownames(meta)[meta$CellType %in% c('T_CD8_CytT')]
T_CD8_CytT  <- sample(x=T_CD8_CytT ,size = 1000,replace = F)
ad_T_CD8_CytT <- ad[rownames(ad) %in% T_CD8_CytT,]
ad_T_CD8_CytT <- as.matrix(ad_T_CD8_CytT)
ad_T_CD8_CytT <- t(ad_T_CD8_CytT)
ad_T_CD8_CytT <- CreateSeuratObject(ad_T_CD8_CytT,project = 'Healthy T_CD8_CytT')
saveRDS(ad_T_CD8_CytT,file = 'HealthyT_CD8_CytT.rds')
rm(ad_T_CD8_CytT)

set.seed(111)
B_cells <- rownames(meta)[meta$CellType %in% c('B_cells')]
B_cells <- sample(B_cells,size = 1000,replace = T)
ad_B_cells <- ad[rownames(ad) %in% B_cells,]
ad_B_cells <- as.matrix(ad_B_cells)
ad_B_cells <- t(ad_B_cells)
ad_B_cells <- CreateSeuratObject(ad_B_cells,project = 'Healthy B_cells')
saveRDS(ad_B_cells,file = 'HealthyB_cells.rds')
rm(ad_B_cells)

Monocyte <- rownames(meta)[meta$CellType %in% c('Monocyte')]
length(Monocyte)
set.seed(111)
Monocyte <- sample(x=Monocyte,size = 1000,replace = F)
ad_Monocyte <- ad[rownames(ad) %in% Monocyte,]
ad_Monocyte <- as.matrix(ad_Monocyte)
ad_Monocyte <- t(ad_Monocyte)
ad_Monocyte <- CreateSeuratObject(ad_Monocyte,project = 'Healthy Monocyte')
saveRDS(ad_Monocyte,file = 'HealthyMonocyte.rds')
rm(ad_Monocyte)

DC <- rownames(meta)[meta$CellType %in% c('DC_1','DC_2')]
set.seed(111)
DC <- sample(x=DC,size = 1000,replace = F)
ad_DC <- ad[rownames(ad) %in% DC,]
ad_DC <- as.matrix(ad_DC)
ad_DC <- t(ad_DC)
ad_DC <- CreateSeuratObject(ad_DC,project = 'Healthy DC')
saveRDS(ad_DC,file = 'HealthyDC.rds')
rm(ad_DC)


Healthymacrophage <- readRDS("Z:/Cailab/Qian_writting/pyroptosis_covid/healthy_lung_upper_airway/Healthymacrophage.rds")
#HealthyFibroblast <- readRDS("Z:/Cailab/Qian_writting/pyroptosis_covid/healthy_lung_upper_airway/HealthyFibroblast.rds")
#A <- merge(Healthymacrophage,HealthyFibroblast)
A <- Healthymacrophage
rm(Healthymacrophage,HealthyFibroblast)
HealthyNK <- readRDS("Z:/Cailab/Qian_writting/pyroptosis_covid/healthy_lung_upper_airway/HealthyNK.rds")
A <- merge(A,HealthyNK )
rm(HealthyNK)
HealthyT_CD4 <- readRDS("Z:/Cailab/Qian_writting/pyroptosis_covid/healthy_lung_upper_airway/HealthyT_CD4.rds")
A <- merge(A,HealthyT_CD4)
rm(HealthyT_CD4)
HealthyT_CD8_CytT <- readRDS("Z:/Cailab/Qian_writting/pyroptosis_covid/healthy_lung_upper_airway/HealthyT_CD8_CytT.rds")
A <- merge(A,HealthyT_CD8_CytT)
rm(HealthyT_CD8_CytT)
HealthyB_cells <- readRDS("Z:/Cailab/Qian_writting/pyroptosis_covid/healthy_lung_upper_airway/HealthyB_cells.rds")
A <- merge(A,HealthyB_cells)
rm(HealthyB_cells)
HealthyMonocyte <- readRDS("Z:/Cailab/Qian_writting/pyroptosis_covid/healthy_lung_upper_airway/HealthyMonocyte.rds")
A <- merge(A,HealthyMonocyte)
rm(HealthyMonocyte)
HealthyDC <- readRDS("Z:/Cailab/Qian_writting/pyroptosis_covid/healthy_lung_upper_airway/HealthyDC.rds")
A <- merge(A,HealthyDC)
rm(HealthyDC)

A <- NormalizeData(A)
A <- FindVariableFeatures(A)
A <- ScaleData(A)
A <- RunPCA(A)
A <- RunUMAP(A,dims=1:50)
DimPlot(A)
table(A$orig.ident)
saveRDS(A,file = 'healthy_lung.rds')
