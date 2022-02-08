library(anndata)
library(Seurat)
library(patchwork)

setwd("Z:/Cailab/Qian_writting/pyroptosis_covid/healthy_lung_upper_airway")
ad <- read_h5ad('madissoon19_lung.processed.h5ad')
meta <- ad$obs
table(meta$CellType)
ad <- ad$X

set.seed(111)
macrophage <- rownames(meta)[meta$CellType %in% c('Macrophage_MARCOneg','Macrophage_MARCOpos')]
macrophage <- sample(x = macrophage, size = 800, replace = F)
ad_macrophage <- ad[rownames(ad) %in% macrophage,]
ad_macrophage <- as.matrix(ad_macrophage)
ad_macrophage <- t(ad_macrophage)
ad_macrophage <- CreateSeuratObject(ad_macrophage,project = 'Macrophage')
saveRDS(ad_macrophage,file = 'Healthymacrophage.rds')

rm(ad_macrophage)
set.seed(111)
NK <- rownames(meta)[meta$CellType %in% c('NK')]
NK <- sample(x = NK, size = 800, replace = F)
ad_NK <- ad[rownames(ad) %in% NK,]
ad_NK <- as.matrix(ad_NK)
ad_NK <- t(ad_NK)
ad_NK <- CreateSeuratObject(ad_NK,project = 'NK cell')
saveRDS(ad_NK,file = 'HealthyNK.rds')
dim(ad_NK)
rm(ad_NK)

set.seed(111)
T_CD4 <- rownames(meta)[meta$CellType %in% c('T_CD4')]
T_CD4 <- sample(x = T_CD4, size = 800, replace = F)
ad_T_CD4 <- ad[rownames(ad) %in% T_CD4,]
ad_T_CD4 <- as.matrix(ad_T_CD4)
ad_T_CD4 <- t(ad_T_CD4)
ad_T_CD4 <- CreateSeuratObject(ad_T_CD4,project = 'CD4+ T-cell')
saveRDS(ad_T_CD4,file = 'HealthyT_CD4.rds')
rm(ad_T_CD4)

set.seed(111)
T_CD8_CytT <- rownames(meta)[meta$CellType %in% c('T_CD8_CytT')]
T_CD8_CytT  <- sample(x=T_CD8_CytT ,size = 800,replace = F)
ad_T_CD8_CytT <- ad[rownames(ad) %in% T_CD8_CytT,]
ad_T_CD8_CytT <- as.matrix(ad_T_CD8_CytT)
ad_T_CD8_CytT <- t(ad_T_CD8_CytT)
ad_T_CD8_CytT <- CreateSeuratObject(ad_T_CD8_CytT,project = 'CD8+ T-cell')
saveRDS(ad_T_CD8_CytT,file = 'HealthyT_CD8_CytT.rds')
rm(ad_T_CD8_CytT)

set.seed(111)
B_cells <- rownames(meta)[meta$CellType %in% c('B_cells')]
B_cells <- sample(B_cells,size = 800,replace = T)
ad_B_cells <- ad[rownames(ad) %in% B_cells,]
ad_B_cells <- as.matrix(ad_B_cells)
ad_B_cells <- t(ad_B_cells)
ad_B_cells <- CreateSeuratObject(ad_B_cells,project = 'B cell')
saveRDS(ad_B_cells,file = 'HealthyB_cells.rds')
rm(ad_B_cells)

Healthymacrophage <- readRDS("Z:/Cailab/Qian_writting/pyroptosis_covid/healthy_lung_upper_airway/Healthymacrophage.rds")
HealthyNK <- readRDS("Z:/Cailab/Qian_writting/pyroptosis_covid/healthy_lung_upper_airway/HealthyNK.rds")
A <- merge(Healthymacrophage,HealthyNK )
rm(Healthymacrophage, HealthyNK)
HealthyT_CD4 <- readRDS("Z:/Cailab/Qian_writting/pyroptosis_covid/healthy_lung_upper_airway/HealthyT_CD4.rds")
A <- merge(A,HealthyT_CD4)
rm(HealthyT_CD4)
HealthyT_CD8_CytT <- readRDS("Z:/Cailab/Qian_writting/pyroptosis_covid/healthy_lung_upper_airway/HealthyT_CD8_CytT.rds")
A <- merge(A,HealthyT_CD8_CytT)
rm(HealthyT_CD8_CytT)
HealthyB_cells <- readRDS("Z:/Cailab/Qian_writting/pyroptosis_covid/healthy_lung_upper_airway/HealthyB_cells.rds")
A <- merge(A,HealthyB_cells)
rm(HealthyB_cells)

A <- NormalizeData(A)
A <- FindVariableFeatures(A)
A <- ScaleData(A)
A <- RunPCA(A)
A <- RunUMAP(A,dims=1:50)
DimPlot(A,label = T,repel = T)
table(A$orig.ident)
saveRDS(A,file = 'Z:/Cailab/Qian_writting/pyroptosis_covid/Integrate/Data/healthy_lung2.rds')


