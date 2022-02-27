library(Seurat)
library(patchwork)
library(ggplot2)
library(harmony)
library(ggpubr)
library(dplyr)
library(anndata)
setwd("Z:/Cailab/Qian_writting/pyroptosis_covid/reply_to_reviewer")
BALF <- readRDS('Z:/Cailab/Qian_writting/pyroptosis_covid/Integrate/Data/BALF.rds')
table(BALF$orig.ident)
table(BALF$condition)
table(BALF$celltype)
SF2_BALF <- DimPlot(BALF, group.by = 'orig.ident')+ggtitle('BALF1')

healthy_lung <- readRDS("Z:/Cailab/Qian_writting/pyroptosis_covid/Integrate/Data/healthy_lung2.rds")
setwd("Z:/Cailab/Qian_writting/pyroptosis_covid/healthy_lung_upper_airway")
ad <- read_h5ad('madissoon19_lung.processed.h5ad')
meta <- ad$obs
rm(ad)
healthy_lung$patientID <- meta$Donor[match(colnames(healthy_lung),rownames(meta))]
healthy_lung$celltype <- Idents(healthy_lung)
SF2_healthy_lung <- DimPlot(healthy_lung, group.by = 'patientID')+ggtitle('CTRL')
rm(meta)

blish <- readRDS('Z:/Cailab/Qian_writting/pyroptosis_covid/Integrate/Data/blish2.rds')
meta <- blish@meta.data
table(blish$Ventilated)
blish$celltype <- blish$cell.type
table(blish$celltype)
blish <- NormalizeData(blish)
blish <- FindVariableFeatures(blish)
blish <- ScaleData(blish)
blish <- RunPCA(blish)
blish <- RunUMAP(blish,dims = 1:20)
SF2_PBMC <- DimPlot(blish, group.by = 'orig.ident')+ggtitle('PBMC')

BALF2 <- readRDS('Z:/Cailab/Qian_writting/pyroptosis_covid/reply_to_reviewer/Data/BALF2.rds')
meta <- read.csv('Z:/Cailab/Qian_writting/pyroptosis_covid/reply_to_reviewer/suplementaryData1/2080-Allcells.meta.data.patientID_added.csv')
table(meta$PatientNumber)
rownames(meta) <- meta$Cell
dim(BALF2)
BALF2 <- BALF2[,colnames(BALF2) %in% meta$Cell[meta$CellType %in% c('B_cell', 'CD4_Tcell','CD8_Tcell','cDC',
                                                                                'Md_macrophage','Monocyte','Neutrophil','NK')]]
dim(BALF2)
table(BALF2$CellType)
SF2_BALF2 <- DimPlot(BALF2, group.by = 'orig.ident')+ggtitle('BALF2')
png('Figure/Supplementary Fig. S1.png',res = 300,height = 2000,width = 3000 )
SF2_BALF+SF2_BALF2+SF2_healthy_lung+SF2_PBMC
dev.off()
