library(Seurat)
library(patchwork)
library(ggplot2)
library(readxl)
library(dplyr)
setwd("Z:/Cailab/Qian_writting/pyroptosis_covid/reply_to_reviewer/suplementaryData1")
dataset1 <- readRDS("2069-Allcells.counts.rds")
meta <- read.csv('2080-Allcells.meta.data.patientID_added.csv')
rownames(meta) <- meta$X
table(meta$Condition)
# we only use cells from COVID-19 positive patients
# therefore, we only have two conditions in this BALF2 sample: moderate and severe
# BAL009, BAL037 are the two patients from the moderate group
# the rest 20 patients are from the severe group

meta <- meta[meta$Condition %in% c('Moderate','Severe'),]
# fiter the cell type

table(meta$Condition)
dataset1 <- dataset1[,colnames(dataset1) %in% meta$Cell]
dim(dataset1)

d1 <- CreateSeuratObject(counts = dataset1,meta.data = meta)
d1 <- NormalizeData(object = d1)
d1 <- FindVariableFeatures(d1)
d1 <- ScaleData(d1)
d1 <- RunPCA(d1)
d1 <- RunUMAP(d1,dims = 1:50)
DimPlot(d1)
DimPlot(d1,label = T,repel = T,group.by = 'Condition')
Idents(d1) <- d1$CellType
d1$CellType
# we only use the eight cell types to keep consistency
d1 <- subset(d1,idents = c('Alveolar_macrophage', 'B_cell',       
                           'CD4_Tcell','CD8_Tcell','cDC',
                           'Md_macrophage','Monocyte','Neutrophil',
                           'Neutrophil','NK'))
DimPlot(d1,label = T,repel = T,group.by = 'Condition')
DimPlot(d1,label = T,repel = T,group.by = 'CellType')

saveRDS(d1,'Z:/Cailab/Qian_writting/pyroptosis_covid/reply_to_reviewer/Data/BALF2.rds')
