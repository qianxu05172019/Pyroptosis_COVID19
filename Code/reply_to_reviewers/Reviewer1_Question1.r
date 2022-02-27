library(Seurat)
library(patchwork)
library(ggplot2)
library(harmony)
library(ggpubr)

BALF <- readRDS('Z:/Cailab/Qian_writting/pyroptosis_covid/Integrate/Data/BALF.rds')
healthy_lung <- readRDS("Z:/Cailab/Qian_writting/pyroptosis_covid/Integrate/Data/healthy_lung2.rds")
healthy_lung$celltype <- Idents(healthy_lung)

blish <- readRDS('Z:/Cailab/Qian_writting/pyroptosis_covid/Integrate/Data/blish2.rds')
blish$celltype <- blish$cell.type

A <- readRDS("Z:/Cailab/Qian_writting/pyroptosis_covid/Integrate/Data/Integrated.rds")
################ P1A1 ######################################
PBMC <- subset(A,cells = colnames(A)[colnames(A) %in% colnames(blish)])
S1PBMC <- DimPlot(PBMC,group.by = 'orig.ident',cols = c("#a6cee3","#1f78b4", "#b2df8a",
                                              '#33a02c','#fb9a99','#e31a1c',
                                              '#fdbf6f','#ff7f00','#cab2d6',
                                              '#6a3d9a','#ffff99','#b15928',
                                              '#003c30','#80cdc1'))+xlab(' ')+ylab(' ')+
  ggtitle('PBMC')+theme(panel.grid = element_blank())+
  theme(axis.text = element_blank())+
  theme(axis.ticks = element_blank())+
  theme(panel.border = element_blank())+
  theme(axis.line = element_blank())

################ P1A3 ######################################
BALF <- subset(A,cells = colnames(A)[colnames(A) %in% colnames(BALF)])
S1BALF <- DimPlot(BALF,group.by = 'orig.ident',cols = c('#1b9e77','#d95f02','#7570b3'
                                              ,'#e7298a','#66a61e','#e6ab02'))+
  xlab(' ')+ylab(' ')+ggtitle('BALF')+theme(panel.grid = element_blank())+
  theme(axis.text = element_blank())+
  theme(axis.ticks = element_blank())+
  theme(panel.border = element_blank())+
  theme(axis.line = element_blank())

################ P1A4 ######################################
healthy_lung <- subset(A,cells = colnames(A)[colnames(A) %in% colnames(healthy_lung)])
ad <- read_h5ad('Z:/Cailab/Qian_writting/pyroptosis_covid/healthy_lung_upper_airway/madissoon19_lung.processed.h5ad')
meta <- ad$obs
table(meta$CellType)
ad <- ad$X
healthy_lung$Donor <- meta$Donor[match(colnames(healthy_lung),rownames(meta))]
table(healthy_lung$Donor)

S1HealthyLung <- DimPlot(healthy_lung,group.by = 'Donor',
        cols = c( '#1b9e77','#d95f02','#7570b3'
                  ,'#e7298a','#66a61e','#e6ab02'))+xlab(' ')+ylab(' ')+
  ggtitle('Healthy lung')+theme(panel.grid = element_blank())+
  theme(axis.text = element_blank())+
  theme(axis.ticks = element_blank())+
  theme(panel.border = element_blank())+
  theme(axis.line = element_blank())


png('Z:/Cailab/Qian_writting/pyroptosis_covid/Integrate/Figure/FigS1.png',res = 700,height = 3000,width = 9000)
S1PBMC+S1BALF+S1HealthyLung
dev.off()

