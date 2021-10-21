library(Seurat)
library(patchwork)
library(ggplot2)
library(harmony)
library(ggpubr)
setwd("Z:/Cailab/Qian_writting/pyroptosis_covid/Integrate")

BALF <- readRDS('Z:/Cailab/Qian_writting/pyroptosis_covid/Integrate/Data/BALF.rds')
table(BALF$orig.ident)
table(BALF$condition)
table(BALF$celltype)


healthy_lung <- readRDS("Z:/Cailab/Qian_writting/pyroptosis_covid/Integrate/Data/healthy_lung2.rds")
table(healthy_lung$orig.ident)
table(healthy_lung$orig.ident)
healthy_lung$celltype <- Idents(healthy_lung)

blish <- readRDS('Z:/Cailab/Qian_writting/pyroptosis_covid/Integrate/Data/blish2.rds')
table(blish$Ventilated)
blish$celltype <- blish$cell.type
table(blish$celltype)


################################# load data #############################################################
A <- merge(BALF,healthy_lung)
A <- merge(A,blish)
dim(A)
table(A$celltype)

A <- NormalizeData(A)
A <- FindVariableFeatures(A)
A <- ScaleData(A)
A <- RunPCA(A)
A <- RunUMAP(A,dims = 1:20)


A$Experiment <- NA
A$Experiment[A$orig.ident %in% c('COV002','COV004','COV007','COV012','COV034','COV036')] <- 'BALF'
A$Experiment[A$orig.ident %in% c("covid_555_1","covid_555_2","covid_556",
                                 "covid_557" , "covid_558" , "covid_559" ,
                                 "covid_560" , "covid_561",  "HIP002" , 
                                 "HIP015" , "HIP023","HIP043",
                                 "HIP044","HIP045")] <- 'PBMC'
A$Experiment[A$orig.ident %in% c("Macrophage","NK cell", "CD4+ T cell" , "CD8+ T cell","B cell",
                                 "Monocyte" ,"cDC")] <- 'Healthy lung'
table(A$Experiment)
table(A$celltype)
A$celltype[A$celltype %in% c('NK','NK_cell','NK cell')] <- 'NK'
A$celltype[A$celltype %in% c('Neutrophils','Neutrophil')] <- 'Neutrophil'
A$celltype[A$celltype %in% c('B cell','B_cell')] <- 'B cell'
A$harmony_condition <- NA
A$harmony_condition[A$orig.ident %in% c('COV007','COV034','COV036')] <- 'severe'
A$harmony_condition[A$orig.ident %in% c('COV002','COV004','COV012')] <- 'moderate'
A$harmony_condition[A$Ventilated %in% c('Healthy')] <- 'Healthy'
A$harmony_condition[A$Ventilated %in% c('NonVent')] <- 'moderate'
A$harmony_condition[A$Ventilated %in% c('Vent')] <- 'severe'
A$harmony_condition[A$orig.ident %in% c("Macrophage",
                                        "NK cell","CD4+ T-cell",
                                        "CD8+ T-cell" ,"B cell" , 
                                        "cDC" )] <- 'Healthy'

A <- RunHarmony(A,group.by.vars = 'Experiment')
A <- RunUMAP(A,reduction = 'harmony',dims = 1:20)

setwd("Z:/Cailab/Qian_writting/pyroptosis_covid/Integrate")
A <- readRDS("Z:/Cailab/Qian_writting/pyroptosis_covid/Integrate/Data/Integrated.rds")
################ P1A1 ######################################
P1A1 <- DimPlot(A,group.by = 'celltype',label = T,repel = T,
                cols = c('#d53e4f','#f46d43','#fdae61','#fee08b','#e6f598',
                         '#abdda4','#66c2a5','#3288bd'))+xlab(' ')+ylab(' ')+
  ggtitle('Integrated')
P1B1 <- DimPlot(A,group.by = 'harmony_condition',cols = c("#1b9e77","#7570b3", "#d95f02"))+xlab(' ')+ylab(' ')+ggtitle('Integrated')
P1B1
################ P1A2 ######################################
PBMC <- subset(A,cells = colnames(A)[colnames(A) %in% colnames(blish)])
P1A2 <- DimPlot(PBMC,group.by = 'celltype',label = T,repel = T,
                cols = c('#d53e4f','#f46d43' ,'#fee08b' ,'#abdda4','#66c2a5','#3288bd'))+xlab(' ')+ylab(' ')+ggtitle('PBMC')+NoLegend()+theme(panel.grid = element_blank())+
  theme(axis.text = element_blank())+
  theme(axis.ticks = element_blank())+
  theme(panel.border = element_blank())+
  theme(axis.line = element_blank())

P1B2 <- DimPlot(PBMC,group.by = 'harmony_condition',cols = c("#1b9e77","#7570b3", "#d95f02"))+NoLegend()+xlab(' ')+ylab(' ')+
  ggtitle('PBMC')+NoLegend()+theme(panel.grid = element_blank())+
  theme(axis.text = element_blank())+
  theme(axis.ticks = element_blank())+
  theme(panel.border = element_blank())+
  theme(axis.line = element_blank())
P1B2
################ P1A3 ######################################
BALF <- subset(A,cells = colnames(A)[colnames(A) %in% colnames(BALF)])
P1A3 <- DimPlot(BALF,group.by = 'celltype',label = T,repel = T,
                cols = c(        '#f46d43','#fdae61','#fee08b',
                                 '#e6f598',      '#66c2a5','#3288bd'))+
                  xlab(' ')+ylab(' ')+ggtitle('BALF')+NoLegend()+theme(panel.grid = element_blank())+
  theme(axis.text = element_blank())+
  theme(axis.ticks = element_blank())+
  theme(panel.border = element_blank())+
  theme(axis.line = element_blank())


P1B3 <- DimPlot(BALF,group.by = 'harmony_condition',cols = c( "#7570b3", "#d95f02"))+
  xlab(' ')+ylab(' ')+ggtitle('BALF')+theme(panel.grid = element_blank())+NoLegend()+
  theme(axis.text = element_blank())+
  theme(axis.ticks = element_blank())+
  theme(panel.border = element_blank())+
  theme(axis.line = element_blank())
P1B3 
################ P1A4 ######################################
healthy_lung <- subset(A,cells = colnames(A)[colnames(A) %in% colnames(healthy_lung)])
P1A4 <- DimPlot(healthy_lung,group.by = 'celltype',label = T,repel = T,
                cols = c('#d53e4f','#f46d43','#fdae61', '#e6f598','#3288bd'))+xlab(' ')+ylab(' ')+NoLegend()+
  ggtitle('Healthy lung')+NoLegend()+theme(panel.grid = element_blank())+
  theme(axis.text = element_blank())+
  theme(axis.ticks = element_blank())+
  theme(panel.border = element_blank())+
  theme(axis.line = element_blank())

P1B4 <- DimPlot(healthy_lung,group.by = 'harmony_condition',
                cols = c("#1b9e77"))+xlab(' ')+ylab(' ')+
  ggtitle('Healthy lung')+NoLegend()+theme(panel.grid = element_blank())+
  theme(axis.text = element_blank())+
  theme(axis.ticks = element_blank())+
  theme(panel.border = element_blank())+
  theme(axis.line = element_blank())
P1B4

################ P1B ######################################
plots <- VlnPlot(A, features = c("CASP1",'CASP4','PYCARD')
                 ,split.by = "harmony_condition",group.by = 'celltype',
                 pt.size = 0, combine = F,cols = c("#1b9e77","#7570b3", "#d95f02"))
P1C1 <- wrap_plots(plots = plots, ncol = 3)
wrap_plots(plots = plots, ncol = 3)


plots <- VlnPlot(A, features = c('IL1B','NINJ1','GSDMD'),split.by = "harmony_condition",group.by = 'celltype',
pt.size = 0, combine = F,cols = c( "#1b9e77","#7570b3", "#d95f02"))
P1C2 <- wrap_plots(plots = plots, ncol = 3)

plots <- VlnPlot(A, features = c('TLR2','NLRC5','PPIA'),split.by = "harmony_condition",group.by = 'celltype',
                 pt.size = 0, combine = F,cols = c( "#1b9e77","#7570b3", "#d95f02"))
P1C3 <- wrap_plots(plots = plots, ncol = 3)
P1C3

plots <- VlnPlot(A, features = c('TLR2','NLRC5','PPIA'),split.by = "harmony_condition",group.by = 'celltype',
                 pt.size = 0, combine = F,cols = c( "#1b9e77","#7570b3", "#d95f02"))


png('Figure/Fig1.png',res = 600,height = 11000,width = 9000)
(P1A1|P1A2|P1A3|P1A4) /
  (P1B1|P1B2|P1B3|P1B4)/
P1C1/P1C2/P1C3
dev.off()

#saveRDS(A,'Data/Integrated.rds')
